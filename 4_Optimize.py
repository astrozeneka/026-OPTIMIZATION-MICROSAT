
import argparse
import os
import subprocess
import numpy as np
from numpy.random import choice as np_choice

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", type=str, help="Path to the structure input file") # Generated with GenAlex
parser.add_argument("-m", "--markers", type=int, help="The quantity of markers to be used") # Will be varied from m to M
parser.add_argument("-s", "--slug", type=int, help="Replicate Slug id (eg. 1, 2, 3, 4)")
parser.add_argument("-e", "--epochs", type=int, help="Number of epochs")

# Newly added arguments for script
parser.add_argument("-l", "--loci", type=int, default=22, help="The total number of loci (same parameter as structure software)")
parser.add_argument("-i", "--individual", type=int, default=50,  help="The total number of individual (same parameter as we used in structure software)")
parser.add_argument("-b", "--binary", type=str, default="C:\\Program Files\\R\\R-4.2.2\\bin\\Rscript.exe", help="The Rscript binary executable location")
parser.add_argument("-a", "--ants", type=int, default=50, help="Ants number parameter")

# Here, the initial state of the pheromones can be set by the user
parser.add_argument('-p', '--pheromones', type=str, default=None, help="The .npy file specifying the initial state of the pheromones")
parser.add_argument('-n', '--noise', type=float, default=0, help="The noise applied to the initial pheromone (for Primed run only)")
args = parser.parse_args()

#TRUNCATED_INDIVIDUALS = 50
TRUNCATED_INDIVIDUALS = args.individual
#TRUNCATED_MARKERS = 22
TRUNCATED_MARKERS = args.loci
#RSCRIPT_PATH = "C:\\Program Files\\R\\R-4.2.2\\bin\\Rscript.exe" # Make in .env file
RSCRIPT_PATH = args.binary

GLOBAL_REFERENCE_MATRIX = None # A global variable

def get_loci_list(file_path):
    structure_data = open(file_path).read().split("\n")
    structure_data = [a.split() for a in structure_data]
    structure_data = structure_data[:1+TRUNCATED_INDIVIDUALS]
    structure_data[0] = structure_data[0][:TRUNCATED_MARKERS]
    structure_data[1:] = [a[:3 + TRUNCATED_MARKERS * 2] for a in structure_data[1:]]

    loci_list = structure_data[0]
    return loci_list

def get_matrix_from_structure(file_path, marker_array_mask):
    # 1. Load and truncate the data
    structure_data = open(file_path).read().split("\n")
    structure_data = [a.split() for a in structure_data]
    structure_data = structure_data[:1+TRUNCATED_INDIVIDUALS]
    structure_data[0] = structure_data[0][:TRUNCATED_MARKERS]

    individual_list = [a[0] for a in structure_data[1:]]

    # Apply the array mask [[i*2, i*2+1] for i,m in enumerate(marker_array_mask) if m]
    column_list = [0, 1, 2] + sum([[3 + i*2, 3 + i*2+1] for i,m in enumerate(marker_array_mask) if m], [])

    structure_data[0] = [s for s,b in zip(structure_data[0], marker_array_mask) if b]
    structure_data[1:] = [[b for i,b in enumerate(a) if i in column_list] for a in structure_data[1:]]
    slug = str(hex(abs(hash(",".join(structure_data[0]))))).replace("0x", "")
    with open(f"tmp/{slug}.stru", "w") as f:
        f.write("\n".join(["\t".join(a) for a in structure_data]))
    print("Truncated data written to file")

    # 1.b Calculate distance matrix of this first
    n_indiv = len(individual_list)
    n_loci = len(structure_data[0])
    env = os.environ.copy()
    env["R_LIBS"] = "r_lib"
    subprocess.run([RSCRIPT_PATH, "Generate_individual_matrix.R", slug, str(n_indiv), str(n_loci)], env=env)
    command_outfile = open(f"tmp/{slug}.out").read().strip().split("\n")
    command_outfile = [a.split("\t") for a in command_outfile]
    # 2. Initialize the ACO distances
    distance_matrix = np.zeros((n_indiv, n_indiv))
    for row in command_outfile[1:]:
        ind1_index = individual_list.index(row[1])
        ind2_index = individual_list.index(row[2])
        row[3] = row[3].replace('H', '') # Unexpectd error
        distance_matrix[ind1_index,ind2_index] = float(row[3])
    os.remove(f"tmp/{slug}.out")
    os.remove(f"tmp/{slug}.stru")
    return distance_matrix # only the distance matrix will be found

class ACO:
    def __init__(self, loci_list, ants, n_best, iterations, decay=0.90, alpha=0.7, beta=1):
        # Contrary to the shortest path problem of NxN
        # This loci selection will take only Nx2 (much more quicker to optimize)
        # No need to use distances
        self.ants = ants
        self.n_best = n_best
        self.iterations = iterations
        self.decay = decay
        self.alpha = alpha # Attractiveness
        self.beta = beta # Visibility

        # Computed values (this value will be modified for the initial proportion)
        # No need to put complex algorithm for managing exact number of microsatellites
        # The report should include a explaining graph and p-value table
        # For supporting the efficience of the algorithm
        if args.pheromones is not None:
            self.pheromones = np.load(args.pheromones)
            self.pheromones = self.pheromones + np.random.normal(loc=0, scale=args.noise, size=self.pheromones.shape) + 0.2
            print(f"Phenotype data loaded from {args.pheromones}")
        else:
            self.pheromones = np.ones((len(loci_list), 2)) / len(loci_list) # Each loci are considered as a node
            if args.markers < 10:
                alpha = (14 / (args.markers))
                self.pheromones[:, 1] /= alpha
            print(f"Use default parameter for pheromones")
        self.all_inds = range(len(loci_list)) # Each loci are considered as a node

    def run(self):
        most_accurate_mset = None
        all_time_most_accurate_set = ("placeholder", np.inf)
        for i in range(self.iterations):
            print(f"Epoch : {i}")
            all_msets = self.get_all_msets()
            self.spread_pheromone(all_msets, self.n_best, most_accurate_mset)
            most_accurate_mset = min(all_msets, key=lambda x:x[1])
            print(most_accurate_mset)
            if most_accurate_mset[1] < all_time_most_accurate_set[1]:
                all_time_most_accurate_set = most_accurate_mset
            self.pheromones = self.pheromones * self.decay
        return all_time_most_accurate_set

    def spread_pheromone(self, all_mset, n_best, most_accurate_mset):
        sorted_msets = sorted(all_mset, key=lambda x:x[1])
        for mset, dist in sorted_msets[:n_best]:
            for loci, boole in enumerate(mset):
                self.pheromones[loci, boole] += 1.0 / dist  # This formula can be updated


    def get_all_msets(self):
        global GLOBAL_REFERENCE_MATRIX
        all_msets = []
        for i in range(self.ants):
            mset = self.gen_mset()
            all_msets.append((mset, self.compute_euclidian_distance(mset, GLOBAL_REFERENCE_MATRIX)))
        return all_msets

    def compute_euclidian_distance(self, mset, ref_matrix):
        marker_mask = [{0:False, 1:True}[a] for a in mset]
        matrix = get_matrix_from_structure(args.file, marker_mask)
        mat = (matrix / matrix.mean()) - (ref_matrix / ref_matrix.mean())

        dist = np.absolute(mat).mean() # Unit in percentage
        return dist # The objective of the function is to reduce
        # return np.linalg.norm(matrix - ref_matrix) # The aim is to minimize this value

    def gen_mset(self): # No need to put start
        mset = [] # no need to put visited and no need to manage start
        _phes = self.pheromones
        for i in self.all_inds: # TODO: Constraints on the number of loci will be put here
            next = self.pick_next(_phes[i], i)
            mset.append(next)
        if sum(mset) > args.markers or sum(mset) == 0: # Maximum of 3 marker will be selected
            return self.gen_mset() # Recursive infinite loop
        #if sum(mset) != args.markers: # Running this algorithm is not possible due to the exponential number of possibilities
        #    return self.gen_mset()
        return mset # Much more easy than the shortest path problem

    def pick_next(self, pheromone, index): # No need to manage distance, no need to manage visited
        # No need to use beta also
        possibilities = range(pheromone.shape[0])
        row = pheromone ** self.alpha # A simplified version
        norm_row = row / row.sum()
        choice = np_choice(possibilities, 1, p=norm_row)[0]
        return choice


def aco_run_on_markers(file_path):
    loci_list = get_loci_list(args.file)
    instance = ACO(loci_list, args.ants, 1, args.epochs)
    optimal_combination, distance = instance.run()
    # The data to be saved will be instance.pheromones$
    np.save(f"data/aco-phe/phe-{args.markers}-{args.slug}.npy", instance.pheromones)

if __name__ == '__main__':
    loci_list = get_loci_list(args.file)
    GLOBAL_REFERENCE_MATRIX = get_matrix_from_structure(args.file, [True]*len(loci_list))
    # s_matrix = get_matrix_from_structure(args.file, [True, True, True, True, False])
    print("Done")
    # Benchmarking
    aco_run_on_markers(args.file)

    print()
