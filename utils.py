
import os.path
import argparse
import glob
import subprocess

import numpy as np
import pandas as pd
RSCRIPT_PATH = "C:\\Program Files\\R\\R-4.2.2\\bin\\Rscript.exe" # Make in .env file
env = os.environ.copy()
env["R_LIBS"] = "C:\\Users\\KU\\Desktop\\AGB\\Ryan\\extra-projects\\MARKERS-OPTIMIZATION-ACO\\src\\3_Truncated_Population\\r_lib"

color_palette = ["#2BBCD9", "#E445C6", "#FFA07A", "#F1C232", "#F6647C", "#85506E"]

missing_data = -9

def get_structure_data(file_path):
    structure_data = open(file_path).read().split("\n")
    structure_data = [a.split() for a in structure_data]
    structure_data = structure_data[:]
    structure_data[0] = structure_data[0]
    # structure_data[1:] = [a[:3] for a in structure_data[1:]]
    structure_data[1:] = [list(zip(r[3:][0::2],r[3:][1::2])) for r in structure_data[1:]]
    structure_data[1:] = [[(int(c), int(d)) for c,d in r] for r in structure_data[1:]]
    return structure_data

def get_structure_data_with_pop(file_path):
    structure_data = open(file_path).read().split("\n")
    structure_data = [a.split() for a in structure_data]
    structure_data = structure_data[:]
    structure_data[0] = structure_data[0]
    # structure_data[1:] = [a[:3] for a in structure_data[1:]]

    pop_data = [r[1:3] for r in structure_data[1:]]
    structure_data[1:] = [ list(zip(r[3:][0::2],r[3:][1::2])) for r in structure_data[1:]]
    structure_data[1:] = [[(int(c), int(d)) for c,d in r] for r in structure_data[1:]]
    output = [p+r for p,r in zip(pop_data, structure_data[1:])]
    return output

def get_loci_list(file_path):
    structure_data = get_structure_data(file_path)
    loci_list = structure_data[0]
    return loci_list

def calculate_he(stru_file, mask):
    structure_data = get_structure_data(stru_file)
    loci_list = get_loci_list(stru_file)
    loci_list = [l for l, m in zip(loci_list, mask) if m]

    # Here the mask will be applied
    structure_data[1:] = [[c for c, m in zip(r, mask) if m] for r in structure_data[1:]]
    structure_data[0] = loci_list

    # Using the 2pq formula
    sum_lc = 0
    for l, locus in enumerate(loci_list):
        genotypes = [r[l] for r in structure_data[1:]]
        genotypes_flat = list(sum(genotypes, ()))
        allele_list = [a for a in list(set(sum(genotypes, ()))) if a != missing_data]
        sum_af = 0
        for i, allele in enumerate(allele_list):
            p = genotypes_flat.count(allele) / len(genotypes_flat)
            sum_af += p*p
            # print(f"{allele} : {p}")
        sum_lc += sum_af
    sum_lc /= len(loci_list)
    output = 1 - sum_lc
    return output

def calculate_ho(stru_file, mask):

    structure_data = get_structure_data(stru_file)
    loci_list = get_loci_list(stru_file)
    loci_list = [l for l, m in zip(loci_list, mask) if m]

    # Here the mask will be applied
    structure_data[1:] = [[c for c, m in zip(r, mask) if m] for r in structure_data[1:]]
    structure_data[0] = loci_list

    loci_indivs = [a for a in sum(structure_data[1:], []) if a.count(missing_data) == 0]
    total_hz = len([a for a in loci_indivs if a[0] != a[1]])
    output = total_hz / len(loci_indivs)
    return output


def get_indiv_list(file_path):
    structure_data = open(file_path).read().split("\n")
    structure_data = [a.split() for a in structure_data]
    return [a[0] for a in structure_data[1:]]

def get_pop_data(file_path):
    structure_data = open(file_path).read().split("\n")
    structure_data = [a.split() for a in structure_data]
    return [a[1:3] for a in structure_data[1:]]

def build_tmp_file(file_path, mask):
    structure_data = get_structure_data(file_path)
    indiv_list = get_indiv_list(file_path)
    pop_data = get_pop_data(file_path)
    pop_data = ["\t".join(a) for a in pop_data]

    loci_list = get_loci_list(file_path)
    loci_list = [l for l, m in zip(loci_list, mask) if m]

    # Here the mask will be applied
    structure_data[1:] = [[c for c, m in zip(r, mask) if m] for r in structure_data[1:]]
    structure_data[0] = loci_list

    # Write it to file
    stru_header = "\t".join(structure_data[0])
    stru_content = ["\t".join(["\t".join([str(d) for d in c]) for c in indiv]) for indiv in structure_data[1:]]
    stru_content = ["\t".join(a) for a in zip(indiv_list, pop_data, stru_content)]
    stru_content = "\n".join(stru_content)
    file_content = stru_header + "\n" + stru_content

    slug = str(hex(abs(hash(",".join(structure_data[0]))))).replace("0x", "")
    open(f"tmp/{slug}.stru", "w").write(file_content)
    return f"tmp/{slug}.stru"

def calculate_Fst(file_path, mask):

    # Create the mask file first (same as in the 2_matrix_consecutive file)
    file_path = build_tmp_file(file_path, mask)
    n_loci = str(mask.count(True))

    # Run the process
    res = subprocess.run([RSCRIPT_PATH, "R/4_Fst_test_2.R", file_path, n_loci], env=env, capture_output=True, text=True)
    output = res.stdout
    output = output.split("\n\n")[1:3]
    output = [a.replace("[1] ", "") for a in output]
    output = [a.split("\n") for a in output]
    output = {a[0][1:]: float(a[1]) for a in output}
    return output

def get_matrix(file_path, mask):
    # The same function as in the 2_marker_consecutive function
    structure_data = open(file_path).read().split("\n")
    structure_data = [a.split() for a in structure_data]

    individual_list = [a[0] for a in structure_data[1:]]
    column_list = [0, 1, 2] + sum([[3 + i*2, 3 + i*2+1] for i,m in enumerate(mask) if m], [])
    structure_data[0] = [s for s, b in zip(structure_data[0], mask) if b]
    structure_data[1:] = [[b for i, b in enumerate(a) if i in column_list] for a in structure_data[1:]]
    slug = str(hex(abs(hash(",".join(structure_data[0]))))).replace("0x", "")
    with open(f"tmp/{slug}.stru", "w") as f:
        f.write("\n".join(["\t".join(a) for a in structure_data]))
    print("Truncated data written to file")
    n_indiv = len(individual_list)
    n_loci = len(structure_data[0])
    subprocess.run([RSCRIPT_PATH, "Generate_individual_matrix.R", slug, str(n_indiv), str(n_loci)])
    command_outfile = open(f"tmp/{slug}.out").read().strip().split("\n")
    command_outfile = [a.split("\t") for a in command_outfile]
    distance_matrix = np.zeros((n_indiv, n_indiv))
    for row in command_outfile[1:]:
        ind1_index = individual_list.index(row[1])
        ind2_index = individual_list.index(row[2])
        row[3] = row[3].replace('H', '')  # Unexpectd error
        distance_matrix[ind1_index, ind2_index] = float(row[3])
    os.remove(f"tmp/{slug}.out")
    os.remove(f"tmp/{slug}.stru")

    return distance_matrix

def compute_euclidian_distance(A, B):
    return np.linalg.norm(B - A) # The aim is to minimize this value


FULL_MATRIX = None # used as a singleton
def calculate_similarity(file_path, mask): # DON'T USE THIS FUNCTION
    global FULL_MATRIX
    if FULL_MATRIX is None:
        FULL_MATRIX = get_matrix(file_path, [True]*len(get_loci_list(file_path)))
    MAT = get_matrix(file_path, mask)
    return compute_euclidian_distance(MAT, FULL_MATRIX)

#def compute_euclidian_distance(mset, ref_matrix):
#    marker_mask = [{0:False, 1:True}[a] for a in mset]
#    matrix = get_matrix_from_structure(args.file, marker_mask)
#    return dist # The objective of the function is to reduce

GLOBAL_FILE_PATH = "data/gallus-gallus.str" # IMPORTANT, UPDATE IT ACCORDING TO THE
GLOBAL_FILE_PATH = "data/genotype-data/Naemorhedus-griseus.str"
FULL_MATRIX = None
def calculate_average_genetic_distance(file_path, mask):
    global FULL_MATRIX
    global GLOBAL_FILE_PATH

    MAT = get_matrix(file_path, mask)
    if FULL_MATRIX is None:
        FULL_MATRIX = get_matrix(GLOBAL_FILE_PATH, [True]*28) # Ano
    mat = (MAT / MAT.mean()) - (FULL_MATRIX / FULL_MATRIX.mean()) # Should change the GLOBAL_FILE_PATH
    return np.absolute(mat).mean()
    #print()
    #mat = (matrix / matrix.mean()) - (ref_matrix / ref_matrix.mean())
    #dist = np.absolute(mat).mean() # Unit in percentage