from utils import *
import random

def build_tmp_file_im(file_path, marker_mask, individual_mask):
    structure_data = get_structure_data(file_path)
    indiv_list = get_indiv_list(file_path)


    # Population data (+ applying the mask)
    pop_data = get_pop_data(file_path)
    pop_data = [p for p, m in zip(pop_data, individual_mask) if m]

    # Loci data (+ apply the mask)
    loci_list = get_loci_list(file_path)
    loci_list = [l for l, m in zip(loci_list, marker_mask) if m]

    # Indiv list (+ apply the mask)
    indiv_list = [i for i,m in zip(indiv_list, individual_mask) if m]

    # Genotype data (+individual mask + marker mask)
    genotypes = structure_data[1:]
    genotypes = [g for g,m in zip(genotypes, individual_mask) if m]
    genotypes = [[c for c,m in zip(i, marker_mask) if m] for i in genotypes]
    genotypes = [[f"{c[0]}\t{c[1]}" for c in i] for i in genotypes]

    stru_header = "\t".join(loci_list)
    stru_content = [i + "\t" + "\t".join(p) + "\t" + "\t".join(g) for i,p,g in zip(indiv_list, pop_data, genotypes)]

    file_content = stru_header + "\n" + "\n".join(stru_content)
    slug = str(hex(abs(hash(",".join(structure_data[0]))))).replace("0x", "")
    open(f"tmp/{slug}.stru", "w", newline='\n').write(file_content)
    return f"tmp/{slug}.stru"



def calculate_Fst_im(file_path, marker_mask, individual_mask):
    file_path = build_tmp_file_im(file_path, marker_mask, individual_mask)
    n_loci = str(marker_mask.count(True))
    n_indiv = str(individual_mask.count(True))
    # Run the process
    res = subprocess.run([RSCRIPT_PATH, "R/4_Fst_test_2.R", file_path, n_loci, n_indiv], env=env, capture_output=True, text=True)
    output = res.stdout
    output = output.split("\n\n")[1:3]
    output = [a.replace("[1] ", "") for a in output]
    output = [a.split("\n") for a in output]
    output = {a[0][1:]: float(a[1]) for a in output}
    return output

def calculate_Fst_fst_im(file_path, marker_mask, individual_mask):
    return calculate_Fst_im(file_path, marker_mask, individual_mask)["FST"]

def calculate_Fst_fis_im(file_path, marker_mask, individual_mask):
    return calculate_Fst_im(file_path, marker_mask, individual_mask)["FIS"]

def calculate_similarity_im(file_path, marker_mask, individual_mask):
    tmp_file_path = build_tmp_file_im(file_path, marker_mask, individual_mask)
    n_loci = marker_mask.count(True)
    agd = calculate_average_genetic_distance(tmp_file_path, [True]*n_loci)
    return agd

def calculate_ho_im(file_path, marker_mask, individual_mask):
    tmp_file_path = build_tmp_file_im(file_path, marker_mask, individual_mask)
    n_loci = marker_mask.count(True)
    ho = calculate_ho(tmp_file_path, [True]*n_loci)
    return ho

def calculate_he_im(file_path, marker_mask, individual_mask):
    tmp_file_path = build_tmp_file_im(file_path, marker_mask, individual_mask)
    n_loci = marker_mask.count(True)
    he = calculate_he(tmp_file_path, [True]*n_loci)
    return he

def select_set_random(loci_list, n):
    random.seed()
    random_set = random.sample(loci_list, int(n))
    random_mask = [l in random_set for l in loci_list]
    return random_mask

def select_set_pic(loci_list, n):
    # Read the pic sorted loci
    pic_tsv = open("data/PIC.tsv").read().replace('Overall', '').replace('mean', '')
    pic_tsv = pic_tsv.replace('NA', '0.00')
    pic_tsv = pic_tsv.strip().split("\n")
    pic_tsv = [a.split("\t") for a in pic_tsv]

    ordered = [(a,float(b)) for a,b in pic_tsv if len(a) != 0]
    ordered.sort(key=lambda x:x[1], reverse=True)
    ordered = ordered[:n]
    ordered = [a for a,b in ordered]

    output = [l in ordered for l in loci_list]
    return output

def select_set_pic_aco(loci_list, n):
    pic_aco_data = pd.read_csv("data/marker-set/aco-pic.csv", index_col=0)
    pic_aco_data = dict(pic_aco_data[str(n)])
    pic_aco_data = [(m,pic_aco_data[m]) for m in pic_aco_data]
    pic_aco_data.sort(key=lambda x: x[1], reverse=True)
    marker_set = [a for a,b in pic_aco_data[:n]]
    marker_mask = [l in marker_set for l in loci_list]
    return marker_mask
def select_set_aco(loci_list, n):
    aco_data = pd.read_csv("data/marker-set/aco.csv", index_col=0)
    aco_data = dict(aco_data[str(n)])
    aco_data = [(m, aco_data[m]) for m in aco_data]
    aco_data.sort(key=lambda x:x[1], reverse=True)
    marker_set = [a for a,b in aco_data[:n]]
    marker_mask = [l in marker_set for l in loci_list]
    return marker_mask

def calculate_ne(file_path, mask):
    structure_data = get_structure_data(file_path)
    loci_list = get_loci_list(file_path)
    loci_list = [l for l, m in zip(loci_list, mask) if m]

    # Here the mask will be applied
    structure_data[1:] = [[c for c, m in zip(r, mask) if m] for r in structure_data[1:]]
    structure_data[0] = loci_list

    # Only one allele expected
    if len(loci_list) != 1:
        raise "One allele expected"

    # Here the formula will be applied
    # https://cropgenebank.sgrp.cgiar.org/images/file/learning_space/molecular_markers/volume2/04_Measures.pdf
    output = 0
    locus = loci_list[0]
    genotypes = [r[0] for r in structure_data[1:]]
    genotypes_flat = list(sum(genotypes, ()))
    allele_list = [a for a in list(set(sum(genotypes, ()))) if a != missing_data]
    for i, allele in enumerate(allele_list):
        p = genotypes_flat.count(allele) / len(genotypes_flat)
        output += p*p
    output = 1 / output
    return output




def calculate_ne_im(file_path, marker_mask, indiv_mask):
    tmp_file_path = build_tmp_file_im(file_path, marker_mask, indiv_mask)
    n_loci = marker_mask.count(True)
    ne = calculate_ne(tmp_file_path, [True]*n_loci)
    return ne

def calculate_na(file_path, mask):
    structure_data = get_structure_data(file_path)
    loci_list = get_loci_list(file_path)
    loci_list = [l for l, m in zip(loci_list, mask) if m]

    # Here the mask will be applied
    structure_data[1:] = [[c for c, m in zip(r, mask) if m] for r in structure_data[1:]]
    structure_data[0] = loci_list

    # Only one allele expected
    if len(loci_list) != 1:
        raise "One allele expected"

    plain_alleles = sum([list(a[0]) for a in structure_data[1:]], [])
    plain_alleles = [a for a in plain_alleles if a != -9]
    alleles = set(plain_alleles)
    na = len(alleles)
    return na

def calculate_na_im(file_path, marker_mask, indiv_mask):
    tmp_file_path = build_tmp_file_im(file_path, marker_mask, indiv_mask)
    n_loci = marker_mask.count(True)
    na = calculate_na(tmp_file_path, [True]*n_loci)
    return na

def calculate_matrix_im(file_path, marker_mask, individual_mask):
    file_path = build_tmp_file_im(file_path, [True]*28, individual_mask)
    return get_matrix(file_path, marker_mask)

SELECTION_METHODS = {
    "pic-aco": select_set_pic_aco,
    "random": select_set_random,
    "pic": select_set_pic,
    "aco": select_set_aco
}