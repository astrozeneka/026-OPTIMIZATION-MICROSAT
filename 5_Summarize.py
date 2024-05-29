
from glob import glob
import numpy as np
import pandas as pd

from utils import get_loci_list, get_indiv_list
from utils_im import calculate_similarity_im


import argparse
argparser = argparse.ArgumentParser()
argparser.add_argument('--input', help='Input file as a distance matrix')
argparser.add_argument('--slug', type=str, default='5')
args = argparser.parse_args()

if __name__ == '__main__':
    loci_list = get_loci_list(args.input)
    indiv_list = get_indiv_list(args.input)
    phe_files = glob(f"data/aco-phe/phe-*-{args.slug}.npy")

    error_df = pd.DataFrame(columns=["Error", "Loci"], index=(range(2, len(loci_list))))

    optimized_list = {a: [] for a in range(2, len(loci_list))}
    for file in phe_files:
        index = int(file.replace("data/aco-phe\\phe-", "").replace(f"-{args.slug}.npy", ""))
        loci_mask = [a[1]>a[0] for a in phe_data]
        optimized_loci = [l for m, l in zip(loci_mask, loci_list) if m]
        optimized_list[index] = optimized_loci

        # 1. Calculate the AGD Loss
        agd_loss = calculate_similarity_im(args.input, loci_mask, [True]*len(indiv_list))
        error_df["Error"][index] = agd_loss
        error_df["Loci"][index] = ", ".join(optimized_loci)
    print()

    outfile = args.input.replace(".str", "") + "_report.csv"
    error_df.to_csv(outfile, index=True)



    print()