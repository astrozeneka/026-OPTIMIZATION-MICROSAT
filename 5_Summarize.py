
from glob import glob
import numpy as np
import pandas as pd

from utils import get_loci_list, get_indiv_list
from utils_im import calculate_similarity_im

file_path = "data/EMA_pop_clean.str"
if __name__ == '__main__':
    GLOBAL_FILE_PATH=file_path

    loci_list = get_loci_list(file_path)
    indiv_list = get_indiv_list(file_path)
    phe_files = glob("data/aco-phe/phe-*-3.npy")
    error_df = pd.DataFrame(columns=["Error", "Loci"], index=(range(2, len(loci_list))))

    optimized_list = {a: [] for a in range(2, len(loci_list))}
    for file in phe_files:
        index = int(file.replace("data/aco-phe\\phe-", "").replace("-3.npy", ""))
        phe_data = np.load(file)
        loci_mask = [a[1]>a[0] for a in phe_data]
        optimized_loci = [l for m, l in zip(loci_mask, loci_list) if m]
        optimized_list[index] = optimized_loci

        # 1. Calculate the AGD Loss
        agd_loss = calculate_similarity_im(file_path, loci_mask, [True]*len(indiv_list))
        error_df["Error"][index] = agd_loss
        error_df["Loci"][index] = ", ".join(optimized_loci)
    print()
    error_df.to_csv("data/report.csv", index=True)



    print()