from utils import *
import pandas as pd
from vars import *

import argparse
argparser = argparse.ArgumentParser()
argparser.add_argument('--input', help='Input file as a distance matrix')
argparser.add_argument('--pic', type=str, help='PIC file')
args = argparser.parse_args()

def generate_phe_map(loci_sorted, n):
    output = np.ones((len(loci_sorted), 2))

    if n == len(loci_sorted):
        output[:, 1] = 1
    if n == len(loci_sorted) // 2:
        output[:, 1] = np.linspace(0, 1, len(loci_sorted))
    if n > len(loci_sorted) // 2:
        # Should be update
        a = (n - (len(loci_sorted) // 2)) / (len(loci_sorted) // 2)
        output[:, 1] = np.linspace(a, 1, len(loci_sorted))
    if n < len(loci_sorted) // 2:
        b = n / (len(loci_sorted) // 2)
        output[:, 1] = np.linspace(0, b, len(loci_sorted))
        print()
    output[:, 0] = 1 - output[:, 1]
    df = pd.DataFrame(output, index=reversed(loci_sorted), columns=[0, 1])

    # Now, the output, should be aligned with the real order value of microsatellites
    # Not ordered like in this picture
    structure_data= get_structure_data(args.input)
    loci_list = structure_data[0]

    output = [(df[0][a], df[1][a]) for a in loci_list]
    output = np.array(output)
    np.save(f"data/aco-phe/ini-{n}.npy", output)
    return output

if __name__ == '__main__':
    loci_sorted = open("data/4_thai_cat_v2/PIC.tsv").read().strip().split("\n")[1:][:-1]
    loci_sorted = [a.replace('NA', '0').split("\t") for a in loci_sorted]
    loci_sorted = [(a, float(b)) for a,b in loci_sorted]
    loci_sorted.sort(key=lambda x:x[1], reverse=True)
    loci_sorted = [a for a,b in loci_sorted]

    generate_phe_map(loci_sorted, 14) # TODO, should be updated
    for i in range(2, 14):
        generate_phe_map(loci_sorted, i)

    print("Done")