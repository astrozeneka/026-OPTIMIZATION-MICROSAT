import numpy as np
import argparse

# User argument parameter instead "data/EMA_pop.str"
argparser = argparse.ArgumentParser()
argparser.add_argument('--input', help='Input file as a distance matrix')
args = argparser.parse_args()

if __name__ == '__main__':
    data = open(args.input).read().strip().split("\n")
    data = [a.split(" ") for a in data]

    for i,row in enumerate(data):
        for j, cell in enumerate(row):
            if "." in cell:
                try:
                    data[i][j] = str(int(np.round(float(cell), 0)))
                except:
                    pass

    open(args.input.replace(".str", "") + "_clean.str", "w").write("\n".join([" ".join(a) for a in data]))
    print()