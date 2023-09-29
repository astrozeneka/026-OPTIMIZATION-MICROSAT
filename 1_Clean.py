import numpy as np
FILE_PATH="data/EMA_pop.str"
if __name__ == '__main__':
    data = open("data/EMA_pop.str").read().strip().split("\n")
    data = [a.split(" ") for a in data]

    for i,row in enumerate(data):
        for j, cell in enumerate(row):
            if "." in cell:
                try:
                    data[i][j] = str(int(np.round(float(cell), 0)))
                except:
                    pass

    open(FILE_PATH.replace(".str", "") + "_clean.str", "w").write("\n".join([" ".join(a) for a in data]))
    print()