from collections import Counter



if __name__ == '__main__':
    data = open("data/EMA_pop_clean.str").read().strip().split("\n")
    data = [a.split() for a in data]
    loci_list = data[0]
    data = [a[3:] for a in data[1:]]
    data = [{l:(a[j*2], a[j*2+1]) for j,l in enumerate(loci_list)} for a in data]

    for l in loci_list:
        print("")
        print(l)

        column = [a[l] for a in data]
        genotypes = list(set(column))
        frequency = [
            [genotype,
                len([a for a in column if a == genotype]) / len(column)]
            for genotype in genotypes
        ]
        frequency.sort(key=lambda x:x[1], reverse=True)
        print("\n".join([f"{','.join(a)},{b}" for a,b in frequency]))

    print()
