options(
  repos = c(
    zkamvar = "https://zkamvar.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
)

install.packages("vctrs") # Maybe need to run this before adegenet
install.packages("adegenet")
install.packages("poppr")
install.packages("ape")
install.packages("adegenet")
install.packages("reshape2")



library("ape")
library("adegenet")
library("reshape2")
# No polysat

setwd("C:\\Users\\KU\\Desktop\\AGB\\Ryan\\projects\\061-MS-OPTIMIZATION")
install.packages("remotes")
remotes::install_github("nikostourvas/PopGenUtils", force=TRUE)
library("PopGenUtils")

genotype_data = "data/3_thai_cat/thai_cat_clean.str"
tsv_output = "data/3_thai_cat/PIC.tsv"
n_ind = 184 # TODO: this should be loaded automatically using python utilities
n_loc = 15 # TODO: this should be loaded automatically


obj = read.structure(genotype_data, n.ind=n_ind, n.loc=n_loc, onerowperind=TRUE, col.lab=1, col.pop=2, row.marknames=1, NA.char="-9", ask=FALSE)
pic = pic_calc(obj, ploidy=2, bypop=FALSE)
write.table(pic, file=tsv_output, sep="\t", row.names = TRUE, quote=FALSE)

