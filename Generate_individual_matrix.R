

#list.of.packages = c("ggplot2", "ape")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]


#install.packages("ape", repos='http://cran.us.r-project.org', lib='Rlib')
#install.packages("ggplot2", repos='http://cran.us.r-project.org', lib='Rlib')
#install.packages("adegenet", repos='http://cran.us.r-project.org', lib='Rlib')
#install.packages("eimpute", repos='http://cran.us.r-project.org', lib='Rlib')
library("ape")
library("adegenet") # Maybe don't need to make pcoa and ggplot
library("eimpute")
library("reshape2")

args = commandArgs(trailingOnly = TRUE)
slug = args[1]
n_indiv = as.integer(args[2])
n_loci = as.integer(args[3])

input_file = paste("tmp/", slug, ".stru", sep="")
output_file = paste("tmp/", slug, ".out", sep="")
print(input_file)
obj = read.structure(input_file, n.ind=n_indiv, n.loc=n_loci, onerowperind=TRUE, col.lab=1, col.pop=2, row.marknames=1, NA.char="-9", ask=FALSE)
matrix = dist(obj, method="euclidian", diag=TRUE)
matrix[is.na(matrix)] <- 0.1 # TODO: This can be a potential cause of a bug

# Write it to the output_file
df = melt(as.matrix(matrix), varnames = c("ind1", "ind2"))
write.table(df, file=output_file, sep="\t", quote=FALSE, na="", row.names=TRUE, col.names=TRUE)
