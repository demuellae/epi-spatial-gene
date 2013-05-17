library(biomaRt)
eurexpress = useMart("Eurexpress Biomart", dataset="template")
tissuetab = getBM(attributes = c("emap_id","emap_term"), mart=eurexpress)
intM = read.table(file="matrix.csv", header=TRUE, row.names=2, check.names=FALSE, sep=",")

