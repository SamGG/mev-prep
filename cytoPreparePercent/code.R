# Test
file = "TDMS_format_sample.txt" ; sep = "\t" ; nrows = 20
file = "TDMS_0x1.txt"
file = "TDMS_2_0.txt"
file = "TDMS_2x3.txt"




library(RColorBrewer)
col = brewer.pal(11, "RdYlBu")
heatmap(tdms$exprs, Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)
heatmap(tdms$exprs[1:9,], Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)

write.tdms(tdms, "toto.txt")

cc = tdms.t(tdms)
write.tdms(cc, "titi.txt")

cc = tdms.summarize(tdms, group = tdms$pData$Strain, probs = c(0.5))
write.tdms(cc, "titi.txt")

cc = tdms.summarize(tdms, group = tdms$pData$Strain, probs = c(0.05, 0.5, 0.95))
write.tdms(cc, "titi.txt")
cc$pData$Strain = cc$pData$id
write.tdms(cc, "titi.txt")

dim(cc$exprs)

dim(s)
tdms = read.tdms(file)
write.tdms(tdms, "toto.txt")

aa = tdms.center(tdms)
aa = tdms.center(tdms, meth="zz")
aa = tdms.center(tdms, meth="percent")
heatmap(aa$exprs, Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)
heatmap(aa$exprs[1:9,], Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)
round(aa$exprs[1:9,], 3)
bb = tdms.scale(aa, meth="percent")
heatmap(bb$exprs, Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)
heatmap(bb$exprs[1:9,], Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)
round(bb$exprs[1:9,], 3)




library(Biobase)
# Read a TDMS and return an eSet
read.tdms.eset = function(file, ...) {
  tdms.list = read.tdms(file, ...)
  tdms.pData = new("AnnotatedDataFrame", data=tdms.list$pData)
  tdms.fData = new("AnnotatedDataFrame", data=tdms.list$fData)
  ExpressionSet(assayData = tdms.list$exprs,
                phenoData = tdms.pData, featureData = tdms.fData)
}
