# Test
file = "TDMS_format_sample.txt" ; sep = "\t" ; nrows = 20
file = "TDMS_0x1.txt"
file = "TDMS_2_0.txt"
file = "TDMS_2_3.txt"


tdms = read.tdms(file)


library(RColorBrewer)
col = brewer.pal(11, "RdYlBu")
heatmap(tdms$exprs, Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)
heatmap(tdms$exprs[1:9,], Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)

write.tdms(tdms, "toto.txt")

# Transpose
cc = tdms.t(tdms)
write.tdms(cc, "titi.txt")

# Summarize
# median alone
cc = tdms.summarize(tdms, group = tdms$pData$Strain, probs = c(0.5))
cc$pData
write.tdms(cc, "titi.txt")
# median plus 5% quantile
cc = tdms.summarize(tdms, group = tdms$pData$Strain, probs = c(0.05, 0.5, 0.95))
cc$pData
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
