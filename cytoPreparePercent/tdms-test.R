# Test
file = "TDMS_format_sample.txt" ; sep = "\t" ; nrows = 20
file = "TDMS_0x1.txt"
file = "TDMS_2_0.txt"
file = "TDMS_2_3.txt"


file = "TDMS_ultra.txt"
file = "panel_plus_ultra.txt"

# File base
filen = gsub("\\..{3,4}$", "", file)

# Read a file
tdms = read.tdms(file)
lapply(tdms, dim)
tdms.head(tdms)

# Direct output
write.tdms(tdms, sprintf("%s-copy.txt", filen))

# Check color ramp
library(RColorBrewer)
col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")), space = "Lab")(100)
image(matrix(seq(0,1,length.out = 100),nr=100), col=col)

# Heatmap view
heatmap(tdms$exprs, Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)

# Transpose
cc = tdms.t(tdms)
heatmap(cc$exprs, Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)
write.tdms(cc, sprintf("%s-transpose.txt", filen))

# Transform
cc = tdms.transform(tdms, func = "none")
head(cc$exprs)
heatmap(cc$exprs, Rowv = NA, Colv = NA, scale = "none", col = col, revC = TRUE)

# Check param is negative
cc = tdms.transform(tdms, func = "none", param = -1) # ignore param if "none"
cc = tdms.transform(tdms, func = "log2", param = -1) # correct param
cc = tdms.transform(tdms, func = "asinh", param = -1) # correct param

# Check exprs contains negative

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


tdms = read.tdms("ex_panel.txt")

cc = tdms.transform(tdms, func = "log2", param = -1) # correct param
aa = tdms.center(cc, meth="median")
bb = tdms.scale(aa, method = "sd", group = "Strain")
round(aa$exprs[14:16,],3)
round(bb$exprs[14:16,],3)
round(tdms$exprs[14:16,],3)
dd = tdms.center.ref(bb, method = "median", grouping = "Strain", group = "S")
round(dd$exprs[14:16,],3)

dd = tdms.center.ref(bb, method = "none", grouping = "Strain", group = "S")



cc$exprs
View(cc$exprs)
