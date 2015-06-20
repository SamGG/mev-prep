m.input = 0

transf.meth = ""
trans.log2 = 0
trans.asinh = 1

center.to = ""
center.meth = ""
# TODO: define ref.col.id
# TODO: use center.trim
# TODO: typical uses with advantages
# TODO: add reference to workflow: reference, anonymous, illustrative, ignored/removed
# TODO: input TSV with annotation
# TODO: output TSV with annotation
# TODO: minimal # samples to form a group
# TODO: add alerts
# TODO: level of the reference group for final centering
# TODO: centering: median of group medians

# Remove reference if asked
m.input = m.input[ , -ref.col.id, drop = F]

# Transform

if (transf.meth == "none") {
  m.trans = m.input
} else (transf.meth == "log2") {
  m.trans = log2( m.input + trans.log2 )
} else (transf.meth == "asinh") {
  m.trans = asinh( m.input / trans.asinh )
  # TODO: scale to mimic log2 at 50%?
} else {
  stop(paste("Unknown transf.meth:", transf.meth))
}

# Centering
# TODO: centering to the reference at the final step, insuring a real zero?
# TODO: median of group median

center.to = switch(reference,
                   reference = "line.but.ref",
                   anonymous = "line",
                   illustrative = "line.but.ref",
                   ignored = "line")

if (center.to == "ref.group") {  # by line
  m.ref = m.trans[ , ref.col.id, drop=F]
} else if (center.to == "line") {  # center
  m.ref = m.trans[ , ]
} else if (center.to == "line.but.ref") {  # center
  m.ref = m.trans[ , -ref.col.id, drop=F]
} else {
  stop(paste("Unknown center.to:", center.to))
}

if (center.meth == "median") {
  centers = apply( m.ref, 1, median, na.rm=T)
} else (center.meth == "mean") {  # center.meth.mean
  centers = apply( m.ref, 1, mean, na.rm=T, trim=center.trim)
} else {
  stop(paste("Unknown center.meth:", center.meth))
}

m.center = m.trans - centers  # recycling by column is nice

# Standardize
# TODO?: pre-trimming
# TODO: winsoring?

# min.max
# Trim, or percentile
m.min = apply(m.center, 1, min, na.rm=T)
m.max = apply(m.center, 1, max, na.rm=T)
sweep(m.center - m.min, 1, 2 / (m.max - m.min), FUN = "*" )

# Percentile
# Lower and higher percentiles suffer also from instability
m.qmin = apply(m.center, 1, quantile, probs = p.percentile, na.rm=T)
m.qmax = apply(m.center, 1, quantile, probs = 1 - p.percentile, na.rm=T)
sweep(m.center - m.qmin, 1, 2 / (m.qmax - m.qmin), FUN = "*" )

# sd.line == 1
# Over line scaling mix intra and inter group variance
sd.line = sd
getAnywhere("sd")
# sd.group == 1 or 1/2.355
# Scaling versus average intra group variance

# MAD
# http://en.wikipedia.org/wiki/Trimmed_estimator

# IQR
# http://en.wikipedia.org/wiki/Interquartile_range
# Robust standard deviation could be estimated with IQR * 1.349

# Trimming = removing
# Winsoring = clipping
# http://en.wikipedia.org/wiki/Winsorising
# http://www.personality-project.org/r/psych/help/winsor.html



qnorm(.25)
pnorm(2) - pnorm(-2)
pnorm(-2) * 2
pnorm(-.6745) * 2


# For more applications of robust procedure, see RR Wilcox package
# http://www.personality-project.org/r/psych/help/winsor.html

# Scale by line
# https://github.com/dgrapov/DeviumWeb/blob/master/R/Devium%20data%20transformation.r
# could add level and vast scaling http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1534033/table/T1/

# Color scaling

# Inform

# color scale

# hEATMAP
# pheatmap
# NMF
# heatmap_2
# https://github.com/dgrapov/devium/blob/master/R/Devium%20Cluster%20Analysis.r

# How Do I Draw A Heatmap In R With Both A Color Key And Multiple Color Side Bars?
# https://www.biostars.org/p/18211/
# https://gist.github.com/nachocab/3853004

# How to assign your color scale on raw data in heatmap.2()
# http://stackoverflow.com/questions/20535635/how-to-assign-your-color-scale-on-raw-data-in-heatmap-2



# Trimming from mean source code
# getAnywhere(mean.default)
n <- length(x)
if (trim > 0 && n) {
  if (is.complex(x)) 
    stop("trimmed means are not defined for complex data")
  if (any(is.na(x))) 
    return(NA_real_)
  if (trim >= 0.5) 
    return(stats::median(x, na.rm = FALSE))
  lo <- floor(n * trim) + 1
  hi <- n + 1 - lo
  x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
}

# Colorbar
# http://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette



distEisen <- function(x, use = "pairwise.complete.obs") {
  co.x <- cor(x, use = use)
  dist.co.x <- 1 - co.x
  return(as.dist(dist.co.x))
}


# http://ma-bioinformatics.googlecode.com/svn/trunk/R/