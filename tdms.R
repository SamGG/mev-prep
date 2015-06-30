# Read a TDMS file
# A TDMS is a tabulated text file carrying numerical data as well as
# annotation about objects on row and objects on columns.
# For micro-arrays, columns are usually samples, rows are genes.
# For cytometry, rows are usually populations and data are percentages.
# Default values correspond to no extra annotation: columns are named,
# rows are also named, and data begins at 2nd line and 2nd column.
read.tdms = function(file, sep = "\t", rd1 = NA, cd1 = NA, nrows = 20, ...) {
  # TODO: guess separator and deduce decimal
  # First read of top lines in order to guess header end
  hdr = read.table(file, sep = sep, nrows = nrows, header = FALSE, as.is = T, comment.char = "", ...)
  # Automatic detection: select cells from lower right corner up to an
  # increasing upper left corner and try to convert those columns to numerical
  # values; memorize the highest upper left corner.
  if (is.na(rd1) | is.na(cd1)) {
    rd1.auto = cd1.auto = NA
    for (j in seq(ncol(hdr)-6, 2)) {
      for (i in seq(nrow(hdr)-10, 2)) {
        r = sapply(hdr[i:nrow(hdr),j:ncol(hdr)],
                   function(x) is.numeric(type.convert(x)))
        if (all(r==TRUE)) {
          rd1.auto = i ; cd1.auto = j
        }
        # cat(i, j, r, "\n")
      }
    }
    message(sprintf("Automatic detection: (%d, %d)", rd1.auto, cd1.auto))
    if (is.na(rd1.auto))
      stop("Automatic detection of row data start failed!")
    if (is.na(rd1)) rd1 = rd1.auto
    if (rd1.auto > rd1)
      stop(sprintf("Row data start is greater than given value: %d > %d!"), rd1.auto, rd1)
    if (is.na(cd1.auto))
      stop("Automatic detection of column data start failed!")
    if (is.na(cd1)) cd1 = cd1.auto
    if (cd1.auto > cd1)
      stop(sprintf("Column data start is greater than given value: %d > %d!"), cd1.auto, cd1)
  }
  if (rd1 < 2 | cd1 < 2)
    stop( sprintf("Incorrect start of data: (row,col) = (%d,%d) < (2,2)", rd1, cd1))
  # TODO: skip some columns such as counts, percent, cluster
  # Extract annotations on columns
  ann.col = as.data.frame(t(hdr[seq(1,rd1-1), seq(cd1, ncol(hdr))]))
  colnames(ann.col) = c("Id", hdr[seq(2,rd1-1), cd1-1])
  # Read full file and extract row annotations and data
  dat = read.table(file, sep = sep, skip = rd1-1, header = FALSE, comment.char = "", ...)
  ann.row = as.data.frame(dat[, seq(1, cd1-1)])
  colnames(ann.row) = hdr[1, seq(1, cd1-1)]
  dat = as.matrix(dat[, -seq(1, cd1-1)])
  # if enhanced output as eSet
  if (any(duplicated(ann.row[,1]))) {
    message("row id are duplicated")
  } else {
    rownames(dat) = ann.row[,1]
    rownames(ann.row) = ann.row[,1]
  }
  if (any(duplicated(ann.col[,1]))) {
    message("col id are duplicated")
  } else {
    colnames(dat) = ann.col[,1]
    rownames(ann.col) = ann.col[,1]
  }
  list(exprs = dat, pData = ann.col, fData = ann.row)
}

# Write TDMS file
write.tdms = function(tdms, file, ...) {
  # TODO: string for NA
  # Define the shifted area
  rd0 = ncol(tdms$pData)
  cd0 = ncol(tdms$fData)
  # Define the complete area
  nr = nrow(tdms$exprs) + rd0
  nc = ncol(tdms$exprs) + cd0
  # Setup a string matrix
  s = matrix('', nrow = nr, ncol = nc)
  s[seq(nrow(tdms$exprs))+rd0, seq(ncol(tdms$exprs))+cd0] = tdms$exprs
  # TODO: squeeze digits
  # Add annotations
  s[seq(ncol(tdms$pData)), seq(nrow(tdms$pData))+cd0] = t(tdms$pData)
  s[seq(nrow(tdms$fData))+rd0, seq(ncol(tdms$fData))] = as.matrix(tdms$fData)
  # Add names
  s[seq(ncol(tdms$pData)), cd0] = colnames(tdms$pData)
  s[1, seq(ncol(tdms$fData))] = colnames(tdms$fData)
  # Write to disk
  write.table(s, file, sep = "\t", dec = ".",
              quote = FALSE, row.names = FALSE, col.names = FALSE, ...)
}

# head
tdms.head <- function(tdms, digits = 3, rows = 6, cols = 7) {
  .exprs = tdms$exprs
  rows = seq(min(rows, nrow(.exprs)))
  cols = seq(min(cols, ncol(.exprs)))
  list( exprs = head(round(tdms$exprs[rows, cols], digits)),
        pData = tdms$pData[cols, ],
        fData = tdms$fData[rows, ]
  )
}

# transforming
tdms.transform = function(tdms, func = "none", param = NA) {
  if (func == "none") {
    # TODO: allow negative values?
  } else if (func == "log2") {
    if (is.na(param)) param = 0  # default to 0
    #if (param < 0) param = 0
    tdms$exprs[tdms$exprs <= 0] = NA  # avoid NaN, prefill with NA
    tdms$exprs = log2( tdms$exprs + param )
  } else if (func == "asinh") {
    if (is.na(param)) param = 1  # default to 1
    if (param < 0) param = 1  # avoid opposite
    tdms$exprs[tdms$exprs <= 0] = NA  # avoid NaN, prefill with NA
    tdms$exprs = asinh( tdms$exprs / param )
  } else {
    warning(sprintf("Unknown transformation \"%s\", data are unchanged.", func))
  }
  tdms
}

# centering
tdms.center = function(tdms, method="median", percent = 0.05) {
  # TODO: set and check trimming parameter: if <=1 then percent else count
  # TODO: percent as parameter
  if (!method %in% strsplit("none:mean:median:trimmed:minmax", ":")[[1]]) {
    warning(sprintf("Method '%s' is unknown, switching to median", method))
    method = "median"
  }
  # Check percent
  if (percent < 0)  percent = -percent  # invert negative value
  if (percent > 50) percent = 0.5  # limit up 50%
  if (percent > 1) percent = percent / 100  # convert to 0..1
  if (percent > 0.5) percent = 0.5  # limit up 50%
  # Compute centers
  .centers =
    switch(method,
           none    = rep(0, ncol(tdms$exprs)),
           mean    = apply(tdms$exprs, 1, mean, na.rm = TRUE),
           median  = apply(tdms$exprs, 1, median, na.rm = TRUE),
           trimmed = apply(tdms$exprs, 1, mean, trim = percent, na.rm = TRUE),
           minmax  = apply(tdms$exprs, 1, function(x)
             sum(quantile(x, probs = c(percent, 1- percent), na.rm = TRUE))/2)
           )
  # Do centering
  tdms$exprs = sweep(tdms$exprs, 1, .centers, FUN = "-")
  tdms
}

# Split the matrix using the given factor
# Build a matrix list, one per level
# A reference level could be included, excluded or the only level of interest
tdms.split.matrix.1 = function(m, f, filter = "none", group = "ref", min.size = 1) {
  print("filter")
  print(filter)
  # Prepare a list
  lf <- levels(f)
  y <- vector("list", length(lf))
  names(y) <- lf
  # Split the indices of columns
  ind <- split(seq(ncol(m)), f)
  # Split the matrix according the indices
  for (k in lf)
    y[[k]] <- m[ ,ind[[k]], drop = FALSE]
  # Remove reference
  if (filter == "excl.ref") {
    ind.ref = match(tolower(group), tolower(lf))
    if (all(!is.na(ind.ref))) y = y[-ind.ref] else
      warning(sprintf("Reference '%s' is unknown.", paste(group[is.na(ind.ref)])))
  # Keep reference alone
  } else if (filter == "only.ref") {
    ind.ref = match(tolower(group), tolower(lf))
    if (all(!is.na(ind.ref))) y = y[ind.ref] else
      warning(sprintf("Reference '%s' is unknown.", paste(group[is.na(ind.ref)])))
  # None
  } else if (filter != "none") {
    warning(sprintf("Filter '%s' is unknown.", filter))
    # TODO: return NULL?
  }
  # Return, debug lapply(y, head)
  y
}

# Split the matrix using the given factor
# Build a matrix list, one per level
# A reference level could be included, excluded or the only level of interest
tdms.split.matrix = function(m, grouping, filter = "none", group = "ref", min.size = 1) {
  # Split the indices of columns
  ind <- split(seq(ncol(m)), grouping)
  lgrouping <- names(ind)
  # Filter
  if (filter == "excl" || filter == "only") {
    # Check group name
    ind.ref = match(tolower(group), tolower(lgrouping))
    if (any(is.na(ind.ref)))
      stop(sprintf("Reference '%s' is unknown.", paste(group[is.na(ind.ref)])))
    # Set up filter
    if (filter == "excl") ind = ind[-ind.ref]  # Remove reference
    if (filter == "only") ind = ind[ ind.ref]  # Keep reference only
    # None
  } else if (filter != "none") {
    warning(sprintf("Filter '%s' is unknown.", filter))
    # TODO: return NULL?
  }
  # TODO: Filter according min.size
  # Prepare a list
  y <- vector("list", length(ind))
  names(y) <- names(ind)
  # Split the matrix according the indices
  for (k in names(ind))
    y[[k]] <- m[ ,ind[[k]], drop = FALSE]
  # Return, debug lapply(y, head)
  y
}

# scaling
# if percentile, group is irrelevant
# if sd or mad, group is used as follows
# if group == NA, global scale is computed
# if group != NA, per group is computed
tdms.scale = function(tdms, method = "sd", percent = 0.05, group = NA, filter = "none", min.size = 3) {
  if (!method %in% strsplit("none:sd:mad:minmax", ":")[[1]]) {
    warning(sprintf("Method '%s' is unknown, switching to sd", method))
    method = "sd"
  }

  # TODO: group defined, method=none, error
  # Check grouping
  if (tolower(group) == "none" || is.na(group[1])) {
    group = NA
  } else {
    # named factor
    if (is.character(group) && group[1] %in% colnames(tdms$pData))
      group = factor(tdms$pData[,group])
    # Unrecognized grouping
    if (!is.factor(group) && length(group) != nrow(tdms$pData))
      stop("Grouping is unrecognized")
  }
  if (!is.na(group) && (method %in% c("minmax", "none"))) {
    warning(sprintf("Method '%s' ignores grouping", method))
    group = NA
  }

  # Compute scales, no grouping
  if (length(group) == 1 && is.na(group)) {
    # TODO: set and check trimming parameter: if <=1 then percent else count
    # TODO: trimming, windsoring
    # Check percentile
    if (method == "minmax") {
      if (percent <= 0) percent = 0
      if (percent >= 0.25) percent = 0.25
      probs = c(percent, 1-percent)
    }
    # Compute scaling
    .scales =
      switch(method,
             none   = rep(1, ncol(tdms$exprs)),
             sd     = apply(tdms$exprs, 1, sd, na.rm = TRUE)/2,
             mad    = apply(tdms$exprs, 1, mad, na.rm = TRUE)/2,
             minmax = apply(tdms$exprs, 1, function(x)
               diff(quantile(x, probs = probs, na.rm = TRUE))/2/2.5)
      )
  # Compute scales, grouping
  } else {
    # Split matrix by group
    l = tdms.split.matrix(m = tdms$exprs, grouping = group, filter = filter, min.size = min.size)
    .scales =
      switch(method,
             sd    = lapply(l, apply, 1, sd, na.rm = TRUE),
             mad   = lapply(l, apply, 1, mad, na.rm = TRUE)
      )
    if (is.null(.scales))
      stop(sprintf("Method '%s' is unknown in per group mode", method))
    .scales = do.call("cbind", .scales)  # list to matrix
    .scales = apply(.scales, 1, median)  # extract median values
    .scales = .scales / 2
  }
  # Ready to apply scaling
  # Avoid scaling == 0, switch to 5% of all scales
  .min = quantile(.scales, probs = 0.05, na.rm = TRUE)
  .scales[.scales == 0] = .min
  # Do scaling
  tdms$exprs = sweep(tdms$exprs, 1, .scales, FUN = "/")
  tdms
}

# Sort within group
# TODO: to be continued
tdms.sort = function(tdms, group, filter = "none", min.size = 3) {
  # named factor
  if (is.character(group) & group[1] %in% colnames(tdms$pData))
    group = tdms$pData[group]
  # Unrecognized grouping
  if (!is.factor(group) & length(group) != nrow(tdms$pData))
    stop("Grouping is unrecognized")
  # Split matrix by group
  l = tdms.split.matrix(m = tdms$exprs, grouping = group, filter = filter, min.size = min.size)
  # Impute: replace NA with median, optional
  .scales = lapply(l, apply(1, sd, na.rm = TRUE))
  # Sort
  if (is.null(.scales))
    stop(sprintf("Method '%s' is unknown in per group mode", method))
  .scales = do.call("cbind", .scales)  # list to matrix
  .scales = apply(.scales, 1, median)  # extract median values
  tdms
}


# Sort by factor
tdms.sort.byfactor = function(tdms, group) {
  # named factor
  if (tolower(group) == "none") return(tdms)
  if (is.character(group) & group[1] %in% colnames(tdms$pData))
    group = factor(tdms$pData[,group])
  # Unrecognized grouping
  if (!is.factor(group) & length(group) != nrow(tdms$pData))
    stop("Grouping is unrecognized")
  # Re-factor, keeping original order with help of unique()
  ff = factor(group, unique(group))
  oo = order(ff)  # sort
  # Apply the new order
  tdms$exprs = tdms$exprs[, oo]
  tdms$pData = tdms$pData[oo, ]
  tdms
}


# Final centering of reference group
tdms.center.ref = function(tdms, method = "median", grouping = NA, group = NA) {
  if (!method %in% strsplit("median:mean", ":")[[1]]) {
    warning(sprintf("Method '%s' is unknown, switching to median", method))
    method = "median"
  }

  # TODO: group defined, method=none, error
  # Check grouping
  if (tolower(grouping) == "none" || is.na(grouping[1])) {
    grouping = NA
  } else {
    # named factor
    if (is.character(grouping) && grouping[1] %in% colnames(tdms$pData))
      grouping = factor(tdms$pData[,grouping])
    # Unrecognized grouping
    if (!is.factor(grouping) && length(grouping) != nrow(tdms$pData))
      stop("Grouping is unrecognized")
  }
  print(sprintf("Meth %s Group %s\n", method, group)); print(grouping)

  # Compute scales, grouping
  # Split matrix by group
  l = tdms.split.matrix(m = tdms$exprs, grouping = grouping, filter = "only", group = group, min.size = min.size)
  if (!group %in% names(l)) {
    warning(sprintf("%s is not found", group))
    return(NULL)
  }
    .centers =
      switch(method,
             median = lapply(l, apply, 1, median, na.rm = TRUE),
             mean   = lapply(l, apply, 1, mean, na.rm = TRUE)
      )
    print(.centers)
    if (is.null(.centers))
      stop(sprintf("Method '%s' is unknown in per group mode", method))
    .centers = do.call("cbind", .centers)  # list to matrix
    .centers = apply(.centers, 1, median)  # extract median values
  # Do zeroing
  tdms$exprs = sweep(tdms$exprs, 1, .centers, FUN = "-")
  tdms
}

# Summarize
tdms.summarize = function(tdms, group = "none", probs = c(0.5), filter = "none", min.size = 1, inter.leave = F) {
  # named factor
  if (tolower(group) == "none") return(tdms)
  if (is.character(group) && group[1] %in% colnames(tdms$pData))
    group = factor(tdms$pData[,group], unique(tdms$pData[,group]))
  # Unrecognized grouping
  if (!is.factor(group) || length(group) != nrow(tdms$pData))
    stop("Grouping is unrecognized")
  # TODO: check min.size vs length(probs)
  # Re-factor, keeping original order with help of unique()
  ff = factor(group, unique(group))
  oo = order(ff)
  # Apply the new order
  tdms$exprs = tdms$exprs[, oo]
  tdms$pData = tdms$pData[oo,]
  # Split matrix by group
  l = tdms.split.matrix(m = tdms$exprs, grouping = group, filter = filter, min.size = min.size)
  # Summarize exprs
  .exprs = lapply(l, apply, 1, quantile, probs = probs, na.rm = TRUE)
  if (length(probs) > 1) .exprs = lapply(.exprs, t)  # list to matrix
  .exprs = do.call("cbind", .exprs)  # list to matrix
  dd = data.frame(id=rep(names(l), each=length(probs)))
  # TODO: add group  # dd$group = dd$id
  pp = sapply(rep(probs*100, times=length(l)), function(x) sprintf("%02d%%", x))
  #rownames(dd) = paste(dd$id, colnames(.exprs), sep="_")
  colnames(.exprs) = rownames(dd) = dd[,1] = paste(dd$id, pp, sep="_")
  # Modify the object
  tdms$exprs = .exprs
  tdms$pData = dd
  # Interleave
  if (inter.leave) {
    .cgrp = length(unique(group))
    .cprb = length(probs)
    oo = rep(seq(.cprb), each = .cgrp) + (seq(.cgrp)-1) * .cprb
    tdms$exprs = .exprs[ , oo]
    tdms$pData = dd[oo,]
  }
  tdms
}

# Transpose
tdms.t = function(tdms) {
  tdms$exprs = t(tdms$exprs)
  swap = tdms$pData
  tdms$pData = tdms$fData
  tdms$fData = swap
  tdms
}

# Remove uninformative rows and columns
tdms.remove = function(tdms, margins, emptyness = 0.20, unicity = 1) {
  for (margin in margins) {
    # Check margin
    margin = trunc(margin)
    if (margin<1 || margin>2) warning(sprintf("Incorrect margin: %d", margin))
    # Init indices
    idx = NULL
    # Look for emptyness
    if (!is.na(emptyness)) {
      res = apply(tdms$exprs, margin, function(x) sum(is.na(x)))
      if (emptyness < 1) {  # as a ratio
        if (emptyness > 0.5) message(sprintf("Emptyness criteria is not stringent: %f", emptyness))
        idx = which(res / dim(tdms$exprs)[margin] >= emptyness)
      } else if (emptyness >= 1) {  # as a count
        emptyness = trunc(emptyness)
        if (emptyness < 5) message(sprintf("Emptyness criteria is very stringent: %d", emptyness))
        idx = which(res >= emptyness)
      } else {
        warning(sprintf("Negative emptyness is ignored: %f", emptyness))
      }
    }
    # Look for unicity
    if (!is.na(unicity)) {
      res = apply(tdms$exprs, margin, function(x) length(unique(na.omit(x))))
      if (unicity < 1) {
        if (unicity < 0.5) message(sprintf("Unicity criteria is not very stringent: %f", unicity))
        idx = c(idx, which(res / dim(tdms$exprs)[margin] <= unicity))
      } else if (unicity >= 1) {
        unicity = trunc(unicity)
        if (unicity < 5) message(sprintf("Unicity criteria is not very stringent: %d", unicity))
        idx = c(idx, which(res <= unicity))
      } else {
        warning(sprintf("Negative unicity is ignored: %f", unicity))
      }
    }
    # Remove indices
    if (length(idx)) {
      idx = unique(idx)  # idx may have matched both criteria
      if (margin == 1) {
        tdms$fData = tdms$fData[-idx, , drop = FALSE]
        tdms$exprs = tdms$exprs[-idx, , drop = FALSE]
      } else {
        tdms$pData = tdms$pData[-idx, , drop = FALSE]
        tdms$exprs = tdms$exprs[ ,-idx, drop = FALSE]
      }
    }
  }
  tdms
}

# TODO: add clustering function
