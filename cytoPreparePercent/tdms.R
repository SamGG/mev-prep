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
  hdr = read.table(file, sep = sep, nrows = nrows, header = FALSE, as.is = T)
  # Automatic detection: select cells from lower right corner up to an
  # increasing upper left corner and try to convert those columns to numerical
  # values; memorize the highest upper left corner.
  if (is.na(rd1) | is.na(cd1)) {
    for (i in seq(nrow(hdr)-10, 2)) {
      for (j in seq(ncol(hdr)-6, 2)) {
        r = sapply(hdr[i:nrow(hdr),j:ncol(hdr)], 
                   function(x) is.numeric(type.convert(x)))
        if (all(r==TRUE)) {
          rd1 = i ; cd1 = j
        }
        # cat(i, j, r, "\n")
      }
    }
    if (is.na(rd1))
      stop("Automatic detection of data start failed!")
  }
  if (rd1 < 2 | cd1 < 2) {
    stop( sprintf("Incorrect start of data: row/col = %d,%d", rd1, cd1))
  }
  # Extract annotations on columns
  ann.col = as.data.frame(t(hdr[seq(1,rd1-1), seq(cd1, ncol(hdr))]))
  colnames(ann.col) = c("Id", hdr[seq(2,rd1-1), cd1-1])
  # Read full file and extract row annotations and data
  dat = read.table(file, sep = sep, skip = rd1-1, header = FALSE)
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
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}


# transforming
tdms.transform = function(tdms, func = "none", param = NA) {
  if (func == "none") {
    # TODO: allow negative values?
  } else if (func == "log2") {
    if (is.na(param)) {
      param = 0
    } else if (param < 0) {
      warning("param can't be negative")
      param = 0
    }
    if (any(tdms$exprs <= 0)) {
      warning("negative values will lead to NA")
      tdms$exprs[tdms$exprs <= 0] = NA
    }
    tdms$exprs = log2( tdms$exprs + param )
  } else if (func == "asinh") {
    if (is.na(param)) {
      param = 1
    } else if (param < 0) {
      warning("param can't be negative")
      param = 1
    }
    if (any(tdms$exprs <= 0)) {
      warning("negative values will lead to NA")
      tdms$exprs[tdms$exprs <= 0] = NA
    }
    tdms$exprs = asinh( tdms$exprs / param )
  } else {
    warning(sprintf("Unknown transformation %s, data are unchanged.", func))
  }
  tdms
}

# centering
tdms.center = function(tdms, method="median") {
  # TODO: set and check trimming parameter: if <=1 then percent else count
  # TODO: percent as paramater
  if (!method %in% strsplit("mean:median:trimmed:minmax", ":")[[1]]) {
    warning(sprintf("Method '%s' is unknown, switching to median", method))
    method = "median"
  }
  # Compute centers
  .centers = 
    switch(method,
           mean    = apply(tdms$exprs, 1, mean, na.rm = TRUE),
           median  = apply(tdms$exprs, 1, median, na.rm = TRUE),
           trimmed = apply(tdms$exprs, 1, mean, trim = 0.05, na.rm = TRUE),
           minmax = apply(tdms$exprs, 1, function(x) 
             sum(quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))/2)
           )
  # Do centering
  tdms$exprs = sweep(tdms$exprs, 1, .centers, FUN = "-")
  tdms
}

tdms.split.matrix = function(m, f, filter = "none", min.size = 1) {
  lf <- levels(f)
  y <- vector("list", length(lf))
  names(y) <- lf
  ind <- split(seq(ncol(m)), f)
  for (k in lf) y[[k]] <- m[ ,ind[[k]], drop = FALSE]
  # Remove reference
  if (filter == "rm.ref") {
    ind.ref = match(tolower("ref"), tolower(lf))
    if (!is.na(ind.ref)) y = y[-ind.ref]
  # Keep reference alone
  } else if (filter == "keep.ref") {
    ind.ref = match(tolower("ref"), tolower(lf))
    if (!is.na(ind.ref)) y = y[ind.ref]
  # None
  } else if (filter != "none") {
    warning(sprintf("Filter '%s' is unknown.", filter))
  }
  # Return, debug lapply(y, head)
  y
}

# scaling
# if percentile, group is irrelevant
# if sd or mad, group is used as follows
# if group == NA, global scale is computed
# if group != NA, per group is computed
tdms.scale = function(tdms, method = "sd", percent = 0.05, group = NA, filter = "none", min.size = 3) {
  if (!method %in% strsplit("sd:mad:minmax", ":")[[1]]) {
    warning(sprintf("Method '%s' is unknown, switching to sd", method))
    method = "sd"
  }
  # Compute scales, no grouping
  if (is.na(group)) {
    # TODO: set and check trimming parameter: if <=1 then percent else count
    # TODO: trimming
    # TODO: winsoring
    if (method == "minmax") {
      if (percent <= 0) percent = 0
      if (percent >= 0.25) percent = 0.25
      probs = c(percent, 1-percent)
    }
    .scales = 
      switch(method,
             sd    = apply(tdms$exprs, 1, sd, na.rm = TRUE),
             mad   = apply(tdms$exprs, 1, mad, na.rm = TRUE),
             minmax = apply(tdms$exprs, 1, function(x) 
               -diff(quantile(x, probs = probs, na.rm = TRUE))/2)
      )
  # Compute scales, grouping
  } else if (is.factor(group)) {
    # Split matrix by group
    l = tdms.split.matrix(m = tdms$exprs, f = group, filter = filter, min.size = min.size)
    .scales = 
      switch(method,
             sd    = lapply(l, apply(1, sd, na.rm = TRUE)),
             mad   = lapply(l, apply(1, mad, na.rm = TRUE))
      )
    if (is.null(.scales))
      stop(sprintf("Method '%s' is unknown in per group mode", method))
    .scales = do.call("cbind", .scales)  # list to matrix
    .scales = apply(.scales, 1, median)  # extract median values
  # Unrecognized grouping
  } else {
    stop("Grouping is unrecognized")
  }
  # Avoid scaling == 0, switch to 5% of all scales
  .min = quantile(.scales, probs = 0.05, na.rm = TRUE)
  .scales[.scales == 0] = .min
  # Do scaling
  tdms$exprs = sweep(tdms$exprs, 1, .scales, FUN = "/")
  tdms
}

# Sort within group
tdms.sort = function(tdms, group, filter = "none", min.size = 3) {
  # Unrecognized grouping
  if (!is.factor(group))
    stop("Grouping is unrecognized")
  # Split matrix by group
  l = tdms.split.matrix(m = tdms$exprs, f = group, filter = filter, min.size = min.size)
  # Impute: replace NA with median, optional
  .scales = lapply(l, apply(1, sd, na.rm = TRUE))
  # Sort
  if (is.null(.scales))
    stop(sprintf("Method '%s' is unknown in per group mode", method))
  .scales = do.call("cbind", .scales)  # list to matrix
  .scales = apply(.scales, 1, median)  # extract median values
  tdms
}

# Final centering of reference group
tdms.center.ref = function(tdms) {
  tdms
}

# Summarize
tdms.summarize = function(tdms, group, probs = c(0.5), filter = "none", min.size = 1) {
  # Unrecognized grouping
  if (!is.factor(group))
    stop("Grouping is unrecognized")
  # Split matrix by group
  l = tdms.split.matrix(m = tdms$exprs, f = group, filter = filter, min.size = min.size)
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
