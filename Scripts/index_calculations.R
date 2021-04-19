

### function that calls the index
calc_index <- function(behaviours = NULL, 
                       valence = NULL, 
                       interaction.list = NULL, 
                       dyads = NULL, 
                       type = c('DSI', 'CRI')) {
  if(type == 'DSI'){xx.res = calc_DSI(behaviours = behaviours, dyads)}
  if(type == 'CRI'){xx.res = calc_CRI(behaviours = behaviours, dyads, valence = valence)}
  return(xx.res)
}




### calculate DSI
calc_DSI <- function(behaviours, dyads) {
  xx <- do.call(cbind, lapply(1:ncol(behaviours), function(x) {
    aggregate(behaviours[,x], by = list(dyads), sum)[,2] / 
      mean(behaviours[,x], na.rm = T)
  }))
  xx <- rowMeans(xx, na.rm = T)
  names(xx) = sort(unique(dyads))
  xx = xx[match(dyads, names(xx))]
  return(xx)
}

### calculate CRI
calc_CRI <- function(behaviours, valence, dyads) {
  xx <- lapply(1:ncol(behaviours), function(x) {
    aggregate(behaviours[,x], by = list(dyads), sum)[,2] / 
      mean(behaviours[, x], na.rm = T)
  })
  vals <- unique(valence)
  xx.vals <- do.call(cbind, lapply(vals, function(x) {
    y <- which(valence == x)
    xx.val <- do.call(cbind, xx[y])
    xx.val <- rowMeans(xx.val)
    return(xx.val)
  }))
  xx.pos <- rowMeans(cbind(xx.vals[, vals != "neg"]))
  xx.neg <- rowMeans(cbind(xx.vals[, vals == "neg"]))
  xx <- xx.pos - xx.neg
  names(xx) = sort(unique(dyads))
  xx = xx[match(dyads, names(xx))]
  return(xx)
}
