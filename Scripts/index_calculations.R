

### function that calls the index
calc_index <- function(behaviours = NULL, 
                       valence = NULL, 
                       interaction.list = NULL, 
                       dyads = NULL, 
                       type = c('CSI')) {
  if(type == 'CSI'){xx.res = calc_CSI(behaviours = behaviours, dyads)}
  return(xx.res)
}




### calculate CSI
calc_CSI <- function(behaviours, dyads) {
  xx <- do.call(cbind, lapply(1:ncol(behaviours), function(x) {
    aggregate(behaviours[,x], by = list(dyads), sum)[,2] / 
      mean(behaviours[,x], na.rm = T)
  }))
  xx <- rowMeans(xx, na.rm = T)
  names(xx) = sort(unique(dyads))
  xx = xx[match(dyads, names(xx))]
  return(xx)
}

