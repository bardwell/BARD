# function to get a vector of states for a drawn value frin filtering posterior # 
# segs is vector of locations of segments
# states vector 
get.states <- function( segs , states ){
  
  state.vec <- numeric(n)
  state.vec[ segs[1]:segs[2] ] <- states[1]
  for (i in 2:( length(segs) - 1 ) ){
    state.vec[ ( segs[i]-1 ):segs[i+1] ] <- states[i]
  }

  return(state.vec)
  
}

