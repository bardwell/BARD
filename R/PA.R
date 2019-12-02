# Calculating log P.A(s,t) for segment (s,t) integrated on uniform prior over mu_seq
# done without normal part at front
P.A <- function(s, t, mu_seq, N, S, p){

  # prior value mu_dens
  mu_dens <- 1/( tail(mu_seq,1) - mu_seq[1] )

  # width of rectangle
  mu_wid <- diff(mu_seq)[1]

  vec <- numeric(length(mu_seq))

  # evaluate at each point of grid
  # do this as typically smaller than dimension
  # evaluating log of quantity
  for (k in 1:length(mu_seq)){

    vec[k] <- N*log(1-p) + sum( log( 1 + exp( mu_seq[k] * ( S[(t+1),] - S[s,]  -  mu_seq[k] * (t-s+1)/2 ) + log(p) - log(1-p) ) ) )

  }

  # finding sum of logs -- for numerical instability
  cmax <- max( vec )

  marg.like <- cmax + log( sum( exp( vec - cmax ) ) ) + log(mu_dens) + log(mu_wid)

  return(marg.like)

}
