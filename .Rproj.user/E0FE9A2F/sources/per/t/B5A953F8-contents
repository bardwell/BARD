#' Marginal probability of being abnormal
#'
#' @param res A bard result
#' @param no.draws Number of draws to take from posterior distribution
#' @param plot If true produce plot
#'
#' @return Marginal probability of being abnormal at each point
#'
#' @examples
#' set.seed(1)
#' k_N <- 10
#' p_N <- 0.1
#' k_A <- 15
#' p_A <- 0.3
#' pi_N <- 0.8
#' affected_dim <- 8
#' lower <- 0.3
#' upper <- 0.7
#' simparams <- c(k_N, p_N, k_A, p_A, pi_N, affected_dim, lower, upper)
#' sim <- sim.data(T = 1000, N = 200, simparams)
#'
#' bardparams <- c(k_N, p_N, k_A, p_A, pi_N, affected_dim)
#' mu_seq <- c(seq(-0.7, -0.3, by=0.05), seq(0.3, 0.7, by=0.05))
#' res <- Rbard(sim$data, bardparams, mu_seq, alpha = 1e-4)
#'
#' margprob = marginal.prob(res, no.draws = 100, plot = T)
#'
#' @export
#'
marginal.prob <- function(res, no.draws=1000, plot = F){

  logprobs <- res$weights
  cpt.locs <- res$locations
  types <- res$type
  params <- res$params

  n = length(logprobs)
  prob.states <- matrix(nrow=no.draws, ncol=n)

  for (i in 1:no.draws){

    f <- draw.from.post(logprobs, cpt.locs, types, params)
    segs <- f$draws
    states <- f$states

    # probs of being N or A
    prob.states[i,] <- get.states(segs, states, n)

  }

  margprob = apply(prob.states,2,sum)/no.draws
  if (plot == T){
    plot(margprob, t="l",
         xlab="t", ylab = expression("Pr("~B[t]~"=A)" ))
  }

  return(margprob)

}


#' Multiple simulations for a heatmap plot
#'
#' @param res A bard result
#' @param no.draws Number of draws to take from posterior distribution
#'
#' @return A list of matrices
#'
#' @examples
#' set.seed(1)
#' k_N <- 10
#' p_N <- 0.1
#' k_A <- 15
#' p_A <- 0.3
#' pi_N <- 0.8
#' affected_dim <- 8
#' lower <- 0.3
#' upper <- 0.7
#' simparams <- c(k_N, p_N, k_A, p_A, pi_N, affected_dim, lower, upper)
#' sim <- sim.data(T = 1000, N = 200, simparams)
#'
#' bardparams <- c(k_N, p_N, k_A, p_A, pi_N, affected_dim)
#' mu_seq <- c(seq(-0.7, -0.3, by=0.05), seq(0.3, 0.7, by=0.05))
#' res <- Rbard(sim$data, bardparams, mu_seq, alpha = 1e-4)
#'
#' hs = heatmap.simulations(res, no.draws = 100)
#'
#' @export
#'
heatmap.simulations <- function(res, no.draws=1000)
{

  logprobs <- res$weights
  cpt.locs <- res$locations
  types <- res$type
  params <- res$params
  # results list
  sampled.res = list()

  for (i in 1:no.draws){

    f <- draw.from.post(logprobs, cpt.locs, types, params)
    segs <- f$draws
    states <- f$states

    end = segs[which(states == 1)]
    begin = segs[which(states == 1) + 1]
    tmpmat = matrix(nrow = 2, ncol = length(end))
    tmpmat[1,] = begin
    tmpmat[2,] = end
    sampled.res[[i]] = tmpmat

  }

  return(sampled.res)

}


# draw one sample of chpts and states from the posterior
draw.from.post <- function(logprobs, cpt.locs, types, params){

  n <- length(logprobs)
  l.probs <- logprobs[[n]]
  seg.locs <- cpt.locs[[n]]
  seg.types <- types[[n]]

  ind <- 1:length(l.probs)
  a <- sample( ind , size=1 , replace = TRUE, exp( l.probs ) )

  # state amd location
  curr.state <- seg.types[a]
  t <- seg.locs[a]

  STATES <- curr.state
  DRAW <- t

  k_N = params[1]
  p_N = params[2]
  k_A = params[3]
  p_A = params[4]
  pi_N = params[5]
  pi_A = 1 - pi_N

  while (t > 1){

    if(curr.state==0){

      l.probs <- logprobs[[t]][ types[[t]] == 1 ]
      i <- cpt.locs[[t]][ types[[t]] ==1 ]
      back.dens <- log(pi_N) + dnbinom(t-i-1,k_A,p_A,log=TRUE) - pnbinom(t-i-2,k_A,p_A, lower.tail=FALSE , log.p=TRUE)
      l.probs <- l.probs + back.dens

      ind <- 1:length(l.probs)
      a <- sample( ind , size=1 , replace = TRUE, exp( l.probs ) )

      # state amd location
      curr.state <- 1
      t <- i[a]
      STATES <- c(STATES,curr.state)
      DRAW <- c(DRAW,t)
    }
    else{

      # currently in abnormal state poss transitions to normal or abnormal
      i.N <- cpt.locs[[t]][ types[[t]] == 0 ]
      i.A <- cpt.locs[[t]][ types[[t]] == 1 ]

      l.probsA <- logprobs[[t]][  types[[t]] == 1 ]
      back.densA <- log(pi_A) + dnbinom(t-i.A-1,k_A,p_A,log=TRUE) - pnbinom(t-i.A-2,k_A,p_A, lower.tail=FALSE , log.p=TRUE)
      l.probsA <- l.probsA + back.densA
      # normal points
      l.probsN <- logprobs[[t]][  types[[t]] == 0 ]
      back.densN <- dnbinom( t-i.N-1 , k_N , p_N , log=TRUE ) - pnbinom( t-i.N-2 , k_N , p_N , lower.tail=FALSE , log.p=TRUE )
      l.probsN <- l.probsN + back.densN
      # now normalise
      i.all <- c( i.N , i.A )
      t.all <- c( rep(0,length(i.N)) , rep(1,length(i.A)) )
      l.probs <- c( l.probsN , l.probsA )
      c <- max(l.probs)
      l.probs.norm <- l.probs - ( c + log( sum( exp( l.probs - c ) ) ) )

      ind <- 1:length(l.probs.norm)
      a <- sample( ind , size=1 , replace = TRUE, exp( l.probs.norm ) )

      # state amd location
      curr.state <- t.all[a]
      t <- i.all[a]
      STATES <- c(STATES,curr.state)
      DRAW <- c(DRAW,t)
    }

  }

  DRAW[DRAW==0] <- 1
  newlist <- list("draws"=c(n,DRAW),"states"=STATES)
  return(newlist)

}

get.states <- function(segs, states, n){

  state.vec <- numeric(n)
  state.vec[ segs[1]:segs[2] ] <- states[1]
  for (i in 2:( length(segs) - 1 ) ){
    state.vec[ ( segs[i]-1 ):segs[i+1] ] <- states[i]
  }

  return(state.vec)

}
