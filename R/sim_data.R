#' Simulate data from model
#'
#' @param T Length of series
#' @param N Number of variates
#' @param simparams Vector of length 8 containing  k_N, p_N, k_A, p_A, pi_N, affected_dim, lower and upper
#'
#' @return List containing data (T by N matrix), states and changepoints
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
#' @export
#'
sim.data = function(T = 1000, N = 200, simparams = NA){

  if (length(simparams) != 8){
    stop("Not enough params should be a vector of length 8.")
  }

  k_N = simparams[1]
  p_N = simparams[2]
  k_A = simparams[3]
  p_A = simparams[4]
  pi_N = simparams[5]
  affected_dim = simparams[6]
  lower = simparams[7]
  upper = simparams[8]

  pi_A = 1 - pi_N

  # start in a normal state encoded as 0
  cpts <- rnbinom( 1 , k_N , p_N )
  state <- rep( "N", cpts )

  data <- matrix( rnorm(cpts*N) , nrow = cpts ,ncol = N)
  curr.state <- 0
  n <- cpts

  while (n < T){

    # if current state is normal - draw abnormal length
    if (curr.state == 0){

      # in normal transition to abnormal length ~ NBinom( k_A , p_A )
      seg.length <- rnbinom( 1 , k_A , p_A )
      cpts <- c( cpts , tail(cpts,1)+seg.length )
      state <- c( state , rep( "A" , seg.length ) )

      # grow data seg.length no of rows
      data <- rbind( data, matrix(rnorm(seg.length*N),ncol=N) )

      # draw value of mu for dims with affected means
      mu <- runif(1 ,lower ,upper)
      data[ (tail(cpts,2)[1]+1):tail(cpts,1) ,  1:affected_dim  ] <- data[ (tail(cpts,2)[1]+1):tail(cpts,1) , 1:affected_dim ] + mu

      # current state abnormal update it
      curr.state <- 1

    }

    else{

      # current state is abnormal curr.state == 1
      # with prob pi_N go to Normal
      # with prob pi_A go to Abnormal with diff mu
      # generate a u bernoulli with probs pi_A, pi_N
      u <- rbinom(1, 1, pi_A)

      # u==0 transition to normal

      if (u==0){

        seg.length <- rnbinom(1, k_N, p_N )
        cpts <- c( cpts , tail(cpts,1)+seg.length )
        state <- c( state , rep( "N" , seg.length ) )

        # grow data seg.length no of rows
        data <- rbind( data, matrix(rnorm(seg.length*N),ncol=N) )

        curr.state <- 0

      }
      else{

        # u==1 thus transition to abnormal
        seg.length <- rnbinom( 1 , k_A , p_A )
        cpts <- c( cpts , tail(cpts,1)+seg.length )
        state <- c( state , rep( "A" , seg.length ) )

        # grow data seg.length no of rows
        data <- rbind( data, matrix(rnorm(seg.length*N),ncol=N) )

        # draw mu from prior
        mu <- runif(1, lower, upper)

        # random dimensions to affect ran.affec
        # ran.affec is first one
        ran.affec <- sample(affected_dim:(N-affected_dim),1)
        data[ (tail(cpts,2)[1]+1):tail(cpts,1), ran.affec:(ran.affec+affected_dim)  ] <- data[ (tail(cpts,2)[1]+1):tail(cpts,1) , ran.affec:(ran.affec+affected_dim) ] + mu

        curr.state <- 1

      }

    }

    n <- tail(cpts,1)

  }

  retlist = list("data"=data, "state"=state, "cpts"= cpts)
  return(retlist)

}
