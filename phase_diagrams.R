## Phase diagram explorations
## Juan Rocha
## 210121

library(deSolve)
library(rootSolve)
library(tictoc)
library(tidyverse)

#### Resource model ####
# This is where we define your event function
# Add this directly above your call to ode()
posfun <- function(t, y, parms){
    with(as.list(y), {
        y[which(y<0)] <- 0
        return(y)
    })
}


## Model
resource_norm <- function(t, y, params){
    with(as.list(c(y, params)), {
        ## Normalization steps:
        Y = y/r0
        R = r/r0
        K0 = K/q
        C = c/(r0*q)
        Qij_parts <- d_ij * q
        Qij = t(Qij_parts) / Qij_parts ## the same as lower triangle / upper triangle and viceversa
        ## the next line corrects the case when Zij is NaN induced by division by zero. In that case, Zij should be a matrix full of zeroes = no connections.
        if(all(is.nan(Qij))) Qij = Qij_parts
        
        # set diag to zero
        diag(Qij) <- 0
        D_ij = d_ij / r0
        
        ## System with normalized parameters
        resource <- R * Y * (1 - (Y/K0)) - C * (Y^b / (1 + Y^b))
        outflow <-  (A_ij * D_ij * Qij) %*% Y
        inflow <-  t(A_ij * D_ij) %*% Y 
        # Qij is not multiplied by D_ij by suggestion of ASC.
        
        dX <- resource + (inflow - outflow)
        
        return(list(c(dX)))
    })
}

## set up flux matrix
n <- 5        # number of systems
#delta_ij <- matrix(runif(n^2, min = 0.02, max = 0.05), ncol = n)
d_ij <- matrix(rep(0,n^2), ncol = n) # turn off difussion terms
diag(d_ij) <- 0
A_ij <- matrix(rbinom(n^2, 1, prob = 0.3), n, n)
diag(A_ij) <- 0


params <- list(
    r0 = rep(1,n),                # scaling parameter
    r = rep(1.5, n),              # intrinsic growth rate
    K = rep(10, n),               # Carrying capacity
    c = rep(2, n),                # max predation rate
    q = 1, #runif(n, min = 2, max = 8),        # threshold
    b = 4 ,                       # sharpness of the shift
    d_ij = d_ij,                  # matrix of diffusion terms
    A_ij = A_ij                   # adjacency matrix
)

## set up time steps
times <- seq(from = 0, to = 100, by = 0.01)

## initial conditions
yini <- runif(n, 0.5, 5)

## run the model
out <- ode(
    y = yini, times = times,  func = resource_norm, parms = params,
    method = "bdf", ## see help("ode") for more methods
    events=list(func = posfun, time = times)
)

## steady state:
tic()
stst <- steady(
    y = yini, parms = params, func = resource_norm, method = "runsteady",
    time = c(1,10))
toc()


stst




df <- out %>% as_tibble() %>%
    gather(key = "species", value = "population", 2:last_col())

df %>%
    ggplot(aes(x=time, y=population)) +
    geom_line(aes(color = species), size = 0.25, show.legend = T) +
    scale_y_continuous(limits = c(0,NA))
