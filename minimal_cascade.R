## Minimalistic cascading regime shifts model
## Developed in collaboration with Anne-Sophie Crepin
## by Juan Rocha
## juan.rocha@su.se

library(tidyverse)
library(deSolve)
library(deTestSet)


## model
resource <- function(t, H, param){
    dH <- H * (1-(H/ a_i * A)) - (1/C_i)*((H^alpha) / (1 + H^alpha))
    return(list(dH, a_i))
}

## initial parameters
H0 = 0.1
a_i = 0.5# efficiency of fish in eating Algae?
A = 1 # Initial avalialble algae? This is the equivalent of carrying capacity K in LV.
C_i = 1 # Coral abundance
alpha = 2 # efficiency term of herbivores which depend on coral cover?

# model parameters
param <- c(a_i, A, C_i, alpha)
yini <- c(H0)
times <- seq(0,200, by = 1)

# run
print(system.time(
    out <- ode(func = resource, y = yini, parms = param, times = times)
))

out %>% as_tibble() %>%
    rename(H = `1`) %>%
    ggplot(aes(time, H)) +
    geom_line() +
    theme_minimal()

## testing sensitivity to one parameter
alpha_test <- c(0.9, 1, 1.5, 2, 2.5, 3) 
out <- list()

for (i in seq_along(alpha)){
    # set value of parameter of interest
    alpha  <- alpha_test[i]
    
    # set model
    param <- c(a_i, A, C_i, alpha)
    yini <- c(H0)
    times <- seq(0,200, by = 1)
    
    # run
    out[[i]] <-  ode(func = resource, y = yini, parms = param, times = times)
}

out %>% map(., as.data.frame) %>%
    bind_rows() %>%
    rename(H = `1`, alpha = `2`) %>%
    ggplot(aes(time, H, group = alpha)) + 
    geom_line(aes(color = as.factor(alpha)))


#### Now with a_i
ai_test <- c(0.1, 0.5, 0.9, 1, 1.5, 2, 2.5, 3) 
out <- list()

for (i in seq_along(ai_test)){
    # set value of parameter of interest
    a_i  <- ai_test[i]
    
    # set model
    param <- c(a_i, A, C_i, alpha)
    yini <- c(H0)
    times <- seq(0,200, by = 1)
    
    # run
    out[[i]] <-  ode(func = resource, y = yini, parms = param, times = times)
}

out %>% map(., as.data.frame) %>%
    bind_rows() %>%
    rename(H = `1`, a_i = `2`) %>%
    ggplot(aes(time, H, group = a_i)) + 
    geom_line(aes(color = as.factor(a_i))) +
    theme_minimal()



### Minimal cascade:
## model
resource <- function(t, H, param){
    with(as.list(H), {
        growth <-  H * (1-(H/ a_i * A)) - (1/C_i)*((H^alpha) / (1 + H^alpha))
        
        ## For exmports and imports, I assume delta is 1 and K_L (demand by population) is also 1
        exports <-  t(A_ij) %*% H
        imports <-  A_ij %*% H
        
        dH <- growth  - exports + imports #+ rnorm(n=n, mean = 0, sd = 1)
        dH <- ifelse(dH < 0, 0, dH)
        return(list(dH)) #,  exports, imports
    })
}

## J190415: Don't add exports and imports yet because in the model they are not 
# modeled pair-wise, they are agregated. One would need to dissagregate or store
# the matrix instead of the vector. Where the vector now is the sum of the exports 
# per species/country, while the matrix is the raw values.

## initial parameters
n <- 5 ## number of patches
H0 <- runif(n, 0.1, 1) ##  resources with random numbers

a_i = 2# efficiency of fish in eating Algae?
A = 1 # Initial avalialble algae? This is the equivalent of carrying capacity K in LV.
C_i = 1 # Coral abundance
alpha = 2 # efficiency term of herbivores which depend on coral cover?

A_ij <- matrix(rbinom(n^2, 1, prob = 0.5), ncol = n) * 0.8
diag(A_ij) <- 0


times <- seq(from = 0, to = 100, by = 0.1)
param <- c(a_i, A, C_i, alpha)

## run the model
print(system.time(
    out <- ode(
        y = H0, times = times, parms = params, func = resource## see help("ode") for more methods
    )
))


df <- out %>% as_data_frame() %>%
    gather(key = "patches", value = "population", 2:(n+1))

df
df %>%
    ggplot(aes(x=time, y=population)) +
    geom_line(aes(color = patches), size = 0.25, show.legend = T) +
    theme_light()

image(A_ij)
