## Model of interconnected regime shifts
## If follows partially the equations of Tu et al 2019
## Juan Rocha
## 190412


library(tidyverse)
library(deSolve)
library(rEDM)

## model

resource <- function(t, R, param){
    with(as.list(R), {
        growth <- alpha * R * (R - c_allee)*(1-(R/K))

        ## For exmports and imports, I assume delta is 1 and K_L (demand by population) is also 1
        exports <- delta * (t(C_ij) * t(A_ij)) %*% R
        imports <- delta *(C_ij * A_ij) %*% R

        dR <- growth  - exports + imports + rnorm(n=n, mean = 0, sd = 1)
        return(list(dR, exports, imports))
    })
}




n <- 5 ## number of countries and resources
R <- runif(n, 10, 100) ## 10 resources with random numbers
alpha <- runif(n, 0.01, 0.1) ## growth rate
c_allee <- runif(n, 10, 30) ## Allee parameter: food necessary to feed L_i population in country i
K <- rep(100, times = n) ## carrying capacity

C_ij <- matrix(runif(n^2, min = 0, max = 0.5), ncol = n)
diag(C_ij) <- 0

A_ij <- matrix(rbinom(n^2, 1, prob = 0.25), ncol = n)
diag(A_ij) <- 0

delta <-  1 # resource depletion coefficient (I dont' get why is necessary)

times <- seq(from = 0, to = 100, by = 0.01)
params <- list(alpha = alpha, c_allee = c_allee, K = K, C_ij = C_ij)

print(system.time(
    out <- ode(y = R, times = times, parms = params, func = resource)
))


df <- out %>% as_data_frame() %>%
    gather(key = "species", value = "population", 2:(n+1))


## plot result
df %>%
    ggplot(aes(x=time, y=population)) +
    geom_line(aes(color = species), show.legend = T) +
    theme_light()

## Test of CCM as identification method:
## follow my own code from 04-DataExploration.R

## Embedding dimension:
# emb <- list()
# for (i in 2:dim(out)[2]){
#     emb[[i-1]] <- out[c(1,i)] %>%
#         simplex()
# }

emb <- simplex(out[,c(1,2)])
plot(emb$E, emb$rho, type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)")




ind <- crossing(
    lib_column = colnames(out)[-1],
    target_column = colnames(out)[-1])

ind <- ind %>% filter(lib_column != target_column)
