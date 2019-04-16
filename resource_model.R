## Model of interconnected regime shifts
## If follows partially the equations of Tu et al 2019
## Juan Rocha
## 190412


library(tidyverse)
library(deSolve)
library(deTestSet)
library(rEDM)

## model

resource <- function(t, R, param){
    with(as.list(R), {
        growth <- alpha * R * (R - c_allee)*(1-(R/K))

        ## For exmports and imports, I assume delta is 1 and K_L (demand by population) is also 1
        exports <- delta * (t(C_ij) * t(A_ij)) %*% R
        imports <- delta *(C_ij * A_ij) %*% R

        dR <- growth  - exports + imports + rnorm(n=n, mean = 0, sd = 1)
        return(list(dR)) #,  exports, imports
    })
}

## J190415: Don't add exports and imports yet because in the model they are not modeled pair-wise, they are agregated. One would need to dissagregate or store the matrix instead of the vector. Where the vector now is the sum of the exports per species/country, while the matrix is the raw values.


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


## run the model
print(system.time(
    out <- ode(
      y = R, times = times, parms = params, func = resource,
      method = "lsoda" ## see help("ode") for more methods
    )
))


df <- out %>% as_data_frame() %>%
    gather(key = "species", value = "population", 2:(n+1))


## plot result
quartz(width = 4, height = 4)
df %>%
    ggplot(aes(x=time, y=population)) +
    geom_line(aes(color = species), size = 0.25, show.legend = T) +
    theme_light()

## Test of CCM as identification method:
## follow my own code from 04-DataExploration.R
lib <- c(1,2500)
pred <- c(2501, 10001)
## Embedding dimension:
# emb <- list()
# for (i in 2:dim(out)[2]){
#     emb[[i-1]] <- out[c(1,i)] %>%
#         simplex()
# }
simplex_output <- out %>%
    as.data.frame() %>%
    select(-time) %>%
    map(simplex, lib, pred)

bestE <- simplex_output %>% map_dbl(~ .$E[which.max(.$rho)] )

sps <- c(1:5)

for (i in seq_along(sps)){simplex_output[[i]]$species <- sps[i]}

simplex_output <- simplex_output %>%
  bind_rows()

simplex_output %>%
  ggplot(aes(E, rho)) +
  geom_line(aes(group = species), size = 0.3) +
  labs(x="Embeding dimension (E)", y = "Forecast skill (rho)") +
  theme_light(base_size = 9)

## Prediction decay
prediction_decay <- out %>%
  as.data.frame() %>%
  select(-time) %>%
  map(simplex, lib, pred, E = 2, tp = 1:10)

for (i in seq_along(sps)){prediction_decay[[i]]$species <- sps[i]}

prediction_decay %>%
  bind_rows() %>%
  ggplot(aes(tp, rho)) +
  geom_line(aes(group = species), size = 0.3) +
  labs(x="Time to prediction", y = "Forecast skill (rho)") +
  theme_light(base_size = 9)

## Non-linearity
non_linear <- out %>%
  as.data.frame() %>%
  select(-time) %>%
  as.list() %>%
  map(., .f = s_map, E = 2, lib = lib, pred = pred) # this takes ages.

## does not work with map2 - for now fixing E = 2
## J190416: Note that E is a list or vector on which the function also parellize, so it has to be before the .function declaration.

# print(system.time((s_map(out[[1]], lib, pred, E = 2))))
for (i in seq_along(sps)){non_linear[[i]]$species <- sps[i]}

non_linear %>%
  bind_rows() %>%
  ggplot(aes(theta, rho)) +
  geom_line(aes(group = species), size = 0.3) +
  labs(x="Nonlinearity (theta)", y = "Forecast skill (rho)") +
  theme_light(base_size = 9)

########

ind <- crossing(
    lib_column = colnames(out)[-1],
    target_column = colnames(out)[-1])

ind <- ind %>% filter(lib_column != target_column)

rho_list <- map2(
    .x = ind$lib_column, .y = ind$target_column ,
    .f = ~ ccm(block = out, E = 2,
        lib_column = .x, target_column = .y,
        lib_sizes = seq(10,dim(out)[1], by = 100),
        replace = FALSE, silent = TRUE,
        random_libs = FALSE)
    )
