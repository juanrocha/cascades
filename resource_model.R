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

## J190415: Don't add exports and imports yet because in the model they are not 
# modeled pair-wise, they are agregated. One would need to dissagregate or store
# the matrix instead of the vector. Where the vector now is the sum of the exports 
# per species/country, while the matrix is the raw values.


n <- 5 ## number of countries and resources
R <- runif(n, 10, 100) ## 10 resources with random numbers
alpha <- runif(n, 0.01, 0.1) ## growth rate
c_allee <- runif(n, 10, 30) ## Allee parameter: food necessary to feed L_i population in country i
K <- rep(100, times = n) ## carrying capacity

C_ij <- matrix(runif(n^2, min = 0.2, max = 0.8), ncol = n)
diag(C_ij) <- 0

A_ij <- matrix(rbinom(n^2, 1, prob = 0.5), ncol = n)
diag(A_ij) <- 0

delta <-  1 # resource depletion coefficient (I dont' get why is necessary)

times <- seq(from = 0, to = 100, by = 0.01)
params <- list(alpha = alpha, c_allee = c_allee, K = K, C_ij = C_ij)

## run the model
print(system.time(
    out <- ode(
      y = R, times = times, parms = params, func = resource,
      method = "lsodes" ## see help("ode") for more methods
    )
))


df <- out %>% as_data_frame() %>%
    gather(key = "species", value = "population", 2:(n+1))


## plot result
# quartz(width = 4, height = 4)
df %>%
    ggplot(aes(x=time, y=population)) +
    geom_line(aes(color = species), size = 0.25, show.legend = T) +
    theme_light()

image(t(C_ij * A_ij))


## Test of CCM as identification method:
## follow my own code from 04-DataExploration.R
## J190416: It's taking a lot of time to run the algorithms of EDM. One of the 
# issues is that I have in my synthetic data 10000 calculations for 100 years due
# to the time step of 'times'. Sample that and use a smaller dataset for the causality tests.

df2 <- out[seq(from = 1, to = 10001, by = 100),]
dim(df2)
lib <- c(1,50)
pred <- c(51, 101)
## Embedding dimension:
# emb <- list()
# for (i in 2:dim(out)[2]){
#     emb[[i-1]] <- out[c(1,i)] %>%
#         simplex()
# }

## Normalize the time series:
normalize <- function(x){
  y <- (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
  return(y)
}


df2[,-1] <- apply(df2[,-1], 2, normalize)



simplex_output <- df2 %>%
    as.data.frame() %>%
    select(-time) %>%
    map(simplex, lib, pred, E=seq(from = 1, to = 10, by = 1))

bestE <- simplex_output %>% map_dbl(~ .$E[which.max(.$rho)] )

sps <- c(1:5)

for (i in seq_along(sps)){simplex_output[[i]]$species <- sps[i]}

simplex_output %>%
  bind_rows() %>%
  ggplot(aes(E, rho)) +
  geom_line(aes(group = species), size = 0.3) +
  labs(x="Embeding dimension (E)", y = "Forecast skill (rho)") +
  theme_light(base_size = 9)

## Prediction decay
prediction_decay <- df2 %>%
  as.data.frame() %>%
  select(-time) %>%
  map(simplex, lib, pred, E = 1, tp = seq(1,20,by= 1))

for (i in seq_along(sps)){prediction_decay[[i]]$species <- sps[i]}

prediction_decay %>%
  bind_rows() %>%
  ggplot(aes(tp, rho)) +
  geom_line(aes(group = species), size = 0.3) +
  labs(x="Time to prediction", y = "Forecast skill (rho)") +
  theme_light(base_size = 9)

## Non-linearity
non_linear <- df2 %>%
  as.data.frame() %>%
  select(-time) %>%
  as.list() %>%
  map(., .f = s_map, E = 3, lib = lib, pred = pred) # this takes ages.

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

## test: about 2mins for one relationship. For a small model with 5sp we have 20 relationships.
print(system.time(  
  test_x <- ccm(
    block = df2, E = 10, lib = lib, pred = pred, 
    lib_column = 3, target_column = 4, num_neighbors = 0.5,
    lib_sizes = seq(11, dim(df2)[1], by = 10),
    first_column_time = TRUE,
    replace = TRUE, silent = TRUE,
    random_libs = TRUE, num_samples = 25,
    stats_only = TRUE)
  )
)


### parallel predicting all species

print(system.time(  
  rho_list <- map2(
    .x = ind$lib_column, .y = ind$target_column ,
    .f = ~ ccm(
      block = df2, E = 5, lib = lib, pred = pred, 
      lib_column = .x, target_column = .y,
      lib_sizes = seq(11, dim(df2)[1], by = 10),
      first_column_time = TRUE, #num_neighbors = 0.5,
      replace = TRUE, silent = TRUE,
      random_libs = FALSE, num_samples = 30,
      stats_only = TRUE)
  )
))


## t-test for each relationship:
t_tests <- map(.x = rho_list, safely(
  .f = ~ t.test(.x$rho, alternative = "greater", mu = 0, na.action = na.exclude))
)

# p_vals <- map(t_tests, function(x) x$result$p.value)
t_tests <- transpose(t_tests)
fail <- t_tests$error %>% map_lgl(is_null) %>% unlist() ## should be TRUE for all
fail

ind <- ind %>%
  mutate(
      rho = map_dbl(.x = rho_list, .f = ~ mean(.x$rho, na.rm = TRUE)),
      rho_t = map_dbl(.x = t_tests$result, function(x) x$estimate ),
      p_value = map_dbl(.x = t_tests$result, function(x) x$p.value ),
      detection = ifelse(p_value <0.05 & rho > 0.05, TRUE, FALSE)
  )

ind$p_value < 0.05

### So the million dollar question: Can I recover the adjacency matrix?

ind %>%
  ggplot(., aes(lib_column, target_column)) +
  geom_tile(aes(fill = detection)) + 
  theme_light(base_size = 9)

quartz(width = 4, height =4)
image(A_ij)

round(C_ij * A_ij, 2)


ind %>% select(1,2, detection) %>% spread(key =  lib_column, value = detection)

## It doesn't work.
## I think I have the math of the model wrong. If it's a model of difussion given 
# state variables (species / countries) and what defines the difussion is the adjacency matrix
# then the model should be a PDE not a ODE. In the deSolve tutorial there is a 
# PDE example that finishes in <1sec with 2000 state vars, why mine with 5 takes longer?

save(out, df, df2, rho_list, ind, A_ij, C_ij,file = "resource_model_experiment.RData")


##############3
## Block model
#############3








