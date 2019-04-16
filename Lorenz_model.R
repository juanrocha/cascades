## Example with Lorenz model

library(deSolve)
library(deTestSet)
library(tidyverse)
library(rEDM)

## set up model, following Soetaert 2012

a <--8/3;b<--10;c<-28           # parameters
yini <- c(X = 1, Y = 1, Z = 1)  # Initial conditions

## function:
Lorenz <- function (t, y, parms) { 
    with(as.list(y), {
        dX <- a * X + Y * Z
        dY <- b * (Y - Z)
        dZ <- -X * Y + c * Y - Z 
        list(c(dX, dY, dZ)) })
}

# set up
times <- seq(from = 0, to = 100, by = 0.01)
out <- ode(y = yini, times = times, func = Lorenz,
           parms = NULL)

## plot:
plot(out, lwd = 2)
plot(out[,"X"], out[,"Y"], type = "l", xlab = "X",
     ylab = "Y", main = "butterfly")

## Does CCM really works?

df2 <- out[seq(from = 1, to = 10001, by = 10),]
dim(df2)
lib <- c(1,500)
pred <- c(501, 1001)

## embeding
emb <- df2 %>%
    as.data.frame() %>%
    select(-time) %>%
    map(simplex, lib, pred, E=seq(from = 1, to = 100, by = 10))

bestE <- emb %>% map_dbl(~ .$E[which.max(.$rho)] )

sps <- c('X','Y',"Z")

for (i in seq_along(sps)){emb[[i]]$species <- sps[i]}

emb %>%
    bind_rows() %>%
    ggplot(aes(E, rho)) +
    geom_line(aes(group = species), size = 0.3) +
    labs(x="Embeding dimension (E)", y = "Forecast skill (rho)") +
    theme_light(base_size = 9)


## Prediction decay
prediction_decay <- df2 %>%
    as.data.frame() %>%
    select(-time) %>%
    map(simplex, lib, pred, E = 11, tp = seq(1,50,by= 5))

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
    map(., .f = s_map, E = 11, lib = lib, pred = pred) 

# print(system.time((s_map(out[[1]], lib, pred, E = 2))))
for (i in seq_along(sps)){non_linear[[i]]$species <- sps[i]}

non_linear %>%
    bind_rows() %>%
    ggplot(aes(theta, rho)) +
    geom_line(aes(group = species), size = 0.3) +
    labs(x="Nonlinearity (theta)", y = "Forecast skill (rho)") +
    theme_light(base_size = 9)


#####

ind <- crossing(
    lib_column = colnames(out)[-1],
    target_column = colnames(out)[-1])

ind <- ind %>% filter(lib_column != target_column)

rho_list <- map2(
    .x = ind$lib_column, .y = ind$target_column ,
    .f = ~ ccm(block = df2, E = 11,
               lib_column = .x, target_column = .y,
               lib_sizes = seq(10, dim(df2)[1], by = 10),
               first_column_time = TRUE,
               replace = TRUE, silent = FALSE,
               random_libs = TRUE)
)

## t-test for each relationship:
t_tests <- map(.x = rho_list, safely(
    .f = ~ t.test(.x$rho, alternative = "greater", mu = 0, na.action = na.exclude))
)

# p_vals <- map(t_tests, function(x) x$result$p.value)
t_tests <- transpose(t_tests)
fail <- t_tests$error %>% map_lgl(is_null) %>% unlist()

ind <- ind %>%
    mutate(
        rho = map_dbl(.x = rho_list, .f = ~ mean(.x$rho, na.rm = TRUE)),
        rho_t = map_dbl(.x = t_tests$result, function(x) x$estimate ),
        p_value = map_dbl(.x = t_tests$result, function(x) x$p.value ),
        detection = ifelse(p_value <0.05 & rho > 0.1, TRUE, FALSE)
    )


ind


### J190416: So it does not work. At least not the pair-wise brute force computation. 
#Takes ages and does not produce anything reliable.


block <- df2 %>% make_block()

block_output <- block_lnlp(
    block,
    lib = c(1,500), pred = c(501, 1001),
    columns = c("X", "Y", "Z"), target_column = "Z",
    method = "s-map", theta = 2, stats_only = FALSE,
    first_column_time = TRUE, save_smap_coefficients = TRUE, silent = FALSE
)

smap_coeffs <- block_output$smap_coefficients[[1]]
str(smap_coeffs)


## from their code
predictions <- block_output$model_output[[1]]
t <- predictions$time

plot(t, predictions$obs, type = "l", col = "black", ylab = "x", xlab = "")
lines(t, predictions$pred, lty = 2)
legend("topright", legend = c("observed", "predicted"), lty = c(1, 2), bty = "n")

plot(t, smap_coeffs[, 1], type = "l", col = "red", ylab = "effect of x", xlab = "")
plot(t, smap_coeffs[, 2], type = "l", col = "blue", ylab = "effect of y", xlab = "")
plot(t, smap_coeffs[, 3], type = "l", col = "magenta", ylab = "effect of z",
     xlab = "")


block_output$model_output[[1]] %>% ggplot(aes(obs,pred)) +
    geom_point() + 
    geom_abline(aes(intercept = 0, slope =1))

## So it does work with the block option, but it is not clear how to identify or not causal effects.
## the coeffs show some variables are stronger than others over time, but there is no yes/no test.
## the algorith is good at prediction.

