## Figures draft 200331
## Juan Rocha
## juan.rocha@su.se

library(deSolve)
# library(deTestSet)
library(tidyverse)
library(tsibble)
library(patchwork)

## Avoiding zero values:
# This is where we define your event function
# Add this directly above your call to ode()
posfun <- function(t, y, parms){
    with(as.list(y), {
        y[which(y<0)] <- 0
        return(y)
    })
}

## Pollution model:
pollution <- function(t, y, params){
    with(as.list(c(y, params)), { ## adding the c() eliminates the error of not finding object 'u'
        
        pollutant <- u - (b*y) + v * (y^beta/(Y^beta + y^beta))
        outflow <-  Dij %*% y
        inflow <-  t(Dij) %*% y
        
        dy <- pollutant + (inflow - outflow)
        
        return(list(c(dy)))
    })
}

## set up flux matrix
n <- 23
Dij <- matrix(
    #runif(n^2, min = 0.02, max = 0.05), 
    0,
    ncol = n, nrow = n)
diag(Dij) <- 0

## Parameters
#params <- list(u = 0.5, v = 1.7, Y = 1, beta = 4, Dij = Dij )
params <- list(
    b = 2.2,
    v = 10,
    u = 0.5,
    Y = 2.5,
    beta = 4 ,
    Dij = Dij
)


## The initial conditions are points in the phase space along the x axis
yini <- seq(from = 0, to = 11, by = 0.5) #runif(n, 0, 5)
times <- seq(from = 0, to = 100, by = 0.5)

## run the model
print(system.time(
    out <- ode(
        y = yini, times = times, parms = params, func = pollution,
        method = "bdf", ## see help("ode") for more methods
        events=list(func = posfun, time = times)
    )
))

block <- rEDM::make_block(out,  max_lag = 2, tau = 1)


### this doesn't work
df_x <- block[, -1] %>% as_tibble() %>% 
    select(1,2, seq(3,47,2)) %>%
    pivot_longer(cols = 3:last_col(), values_to = "state", names_to = "x") %>%
    select(-time_1)

df_x1 <- block[, -1] %>% as_tibble() %>% 
    select(1,2,seq(4,48,2)) %>%
    pivot_longer(cols = 3:last_col(), values_to = "lag_1", names_to = "x") %>%
    mutate(x = str_remove(x, "_1")) %>%
    select(-time_1)

df <- df_x %>% left_join(df_x1) 


g2 <- df %>% 
    mutate(x = as.numeric(x)/2 - 0.5) %>%
    ggplot(aes(lag_1, state), group = as.factor(x)) +
    geom_point(aes(color = as.factor(x), alpha = time)) +
    geom_line(aes(color = as.factor(x)), alpha = 0.24, size = 0.5) +
    geom_density2d(aes(color= state), h = 1, size = 0.2)+
    coord_fixed() + labs(x = "x", y = bquote('x'[t+1]) ) + # another option is expression('x'[t+1])
    # facet_wrap(~x)
    scale_color_discrete("Initial conditions", 
                         guide = guide_legend(title.position = "top", nrow = 3)) +
    labs(tag = "B") +
    theme_light(base_size = 7) + 
    theme(legend.position = c(0.5, 0.8), legend.direction = "horizontal",
          legend.key.size  = unit(0.15, "cm"), legend.text = element_text(size = 5))



## small dataframe
df_mini <- tibble( x = seq(0, 11.5, by = 0.5))

## parameters on the working environment
v = 10    # this is like a carrying capacity
u = 0.5
Y = 2.5
beta = 4
Dij = Dij
b = 2.2


## components of the function, with and without diffusion
pollutant <- function(y) { v * (y^beta/(Y^beta + y^beta))}  ## Original model formulated by ASC
pollutant_sedimentation <- function(y) { v * (y^beta/(Y^beta + y^beta)) } ## adding the saturation parameter
management <- function(y){- u + (b*y)} ## now the diffusion term is in the management line.
management_diffusion <- function(y){ - ((0.01 * y) - (0.05 * y)) - u + (b*y)} ## now the diffusion term is in the management line.
## plot
g1 <- ggplot(data = df_mini, aes(x)) +
    #stat_function(fun = pollutant, color = "blue", linetype = 2) +
    stat_function(fun = pollutant_sedimentation, color = "purple") +
    #stat_function(fun = management_diffusion, color = "purple") +
    stat_function(fun = management, color = "dodgerblue") +
    annotate("text",x = c(9,9), y = c(14,9),
             label = c("- u + (s*x)", "v * (x^beta/(X^beta + x^beta))"),
             parse = TRUE, color = c("dodgerblue", "purple")) +
    #geom_hline(yintercept = 0, color = "gray") + 
    labs(x = "x", y = "f(x)", tag = "A") + #coord_fixed() +
    ylim(NA, 15) +
    theme_light(base_size = 8) #+ theme(axis.text = element_blank())


g1 + g2 + plot_layout(ncol = 2, nrow = 1, widths = c(1,1))

ggsave(filename= "fig1_draft.png", device = "png", width = 7, height = 3)


### Resource model:

# Resource models:
## small dataframe
df_mini <- tibble( x = seq(0, 11.5, by = 0.5))
# parameters:
## Parameters: The delta and d parameters are set here by hand to replicate ASC results.
##  Later they are set randomly on the fly.
R <- 1          # Growth rate
K <- 10         # carryng capacity
delta <- 0.05   # share of stock that diffuses to other plots
a <- 4          # maximum level after crossing the threshold
c <- 1.7        # consumption efficiency
d <- 0.02       # total diffiusion from all other plots to i
# Dij <- matrix(c(0, delta, d, 0), nrow = 2, byrow = FALSE)

## components of the function, with and without diffusion
growth <- function(x) {R*x*(1-(x/K))}
growth_diffusion <- function(x) {R*x*(1-(x/K)) + d*x - delta *x}
consumption <- function(x){c *((x^a) / (1 + (x^a) ))}
## plot
g3 <- ggplot(data = df_mini, aes(x)) +
    stat_function(fun = growth, color = "blue") +
    # stat_function(fun = growth_diffusion, color = "orange") +
    stat_function(fun = consumption, color = "purple") +
    annotate(
        "text", x = c(6, 5), y = c(2.6, 1.5),
        label = c("R*x*(1-(x/K))", "c *((x^alpha) / (1 + (x^alpha) ))"),
        parse = TRUE, color = c("blue", "purple")
    ) +
    # geom_hline(yintercept = 0, color = "gray") + 
    labs(x = "x", y = "f(x)", tag = "A") + xlim(0,10) + ylim(NA,3) +
    theme_light(base_size = 8)

## the model
resource <- function(t, y, params){
    with(as.list(y, params), {
        growth <- R * y * (1- y/K)  ## This should recover Lodka-Volterra
        consumption <- c *((y^a) / (1 + (y^a) ))
        
        immigration <-  t(Dij) %*% y
        emmigration <-  Dij %*% y
        
        dX <- growth - consumption + (immigration - emmigration)
      
    })
}

## set up time steps
yini <- seq(from = 0, to = 11, by = 0.5) #runif(n, 0, 5)
times <- seq(from = 0, to = 100, by = 0.5)
## declare parameters
params <- list(a = a, R = R, K = K, c = c, Dij = Dij)

## run the model
print(system.time(
    out <- ode(
        y = yini, times = times, parms = params, func = resource,
        method = "bdf" ## see help("ode") for more methods
    )
))

block <- rEDM::make_block(out,  max_lag = 2, tau = 1)


### this doesn't work
df_x <- block[, -1] %>% as_tibble() %>% 
    select(1,2, seq(3,47,2)) %>%
    pivot_longer(cols = 3:last_col(), values_to = "state", names_to = "x") %>%
    select(-time_1)

df_x1 <- block[, -1] %>% as_tibble() %>% 
    select(1,2,seq(4,48,2)) %>%
    pivot_longer(cols = 3:last_col(), values_to = "lag_1", names_to = "x") %>%
    mutate(x = str_remove(x, "_1")) %>%
    select(-time_1)

df <- df_x %>% left_join(df_x1) 

g4 <- df %>% 
    mutate(x = as.numeric(x)/2 - 0.5) %>%
    ggplot(aes(lag_1, state), group = as.factor(x)) +
    geom_line(aes(color = as.factor(x)), alpha = 0.24, size = 0.5) +
    geom_point(aes(color = as.factor(x), alpha = time)) +
    geom_density2d(aes(color= state), h = 1, size = 0.2)+
    coord_fixed() + labs(x = "x", y = quote('x'[t+1]) ) + # another option is expression('x'[t+1])
    # facet_wrap(~x)
    scale_color_discrete("Initial conditions", 
                         guide = guide_legend(title.position = "top", nrow = 5)) +
    labs(tag = "B") +
    theme_light(base_size = 7) + 
    # theme(legend.position = "bottom")
    theme(legend.position = c(0.35, 0.9), legend.direction = "horizontal",
          legend.key.size  = grid::unit(0.15, "cm"), legend.text = element_text(size = 5))

g3 + g4 + plot_layout(ncol = 2, nrow = 1, widths = c(1,1))



ggsave(filename= "fig2_draft.png", device = "png", width = 7, height = 3)
