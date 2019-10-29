

# ------ Libraries -------
library(rEDM)
library(deSolve)
library(pracma)
library(viridis)

# ------ Utilities ------

l_v_discrete <- function(current_state, r, A){
  next_state <- matrix(NA, nrow = nrow(A), ncol = 1);
  next_state <- current_state*(r-A%*%current_state);
  return(next_state)
}

make_timeseries <- function(dynamics, start_state, steps, r, A){
  count <- 1;
  dim <- dim(A)[1]
  vars <- rep(NA, dim)
  for(i in 1:dim){
    vars[i] <- paste("x", i, sep = "")
  }
  timeseries <- matrix(NA, nrow = nrow(A), ncol = steps+1);
  timeseries[,1] <- start_state;
  while(count <= steps){
    count <- count + 1;
    timeseries[,count] <- dynamics(timeseries[,count-1], r, A);
  }
  timeseries <- t(as.data.frame(timeseries))
  colnames(timeseries) <- vars
  return(timeseries)
}

get_all_ccm <- function(data, mean = TRUE, time, ...){
  results <- list()
  count <- 0
  vars <- colnames(data)
  if(as.logical(time)){vars <- vars[-1]}
  comb <- rep(NA, length(vars)*(length(vars)-1))
  for(lib in vars){
    for(target in vars){
      if(strcmp(lib,target)){next}
      count <- count + 1
      comb[count] <- paste(lib, target, sep = "_")
      ccm_loop <- ccm(data, lib_column = lib, 
                      target_column = target, first_column_time = time, ...)
      if(as.logical(mean)){
        results[[count]] <- ccm_means(ccm_loop)
      }
      else{
        results[[count]] <- ccm_loop
      }
    }
  }
  return(list(comb, results))
}


add_timelags <- function(data, columns, number_of_lags){
  d <- dim(data)[[1]]
  if(!is.numeric(columns)){
    columns <- match(columns, colnames(data))
  }
  columns <- sort(columns, decreasing = TRUE)
  names <- colnames(data)
  for(i in columns){
    lags <- matrix(NA, d, number_of_lags)
    colnames(lags) <- 1:number_of_lags
    for(j in 1:number_of_lags){
      lags[-(1:j),j] <- data[1:(d-j), i]
      colnames(lags)[j] <- paste(names[i], "-", j, "tau", sep= "")
    }
    if(i== dim(data)[[2]]){
      data <- cbind(data[,1:i],lags)
    }
    else{
      data <- cbind(data[,1:i],lags,data[,(i+1):(dim(data)[[2]])])
    }
  }
  return(data)
}


get_best_embedding <- function(data, dim_try, ...){
  simplex_out <- simplex(data, E = dim_try, ...)
  e <- which.max(simplex_out$rho)
  return(e)
}


# ------ Following Sugihara et al: Detecting Causality in Complex Systems -----


# Figure 3

r <- c(3.7, 3.7);
A <- matrix(c(3.7, 0.32, 0.05 ,3.7),2,2);
start <- c(0.2, 0.4)

data_fig3a <- make_timeseries(dynamics = l_v_discrete, start_state = start, steps = 2000, r = r, A = A)

ccm_out_fig3a <- get_all_ccm(data = data_fig3a, time = FALSE, E = 2, lib_sizes = c(seq(0,201,10),seq(250, 3250, by = 250)), num_samples = 40, 
                             random_libs = TRUE, replace = TRUE)

combs_fig3a <- ccm_out_fig3a[[1]]
ccm_res_fig3a <- ccm_out_fig3a[[2]]

par(mfrow = c(1,1))
plot(ccm_res_fig3a[[1]]$lib_size, ccm_res_fig3a[[1]]$rho, col = "red", type = "l", ylim = c(0,1), xlim = c(0, 3250), xlab = "Library Size L", ylab = "Prediction Skill rho")
lines(ccm_res_fig3a[[2]]$lib_size, ccm_res_fig3a[[2]]$rho, col = "blue", type = "l", ylim = c(0,1), xlim = c(0, 3250), xlab = "Library Size L", ylab = "Prediction Skill rho")
legend("bottomright", c("Y predicted using X", "X predicted using Y"), col = c("red", "blue"), lty = 1)
title("Difference in causation strength")

#

r <- c(3.7, 3.7);
A <- matrix(c(3.7, 0.32, 0,3.7),2,2);
start <- c(0.2, 0.4)

data_fig3cd <- make_timeseries(dynamics = l_v_discrete, start_state = start, steps = 2000, r = r, A = A)

par(mfrow = c(1,1))
plot(0:200, data_fig3cd[1:201,1], col = "blue", type = "l", xlab = "Time", ylab = "Value")
lines(0:200, data_fig3cd[1:201, 2], col = "red", type = "l")
legend("topright", c("x1", "x2"), col = c("blue","red"), lty = 1)
title("Timeseries x1 and x2")

# 1000 cross map estimates are plotted using attractor reconstructions with L = 1000 and E = 2

ccm_out_fig3cd <- get_all_ccm(data = data_fig3cd, mean = FALSE, time = FALSE, lib = c(1,1000), pred = c(1001, 2000), E = 2, lib_sizes = 1000,stats_only = FALSE, num_samples = 1)

combs_fig3cd <- ccm_out_fig3cd[[1]]
ccm_res_fig3cd <- ccm_out_fig3cd[[2]]

par(mfrow = c(1,2),oma = c(0, 0, 2, 0))
plot(x = data_fig3cd[1001:2000,1], y = ccm_res_fig3cd[[2]]$model_output[[1]]$pred, xlim = c(0,1), ylim = c(0,1), xlab = "X(t) observed", ylab = "X(t) predicted", col = "steelblue", pch = '.', cex = 2)
title("Causality from X towards Y")
plot(data_fig3cd[1001:2000,2], ccm_res_fig3cd[[1]]$model_output[[1]]$pred, xlim = c(0,1), ylim = c(0,1), xlab = "Y(t) observed", ylab = "Y(t) predicted", col = "firebrick3", pch = '.', cex = 2)
title("No causality from Y towards X")
mtext("Two species Lotka Volterra System with one-sided causality", outer = TRUE, cex = 1.5)


# Figure 4 (5 Species Model)

r <- c(4, 3.1, 2.12, 3.8, 4.1);
A <- matrix(c(4, 0.31, -0.636, 0.111, 0.082, 2, 3.1, -0.636, 0.111, 0.111, 0.4, 0.93, 2.12, -0.131, 0.125, 0,0,0,3.8,0,0,0,0,0,4.1),5,5);
start <- c(0.1,0.2,2/3,0.5,0.2)

data_fig4 <- make_timeseries(dynamics = l_v_discrete, start_state = start, steps = 1000, r = r, A = A)

ccm_out_fig4 <- get_all_ccm(data = data_fig4, time = FALSE, E = 5, lib_sizes = c(seq(10,90,10),seq(100, 800, by = 100)), num_samples = 40, 
                               random_libs = TRUE, replace = TRUE)

combs <- ccm_out_fig4[[1]]
ccm_res <- ccm_out_fig4[[2]]

par(mfrow = c(1,1))
plot(ccm_res[[1]]$lib_size, ccm_res[[1]]$rho, col = "red", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
for(i in c(2,5,6,9,10)){
  lines(ccm_res[[i]]$lib_size, ccm_res[[i]]$rho, col = "red", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
}
for(i in 13:15){
  lines(ccm_res[[i]]$lib_size, ccm_res[[i]]$rho, col = "blue", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
}
for(i in 17:19){
  lines(ccm_res[[i]]$lib_size, ccm_res[[i]]$rho, col = "pink", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
}
for(i in c(3,7,11)){
  lines(ccm_res[[i]]$lib_size, ccm_res[[i]]$rho, col = "cyan", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
}
for(i in c(4,8,12)){
  lines(ccm_res[[i]]$lib_size, ccm_res[[i]]$rho, col = "yellow", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
}
for(i in c(16,20)){
  lines(ccm_res[[i]]$lib_size, ccm_res[[i]]$rho, col = "black", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
}
legend("topright" ,lty = 1, combs[c(1,2,5,6,9,10,13:15,17:19,3,7,11,4,8,12,16,20)], col = c(rep("red",6),rep("blue",3),rep("pink", 3), rep("cyan", 3), rep("yellow", 3), rep("black", 2)))
title("CCM for 5 Species Model (Fig 4)")




# ------ Lotka-Volterra -----------

# Interdependences and growth rate
r <- rep(0.5, 3);
A <- matrix(c(-0.2,0,0,-0.1,-0.4,-0.1,0.2,0,-0.1),3,3)
parameters <- list(r,A);
state <- c(X = 1, Y = 1, Z = 1);

Lotka_Volterra <- function(t, state, parameters){
  D <- state*(r+A%*%state)
  list(c(D[1],D[2],D[3]))
}

times <- seq(0, 100, by = 0.1)
out_lv <- ode(y = state, times = times, func = Lotka_Volterra, parms = parameters)
set.seed(050619)
noise_lv <- matrix(rnorm(n = 1001*3, mean = 0, sd = 0.5), 1001, 3)
out_noise_lv <- out_lv[,2:4] + noise_lv

data <- out_noise_lv

ccm_out_lv_noise <- get_all_ccm(data = data, time = FALSE, E = 3, lib_sizes = c(seq(10,90,10),seq(100, 800, by = 100)), num_samples = 100, 
                                random_libs = TRUE, replace = TRUE)

combs <- ccm_out_lv_noise[[1]]
ccm_res <- ccm_out_lv_noise[[2]]

# x_xmap_y high and converging if y causes x
par(mfrow = c(1,1))
plot(ccm_res[[1]]$lib_size, ccm_res[[1]]$rho, col = "red", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
lines(ccm_res[[2]]$lib_size, ccm_res[[2]]$rho, col = "blue", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
lines(ccm_res[[3]]$lib_size, ccm_res[[3]]$rho, col = "pink", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
lines(ccm_res[[4]]$lib_size, ccm_res[[4]]$rho, col = "green", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
lines(ccm_res[[5]]$lib_size, ccm_res[[5]]$rho, col = "magenta", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
lines(ccm_res[[6]]$lib_size, ccm_res[[6]]$rho, col = "cyan", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))

legend("topright" ,lty = 1, combs, col = c("red","blue","pink","green","magenta", "cyan"))
title("CCM Lotka-Volterra")


# ------ Lorenz System ------------

parameters <- c(s = 10, r = 28, b = 8/3)
state <- c(X = 0, Y = 1, Z = 1)

Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- s * (Y - X)
    dY <- X * (r - Z) - Y
    dZ <- X * Y - b * Z
    list(c(dX, dY, dZ))
  })
}

times <- seq(0, 100, by = 0.1)
out_l <- ode(y = state, times = times, func = Lorenz, parms = parameters)


# Rotating the system should lead to Z having a valid shadow manifold
data_rot <- out_l[,3:4]
data_rot <- t(apply(data_rot, 1, function(x) matrix(c(cos(pi/3),-sin(pi/3),sin(pi/3),cos(pi/3)),2,2)%*%x))
data_rot <- cbind(out_l[,1:2], data_rot)
colnames(data_rot) <- c("time", "X", "Y", "Z")

# choosing which data version to use
data <- data_rot # out_l or data_rot

ccm_out_lorenz <- get_all_ccm(data, time = TRUE, E = 3, lib_sizes = c(seq(10,90,10),seq(100, 800, by = 100)), num_samples = 20, 
                              random_libs = TRUE, replace = TRUE)

combs <- ccm_out_lorenz[[1]]
ccm_res <- ccm_out_lorenz[[2]]

# x_y high and converging if y causes x
par(mfrow = c(1,1))
plot(ccm_res[[1]]$lib_size, ccm_res[[1]]$rho, col = "red", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
lines(ccm_res[[2]]$lib_size, ccm_res[[2]]$rho, col = "blue", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
lines(ccm_res[[3]]$lib_size, ccm_res[[3]]$rho, col = "pink", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
lines(ccm_res[[4]]$lib_size, ccm_res[[4]]$rho, col = "green", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
lines(ccm_res[[5]]$lib_size, ccm_res[[5]]$rho, col = "magenta", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
lines(ccm_res[[6]]$lib_size, ccm_res[[6]]$rho, col = "cyan", xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))

legend("topright" ,lty = 1, combs, col = c("red","blue","pink","green","magenta", "cyan"))
title("CCM for Lorenz Attractor")

#

# ------ Extended Lorenz System ----------

# (to include uncoupled variables)

parameters <- c(s = 10, r = 28, b = 8/3, a = 0.1, c = 0.5)
state <- c(X = 0, Y = 1, Z = 1, U = 0.5, V = 2)

Lorenz2 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- s * (Y - X)
    dY <- X * (r - Z) - Y
    dZ <- X * Y - b * Z
    dU <- a * Z + X
    dV <- c * U + Y
    list(c(dX, dY, dZ, dU, dV))
  })
}

times <- seq(0, 100, by = 0.1)
timeseries_lorenz2 <- ode(y = state, times = times, func = Lorenz2, parms = parameters)
set.seed(050619)
noise <- matrix(rnorm(n = 1001*5, mean = 0, sd = 1), 1001, 5)

noise_lorenz2 <- timeseries_lorenz2
noise_lorenz2[,2:6] <- timeseries_lorenz2[,2:6] + noise

data_lorenz2 <- timeseries_lorenz2

ccm_out_lorenz2 <- get_all_ccm(data = data_lorenz2, time = TRUE,E = 5, lib_sizes = c(seq(10,90,10),seq(100, 800, by = 100)), num_samples = 20, 
                               random_libs = TRUE, replace = TRUE)

combs_lorenz2 <- ccm_out_lorenz2[[1]]
ccm_res_lorenz2 <- ccm_out_lorenz2[[2]]

palette(plasma(4))
par(mfrow = c(1,1))

l <- dim(data_lorenz2)[[2]] - 1
for(v in 1:l){
  plot(ccm_res_lorenz2[[(v-1)*(l-1)+1]]$lib_size, ccm_res_lorenz2[[(v-1)*(l-1)+1]]$rho, col = 1, xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
  for(i in 2:(l-1)){
    lines(ccm_res_lorenz2[[(v-1)*(l-1)+i]]$lib_size, ccm_res_lorenz2[[(v-1)*(l-1)+i]]$rho, col = i, xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 800), ylim = c(0,1))
  }
  legend("topright" ,lty = 1, combs_lorenz2[((v-1)*(l-1)+1):((v-1)*(l-1)+l-1)], col = 1:4)
  title("CCM for Extended Lorenz System")
  readline(prompt="Press [enter] to continue")
}

#

# ------ Couples Rossler-Lorenz System --------

parameters <- c(s = 10, r = 28, b = 8/3, a = 6, C = 2)
state <- c(X = 0, Y = 1, Z = 1, U = 0.5, V = 2, W = 0.5)

rossler_lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- s * (Y - X)
    dY <- X * (r - Z) - Y + C * V^2
    dZ <- X * Y - b * Z 
    dU <- -a * (V + W)
    dV <- a * (U + 0.2 * V)
    dW <- a * (0.2 + W * (U - 5.7))
    list(c(dX, dY, dZ, dU, dV, dW))
  })
}


times <- seq(0, 100, by = 0.1)
timeseries_rossler_lorenz <- ode(y = state, times = times, func = rossler_lorenz, parms = parameters)
set.seed(60619)
noise <- matrix(rnorm(n = 1001*6, mean = 0, sd = 1), 1001, 6)

noise_rossler_lorenz <- timeseries_rossler_lorenz
noise_rossler_lorenz[,2:7] <- timeseries_rossler_lorenz[,2:7] + noise

# choose version of data 
data_rossler_lorenz <- timeseries_rossler_lorenz


ccm_out_rossler_lorenz <- get_all_ccm(data = data_rossler_lorenz, time = TRUE, E = 5, lib_sizes = c(seq(10,90,10),seq(100, 300, by = 100)), num_samples = 30, 
                                      random_libs = TRUE, replace = TRUE)

combs_rossler_lorenz <- ccm_out_rossler_lorenz [[1]]
ccm_res_rossler_lorenz <- ccm_out_rossler_lorenz [[2]]

l <- dim(data_rossler_lorenz)[[2]]-1
palette(plasma(l-1))
par(mfrow = c(1,1))

for(v in 1:l){
  plot(ccm_res_rossler_lorenz[[(v-1)*(l-1)+1]]$lib_size, ccm_res_rossler_lorenz[[(v-1)*(l-1)+1]]$rho, col = 1, xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 300), ylim = c(0,1))
  for(i in 2:(l-1)){
    lines(ccm_res_rossler_lorenz[[(v-1)*(l-1)+i]]$lib_size, ccm_res_rossler_lorenz[[(v-1)*(l-1)+i]]$rho, col = i, xlab = "Lib Size", ylab = "Prediction Skill", type = "l", xlim =c(10, 300), ylim = c(0,1))
  }
  legend("topright" ,lty = 1, combs_rossler_lorenz[((v-1)*(l-1)+1):((v-1)*(l-1)+l-1)], col = 1:(l-1))
  title("CCM for Rossler-Lorenz with noise")
  readline(prompt="Press [enter] to continue")
}



data_rl_timelagged <- add_timelags(data_rossler_lorenz, columns = 2:7, number_of_lags = 3)

u_from_lorenz <- block_lnlp(data_rl_timelagged, columns = 2:13, target_column = "U")
v_from_lorenz <- block_lnlp(data_rl_timelagged, columns = 2:13, target_column = "V")
w_from_lorenz <- block_lnlp(data_rl_timelagged, columns = 2:13, target_column = "W")

x_from_rossler <- block_lnlp(data_rl_timelagged, columns = 14:25, target_column = "X")
y_from_rossler <- block_lnlp(data_rl_timelagged, columns = 14:25, target_column = "Y")
z_from_rossler <- block_lnlp(data_rl_timelagged, columns = 14:25, target_column = "Z")

