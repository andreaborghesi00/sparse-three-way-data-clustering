# install.packages("devtools")
# devtools::install_github("AndreaCappozzo/sparsemixmat") # you might need several library dependencies
# install.packages("tidyr")
# install.packages("dplyr")
# install.packages("reshape2")
library(reshape2)
library(tidyr)
library(dplyr)
library(sparsemixmat)
library(readr)
library(ggplot2)
library(gplots)

# read csv files
madrid_2001 <- read_csv("uni/statistics/seminar_project/csvs_per_year/csvs_per_year/madrid_2001.csv")
madrid_2002 <- read_csv("uni/statistics/seminar_project/csvs_per_year/csvs_per_year/madrid_2002.csv")
madrid_2003 <- read_csv("uni/statistics/seminar_project/csvs_per_year/csvs_per_year/madrid_2003.csv")
madrid_2004 <- read_csv("uni/statistics/seminar_project/csvs_per_year/csvs_per_year/madrid_2004.csv")
madrid_2005 <- read_csv("uni/statistics/seminar_project/csvs_per_year/csvs_per_year/madrid_2005.csv")
madrid_2006 <- read_csv("uni/statistics/seminar_project/csvs_per_year/csvs_per_year/madrid_2006.csv")
madrid_2007 <- read_csv("uni/statistics/seminar_project/csvs_per_year/csvs_per_year/madrid_2007.csv")
madrid_2008 <- read_csv("uni/statistics/seminar_project/csvs_per_year/csvs_per_year/madrid_2008.csv")
madrid_2009 <- read_csv("uni/statistics/seminar_project/csvs_per_year/csvs_per_year/madrid_2009.csv")

# add year
madrid_2001$Year = 2001
madrid_2002$Year = 2002
madrid_2003$Year = 2003
madrid_2004$Year = 2004
madrid_2005$Year = 2005
madrid_2006$Year = 2006
madrid_2007$Year = 2007
madrid_2008$Year = 2008
madrid_2009$Year = 2009

# intersection of columns
madrid_list <- list(
  names(madrid_2001), names(madrid_2002), names(madrid_2003), names(madrid_2004),
  names(madrid_2005), names(madrid_2006), names(madrid_2007), names(madrid_2008),
  names(madrid_2009)
)

common_columns <- Reduce(intersect, madrid_list)
length(common_columns)

# rbind, only common columns
madrid_data <- rbind(
  madrid_2001[common_columns],
  madrid_2002[common_columns],
  madrid_2003[common_columns],
  madrid_2004[common_columns],
  madrid_2005[common_columns],
  madrid_2006[common_columns],
  madrid_2007[common_columns],
  madrid_2008[common_columns],
  madrid_2009[common_columns]
)

# remove date
madrid_data <- madrid_data[,-1]
head(madrid_data)

# transform into long format
madrid_long <- madrid_data %>%
  pivot_longer(
    cols = -c(Year, station),
    names_to = "Attribute",
    values_to = "Value"
  )

# type conversion
madrid_long$Year <- as.factor(madrid_long$Year)
madrid_long$Attribute <- as.factor(madrid_long$Attribute)
madrid_long$station <- as.factor(madrid_long$station)
madrid_long$Value <- as.numeric(madrid_long$Value)

# three-way data format
three_way <- reshape2::acast(
  data      = madrid_long,
  formula   = Attribute ~ Year ~ station,
  value.var = "Value",
  fun.aggregate = mean,
  na.rm = TRUE
)

result$BIC

# deallocate previous data
rm(madrid_long, madrid_data, madrid_2001, madrid_2002, madrid_2003, madrid_2004, madrid_2005, madrid_2006, madrid_2007, madrid_2008, madrid_2009)



# find percentages of missing data per station
# a station might not have been in use/installed for some years, or unequipped with some sensors,
# hence the high presence of null data
stations_with_missing <- apply(
  X      = three_way,
  MARGIN = 3,            # 3 => station dimension
  FUN    = function(x) mean(is.na(x))
)

# after inspection, a few station had over 90% missing data while the remaining had ~55% or less,
# out of 33, only 3 station have full coverage both timewise and sensor-wise
stations_to_keep <- names(stations_with_missing[stations_with_missing < 0.9])
length(stations_to_keep)

# keep only stations with less than 90% missing values
three_way <- three_way[,, stations_to_keep]

# impute missing values with median for each year
# imputing is necessary as the algorithm is not designed to handle missing data
three_way_imputed <- three_way

for(i in 1:dim(three_way_imputed)[3]){
  for(j in 1:dim(three_way_imputed)[2]){
    for(k in 1:dim(three_way_imputed)[1]){
      if(is.na(three_way_imputed[k,j,i])){
        three_way_imputed[k,j,i] <- median(three_way_imputed[k,j,], na.rm = TRUE)
      }
    }
  }
}

#scale the data
three_way_scaled = sparsemixmat::scale_matrix_data(three_way_imputed)

# restore dimnames
dimnames(three_way_scaled) <- dimnames(three_way_imputed)

# fit the model
# reduce the searching space if you're just testing, i suggest K <= 6
result <- sparsemixmat::fit_sparsemixmat(
  data = three_way_scaled,
  K = 2:13,                                   # number of clusters, if provided with an interval it tries each and returns the best, where the best is identified by the mixture providing the lowest BIC
  penalty_omega = c(0, 0.01, 0.05, 0.1),                       # regularization for omega ! does not support custom matrix nor 0!
  penalty_gamma = c(0, 0.01, 0.05, 0.1),                       # regularization for gamma ! does not support custom matrix nor 0!
  penalty_mu = 0.05,                          # regularization for mean
  type_penalty_mu = "lasso",            # type of penalty ! group-lasso sometimes does not work !
  penalize_diag = c(TRUE, TRUE),              # penalize diagonals ! not penalizing is not working !
  control = EM_controls(
    tol = c(1e-7, sqrt(.Machine$double.eps)), # convergence tolerances
    max_iter = rep(1000, 2),                  # max iterations for EM
    type_start = "hc",                        # initialization type ("hc" or "random"), ! random does not work !
    n_subset_start = NULL,                    # subset size for initialization
    n_random_start = 50,                      # random starts for initialization
    step_width_PGD = 1e-4                     # step width for proximal gradient descent
  ),
  verbose = TRUE
)

print(result)

#save result rdata
time<- format(Sys.time(), "%H:%M")
date<- paste0(Sys.Date(),"_", time)
filename <- paste0("~/uni/statistics/seminar_project/result_", date, ".rdata")
save(result, file = filename)

# plot bics with penalty_omega = 0.01 and penalty_gamma = 0.01, just modify p2 and p3 to plot the others
p2 <- 0.1
p3 <- 0.1
selected <- which(result$BIC$penalty_omega == p2 & result$BIC$penalty_gamma == p3)
selected
result$BIC[selected,]
plot(2:15, result$BIC[selected,]$bic, type = "o", xlab = "n. clusters", ylab = "BIC", main = "BIC", sub = paste0("Penalty omega: ", p2, " Penalty gamma: ", p3, " Penalty mu: ", 0.05))
abline(v = which.max(result$BIC[selected,]$bic)+1, col = "red")

# take row names from threeway
chemicals <- rownames(three_way_scaled)
years <- colnames(three_way_scaled)

# plot covariance/precision matrices, modify n to plot other clusters
n <- 7
omega <- result$parameters$omega[,,n]
gamma <- result$parameters$gamma[,,n]
mu <- result$parameters$mu[,,n]

# adjust plot margins
par(mar = c(5, 4, 8, 2))  # SWNE

# Plot with adjusted margins
omega_log <- sign(omega) * log(1 + abs(omega))
heatmap.2(
  omega_log,
  Rowv = NA,
  Colv = NA,
  dendrogram = "none",
  labRow = chemicals,
  labCol = chemicals,
  trace = "none",
  key = TRUE,
  key.title = "Log Intensity",
  main = paste0("Row Precision Matrix #", n, " tau: ", round(result$parameters$tau[n], 2)),
  col = hcl.colors(100, palette = 'Blue-Red 3'),
  key.par = list(mar = c(5, 1, 5, 3))  # Adjusts margins for the key itself
)

gamma_log <- sign(gamma) * log(1 + abs(gamma))
heatmap.2(
  gamma_log,
  Rowv = NA,
  Colv = NA,
  dendrogram = "none",
  labRow = years,
  labCol = years,
  trace = "none",
  key = TRUE,
  key.title = "Log Intensity",
  main = paste0("Column Precision Matrix #", n, " tau: ", round(result$parameters$tau[n], 2)),
  col = hcl.colors(100, palette='Blue-Red 3'),
  key.par = list(mar = c(5, 1, 5, 3))  # Adjusts margins for the key itself
  
)
sum(!result2$parameters$mu == 0)
sum(!result2$parameters$omega == 0)
sum(!result2$parameters$gamma == 0)

mu <- result2$parameters$mu[,,1]
mu_log <- sign(mu) * log(1 + abs(mu))
heatmap.2(
  mu_log,
  Rowv = NA,
  Colv = NA,
  dendrogram = "none",
  labRow = chemicals,
  labCol = years,
  trace = "none",
  key = TRUE,
  key.title = "Log Intensity",
  main = paste0("Mean Matrix #", n, " tau: ", round(result$parameters$tau[n], 2)),
  col = hcl.colors(100, palette='Blue-Red 3'),
  key.par = list(mar = c(5, 1, 5, 3))  # Adjusts margins for the key itself
  
)

# barplot of cluster distribution
counts<- table(result$classification)/length(result$classification)
barplot(counts, main="Cluster distribution", xlab="Cluster", ylab="Count", col=hcl.colors(14, palette='viridis'), ylim=c(0, 1))