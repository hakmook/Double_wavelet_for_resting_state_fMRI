
# This code is to integrate the results from two data analysis approaches,
# i.e., AVG-FC and double-wavelet approach, performed using the Matlab code,
# 'resting_analysis_submitted.m'. Once you apply the Matlab code and save the
# result csv files into one folder (e.g., '~/Your_result_folder/RESULTS'), 
# then you apply this R code to integrate them and compute p-values 
# for the 14 ROI pairs.


.libPaths("~/Dropbox/RLibrary/")
library("data.table")

setwd('~/Your_result_folder/RESULTS') # The folder contains 12 csv files

file_list <- list.files()
allresult <- data.frame()

for (file in file_list){
  dataset <- fread(file)
  allresult<- rbindlist(list(dataset,allresult)) 
}

allresult <- as.data.frame(allresult)
colnames(allresult) <- c("ID", "method", "pair", "mean")

result_pilot_normal <- allresult[ allresult$ID %in%  c(1000, 1001, 1002, 1003, 1004, 1005),]
result_pilot_depress <- allresult[ allresult$ID %in%  c(2000, 2001, 2002, 2003, 2004, 2005),]

result_test <- data.frame()

ztranform <- function(x) {
  0.5 * log((1+x)/(1-x))
}

iztranform <- function(x){
  (exp(2*x) -1 )/  (exp(2*x) +1 )
}

for (i in unique(allresult$pair)){
  temp0 <- result_pilot_normal[ result_pilot_normal$pair == i ,]
  temp1 <- result_pilot_depress[ result_pilot_depress$pair == i ,]
  
  temp0$mean <- ztranform(temp0$mean )
  temp1$mean <- ztranform(temp1$mean )
  
  test_mean <- t.test(temp0$mean[temp0$method == "mean"], temp1$mean[temp1$method == "mean"])
  test_dw_waveform <- t.test(temp0$mean[temp0$method == "dw_waveform"], temp1$mean[temp1$method == "dw_waveform"])
  test_dw_variance <- t.test(temp0$mean[temp0$method == "dw_variance"], temp1$mean[temp1$method == "dw_variance"])
  
  result_temp <- data.frame( pair = i, method = c("mean", "dw_waveform","dw_variance"), 
                             normal =tanh(c(test_mean$estimate[1],
                                      test_dw_waveform$estimate[1] , 
                                      test_dw_variance$estimate[1] )
                             ),
                             depress = tanh(c(test_mean$estimate[2],
                                      test_dw_waveform$estimate[2], 
                                     test_dw_variance$estimate[2])
                             ),
                             
                             diff = tanh(c(test_mean$estimate[1] - test_mean$estimate[2],
                                      test_dw_waveform$estimate[1] - test_dw_waveform$estimate[2], 
                                      test_dw_variance$estimate[1] - test_dw_variance$estimate[2]
                            )),
                            pval = c( test_mean$p.value, test_dw_waveform$p.value, test_dw_variance$p.value   ))
  result_test <- rbind(result_temp, result_test)
}

# with Z transform
result_test <- subset(result_test, method != "dw_variance")

pair_choice <- c(1,8,14,15,16,18,24,25,28,29,30,32,33,35)
subresult <- result_test[ result_test$pair %in%  pair_choice, ]
subresult <- subset(subresult, method != "dw_variance")

dw_data <- subset(subresult, method == "dw_waveform") # dw_data is for our approach
mean_data <- subset(subresult, method == "mean") # mean_data is for AVG-FC
