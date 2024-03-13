library(MASS)
library(psych)
library(dplyr)
library(glmnet)
library(mice)
library(stats)

set.seed(20231105)

# fit a LASSO model to the non missing data
LASSO_model_fitting <- function(sim_data){
  data <- sim_data
  # adaptive lasso 
  #step 1: get the weights
  model_lasso <- glmnet(x = as.matrix(data[, c(2:11)]), 
                        y = as.matrix(data[, 1]),  
                        alpha = 1, 
                        family = "binomial", 
                        lambda = 0, 
                        standardized = FALSE, 
                        intercept = TRUE)
  weights <- as.numeric(coef(model_lasso))[-1]
  
  cv_model <- cv.glmnet(x = as.matrix(data[, c(2:11)]), 
                        y = as.matrix(data[, 1]), 
                        family = "binomial", 
                        alpha = 1,
                        nfold = 10,
                        penalty.factor= 1/abs(weights))
  
  best1_model <- glmnet(x = as.matrix(data[, c(2:11)]), 
                        y = as.matrix(data[, 1]), 
                        alpha=1, 
                        lambda=cv_model$lambda.min,
                        standardize=FALSE, 
                        intercept=TRUE, 
                        penalty.factor=1/abs(weights),
                        family = "binomial")
  
  all_model_data = matrix(nrow = 11, ncol = 1)
  rownames(all_model_data) <- c("intercept", "x1","x2", "x3", "d1", "d2", "d3", "d4", "d5", "d6", "d7")
  all_model_data[1:11, 1] = as.matrix(coef(best1_model))[, 1]
  
  selection_result = matrix(nrow = 10, ncol = 1)
  rownames(selection_result) <- c("x1","x2", "x3", "d1", "d2", "d3", "d4", "d5", "d6", "d7")
  
  for(var in 1:10){
    if (all_model_data[var+1, 1] == 0){
      selection_result[var, 1] = 0
    }
    else{
      selection_result[var, 1] = 1
    }
  }
  
  result = list()
  result$selection_result = selection_result
  result$coef_est = coef(best1_model)
  result$lambda = cv_model$lambda.min
  return(result)
}


get_metrics <- function(sim_selection, sim_est, t_b1, t_b2, t_b3, rep_num){
  selection = as.data.frame(sim_selection)
  colnames(selection) <- c("x1","x2", "x3", "d1", "d2", "d3", "d4", "d5", "d6", "d7")
  # success rate
  success_count = 0
  correct_list = c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
  for (rep in 1:rep_num){
    if (sum(selection[rep,] == correct_list) == 10){
      success_count = success_count + 1
    }
  }
  success_rate = success_count/rep_num
  
  # false positive rate
  fp_count = 0 
  for (rep in 1:rep_num){
    if (sum(selection[rep,1:3]) == 3 & sum(selection[rep,4:10]) > 0){
      fp_count = fp_count + 1 
    }
  }
  fpr = fp_count/rep_num
  
  # false negative rate
  fn_count = 0 
  for (rep in 1:rep_num){
    if (sum(selection[rep,1:3]) < 3){
      fn_count = fn_count + 1 
    }
  }
  fnr = fn_count/rep_num
  
  # selection rate
  selection_rate = cbind(mean(selection$x1), mean(selection$x2), mean(selection$x3))
  
  # false discover rate
  fdr = cbind(mean(selection$d1), mean(selection$d2), mean(selection$d3), mean(selection$d4),
              mean(selection$d5), mean(selection$d6),mean(selection$d7))
  
  # relative bias
  est_data = as.data.frame(sim_est)
  colnames(est_data ) <- c("intercept","x1","x2", "x3", "d1", "d2", "d3", "d4", "d5", "d6", "d7")
  est = replace(est_data, is.na(est_data), 0)
  est$rb_int = abs((est$intercept - (-1))/(-1))
  est$rb_x1 = abs((est$x1 - (t_b1))/(t_b1))
  est$rb_x2 = abs((est$x2 - (t_b2))/(t_b2))
  est$rb_x3 = abs((est$x3 - (t_b3))/(t_b3))
  
  b0_mean = mean(est$intercept)
  b1_mean = mean(est$x1)
  b2_mean = mean(est$x2)
  b3_mean = mean(est$x3)
  
  b0_sd = sd(est$intercept)
  b1_sd = sd(est$x1)
  b2_sd = sd(est$x2)
  b3_sd = sd(est$x3)
  
  rb_int = mean(est$rb_int)
  rb_x1 = mean(est$rb_x1)
  rb_x2 = mean(est$rb_x2)
  rb_x3 = mean(est$rb_x3)
  
  output = cbind(success_rate, fpr, fnr, selection_rate, fdr, b0_mean, b1_mean, b2_mean, b3_mean, 
                 b0_sd, b1_sd, b2_sd, b3_sd, rb_int, rb_x1, rb_x2, rb_x3)
  
  return(output)
  
}

condition = 101 
repli = 500
N = 500
b0 = 1.734 
b1 = 0 
beta1 = 0.3 
beta2 = 0.4 
beta3 = 0.5 
cor5_7 = 0.3 
d1_mis_pre = 0.05
re = 1

# data generation function for MAR with d1
simulate_lasso_data <- function(condition, repli, N, b0, b1, beta1, beta2, beta3, cor5_7, d1_mis_pre){
  start = Sys.time()
  num_of_pred = 10
  
  # set up selection record
  selection_complete = matrix(NA, repli,num_of_pred)
  selection_lw = matrix(NA, repli,num_of_pred)
  selection_meinp = matrix(NA, repli,num_of_pred)
  selection_mi_ave = matrix(NA, repli,num_of_pred)
  selection_mi_thr = matrix(NA, repli,num_of_pred)
  selection_mi_stacked = matrix(NA, repli,num_of_pred)
  selection_mi_al = matrix(NA, repli,num_of_pred)
  
  # set up parameter estimates record
  est_complete = matrix(NA,repli, num_of_pred+1)
  est_lw = matrix(NA, repli, num_of_pred+1)
  est_meinp = matrix(NA, repli, num_of_pred+1)
  est_mi_ave = matrix(NA, repli, num_of_pred+1)
  est_mi_thr = matrix(NA, repli, num_of_pred+1)
  est_mi_stacked = matrix(NA, repli, num_of_pred+1)
  est_mi_al = matrix(NA, repli, num_of_pred+1)
  
  for (re in 1: repli){
    l1 <- sqrt(0.012/0.04)
    l2 <- sqrt(0.012/0.04)
    l3 <- sqrt(0.012/0.09)
    
    common_factor <- rnorm(N, 0, 1)
    x1 <- l1 * common_factor + rnorm(N, 0, sqrt(1-l1^2))
    x2 <- l2 * common_factor + rnorm(N, 0, sqrt(1-l2^2))
    x3 <- l3 * common_factor + rnorm(N, 0, sqrt(1-l3^2))
    y <- -1 + beta1 * x1 + beta2 * x2 + beta3 * x3 + rlogis(N, 0, 1)
    y_binary <- if_else(y > 0, 1, 0)
    
    d1 <- sqrt(0.1) * common_factor + rnorm(N, 0, sqrt(1- 0.1))
    d2 <- sqrt(0.1) * common_factor + rnorm(N, 0, sqrt(1- 0.1))
    d3 <- sqrt(0.1) * common_factor + rnorm(N, 0, sqrt(1- 0.1))
    d4 <- sqrt(0.1) * common_factor + rnorm(N, 0, sqrt(1- 0.1))
    
    d5 <- cor5_7 * x3 + rnorm(N, 0, sqrt(1- cor5_7^2))
    d6 <- cor5_7 * x2 + rnorm(N, 0, sqrt(1- cor5_7^2))
    d7 <- cor5_7 * x1 + rnorm(N, 0, sqrt(1- cor5_7^2))
    
    sim_data <- data.frame(y_binary,x1, x2, x3, d1, d2, d3, d4, d5, d6, d7)
    names(sim_data) = c('y','x1','x2','x3','d1','d2','d3','d4','d5','d6','d7')
    # write Out Data 
    write.table(sim_data,file=paste('/Users/yibinni/Desktop/LASSO All the way/MissingSimulation/complete/',condition,'_d',re,'.dat',sep=''),row.names = TRUE,col.names = FALSE,sep = "\t")
    
    # complete data LASSO
    complete_model =  LASSO_model_fitting(sim_data)
    est_complete[re, ] = as.matrix(complete_model$coef_est)[,1]
    selection_complete[re, ] = complete_model$selection_result
    
    # missing data generation
    miss_x1 <- c()
    u_x1 <- c()
    mult_pred_miss <- c()
    
    # Create Missing Data Patterns for x1 Only
    u_x1 <- runif(N,0,1)
    p_x1 <- exp(b0 + b1*d1)/(1 + exp(b0 + b1*d1))
    
    for(i in 1:N){
      if(p_x1[i] > u_x1[i]){miss_x1[i] <- 1}
      else(miss_x1[i] <- 0)
    }
    
    # create missing data for d1 
    u_d1 = c()
    miss_d1 <- c()
    u_d1 = runif(N, 0, 1)
    
    for(i in 1:N){
      if(u_d1[i] < d1_mis_pre){miss_d1[i] <- 0}
      else(miss_d1[i] <- 1)
    }
    
    # Create Missing Data Patterns for Other Predictors: x2, x3 & d2 - d7 
    miss = c()
    for(k in 1:8){
      u <- runif(N,0,1)
      
      p <- runif(N,.95,1)
      
      for(i in 1:N){
        if(p[i] > u[i]){miss[i] <- 1}
        else(miss[i] <- 0)
      }
      
      mult_pred_miss[[k]] = miss
    }
    
    x1_miss = x1
    x2_miss = x2
    x3_miss = x3
    d1_miss = d1
    d2_miss = d2
    d3_miss = d3
    d4_miss = d4
    d5_miss = d5
    d6_miss = d6
    d7_miss = d7
    
    # Impose missing values on target predictor: x1
    for(i in 1:N){if(miss_x1[i] == 0){x1_miss[i] <- NA}}
    
    # Impose missing values on target predictor: d1
    for(i in 1:N){if(miss_d1[i] == 0){d1_miss[i] <- NA}}
    
    # Impose missing values on all other predictors x2, x3 & d2-d7
    for(k in 1:8){
      single_pred_miss = mult_pred_miss[[k]]
      
      if(k == 1){
        for(i in 1:N){if(single_pred_miss[i] == 0){x2_miss[i] <- NA}}}
      
      if(k == 2){
        for(i in 1:N){if(single_pred_miss[i] == 0){x3_miss[i] <- NA}}}
      
      if(k == 3){
        for(i in 1:N){if(single_pred_miss[i] == 0){d2_miss[i] <- NA}}}
      
      if(k == 4){
        for(i in 1:N){if(single_pred_miss[i] == 0){d3_miss[i] <- NA}}}
      
      if(k == 5){
        for(i in 1:N){if(single_pred_miss[i] == 0){d4_miss[i] <- NA}}}
      
      if(k == 6){
        for(i in 1:N){if(single_pred_miss[i] == 0){d5_miss[i] <- NA}}}
    
      if(k == 7){
        for(i in 1:N){if(single_pred_miss[i] == 0){d6_miss[i] <- NA}}}
      
      if(k == 8){
        for(i in 1:N){if(single_pred_miss[i] == 0){d7_miss[i] <- NA}}}
    }
    
    data <- data.frame(miss_x1, d1)
    model = glm(miss_x1~d1, data = data)
    print(with(summary(model), 1 - deviance/null.deviance))
    miss_data = as.data.frame(cbind(y_binary,x1_miss,x2_miss,x3_miss,d1_miss,d2_miss,d3_miss,d4_miss,d5_miss,d6_miss,d7_miss))
    names(miss_data) = c('y','x1','x2','x3','d1','d2','d3','d4','d5','d6','d7')
    write.table(miss_data,file=paste('/Users/yibinni/Desktop/LASSO All the way/MissingSimulation/missing/',condition,'_d', re,'.dat',sep=''),row.names = TRUE,col.names = FALSE,sep = "\t")
    
    describe(miss_data)
    ## listwise deletion
    lw_miss_data = miss_data[complete.cases(miss_data),]
    lw_model_result = LASSO_model_fitting(lw_miss_data)
    lw_model_selection = lw_model_result$selection_result
    lw_model_coef = lw_model_result$coef_est
    selection_lw[re,] = lw_model_selection
    est_lw[re, ] = as.matrix(lw_model_result$coef_est)[,1]
    
    ## Mean Imputation
    mean_miss_data = miss_data
    for(i in 1:ncol(miss_data)){
      mean_miss_data[is.na(miss_data[,i]), i] <- mean(miss_data[,i], na.rm = TRUE)
    }
    mean_model_result = LASSO_model_fitting(mean_miss_data)
    mean_model_coefs = mean_model_result$coef_est
    mean_model_selection = mean_model_result$selection_result
    selection_meinp[re,] = mean_model_selection
    est_meinp[re,] = as.matrix(mean_model_coefs)[,1]
    
    ## Multiple Imputation
    m = 20
    pred_vars = cbind(miss_data[,2:ncol(miss_data)])
    dep_var = as.data.frame(miss_data[,1])
    names(dep_var) = c('y')
    
    ## Impute
    imp_data = mice(pred_vars, m=20)
    
    # fitting all imputed datasets 
    all_model_selection = matrix(NA, 10, m)
    all_model_coef = matrix(NA, 11, m)
    all_lambda = matrix(NA, 1, m)
    all_imp_data = rbind(cbind(dep_var, complete(imp_data, action=1)),
                         cbind(dep_var, complete(imp_data, action=2)),
                         cbind(dep_var, complete(imp_data, action=3)),
                         cbind(dep_var, complete(imp_data, action=4)),
                         cbind(dep_var, complete(imp_data, action=5)),
                         cbind(dep_var, complete(imp_data, action=6)),
                         cbind(dep_var, complete(imp_data, action=7)),
                         cbind(dep_var, complete(imp_data, action=8)),
                         cbind(dep_var, complete(imp_data, action=9)),
                         cbind(dep_var, complete(imp_data, action=10)),
                         cbind(dep_var, complete(imp_data, action=11)),
                         cbind(dep_var, complete(imp_data, action=12)),
                         cbind(dep_var, complete(imp_data, action=13)),
                         cbind(dep_var, complete(imp_data, action=14)),
                         cbind(dep_var, complete(imp_data, action=15)),
                         cbind(dep_var, complete(imp_data, action=16)),
                         cbind(dep_var, complete(imp_data, action=17)),
                         cbind(dep_var, complete(imp_data, action=18)),
                         cbind(dep_var, complete(imp_data, action=19)),
                         cbind(dep_var, complete(imp_data, action=20))
    )
    
    for(num_imp in 1:m){
      complete_data = cbind(dep_var, complete(imp_data, action=num_imp))
      model_fitting = LASSO_model_fitting(complete_data)
      all_model_coef[1:11, num_imp] = as.matrix(model_fitting$coef_est)[, 1]
      all_model_selection[1:10, num_imp] = model_fitting$selection_result
      all_lambda[1, num_imp] = model_fitting$lambda
    }
    
    # getting averaged lambda and stacked data
    ave_lambda = mean(all_lambda)
    
    # average pooling
    ap_para_result = rowMeans(all_model_coef)
    est_mi_ave[re, 1:11] = matrix(ap_para_result, nrow = 1, ncol = 11)
    ap_para_result_selection = unlist(lapply(ap_para_result[2:11], function(x) ifelse(abs(x)>0, 1, 0)))
    selection_mi_ave[re,1:10] = ap_para_result_selection
    
    # pooling with a threshold
    pool_thre = .8
    pooling_para_result_if_selected = rowMeans(all_model_selection)
    pooling_para_result_selection = ifelse(pooling_para_result_if_selected >= pool_thre, 1, 0)
    selection_mi_thr[re, ] = pooling_para_result_selection
    # print(pooling_para_result_selection)
    # fit a relaxed lasso with selected variable 
    if (sum(pooling_para_result_selection) == 0){
      est_mi_thr[re,] = NA
    }else{
      relaxed_pre_list = c()
      name_list = c('x1','x2','x3','d1','d2','d3','d4','d5','d6','d7')
      for (i in 1:10){
        if (pooling_para_result_selection[i] == 1){
          relaxed_pre_list = append(relaxed_pre_list, name_list[i])
        }
      }
      equ = paste(relaxed_pre_list, collapse = " + ")
      formu = paste("y ~ ", equ)
      print(formu)
      relax_lasso_model <- glm(formu, data = all_imp_data, family = "binomial")
      result = as.matrix(coef(relax_lasso_model))
      
      all_list = c("(Intercept)",'x1','x2','x3','d1','d2','d3','d4','d5','d6','d7')
      rn = rownames(result)
      for(a in 1:length(rn)){
        ind = which(all_list == rn[a])[1]
        est_mi_thr[re,ind] = result[rn[a],1]
      }
      
    }
    
    
    # stacked mi 
    # model fitting with stacked data
    stacked_model_lasso <- glmnet(x = as.matrix(all_imp_data[, c(2:11)]), 
                                  y = all_imp_data$y,  
                                  alpha = 1, 
                                  family = "binomial", 
                                  lambda = 0, 
                                  standardized = FALSE, 
                                  intercept = TRUE)
    stacked_weights <- as.numeric(coef(stacked_model_lasso))[-1]
    
    stacked_best1_model <- glmnet(x = as.matrix(all_imp_data[, c(2:11)]), 
                                  y = all_imp_data$y, 
                                  alpha=1, 
                                  lambda=ave_lambda,
                                  standardize=FALSE, 
                                  intercept=TRUE, 
                                  penalty.factor=1/abs(stacked_weights),
                                  family = "binomial")
    mi_stacked_coef = coef(stacked_best1_model)
    est_mi_stacked[re, 1:11] = as.matrix(mi_stacked_coef)[,1]
    mi_stacked_result_selection = unlist(ifelse(abs(mi_stacked_coef) > 0, 1, 0)[2:11])
    selection_mi_stacked[re,] = mi_stacked_result_selection
    
    # averaged lambda with .8 threshold
    mi_al_model_selection = matrix(NA, 10, m)
    mi_al_model_coef = matrix(NA, 11, m)
    # refit all models with average lambda
    for(num_imp in 1:m){
      complete_data = cbind(dep_var, complete(imp_data, action=num_imp))
      al_model_lasso <- glmnet(x = as.matrix(complete_data[, c(2:11)]), 
                               y = complete_data$y,  
                               alpha = 1, 
                               family = "binomial", 
                               lambda = 0, 
                               standardized = FALSE, 
                               intercept = TRUE)
      al_weights <- as.numeric(coef(al_model_lasso))[-1]
      
      al_best1_model <- glmnet(x = as.matrix(complete_data[, c(2:11)]), 
                               y =complete_data$y, 
                               alpha=1, 
                               lambda=ave_lambda,
                               standardize=FALSE, 
                               intercept=TRUE, 
                               penalty.factor=1/abs(al_weights),
                               family = "binomial")
      # print(al_best1_model$lambda)
      coef_result = as.matrix(coef(al_best1_model))[, 1]
      selection_result = ifelse(abs(coef_result) > 0, 1, 0)
      mi_al_model_coef[1:11, num_imp] = coef_result
      mi_al_model_selection[1:10, num_imp] = selection_result[-1]
    }
    # pooling the model fitting results
    mi_al_pool_thre = .8
    mi_al_selection_prec = rowMeans(mi_al_model_selection)
    mi_al_selection_result = ifelse(mi_al_selection_prec >= mi_al_pool_thre, 1, 0)
    print(mi_al_selection_result)
    selection_mi_al[re,] = mi_al_selection_result
    
    
    # fit a relaxed lasso with selected parameters
    if(sum(mi_al_selection_result) == 0){
      est_mi_al[re,] = NA
    }else{
      relaxed_pre_list_ave = c()
      name_list = c('x1','x2','x3','d1','d2','d3','d4','d5','d6','d7')
      all_list = c("(Intercept)",'x1','x2','x3','d1','d2','d3','d4','d5','d6','d7')
      
      for (i in 1:10){
        if (mi_al_selection_result[i] == 1){
          relaxed_pre_list_ave = append(relaxed_pre_list_ave, name_list[i])
        }
      }
      equ_ave = paste(relaxed_pre_list_ave, collapse = " + ")
      formu_ave = paste("y ~ ", equ_ave)
      print(formu_ave)
      relax_lasso_model <- glm(formu_ave, data = all_imp_data, family = "binomial")
      result = as.matrix(coef(relax_lasso_model))
      rn = rownames(result)
      
      for(a in 1:length(rn)){
        ind = which(all_list == rn[a])[1]
        est_mi_al[re, ind] = result[rn[a],1]
      }
    }
    
    
    print(re)
  }
  end = Sys.time()
  
  output_pre_replication = as.data.frame(cbind(selection_complete, est_complete, 
                                               selection_lw, est_lw,
                                               selection_meinp, est_meinp, 
                                               selection_mi_ave, est_mi_ave,
                                               selection_mi_thr, est_mi_thr, 
                                               selection_mi_stacked,est_mi_stacked,
                                               selection_mi_al, est_mi_al))
  
  write.table(output_pre_replication,file=paste('/Users/yibinni/Desktop/LASSO All the way/MissingSimulation/output2/',condition,'_d', re, 'result_output_pre_rep', '.csv',sep=''),row.names = FALSE,col.names = FALSE,sep = "\t")
  
  complete_data_output = get_metrics(selection_complete, est_complete, beta1, beta2, beta3, repli)
  lw_data_output = get_metrics(selection_lw, est_lw, beta1, beta2, beta3, repli)
  meinp_data_output = get_metrics(selection_meinp, est_meinp, beta1, beta2, beta3, repli)
  mi_ave_data_output = get_metrics(selection_mi_ave, est_mi_ave, beta1, beta2, beta3, repli)
  mi_thr_data_output = get_metrics(selection_mi_thr, est_mi_thr, beta1, beta2, beta3, repli)
  mi_stacked_data_output = get_metrics(selection_mi_stacked, est_mi_stacked, beta1, beta2, beta3, repli)
  mi_al_data_output = get_metrics(selection_mi_al, est_mi_al, beta1, beta2, beta3, repli)
  
  all_output = as.data.frame(rbind(complete_data_output, lw_data_output, meinp_data_output, mi_ave_data_output,
                                   mi_thr_data_output, mi_stacked_data_output, mi_al_data_output))
  colnames(all_output) = c("success_rate", "fpr", "fnr", "x1_sr", "x2_sr", "x3_sr", "d1_sr", "d2_sr", "d3_sr",
                           "d4_sr", "d5_sr", "d6_sr", "d7_sr", "b0_mean", "b1_mean", "b2_mean", "b3_mean", "b0_sd",
                           "b1_sd", "b2_sd", "b3_sd", "rb_b0", "rb_b1", "rb_b2", "rb_b3")
  rownames(all_output) <- c("complete_data", "listwise", "mean_imputation", "MI_averaging", "MI_threshold",
                            "MI_stacked", "MI_averaged_lambda")
  write.table(all_output,file=paste('/Users/yibinni/Desktop/LASSO All the way/MissingSimulation/output2/',condition,'_d', re, 'result_output', '.csv',sep=''),row.names = TRUE,col.names = TRUE,sep = "\t")
  
  return(list(start, end))
}

# read in conditions
conditions = read.csv(file = "/Users/yibinni/Desktop/LASSO All the way/MissingSimulation/conditions.csv", sep = "\t", header = TRUE)
con = 186
for (con in 101:300){
  simulate_lasso_data(conditions$cond[con], conditions$replication_num[con], conditions$sample_sizes[con], conditions$b0[con], 
                      conditions$b1[con], conditions$beta1[con], conditions$beta2[con], conditions$beta3[con], conditions$cor5_7[con], conditions$d1_mis_pre[con])
}
