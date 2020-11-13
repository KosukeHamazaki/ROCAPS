#######################################################################################################################################
###### 2.1.3_Rice_ROCAPS_Rice_Optimization_for_Cost_and_Accuracy_over_Phenotyping_in_early_growth_stages_Real_for_terminal ######
#######################################################################################################################################

###### 1. Settings ######
##### 1.0. Reset workspace ######
#rm(list=ls())



##### 1.1. Setting working directory to the "ROCAPS" directory #####
crop_name <- "rice"
project <- "ROCAPS"
# setwd(paste0("/Users/hamazaki/research/", crop_name, "/Project/", project))   ### for mac OS
# setwd(paste0("/media/hamazaki/d4000953-5e56-40ce-97ea-4cfee57fc91a/research/",
#              crop_name, "/Project/", project))    ### for Ubuntu 1
setwd(paste0("/media/hamazaki/D6B89D51B89D314B/research/",
             crop_name, "/Project/", project))    ### for Ubuntu 2
scriptID <- 2.1



##### 1.2. Import packages #####
require(data.table)
require(MASS)
require(rrBLUP)
require(BGLR)
require(foreach)
require(doParallel)
require(ggplot2)


# dir.script <- paste0("~/research_secret/", crop_name, "/Project/", project, "/")
# dir.script <- paste0("~/research_secret/", crop_name, "/Project/", project, "/")
dir.script <- paste0("~/GitHub/research_secret/", crop_name, "/Project/", project, "/")
source(paste0(dir.script, "1.0_Rice_ROCAPS_data_scaled_imputation.R"))
source(paste0(dir.script, "1.1.0_Rice_ROCAPS_MTM2.R"))
source(paste0(dir.script, "1.1_Rice_ROCAPS_MTM_prediction.R"))
# source(paste0(dir.script, "1.2_Rice_ROCAPS_rrBLUP_prediction.R"))
# source(paste0(dir.script, "1.3_Rice_ROCAPS_RF_prediction.R"))
source(paste0(dir.script, "1.4_Rice_ROCAPS_Gaussian_process_for_BO.R"))
source(paste0(dir.script, "1.5_Rice_ROCAPS_Acquisition_func.R"))
source(paste0(dir.script, "1.6_Rice_ROCAPS_convex_check.R"))
# source(paste0(dir.script, "1.7_Rice_ROCAPS_pforeach.R"))



##### 1.3. Setting some parameters #####
trial.no <- 1

eval.start <- 21
eval.end <- 49
eval.space <- 14
eval.day <- seq(eval.start, eval.end, by = eval.space)
n.group <- 3

seed.range <- 1:10000

dir.mid <- paste0("midstream/", project, "_st=", eval.start, "_end=", eval.end,
                  "_by_", eval.space, "days_with_", n.group, "groups_trial=",
                  trial.no, "/")
file.params <- paste0(dir.mid, "0.2_parameter_informations.RData")
load(file = file.params)

scriptID <- 2.1
trial.no.pred <- 1
target.traits <- "Real"
method <- "MTM"
convergence.type <- "prop"
fitting.scale <- "R2"
Aq.weighted <- FALSE
relax.param <- 0.01
NDCG.prop <- 0.03
n.fold <- 10
n.core <- 4

count2.thres <- 0.5
count.prop <- 0.1
count2.max <- 5

if(method == "MTM"){
  n.iter.pred <- 3
}else{
  if(n.fold != 371){
    n.iter.pred <- 10
  }else{
    n.iter.pred <- 1
  }
}
seeds.pred <- sample(seed.range, n.iter.pred)

if(length(target.traits) < 2){
  dir.mid.ROCAPS.0 <- paste0(dir.mid, scriptID, "_", target.traits, "_", method, "_", convergence.type, "_",
                             fitting.scale, "_Aq.weighted=", Aq.weighted, "_",
                             n.iter.pred, "_", n.fold, "_", trial.no.pred, "/")

  dir.res.0 <- paste0("results/", scriptID, "_", target.traits, "_", method, "_", convergence.type, "_",
                      fitting.scale, "_Aq.weighted=", Aq.weighted, "_",
                      n.iter.pred, "_", n.fold, "_", trial.no.pred, "/")

} else {
  dir.mid.ROCAPS.0 <- paste0(dir.mid, scriptID, "_", target.traits[1], "_herit=",
                             target.traits[2], "_scenario=", target.traits[3], "_",
                             method, "_", convergence.type, "_",
                             fitting.scale, "_Aq.weighted=", Aq.weighted, "_",
                             n.iter.pred, "_", n.fold, "_", trial.no.pred, "/")

  dir.res.0 <- paste0("results/", scriptID, "_", target.traits[1], "_herit=",
                      target.traits[2], "_scenario=", target.traits[3], "_",
                      method, "_", convergence.type, "_",
                      fitting.scale, "_Aq.weighted=", Aq.weighted, "_",
                      n.iter.pred, "_", n.fold, "_", trial.no.pred, "/")
}

dir.create(dir.mid.ROCAPS.0)
dir.create(dir.res.0)

file.params.ROCAPS <- paste0(dir.mid.ROCAPS.0, scriptID, "_", project, "_all_parameters.RData")
save.image(file.params.ROCAPS)



###### 2. Modification of data ######
##### 2.1. Read phenotype, relationship matrix and other files into R #####
file.pheno <- "data/phenotype/0.2_observed_and_simulated_phenotype_for_ROCAPS.RData"
load(file = file.pheno)


trait.names <- colnames(Y.all)
line.names <- rownames(Y.all)
n.lines <- nrow(Y.all)

file.amat <- "data/genotype/0.1_additive_genetic_relationship_matrix.RData"
load(file = file.amat)

file.case.cross <- paste0(dir.mid, "0.2_all_case_with_cross.RData")
load(file = file.case.cross)

file.case.cost <- paste0(dir.mid, "0.2_all_case_with_cost.csv")
all.case.cost <- read.csv(file = file.case.cost, row.names = 1)

costs <- c(all.case.cost[, ncol(all.case.cost)])
all.case <- as.matrix(all.case.cost[, -ncol(all.case.cost)])
n.combinations <- nrow(all.case)

file.group <- "data/extra/0.2_accession_group_information.csv"
group.data <- read.csv(file = file.group)
group.info <- group.data[, 3]
names(group.info) <- line.names

if(length(target.traits) >= 2){
  file.target.sim <- "data/phenotype/0.2_all_simulated_values_of_target_traits.RData"
  load(file = file.target.sim)
}



##### 2.2. Target trait and other traits to improve genomic prediction #####
predictors.trait.all.no <- c(7, 10, 8, 11, 9, 12)
predictors.trait.all <- trait.names[predictors.trait.all.no]
predictors <- Y.all[, predictors.trait.all.no]
n.predictors <- ncol(predictors)


if(length(target.traits) < 2){
  target.trait.nos <- 1:3
  target.trait.names <- colnames(Y.all)[target.trait.nos]

  Y.target <- Y.all[, target.trait.nos]
} else {
  herit.no <- as.numeric(target.traits[2])
  scenario.no <- as.numeric(target.traits[3])

  Y.target <- target.sim.array[, , herit.no, scenario.no, 3]
  target.trait.names <- colnames(Y.target)
}



for(trait.no in 1:length(target.trait.names)){
  target.trait <- target.trait.names[trait.no]
  y.target <- Y.target[, trait.no]
  names(y.target) <- line.names

  dir.mid.ROCAPS <- paste0(dir.mid.ROCAPS.0,
                           scriptID, "_", target.trait, "_",
                           method, "_", convergence.type, "_",
                           fitting.scale, "_Aq.weighted=", Aq.weighted, "_",
                           n.iter.pred, "_", n.fold, "_", trial.no.pred, "/")
  dir.create(dir.mid.ROCAPS)


  ##### 2.3. Modify data to GS format #####
  Y.0 <- as.matrix(y.target)
  colnames(Y.0) <- target.trait

  Y.max <- cbind(predictors, y.target)
  colnames(Y.max) <- c(predictors.trait.all, target.trait)

  Z <- diag(n.lines)
  XF <- NULL
  Y.0.scaled <- data_scaled_imputation(Y = Y.0, XF = XF, Z = Z, imp = FALSE, rm.NA = FALSE)$Y
  Y.max.scaled <- data_scaled_imputation(Y = Y.max, XF = XF, Z = Z, imp = FALSE, rm.NA = FALSE)$Y





  ###### 3. Start of the prediction ######
  ##### 3.1. Perform genomic prediction for the cost-max case and the cost-nothing case #####
  #### 2.1.1. The cost-max case (All traits of Y.mat are used for prediction) ####
  NDCG.k <- round(n.lines * NDCG.prop, 0)


  iter.now <- "0.max"
  dir.names.mid.now <- paste0(dir.mid.ROCAPS, iter.now)
  dir.create(dir.names.mid.now)


  if(method == "MTM"){
    resCov <- list(type = "UN", df0 = ncol(Y.max.scaled), S0 = diag(ncol(as.matrix(Y.max.scaled))))
    scale.maxs <- MTM_prediction_iter(Y = Y.max.scaled, Y.test.no = NULL, XF = NULL, x = NULL, Z = Z, K = K.A,
                                      resCov = resCov, nIter = 2000, burnIn = 500, thin = 5,
                                      tolD = 1e-05, NDCG.k = NDCG.k, n.iter = n.iter.pred,
                                      fitting.scale = fitting.scale, saveAt = dir.names.mid.now,
                                      CV = TRUE, n.fold = n.fold, seeds = seeds.pred, silent = TRUE,
                                      time = TRUE, verbose = FALSE, Parallel = TRUE,  n.core = n.core)
  }

  # if(method == "rrBLUP"){
  #   scale.maxs <- rrBLUP_prediction_iter(Y = Y.scaled, Y.test.no = NULL, XF = NULL, x = NULL, Z = Z2, K = K.A,
  #                                        method = "REML", bounds = c(1e-09, 1e+09), SE = TRUE, return.Hinv = TRUE, NDCG.k = NDCG.k,
  #                                        n.iter = n.iter.pred, fitting.scale = fitting.scale, saveAt = dir.names.mid.now,
  #                                        CV = TRUE, n.fold = n.fold, seeds = seeds.pred, silent = TRUE, time = TRUE)
  # }
  #
  # if(method == "RF"){
  #   scale.maxs <- RF_prediction_iter(Y = Y.scaled, Y.test.no = NULL, XF = NULL, x = NULL, Z = Z2, K = K.A,
  #                                    pca.option = "prop", pca.prop = 0.5, pca.num = 10, NDCG.k = NDCG.k,
  #                                    n.iter = n.iter.pred, fitting.scale = fitting.scale, saveAt = dir.names.mid.now,
  #                                    CV = TRUE, n.fold = n.fold, seeds = seeds.pred, silent = TRUE, time = TRUE)
  # }

  if(fitting.scale == "MSE"){
    scale.maxs <- 1 / scale.maxs
  }
  scale.max <- mean(scale.maxs)

  #### 2.1.2. The cost-max case (No trait of Y.mat is used for prediction, genotype only) ####
  iter.now <- "0.0"
  dir.names.mid.now <- paste0(dir.mid.ROCAPS, iter.now)
  dir.create(dir.names.mid.now)

  if(method == "MTM"){
    resCov <- list(type = "UN", df0 = ncol(Y.0.scaled), S0 = diag(ncol(as.matrix(Y.0.scaled))))
    scale.nos <- MTM_prediction_iter(Y = Y.0.scaled, Y.test.no = NULL, XF = NULL, x = NULL, Z = NULL, K = K.A,
                                     resCov = resCov, nIter = 2000, burnIn = 500, thin = 5,
                                     tolD = 1e-05, NDCG.k = NDCG.k, n.iter = n.iter.pred,
                                     fitting.scale = fitting.scale, saveAt = dir.names.mid.now,
                                     CV = TRUE, n.fold = n.fold, seeds = seeds.pred, silent = FALSE,
                                     time = TRUE, verbose = FALSE, Parallel = TRUE,  n.core = n.core)
  }

  # if(method == "rrBLUP"){
  #   scale.nos <- rrBLUP_prediction_iter(Y = Y.scaled.no, Y.test.no = NULL, XF = NULL, x = NULL, Z = Z2, K = K.A,
  #                                       method = "REML", bounds = c(1e-09, 1e+09), SE = TRUE, return.Hinv = TRUE, NDCG.k = NDCG.k,
  #                                       n.iter = n.iter.pred, fitting.scale = fitting.scale, saveAt = dir.names.mid.now,
  #                                       CV = TRUE, n.fold = n.fold, seeds = seeds.pred, silent = TRUE, time = TRUE)
  # }
  #
  # if(method == "RF"){
  #   scale.nos <- RF_prediction_iter(Y = Y.scaled.no, Y.test.no = NULL, XF = NULL, x = NULL, Z = Z2, K = K.A,
  #                                   pca.option = "prop", pca.prop = 0.5, pca.num = 10, NDCG.k = NDCG.k,
  #                                   n.iter = n.iter.pred, fitting.scale = fitting.scale, saveAt = dir.names.mid.now,
  #                                   CV = TRUE, n.fold = n.fold, seeds = seeds.pred, silent = TRUE, time = TRUE)
  # }

  if(fitting.scale == "MSE"){
    scale.nos <- 1 / scale.nos
  }

  scale.no <- mean(scale.nos)




  ##### 2.2. Bayesian Optimization by using gaussian process #####
  convex.no <- choice.no <- c(1, n.combinations)
  scales.obs <- rep(NA, n.combinations)
  scales.obs[choice.no] <- c(scale.no, scale.max)

  save.GP <- paste0(dir.names.mid.now, "/BGLR_files_GP/")
  dir.create(save.GP)
  GP.res.0 <- GP(y = scales.obs, K = NULL, X = all.case.cross.0, error = FALSE,
                 method = "BRR", nIter = 6000, burnIn = 1000,
                 thin = 5, saveAt = save.GP, verbose = FALSE)

  m.0 <- GP.res.0$mean
  sd.0 <- GP.res.0$sd

  GP_name <- paste0(dir.names.mid.now, "/GP_results.RData")
  save(GP.res.0, file = GP_name)



  ##### 2.3. Calculate acquisition value #####
  Acquisition.res <- Acquisition_func(x = costs, y = scales.obs, selected = convex.no,
                                      m = m.0, sd = sd.0)
  Aq <- Acquisition.res$Acquire
  if(Aq.weighted){
    Aq.weight <- 1 / (sqrt(costs) + 1)
    Aq <- Aq * Aq.weight
  }
  Aq.order <- order(Aq, decreasing = T)
  scale.thres <- Acquisition.res$y.thres
  Aq.scaled <- Aq * max(scales.obs[choice.no]) / max(Aq)
  Aq.maxes <- max(Aq[!(Aq.order %in% choice.no)])

  Aq_name <- paste0(dir.names.mid.now, "/Aq_results.RData")
  save(Acquisition.res, file = Aq_name)


  #### 2.3.1. Plot the midstream of Bayesian Optimization ####
  # xlim <- c(0, max(costs))
  # ylim <- c(0, max(scales.obs[choice.no]) + 0.3)
  # png(paste0(dir.names.mid.now, "/cost_cor_and_Acquisition.png"), width = 800, height = 600)
  # plot(costs[setdiff(choice.no, convex.no)],
  #      scales.obs[setdiff(choice.no, convex.no)], xlim = xlim,
  #      ylim = ylim, pch = 16, cex = 0.8, col = "grey",
  #      xlab = "costs", ylab = "R2 or Acquisition (scaled)",
  #      cex.lab = 1.5, cex.axis = 1.3)
  # points(costs[convex.no], scales.obs[convex.no],
  #        pch = 16, cex = 1.4)
  # points(sort(costs), scale.thres[order(costs)], type = "l", lty = 3)
  # points(costs, Aq.scaled, col = 2)
  # points(costs[Aq.order[min(which(!(Aq.order %in% choice.no)))]],
  #        Aq.scaled[Aq.order[min(which(!(Aq.order %in% choice.no)))]],
  #        pch = 17, cex = 3, col = 2)
  # legend(-1, max(scales.obs[choice.no]) + 0.3,
  #        legend = c("R2 (convex hull)", "R2 (other cases)", "R2 (new convex)",
  #                   "R2 (new non-convex)", "Scaled acquisition", "Scaled acquisition (next)"),
  #        pch = c(16, 16, 16, 4, 1, 17), col = c("black", "grey", "green4", "purple4", "red", "red"), cex = 1)
  # dev.off()


  xrange <- c(0, max(costs))
  yrange <- c(0, max(scales.obs[choice.no]) + 0.3)


  state.ggplot <- rep(NA, n.combinations)
  state.ggplot[sort(convex.no)] <- "convex"
  state.ggplot[is.na(state.ggplot)] <- "nonchecked"
  state.ggplot.2 <- factor(c(state.ggplot, rep("acquisition", n.combinations)),
                           levels = c("convex", "acquisition", "nonchecked"))
  df.ggplot <- data.frame(cost = rep(costs, 2), acc = c(scales.obs, Aq.scaled),
                          state = state.ggplot.2)
  df.ggplot <- df.ggplot[as.character(state.ggplot.2) != "nonchecked", ]


  cols <- c("black", "red")
  pchs <- c(16, 2)
  sizes <- c(4, 1.5)
  alphas <- c(1, 0)
  p.ggplot <- ggplot(NULL)
  p.ggplot <- p.ggplot +
    geom_point(data = df.ggplot,
               aes(x = cost, y = acc,
                   colour = state, shape = state,
                   size = state, alpha = state)) +
    scale_colour_manual(values = cols) +
    scale_shape_manual(values = pchs) +
    scale_size_manual(values = sizes) +
    scale_alpha_manual(values = alphas) +
    xlim(xrange) +
    ylim(yrange)

  p.ggplot <- p.ggplot +
    geom_line(data = df.ggplot[df.ggplot$state == "convex", ],
              aes(x = cost, y = acc), size = 1.1) +
    geom_point(data = df.ggplot[df.ggplot$state == "convex", ],
               aes(x = cost, y = acc), size = 4)

  p.ggplot <- p.ggplot +
    labs(x = "Phenotyping cost", y = "Prediction accuracy") +
    theme(axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.key.height = grid::unit(0.8, units = "cm"))


  png(paste0(dir.names.mid.now, "/cost_correlation.png"),
      width = 800, height = 600)
  print(p.ggplot)
  dev.off()



  if(!all(Aq.order %in% choice.no)){
    df.Aq <- data.frame(cost = costs[Aq.order[min(which(!(Aq.order %in% choice.no)))]],
                        acquisition = Aq.scaled[Aq.order[min(which(!(Aq.order %in% choice.no)))]])
    p.ggplot <- p.ggplot +
      geom_point(data = df.Aq,
                 aes(x = cost, y = acquisition), size = 5, shape = 17, colour = 2)
  }

  alphas <- c(1, 0.4)
  p.ggplot <- p.ggplot +
    scale_alpha_manual(values = alphas)


  png(paste0(dir.names.mid.now, "/cost_cor_and_Acquisition.png"),
      width = 800, height = 600)
  print(p.ggplot)
  dev.off()



  ###### 2.3.2. Calculate the current area ######
  cost.convex <- costs[convex.no]
  scales.convex <- scales.obs[convex.no]

  cost.convex.sorted <- sort(cost.convex)
  scales.convex.sorted <- scales.convex[order(cost.convex)]

  Area <- 0

  for(j in 1:(length(cost.convex.sorted) - 1)) {
    cost.height <- cost.convex.sorted[j + 1] - cost.convex.sorted[j]

    Area.now <- (scales.convex.sorted[j] + scales.convex.sorted[j + 1]) *
      cost.height / 2

    Area <- Area + Area.now
  }
  print(paste0("Area : ", round(Area, 2)))

  Areas <- Area


  ##### 2.4. Iterate 2 steps (Peform genomic prediction and perform gaussian process) #####
  count <- 0
  count2 <- 0


  if(convergence.type == "prop"){
    convergence.condition <- count <= (count.prop * nrow(all.case) - 1)
  }

  if(convergence.type == "Acquisition"){
    convergence.condition <- count2 <= count2.max
  }

  while(convergence.condition){
    count <- count + 1
    print(paste("Iteration No. :", count + 2))

    sel.check <- !(Aq.order %in% choice.no)
    case.cand <- Aq.order[min(which(sel.check))]
    print(paste0("No.", case.cand, " : ", paste0(all.case[case.cand, ], collapse = " ")))

    group.NA.TRUE <- sapply(all.case[case.cand, ], function(x) {
      sel.group.here <- sort(sample(1:n.group, n.group - x))
      group.NA.TRUE.now <- group.data[, 3] %in% sel.group.here
    })
    group.NA.TRUE <- cbind(group.NA.TRUE, rep(FALSE, n.lines))
    Y.scaled.now <- Y.max.scaled
    Y.scaled.now[group.NA.TRUE] <- NA
    all.NA.col <- apply(Y.scaled.now, 2, function(x) all(is.na(x)))
    Y.scaled.now <- Y.scaled.now[, !all.NA.col]


    #### 2.4.1. Perform genomic prediction ####
    iter.now <- count
    dir.names.mid.now <- paste0(dir.mid.ROCAPS, iter.now)
    dir.create(dir.names.mid.now)

    if(method == "MTM"){
      resCov <- list(type = "UN", df0 = ncol(Y.scaled.now), S0 = diag(ncol(as.matrix(Y.scaled.now))))
      scale.nows <- MTM_prediction_iter(Y = Y.scaled.now, Y.test.no = NULL, XF = NULL, x = NULL, Z = Z, K = K.A,
                                        resCov = resCov, nIter = 2000, burnIn = 500, thin = 5,
                                        tolD = 1e-05, NDCG.k = NDCG.k, n.iter = n.iter.pred,
                                        fitting.scale = fitting.scale, saveAt = dir.names.mid.now,
                                        CV = TRUE, n.fold = n.fold, seeds = seeds.pred, silent = TRUE,
                                        time = TRUE, verbose = FALSE, Parallel = TRUE,  n.core = n.core)
    }

    # if(method == "rrBLUP"){
    #   scale.nows <- rrBLUP_prediction_iter(Y = Y.scaled.now, Y.test.no = NULL, XF = NULL, x = NULL, Z = Z2, K = K.A,
    #                                        method = "REML", bounds = c(1e-09, 1e+09), SE = TRUE, return.Hinv = TRUE, NDCG.k = NDCG.k,
    #                                        n.iter = n.iter.pred, fitting.scale = fitting.scale, saveAt = dir.names.mid.now,
    #                                        CV = TRUE, n.fold = n.fold, seeds = seeds.pred, silent = TRUE, time = TRUE)
    # }
    #
    # if(method == "RF"){
    #   scale.nows <- RF_prediction_iter(Y = Y.scaled.now, Y.test.no = NULL, XF = NULL, x = NULL, Z = Z2, K = K.A,
    #                                    pca.option = "prop", pca.prop = 0.5, pca.num = 10, NDCG.k = NDCG.k,
    #                                    n.iter = n.iter.pred, fitting.scale = fitting.scale, saveAt = dir.names.mid.now,
    #                                    CV = TRUE, n.fold = n.fold, seeds = seeds.pred, silent = TRUE, time = TRUE)
    # }

    if(fitting.scale == "MSE"){
      scale.nows <- 1 / scale.nows
    }

    scale.now <- mean(scale.nows)

    if(scale.now <= scale.thres[case.cand]){
      choice.no <- c(choice.no, case.cand)
      scales.obs[case.cand] <- scale.now

      plt.col <- "purple4"
      plt.pch <- 4

      print(paste0(fitting.scale, " = ", round(scale.now, 3), ", not exceeding threshold..."))
    }else{
      choice.no <- c(choice.no, case.cand)
      scales.obs[case.cand] <- scale.now
      convex.no <- c(convex.no, case.cand)
      convex.no <- convex_check(costs, scales.obs, convex.no)
      plt.col <- "green4"
      plt.pch <- 16
      print(paste0(fitting.scale, " = ", round(scale.now, 3), ", exceeding threshold!!!"))
    }



    ##### 2.4.2. Bayesian Optimization by using gaussian process #####
    dir.GP.now <- paste0(dir.names.mid.now, "/BGLR_files_GP/")
    dir.create(dir.GP.now)
    GP.res.now <- GP(y = scales.obs, K = NULL, X = all.case.cross.0, error = FALSE,
                     method = "BRR", nIter = 6000, burnIn = 1000,
                     thin = 5, saveAt = dir.GP.now, verbose = FALSE)

    m.now <- GP.res.now$mean
    sd.now <- GP.res.now$sd

    GP_name <- paste0(dir.names.mid.now, "/GP_results.RData")
    save(GP.res.now, file = GP_name)


    ##### 2.4.3. Calculate acquisition value #####
    Acquisition.res <- Acquisition_func(x = costs, y = scales.obs, selected = convex.no, m = m.now, sd = sd.now)
    Aq <- Acquisition.res$Acquire
    if(Aq.weighted){
      Aq.weight <- 1 / (sqrt(costs) + 1)
      Aq <- Aq * Aq.weight
    }
    Aq.order <- order(Aq, decreasing = T)
    scale.thres <- Acquisition.res$y.thres
    Aq.scaled <- Aq * max(scales.obs[choice.no]) / max(Aq)
    if(all(Aq.order %in% choice.no)){
      Aq.max <- 0
    }else{
      Aq.max <- max((Aq[Aq.order])[!(Aq.order %in% choice.no)])
    }

    Aq.maxes <- c(Aq.maxes, Aq.max)
    count2 <- count2 + (Aq.max <= count2.thres)

    Aq_name <- paste0(dir.names.mid.now, "/Aq_results.RData")
    save(Acquisition.res, file = Aq_name)


    #### 2.4.4. Plot the midstream of Bayesian Optimization ####
    # xlim <- c(0, max(costs))
    # ylim <- c(0, max(scales.obs[choice.no]) + 0.3)
    # png(paste0(dir.names.mid.now, "/cost_cor_and_Acquisition.png"), width = 800, height = 600)
    # plot(costs[setdiff(choice.no, convex.no)],
    #      scales.obs[setdiff(choice.no, convex.no)], xlim = xlim,
    #      ylim = ylim, pch = 16, cex = 0.8, col = "grey",
    #      xlab = "costs", ylab = "R2 or Acquisition (scaled)",
    #      cex.lab = 1.5, cex.axis = 1.3)
    # points(costs[convex.no], scales.obs[convex.no],
    #        pch = 16, cex = 1.4)
    # points(sort(costs), scale.thres[order(costs)], type = "l", lty = 3)
    # points(costs[case.cand], scales.obs[case.cand], pch = plt.pch, col = plt.col, lwd = 2, cex = 2)
    # points(costs, Aq.scaled, col = 2)
    # if(!all(Aq.order %in% choice.no)){
    #   points(costs[Aq.order[min(which(!(Aq.order %in% choice.no)))]],
    #          Aq.scaled[Aq.order[min(which(!(Aq.order %in% choice.no)))]],
    #          pch = 17, cex = 3, col = 2)
    # }
    #
    # legend(-1, max(scales.obs[choice.no] + 0.3),
    #        legend = c("R2 (convex hull)", "R2 (other cases)", "R2 (new convex)",
    #                   "R2 (new non-convex)", "Scaled acquisition", "Scaled acquisition (next)"),
    #        pch = c(16, 16, 16, 4, 1, 17), col = c("black", "grey", "green4", "purple4", "red", "red"), cex = 1)
    # dev.off()


    xrange <- c(0, max(costs))
    yrange <- c(0, max(scales.obs[choice.no]) + 0.3)


    state.ggplot <- rep(NA, n.combinations)
    state.ggplot[sort(setdiff(choice.no, convex.no))] <- "checked"
    state.ggplot[sort(convex.no)] <- "convex"
    state.ggplot[is.na(state.ggplot)] <- "nonchecked"
    state.ggplot[case.cand] <- "check_now"
    state.ggplot.2 <- factor(c(state.ggplot, rep("acquisition", n.combinations)),
                             levels = c("convex", "checked", "check_now",
                                        "acquisition", "nonchecked"))
    df.ggplot <- data.frame(cost = rep(costs, 2), acc = c(scales.obs, Aq.scaled),
                            state = state.ggplot.2)
    df.ggplot <- df.ggplot[as.character(state.ggplot.2) != "nonchecked", ]


    cols <- c("black", "gray60", plt.col, "red")
    pchs <- c(16, 16, plt.pch, 2)
    sizes <- c(4, 1.5, 5, 1.5)
    alphas <- c(1, 0.6, 1, 0)

    if(length(setdiff(choice.no, convex.no)) == 0) {
      cols <- cols[-2]
      pchs <- pchs[-2]
      sizes <- sizes[-2]
      alphas <- alphas[-2]
    }
    p.ggplot <- ggplot(NULL)
    p.ggplot <- p.ggplot +
      geom_point(data = df.ggplot,
                 aes(x = cost, y = acc,
                     colour = state, shape = state,
                     size = state, alpha = state)) +
      scale_colour_manual(values = cols) +
      scale_shape_manual(values = pchs) +
      scale_size_manual(values = sizes) +
      scale_alpha_manual(values = alphas) +
      xlim(xrange) +
      ylim(yrange)

    p.ggplot <- p.ggplot +
      geom_line(data = df.ggplot[df.ggplot$state == "convex", ],
                aes(x = cost, y = acc), size = 1.1) +
      geom_point(data = df.ggplot[df.ggplot$state == "convex", ],
                 aes(x = cost, y = acc), size = 4)

    p.ggplot <- p.ggplot +
      labs(x = "Phenotyping cost", y = "Prediction accuracy") +
      theme(axis.title.x = element_text(size = 24),
            axis.title.y = element_text(size = 24),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 16),
            legend.key.height = grid::unit(0.8, units = "cm"))


    png(paste0(dir.names.mid.now, "/cost_correlation.png"),
        width = 800, height = 600)
    print(p.ggplot)
    dev.off()



    if(!all(Aq.order %in% choice.no)){
      df.Aq <- data.frame(cost = costs[Aq.order[min(which(!(Aq.order %in% choice.no)))]],
                          acquisition = Aq.scaled[Aq.order[min(which(!(Aq.order %in% choice.no)))]])
      p.ggplot <- p.ggplot +
        geom_point(data = df.Aq,
                   aes(x = cost, y = acquisition), size = 5, shape = 17, colour = 2)
    }

    alphas <- c(1, 0.6, 1, 0.4)
    if(length(setdiff(choice.no, convex.no)) == 0) {
      alphas <- alphas[-2]
    }
    p.ggplot <- p.ggplot +
      scale_alpha_manual(values = alphas)


    png(paste0(dir.names.mid.now, "/cost_cor_and_Acquisition.png"),
        width = 800, height = 600)
    print(p.ggplot)
    dev.off()



    #### 2.4.5. Convergence condition ####
    if(convergence.type == "prop"){
      convergence.condition <- count <= (count.prop * nrow(all.case))
    }

    if(convergence.type == "Acquisition"){
      convergence.condition <- count2 <= count2.max
    }


    ###### 2.4.6. Save the midstream results ######
    midstream.res <- cbind(all.case[convex.no, ], costs[convex.no], scales.obs[convex.no])
    colnames(midstream.res)[(ncol(midstream.res) - 1):ncol(midstream.res)] <- c("cost", "R2")
    rownames(midstream.res) <- convex.no
    write.csv(midstream.res, paste0(dir.names.mid.now, "/midstream_result.csv"))

    convex.relax.unsorted <- choice.no[scale.thres[choice.no] <= (scales.obs[choice.no] + relax.param)]
    convex.relax <- convex.relax.unsorted[order(costs[convex.relax.unsorted])]

    convex.in <- convex.relax %in% convex.no
    convex.sel <- match(convex.relax, choice.no)

    midstream.res.detail <- cbind(convex.sel, convex.in, all.case[convex.relax, ], costs[convex.relax], scales.obs[convex.relax])
    colnames(midstream.res.detail)[c(1:2, (ncol(midstream.res.detail) - 1):ncol(midstream.res.detail))] <- c("check.order", "convex_or_not", "cost", "R2")
    rownames(midstream.res.detail) <- convex.relax
    write.csv(midstream.res.detail, paste0(dir.names.mid.now, "/midstream_result_detail.csv"))


    ###### 2.4.7. Calculate the current area ######
    cost.convex <- costs[convex.no]
    scales.convex <- scales.obs[convex.no]

    cost.convex.sorted <- sort(cost.convex)
    scales.convex.sorted <- scales.convex[order(cost.convex)]

    Area <- 0

    for(j in 1:(length(cost.convex.sorted) - 1)) {
      cost.height <- cost.convex.sorted[j + 1] - cost.convex.sorted[j]

      Area.now <- (scales.convex.sorted[j] + scales.convex.sorted[j + 1]) *
        cost.height / 2

      Area <- Area + Area.now
    }
    print(paste0("Area : ", round(Area, 2)))

    Areas <- c(Areas, Area)


    Area.max <- max(Areas)
    Areas.sc <- Areas / Area.max

    df.Area <- data.frame(iteration = (1:length(Areas) + 1),
                          Area = Areas.sc)

    p.Area <- ggplot(data = df.Area, aes(x = iteration,
                                         y = Areas.sc)) +
      geom_point(size = 0.5) +
      geom_line(size = 0.5)
    png(paste0(dir.names.mid.now, "/current_Area_Improvement.png"),
        width = 800, height = 600)
    print(p.Area)
    dev.off()



    ###### 2.4.8. Save the all enviroment ######
    save.image(paste0(dir.names.mid.now, "/midstream_all_info.Rdata"))
  }


  ###### 3. Final results ######
  ##### 3.1. Make and save final results' dataframe #####
  dir.res <- paste0(dir.res.0,
                    scriptID, "_", target.trait, "_",
                    method, "_", convergence.type, "_",
                    fitting.scale, "_Aq.weighted=", Aq.weighted, "_",
                    n.iter.pred, "_", n.fold, "_", trial.no.pred, "/")
  dir.create(dir.res)

  final.res <- cbind(all.case[convex.no, ], costs[convex.no], scales.obs[convex.no])
  colnames(final.res)[(ncol(final.res) - 1):ncol(final.res)] <- c("cost", "R2")
  rownames(final.res) <- convex.no
  write.csv(final.res, paste0(dir.res, "/final_result.csv"))

  convex.relax.unsorted <- choice.no[scale.thres[choice.no] <= (scales.obs[choice.no] + relax.param)]
  convex.relax <- convex.relax.unsorted[order(costs[convex.relax.unsorted])]

  convex.in <- convex.relax %in% convex.no
  convex.sel <- match(convex.relax, choice.no)

  final.res.detail <- cbind(convex.sel, convex.in, all.case[convex.relax, ], costs[convex.relax], scales.obs[convex.relax])
  colnames(final.res.detail)[c(1:2, (ncol(final.res.detail) - 1):ncol(final.res.detail))] <- c("check.order", "convex_or_not", "cost", "R2")
  rownames(final.res.detail) <- convex.relax
  write.csv(final.res.detail, paste0(dir.res, "/final_result_detail.csv"))

  save(choice.no, file = paste0(dir.res, "/choice_no.RData"))


  ##### 3.2. Final Plot #####
  # png(paste0(dir.res, "/final_cost_cor_and_Acquisition.png"), width = 800, height = 600)
  # plot(costs[setdiff(choice.no, convex.no)],
  #      scales.obs[setdiff(choice.no, convex.no)], xlim = xlim,
  #      ylim = ylim, pch = 16, cex = 0.8, col = "grey",
  #      xlab = "costs", ylab = "R2 or Acquisition (scaled)",
  #      cex.lab = 1.5, cex.axis = 1.3)
  # points(costs[convex.no], scales.obs[convex.no],
  #        pch = 16, cex = 1.4)
  # points(sort(costs), scale.thres[order(costs)], type = "l", lty = 3)
  # points(costs[case.cand], scales.obs[case.cand], pch = plt.pch, col = plt.col, lwd = 2, cex = 2)
  # points(costs, Aq.scaled, col = 2)
  # if(!all(Aq.order %in% choice.no)){
  #   points(costs[Aq.order[min(which(!(Aq.order %in% choice.no)))]],
  #          Aq.scaled[Aq.order[min(which(!(Aq.order %in% choice.no)))]],
  #          pch = 17, cex = 3, col = 2)
  # }
  # legend(-1, max(scales.obs[choice.no] + 0.3),
  #        legend = c("R2 (convex hull)", "R2 (other cases)", "R2 (new convex)",
  #                   "R2 (new non-convex)", "Scaled acquisition", "Scaled acquisition (next)"),
  #        pch = c(16, 16, 16, 4, 1, 17), col = c("black", "grey", "green4", "purple4", "red", "red"), cex = 1)
  # dev.off()


  xrange <- c(0, max(costs))
  yrange <- c(0, max(scales.obs[choice.no]) + 0.3)


  state.ggplot <- rep(NA, n.combinations)
  state.ggplot[sort(setdiff(choice.no, convex.no))] <- "checked"
  state.ggplot[sort(convex.no)] <- "convex"
  state.ggplot[is.na(state.ggplot)] <- "nonchecked"
  state.ggplot[case.cand] <- "check_now"
  state.ggplot.2 <- factor(c(state.ggplot, rep("acquisition", n.combinations)),
                           levels = c("convex", "checked", "check_now",
                                      "acquisition", "nonchecked"))
  df.ggplot <- data.frame(cost = rep(costs, 2), acc = c(scales.obs, Aq.scaled),
                          state = state.ggplot.2)
  df.ggplot <- df.ggplot[as.character(state.ggplot.2) != "nonchecked", ]


  cols <- c("black", "gray60", plt.col, "red")
  pchs <- c(16, 16, plt.pch, 2)
  sizes <- c(4, 1.5, 5, 1.5)
  alphas <- c(1, 0.6, 1, 0)
  p.ggplot <- ggplot(NULL)
  if(length(setdiff(choice.no, convex.no)) == 0) {
    cols <- cols[-2]
    pchs <- pchs[-2]
    sizes <- sizes[-2]
    alphas <- alphas[-2]
  }
  p.ggplot <- p.ggplot +
    geom_point(data = df.ggplot,
               aes(x = cost, y = acc,
                   colour = state, shape = state,
                   size = state, alpha = state)) +
    scale_colour_manual(values = cols) +
    scale_shape_manual(values = pchs) +
    scale_size_manual(values = sizes) +
    scale_alpha_manual(values = alphas) +
    xlim(xrange) +
    ylim(yrange)

  p.ggplot <- p.ggplot +
    geom_line(data = df.ggplot[df.ggplot$state == "convex", ],
              aes(x = cost, y = acc), size = 1.1) +
    geom_point(data = df.ggplot[df.ggplot$state == "convex", ],
               aes(x = cost, y = acc), size = 4)

  p.ggplot <- p.ggplot +
    labs(x = "Phenotyping cost", y = "Prediction accuracy") +
    theme(axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.key.height = grid::unit(0.8, units = "cm"))


  png(paste0(dir.res, "/final_cost_and_correlation.png"),
      width = 800, height = 600)
  print(p.ggplot)
  dev.off()



  if(!all(Aq.order %in% choice.no)){
    df.Aq <- data.frame(cost = costs[Aq.order[min(which(!(Aq.order %in% choice.no)))]],
                        acquisition = Aq.scaled[Aq.order[min(which(!(Aq.order %in% choice.no)))]])
    p.ggplot <- p.ggplot +
      geom_point(data = df.Aq,
                 aes(x = cost, y = acquisition), size = 5, shape = 17, colour = 2)
  }

  alphas <- c(1, 0.6, 1, 0.4)
  if(length(setdiff(choice.no, convex.no)) == 0) {
    alphas <- alphas[-2]
  }
  p.ggplot <- p.ggplot +
    scale_alpha_manual(values = alphas)


  png(paste0(dir.res, "/final_cost_cor_and_Acquisition.png"),
      width = 800, height = 600)
  print(p.ggplot)
  dev.off()


  Area.min <- min(Areas)
  Area.max <- max(Areas)
  Areas.sc <- (Areas - Area.min) / (Area.max - Area.min)

  df.Area <- data.frame(iteration = (1:length(Areas) + 1),
                        Area = Areas.sc)

  p.Area <- ggplot(data = df.Area, aes(x = iteration,
                                       y = Areas.sc)) +
    geom_point(size = 0.5) +
    geom_line(size = 0.5) +
    labs(x = "Iteration", y = "Area Improvement") +
    theme(axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 16),
          legend.key.height = grid::unit(0.8, units = "cm"))



  png(paste0(dir.res, "/final_Area_Improvement.png"),
      width = 800, height = 600)
  print(p.Area)
  dev.off()
}
