##########################################################################################
######  Title: 3.1_Rice_ROCAPS_Summary_the_results_of_ROCAPS                        ######
######  Author: Kosuke Hamazaki (hamazaki@ut-biomet.org)                            ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2019/12/11 (Created), 2019/12/11 (Last Updated)                       ######
##########################################################################################





###### 1. Settings ######
##### 1.0. Reset workspace ######
# rm(list=ls())



##### 1.1. Setting working directory to the "ROCAPS" directory #####
cropName <- "rice"
project <- "ROCAPS"
os <- osVersion

isRproject <- function(path = getwd()) {
  files <- list.files(path)

  if (length(grep(".Rproj", files)) >= 1) {
    out <- TRUE
  } else {
    out <-  FALSE
  }
  return(out)
}

if (!isRproject()) {
  if (os == "macOS Mojave 10.14.5") {
    dirResearchBase <- "/Users/hamazaki/research/"  ### for mac OS
    dirScriptBase <- "~/research_secret/"
  } else if (os == "Ubuntu 16.04.6 LTS") {
    dirResearchBase <- "/media/hamazaki/d4000953-5e56-40ce-97ea-4cfee57fc91a/research/"   ### for Ubunutu 1
    dirScriptBase <- "~/research_secret/"
  } else if (os == "Ubuntu 18.04.3 LTS") {
    dirResearchBase <- "/media/hamazaki/D6B89D51B89D314B/research/"     ### for Ubuntu 2
    dirScriptBase <- "~/GitHub/research_secret/"
  } else {
    stop("Which type of work-station do you use?")
  }

  setwd(paste0(dirResearchBase, cropName, "/Project/", project))
}




##### 1.2. Setting some parameters #####
load("midstream/ROCAPS_st=21_end=49_by_14days_with_3groups_trial=1/0.2_parameter_informations.RData")

scriptIDBefore <- "2.1"


targetTraitType <- "Sim"
nHerit <- 2
nScenario <- 4

for (heritNow in 1:nHerit) {
  for (scenarioNow in 1:nScenario) {
    # heritNow <- 1
    # scenarioNow <- 1
    print(paste0("herit: ", heritNow, " scenario: ", scenarioNow))

    method <- "MTM"
    convergence.type <- "prop"
    fitting.scale <- "R2"
    Aq.weighted <- FALSE
    n.iter.pred <- 3
    n.fold <- 10
    trial.no.pred <- 1


    dirMid <- dir.mid
    dirMidROCAPSBase <- paste0(dirMid, scriptIDBefore, "_", targetTraitType, "_herit=",
                               heritNow, "_scenario=", scenarioNow, "_",
                               method, "_", convergence.type, "_",
                               fitting.scale, "_Aq.weighted=", Aq.weighted, "_",
                               n.iter.pred, "_", n.fold, "_", trial.no.pred, "/")

    file.params.ROCAPS <- paste0(dirMidROCAPSBase, scriptIDBefore, "_", project, "_all_parameters.RData")
    load(file.params.ROCAPS)

    scriptID <- "3.1"


    dirRes <- paste0("results/ROCAPS_st=21_end=49_by_14days_with_3groups_trial=1/")
    # dir.create(dirRes)
    dirPredResBase <- paste0(dirRes, "predictionResults/")
    # dir.create(dirPredResBase)
    dirPredResEachScenario <- paste0(dirPredResBase, scriptID,
                                     "_", targetTraitType, "_herit=",
                                     heritNow, "_scenario=", scenarioNow, "/")
    dir.create(dirPredResEachScenario)


    nIterSVM <- 10
    nFoldSVM <- 10


    methodClassNames <- c("SVM", "RF")
    ssNames <- c("accuracy", "recall",
                 "precision", "fmeasure")
    predictorClassNames <- c("casecost", "useR2")

    nMethodClass <- length(methodClassNames)
    nSS <- length(ssNames)
    nPredictorClass <- length(predictorClassNames)




    ##### 1.3. Import packages #####
    require(data.table)
    require(RAINBOWR)
    require(ggplot2)
    require(corrplot)
    require(mgcv)
    require(lattice)
    require(colorspace)
    require(kernlab)
    require(ranger)
    library(dplyr)
    library(tidyr)
    library(viridis)
    require(cluster)



    dirScript <- paste0(dirScriptBase, cropName, "/Project/", project, "/")

    # source(paste0(dirScript, "dateLastUpdated"))
    # source(paste0(dirScript, "scriptName2"))



    ##### 1.4. Project options #####
    options(stringAsFactors = FALSE)




    ###### 2. Modification of data & Some preparations ######
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


    if(targetTraitType == "Real"){
      target.trait.nos <- 1:3
      target.trait.names <- colnames(Y.all)[target.trait.nos]

      Y.target <- Y.all[, target.trait.nos]
    } else {
      Y.target <- target.sim.array[, , heritNow, scenarioNow, 3]
      target.trait.names <- colnames(Y.target)
    }



    ##### 2.3. Read simulation results into R & some calculation #####
    #### 2.3.1. Prepare empty boxes ####
    resDetails <- NULL
    resAllChecks <- NULL
    resDetailsTrunc <- NULL


    nTrait <- 10
    corMatArray <- array(NA, dim = c(8, 8, nTrait))
    corMatTruncArray <- array(NA, dim = c(8, 8, nTrait))
    ssMeansArray <- array(NA, dim = c(2, 4, 2, nTrait),
                          dimnames = list(method = methodClassNames,
                                          ss = ssNames,
                                          predictors = predictorClassNames,
                                          trait = target.trait.names[1:nTrait]))

    resSelectedArray <- array(NA, dim = c(4, 8, nTrait))

    for (traitNo in 1:nTrait) {
      print(traitNo)

      #### 2.3.2. Read simulation results into R ####
      target.trait.name <- target.trait.names[traitNo]
      dirMidROCAPSNow <- paste0(dirMidROCAPSBase,
                                scriptIDBefore, "_", target.trait.name, "_",
                                method, "_", convergence.type, "_",
                                fitting.scale, "_Aq.weighted=", Aq.weighted, "_",
                                n.iter.pred, "_", n.fold, "_", trial.no.pred, "/")
      maxIter <- max(as.numeric(list.files(dirMidROCAPSNow)), na.rm = TRUE)

      dirPredResEachTrait <- paste0(dirPredResEachScenario, target.trait.name, "/")
      dir.create(dirPredResEachTrait)

      seedsPredClass <- sample(seed.range, nIterSVM)
      fileNameSeedsPredClass <- paste0(dirPredResEachTrait, "seeds_for_classification.csv")
      # write.csv(seedsPredClass, file = fileNameSeedsPredClass, quote = F)


      fileNameRes <- paste0(dirMidROCAPSNow, maxIter, "/midstream_result.csv")
      fileNameResDetail <- paste0(dirMidROCAPSNow, maxIter, "/midstream_result_detail.csv")
      fileNameAllData <- paste0(dirMidROCAPSNow, maxIter, "/midstream_all_info.Rdata")

      resNow <- read.csv(file = fileNameRes, header = TRUE, row.names = 1)
      resDetailNow <- read.csv(file = fileNameResDetail, header = TRUE, row.names = 1)
      resDetailTruncNow <- resDetailNow[1:which.max(resDetailNow$R2), ]

      load(file = fileNameAllData)



      #### 2.3.3. Construct the results for all the checked points ####
      choiceNoSorted <- sort(choice.no)
      scalesChoice <- scales.obs[choiceNoSorted]
      convexRelaxOrNot <- choiceNoSorted %in% convex.relax.unsorted
      checkOrdChice <- order(choice.no)

      resAllCheck <- cbind(checkOrder = checkOrdChice, convexRelax = as.numeric(convexRelaxOrNot),
                           all.case.cost[choiceNoSorted, ], R2 = scalesChoice)
      fileNameResAllCheck <- paste0(dirMidROCAPSNow, maxIter, "/midstream_result_all.csv")

      fileExistAllData <- "midstream_result_all.csv" %in% list.files(paste0(dirMidROCAPSNow, maxIter))
      if (!fileExistAllData) {
        # write.csv(x = resAllCheck, file = fileNameResAllCheck,
                  # quote = FALSE)
      }



      #### 2.3.4. Calculate correlation matrix ####
      corMatNow <- cor(resDetailNow[, -c(1:2)])
      corMatTruncNow <- cor(resDetailTruncNow[, -c(1:2)])



      #### 2.3.5. Perform classification with SVM and Random Forest ####
      # resAllCheckForSVM <- resAllCheck[, -1]
      # resAllCheckForSVM$convexRelax <- factor(resAllCheckForSVM$convexRelax, levels = c(0, 1))
      #
      #
      # accSVMs <- accRFs <- recallSVMs <- recallRFs <-
      #   precisionSVMs <- precisionRFs <-
      #   fmeasureSVMs <- fmeasureRFs <-
      #   accSVMsUseR2 <- accRFsUseR2 <- recallSVMsUseR2 <-
      #   recallRFsUseR2 <- precisionSVMsUseR2 <- precisionRFsUseR2 <-
      #   fmeasureSVMsUseR2 <- fmeasureRFsUseR2 <- rep(NA, nIterSVM)
      #
      # for (iterNo in 1:nIterSVM) {
      #   dirPredResEachIter <- paste0(dirPredResEachTrait, "iteration_", iterNo, "/")
      #   dir.create(dirPredResEachIter)
      #
      #   set.seed(seed = seedsPredClass[iterNo])
      #   idCV <- sample(nrow(resAllCheck)) %% nFoldSVM
      #   idCV[idCV == 0] <- nFoldSVM
      #
      #   predSVMs <- predRFs <-
      #     predSVMsUseR2 <- predRFsUseR2 <- rep(NA, nrow(resAllCheck))
      #   for (foldNo in 1:nFoldSVM) {
      #     # dirPredResEachFold <- paste0(dirPredResEachIter, "fold_", foldNo, "/")
      #     # dir.create(dirPredResEachFold)
      #
      #     resAllCheckTrain <- resAllCheckForSVM[idCV != foldNo, -9]
      #     resAllCheckTest <- resAllCheckForSVM[idCV == foldNo, -9]
      #
      #     svmRes <- try(ksvm(convexRelax ~ ., data = resAllCheckTrain))
      #
      #     if(class(svmRes) != "try-error") {
      #       predSVM <- predict(svmRes, resAllCheckTest)
      #     } else {
      #       predSVM <- factor(rep(1, nrow(resAllCheckTest)), levels = c(0, 1))
      #     }
      #     predSVMs[idCV == foldNo] <- predSVM
      #
      #     rfRes <- ranger(convexRelax ~ ., data = resAllCheckTrain)
      #     predRF <- predict(rfRes, resAllCheckTest)$predictions
      #     predRFs[idCV == foldNo] <- predRF
      #
      #
      #     resAllCheckTrainUseR2 <- resAllCheckForSVM[idCV != foldNo, ]
      #     resAllCheckTestUseR2 <- resAllCheckForSVM[idCV == foldNo, ]
      #
      #     svmUseR2Res <- try(ksvm(convexRelax ~ ., data = resAllCheckTrainUseR2))
      #     if(class(svmUseR2Res) != "try-error") {
      #       predSVMUseR2 <- predict(svmUseR2Res, resAllCheckTestUseR2)
      #     } else {
      #       predSVMUseR2 <- factor(rep(1, nrow(resAllCheckTestUseR2)), levels = c(0, 1))
      #     }
      #     predSVMsUseR2[idCV == foldNo] <- predSVMUseR2
      #
      #     rfUseR2Res <- ranger(convexRelax ~ ., data = resAllCheckTrainUseR2)
      #     predRFUseR2 <- predict(rfUseR2Res, resAllCheckTestUseR2)$predictions
      #     predRFsUseR2[idCV == foldNo] <- predRFUseR2
      #   }
      #
      #   predSVMs <- predSVMs - 1
      #   predRFs <- predRFs - 1
      #   predSVMsUseR2 <- predSVMsUseR2 - 1
      #   predRFsUseR2 <- predRFsUseR2 - 1
      #
      #   tabSVMRes <- table(factor(predSVMs, levels = c(0, 1)),
      #                      resAllCheckForSVM$convexRelax)
      #   accSVM <- sum(diag(tabSVMRes)) / sum(tabSVMRes)
      #   recallSVM <- tabSVMRes[2, 2] / sum(tabSVMRes[, 2])
      #   precisionSVM <- tabSVMRes[2, 2] / sum(tabSVMRes[2, ])
      #   fmeasureSVM <- 2 * (recallSVM * precisionSVM) / (recallSVM + precisionSVM)
      #
      #   accSVMs[iterNo] <- accSVM
      #   recallSVMs[iterNo] <- recallSVM
      #   precisionSVMs[iterNo] <- precisionSVM
      #   fmeasureSVMs[iterNo] <- fmeasureSVM
      #
      #
      #   tabRFRes <- table(factor(predRFs, levels = c(0, 1)),
      #                     resAllCheckForSVM$convexRelax)
      #   accRF <- sum(diag(tabRFRes)) / sum(tabRFRes)
      #   recallRF <- tabRFRes[2, 2] / sum(tabRFRes[, 2])
      #   precisionRF <- tabRFRes[2, 2] / sum(tabRFRes[2, ])
      #   fmeasureRF <- 2 * (recallRF * precisionRF) / (recallRF + precisionRF)
      #
      #   accRFs[iterNo] <- accRF
      #   recallRFs[iterNo] <- recallRF
      #   precisionRFs[iterNo] <- precisionRF
      #   fmeasureRFs[iterNo] <- fmeasureRF
      #
      #
      #   tabSVMUseR2Res <- table(factor(predSVMsUseR2, levels = c(0, 1)),
      #                           resAllCheckForSVM$convexRelax)
      #   accSVMUseR2 <- sum(diag(tabSVMUseR2Res)) / sum(tabSVMUseR2Res)
      #   recallSVMUseR2 <- tabSVMUseR2Res[2, 2] / sum(tabSVMUseR2Res[, 2])
      #   precisionSVMUseR2 <- tabSVMUseR2Res[2, 2] / sum(tabSVMUseR2Res[2, ])
      #   fmeasureSVMUseR2 <- 2 * (recallSVMUseR2 * precisionSVMUseR2) /
      #     (recallSVMUseR2 + precisionSVMUseR2)
      #
      #   accSVMsUseR2[iterNo] <- accSVMUseR2
      #   recallSVMsUseR2[iterNo] <- recallSVMUseR2
      #   precisionSVMsUseR2[iterNo] <- precisionSVMUseR2
      #   fmeasureSVMsUseR2[iterNo] <- fmeasureSVMUseR2
      #
      #
      #   tabRFUseR2Res <- table(factor(predRFsUseR2, levels = c(0, 1)),
      #                          resAllCheckForSVM$convexRelax)
      #   accRFUseR2 <- sum(diag(tabRFUseR2Res)) / sum(tabRFUseR2Res)
      #   recallRFUseR2 <- tabRFUseR2Res[2, 2] / sum(tabRFUseR2Res[, 2])
      #   precisionRFUseR2 <- tabRFUseR2Res[2, 2] / sum(tabRFUseR2Res[2, ])
      #   fmeasureRFUseR2 <- 2 * (recallRFUseR2 * precisionRFUseR2) /
      #     (recallRFUseR2 + precisionRFUseR2)
      #
      #   accRFsUseR2[iterNo] <- accRFUseR2
      #   recallRFsUseR2[iterNo] <- recallRFUseR2
      #   precisionRFsUseR2[iterNo] <- precisionRFUseR2
      #   fmeasureRFsUseR2[iterNo] <- fmeasureRFUseR2
      #
      #
      #   predResDf <- data.frame(obs = as.numeric(resAllCheckForSVM$convexRelax) - 1,
      #                           predSVM = predSVMs,
      #                           predRF = predRFs,
      #                           predSVMUseR2 = predSVMsUseR2,
      #                           predRFsUseR2 = predRFsUseR2)
      #   fileNamePredRes <- paste0(dirPredResEachIter, "predictedValues.csv")
      #   # write.csv(predResDf, file = fileNamePredRes, quote = F)
      # }
      #
      # accSVMsMean <- mean(accSVMs)
      # accRFsMean <- mean(accRFs)
      # recallSVMsMean <- mean(recallSVMs)
      # recallRFsMean <- mean(recallRFs)
      # precisionSVMsMean <- mean(precisionSVMs)
      # precisionRFsMean <- mean(precisionRFs)
      # fmeasureSVMsMean <- mean(fmeasureSVMs)
      # fmeasureRFsMean <- mean(fmeasureRFs)
      #
      # accSVMsMeanUseR2 <- mean(accSVMsUseR2)
      # accRFsMeanUseR2 <- mean(accRFsUseR2)
      # recallSVMsMeanUseR2 <- mean(recallSVMsUseR2)
      # recallRFsMeanUseR2 <- mean(recallRFsUseR2)
      # precisionSVMsMeanUseR2 <- mean(precisionSVMsUseR2)
      # precisionRFsMeanUseR2 <- mean(precisionRFsUseR2)
      # fmeasureSVMsMeanUseR2 <- mean(fmeasureSVMsUseR2)
      # fmeasureRFsMeanUseR2 <- mean(fmeasureRFsUseR2)
      #
      # ssSVMMean <- c(accSVMsMean, recallSVMsMean,
      #                precisionSVMsMean, fmeasureSVMsMean)
      # ssRFMean <- c(accRFsMean, recallRFsMean,
      #               precisionRFsMean, fmeasureRFsMean)
      # ssSVMMeanUseR2 <- c(accSVMsMeanUseR2, recallSVMsMeanUseR2,
      #                     precisionSVMsMeanUseR2, fmeasureSVMsMeanUseR2)
      # ssRFMeanUseR2 <- c(accRFsMeanUseR2, recallRFsMeanUseR2,
      #                    precisionRFsMeanUseR2, fmeasureRFsMeanUseR2)
      #
      # ssMeanArray <- array(NA, dim = c(2, 4, 2),
      #                      dimnames = list(method = methodClassNames,
      #                                      ss = ssNames,
      #                                      predictors = predictorClassNames))
      # ssMeanArray[1, , 1] <- ssSVMMean
      # ssMeanArray[2, , 1] <- ssRFMean
      # ssMeanArray[1, , 2] <- ssSVMMeanUseR2
      # ssMeanArray[2, , 2] <- ssRFMeanUseR2



      #### 2.3.6. Results for the best points ####
      diffImprvement <- (resNow$R2[-1] - resNow$R2[1]) / resNow$cost[-1]

      if (any(diffImprvement > 0)) {
        maxEfficiencyNo <- which.max(diffImprvement) + 1
      } else {
        maxEfficiencyNo <- 1
      }
      maxR2No <- which.max(resNow$R2)

      selectNo <- c(1, maxEfficiencyNo, maxR2No, nrow(resNow))

      resSelected <- as.matrix(resNow)[selectNo, ]



      #### 2.3.7. Save the results needed in the following steps ####
      resDetails <- rbind(resDetails, resDetailNow)
      resDetailsTrunc <- rbind(resDetailsTrunc, resDetailTruncNow)
      resAllChecks <- rbind(resAllChecks, resAllCheck)
      corMatArray[, , traitNo] <- corMatNow
      corMatTruncArray[, , traitNo] <- corMatTruncNow
      # ssMeansArray[, , , traitNo] <- ssMeanArray
      resSelectedArray[, , traitNo] <- resSelected
    }




    ###### 3. Sumary and plot the results ######
    ##### 3.1. Reorder the combined results over multiple simulations #####
    scriptID <- 3.1

    resDetailsOrd <- resDetails[order(resDetails$cost, resDetails$R2), ]
    resDetailsTruncOrd <- resDetailsTrunc[order(resDetailsTrunc$cost, resDetailsTrunc$R2), ]
    resAllChecksOrd <- resAllChecks[order(resAllChecks$cost, resAllChecks$R2), ]



    ##### 3.2. Corrlearion plots (using "corrplot" package) #####
    corMatMean <- apply(corMatArray, c(1, 2), mean, na.rm = TRUE)
    rownames(corMatMean) <- colnames(corMatMean) <- rownames(corMatNow)

    dirCorrPlot <- paste0(dirRes, scriptID, "_corrplot/")
    # dir.create(dirCorrPlot)
    fileNameCorrPlot <- paste0(dirCorrPlot, scriptID, "_",
                               targetTraitType, "_herit=",
                               heritNow, "_scenario=", scenarioNow,
                               "_corrplot.pdf")
    # pdf(file = fileNameCorrPlot, width = 10, height = 12)
    # corrplot(corMatMean)
    # dev.off()


    corMatTruncMean <- apply(corMatTruncArray, c(1, 2), mean, na.rm = TRUE)
    rownames(corMatTruncMean) <- colnames(corMatTruncMean) <- rownames(corMatTruncNow)

    fileNameCorrPlotTrunc <- paste0(dirCorrPlot, scriptID, "_",
                                    targetTraitType, "_herit=",
                                    heritNow, "_scenario=", scenarioNow,
                                    "_corrplot_truncasion.pdf")
    # pdf(file = fileNameCorrPlotTrunc, width = 10, height = 12)
    # corrplot(corMatTruncMean)
    # dev.off()




    ##### 3.3. Density plots for the selected traits #####
    #### 3.3.1. Density plots (levelplot, using "lattice" package) ####
    newData <- data.frame(cost = seq(min(costs), max(costs), 0.1))
    lenNewData <- nrow(newData)
    phenoCovStart <- 3
    phenoCovEnd <- 8
    phenoCovRange <- phenoCovStart:phenoCovEnd
    nPhenoCov <- phenoCovEnd - phenoCovStart + 1
    predictedMat <- predictedMatTrunc <- matrix(NA, nPhenoCov, lenNewData)

    for (phenoCovNo in phenoCovStart:phenoCovEnd) {
      gamModel <- gam(formula = resDetailsOrd[, phenoCovNo] ~ s(cost),
                      data = resDetailsOrd)
      gamModelTrunc <- gam(formula = resDetailsTruncOrd[, phenoCovNo] ~ s(cost),
                           data = resDetailsTruncOrd)
      predictedValue <- predict(gamModel, newData)
      predictedValueTrunc <- predict(gamModelTrunc, newData)

      predictedMat[phenoCovNo - 2, ] <- predictedValue
      predictedMatTrunc[phenoCovNo - 2, ] <- predictedValueTrunc
    }
    predictedMatModi <- predictedMat
    predictedMatModi[predictedMat < 0] <- 0
    predictedMatModi[predictedMat > 3] <- 3

    predictedMatTruncModi <- predictedMatTrunc
    predictedMatTruncModi[predictedMatTrunc < 0] <- 0
    predictedMatTruncModi[predictedMatTrunc > 3] <- 3

    rownames(predictedMatModi) <- rownames(predictedMatTruncModi) <-
      colnames(all.case)
    colnames(predictedMatModi) <- colnames(predictedMatTruncModi) <- NULL

    fileNameMyColor <- "midstream/ROCAPS_st=21_end=49_by_14days_with_3groups_trial=1/3.1_myColor_for_levelplot.RData"
    load(fileNameMyColor)

    # myColForLevelPlot <- choose_palette()
    # save(myColForLevelPlot, file = fileNameMyColor)

    dirDensityRes <- paste0(dirRes, scriptID, "_densityPlots/")
    # dir.create(dirDensityRes)


    fileNameLevelPlot <- paste0(dirDensityRes, scriptID, "_",
                                targetTraitType, "_herit=",
                                heritNow, "_scenario=", scenarioNow,
                                "_levelplot.pdf")
    pdf(file = fileNameLevelPlot, width = 14, height = 10)
    levelplot(x = t(predictedMatModi[6:1, ]), aspect = "fill",
              col.regions = myColForLevelPlot, xlab = "Cost", ylab = "Plant heights")
    dev.off()


    fileNameLevelPlotTrunc <- paste0(dirDensityRes, scriptID, "_",
                                     targetTraitType, "_herit=",
                                     heritNow, "_scenario=", scenarioNow,
                                     "_levelplot_truncasion.pdf")

    pdf(file = fileNameLevelPlotTrunc, width = 14, height = 10)
    levelplot(x = t(predictedMatTruncModi[6:1, ]), aspect = "fill",
              col.regions = myColForLevelPlot, xlab = "Cost", ylab = "Plant heights")
    dev.off()




    #### 3.3.2. Density plots (ggplot, using "ggplot" package) ####
    resDetailsDens <- data.frame(value = unlist(apply(resDetailsOrd[, phenoCovRange], 2,
                                                      function(x) rep(resDetailsOrd[, 9], x))),
                                 covTrait = rep(colnames(resDetailsOrd[, phenoCovRange]),
                                                apply(resDetailsOrd[, phenoCovRange], 2, sum)))

    pDensity <- ggplot(data = resDetailsDens,
                       aes(x = value, group = covTrait, fill = covTrait)) +
      geom_density(adjust = 1.5, alpha = .4)


    fileNameDensityPlot <- paste0(dirDensityRes, scriptID, "_",
                                  targetTraitType, "_herit=",
                                  heritNow, "_scenario=", scenarioNow,
                                  "_densityplot.pdf")
    # pdf(fileNameDensityPlot, width = 12, height = 10)
    # print(pDensity)
    # dev.off()


    resDetailsTruncDens <- data.frame(value = unlist(apply(resDetailsTruncOrd[, phenoCovRange], 2,
                                                           function(x) rep(resDetailsTruncOrd[, 9], x))),
                                      covTrait = rep(colnames(resDetailsTruncOrd[, phenoCovRange]),
                                                     apply(resDetailsTruncOrd[, phenoCovRange], 2, sum)))


    pDensityTrunc <- ggplot(data = resDetailsTruncDens,
                            aes(x = value, group = covTrait, fill = covTrait)) +
      geom_density(adjust = 1.5, alpha = .4)

    fileNameDensityPlotTrunc <- paste0(dirDensityRes, scriptID, "_",
                                       targetTraitType, "_herit=",
                                       heritNow, "_scenario=", scenarioNow,
                                       "_densityplot_truncasion.pdf")
    # pdf(fileNameDensityPlotTrunc, width = 12, height = 10)
    # print(pDensityTrunc)
    # dev.off()


    ##### 3.4. Perform PCA, k-medoids clustering and compare the results  #####
    resAllChecksCase <- resAllChecksOrd[, phenoCovRange]

    resAllChecksCaseString <- apply(resAllChecksCase, 1, paste0, collapse = "")
    resAllChecksGrp0 <- as.numeric(factor(resAllChecksCaseString))
    resAllChecksGrp <- match(resAllChecksGrp0, unique(resAllChecksGrp0))
    resAllChecksGrpUnduplicated <- !duplicated(resAllChecksGrp)

    resAllChecksMean <-
      apply(X = resAllChecksOrd, 2, function(resAllChecksEachcol) {
        tapply(X = resAllChecksEachcol,
               INDEX = resAllChecksGrp,
               FUN = mean)
      })



    pcaRes <- prcomp(resAllChecksMean[, phenoCovRange])
    # pcaRes <- prcomp(all.case)
    plot(pcaRes$x[, 1], pcaRes$x[, 2])


    nPlotOrd <- 100
    plotNow <- resAllChecksMean[, 1] < nPlotOrd

    pamRes <- pam(x = resAllChecksMean[, phenoCovRange], k = 6)
    grpPam <- pamRes$clustering

    dirPCARes <- paste0(dirRes, scriptID, "_pcaResults/")
    # dir.create(dirPCARes)

    fileNamePCARes <- paste0(dirPCARes, scriptID, "_",
                             targetTraitType, "_herit=",
                             heritNow, "_scenario=", scenarioNow,
                             "_pca_results_with_", nPlotOrd, "_points.pdf")

    # pdf(fileNamePCARes, width = 14, height = 10)
    # plot(pcaRes$x[, 1], pcaRes$x[, 2], col = grpPam,
    #      cex = resAllChecksMean[, 2] + 0.5, type = "n")
    # text(pcaRes$x[plotNow, 1], pcaRes$x[plotNow, 2],
    #      labels = round(resAllChecksMean[plotNow, 1], 0), col = grpPam[plotNow],
    #      cex = resAllChecksMean[plotNow, 2] + 0.5)
    # dev.off()


    plot(pcaRes$x[, 1], pcaRes$x[, 2], col = grpPam,
         cex = resAllChecksMean[, 2] + 0.5, type = "p")

    plot(pcaRes$x[, 1], pcaRes$x[, 2], col = ifelse(resAllChecksMean[, 2] >= 0.5, 2, 1),
         pch = 16, type = "p")
    plot(pcaRes$x[, 3], pcaRes$x[, 4], col = ifelse(resAllChecksMean[, 2] >= 0.5, 2, 1),
         pch = 16, type = "p")
    plot(pcaRes$x[, 5], pcaRes$x[, 6], col = ifelse(resAllChecksMean[, 2] >= 0.5, 2, 1),
         pch = 16, type = "p")




    ##### 3.5. Plot and save the results for the best points #####
    rownames(resSelectedArray) <- c("GP", "Best_Eff", "Best_R2", "Max_cost")
    colnames(resSelectedArray) <- colnames(resNow)

    xRange <- range(costs)
    yRange <- c(0, max(resSelectedArray[, 8, ]) + 0.05)
    p.ggplotSc <- ggplot(NULL) +
      xlim(xRange) +
      ylim(yRange) +
      labs(x = "Phenotyping cost", y = "Prediction accuracy") +
      theme(axis.title.x = element_text(size = 24),
            axis.title.y = element_text(size = 24),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 16),
            legend.key.height = grid::unit(0.8, units = "cm"))

    for (traitNo in 1:nTrait) {
      resScPlotDf <- data.frame(resSelectedArray[, 7:8, traitNo])
      p.ggplotSc <- p.ggplotSc +
        geom_point(data = resScPlotDf,
                   aes(x = cost, y = R2),
                   size = c(2, 3.5, 3.5, 2),
                   col = c("gray20", "red", "blue", "gray20"),
                   alpha = 0.8) +
        geom_line(data = resScPlotDf,
                  aes(x = cost, y = R2),
                  size = 1.1,
                  col = "gray50",
                  alpha = 0.7)
    }


    dirSelected <- paste0(dirRes, scriptID, "_", "selectedPoints/")
    # dir.create(dirSelected)


    fileNameSelectedScatter <- paste0(dirSelected, scriptID, "_",
                                      targetTraitType, "_herit=",
                                      heritNow, "_scenario=", scenarioNow,
                                      "_scatter_plots_for_best_points.png")

    png(fileNameSelectedScatter,
        width = 800, height = 600)
    print(p.ggplotSc)
    dev.off()



    resSelectedMean <- round(apply(resSelectedArray[, 7:8, ], c(1, 2), mean), 4)

    fileNameSelectedMean <- paste0(dirSelected, scriptID, "_",
                                   targetTraitType, "_herit=",
                                   heritNow, "_scenario=", scenarioNow,
                                   "_improvement_for_best_points.csv")
    # write.csv(resSelectedMean, file = fileNameSelectedMean, quote = F)




    ##### 3.6. Summary and save the results for the classification #####
    ssArrayRes <- apply(ssMeansArray, c(1, 2, 3), mean, na.rm = TRUE)

    dirClassificationRes <- paste0(dirRes, scriptID, "_classificationResults/")
    # dir.create(dirClassificationRes)

    fileNameClassSS <- paste0(dirClassificationRes, scriptID, "_",
                              targetTraitType, "_herit=",
                              heritNow, "_scenario=", scenarioNow,
                              "_classification_results.RData")
    # save(ssArrayRes, file = fileNameClassSS)
  }
}



#### 3.6.1. Plot the results for the classification (using all the scenarios) ####
allClassResExist <- length(list.files(dirClassificationRes)) == 8


if (allClassResExist) {
  allClassDf <- NULL

  for (heritNo in 1:nHerit) {
    for (scenarioNo in 1:nScenario) {
      fileNameClassSSNow <- paste0(dirClassificationRes, scriptID, "_",
                                   targetTraitType, "_herit=",
                                   heritNo, "_scenario=", scenarioNo,
                                   "_classification_results.RData")

      load(fileNameClassSSNow)
      ssVecRes <- c(ssArrayRes)
      methodClassDfNames <- factor(rep(methodClassNames, nSS * nPredictorClass),
                                   labels = methodClassNames)
      ssClassDfNames <- factor(rep(rep(ssNames, each = nMethodClass), nPredictorClass),
                               labels = ssNames)
      predictorClassDfNames <- factor(rep(predictorClassNames, each = nMethodClass * nSS),
                                      labels = predictorClassNames)

      scenarioDfNames <- rep(paste0(heritNo, "-", scenarioNo),
                             nMethodClass * nSS * nPredictorClass)

      dfNow <- data.frame(value = ssVecRes,
                          method = methodClassDfNames,
                          statistics = ssClassDfNames,
                          predictors = predictorClassDfNames,
                          scenario = scenarioDfNames)

      allClassDf <- rbind(allClassDf, dfNow)
    }
  }


  dim(allClassDf)

  allClassCaseCostDf <- allClassDf[allClassDf$predictors == "casecost", ]
  allClassUseR2Df <- allClassDf[allClassDf$predictors == "useR2", ]

  allClassCaseCostSVMDf <- allClassCaseCostDf[allClassCaseCostDf$method == "SVM", ]
  allClassCaseCostRFDf <- allClassCaseCostDf[allClassCaseCostDf$method == "RF", ]
  allClassUseR2SVMDf <- allClassUseR2Df[allClassUseR2Df$method == "SVM", ]
  allClassUseR2RFDf <- allClassUseR2Df[allClassUseR2Df$method == "RF", ]

  p.ClassCaseCostSVM <- ggplot(data = allClassCaseCostSVMDf, aes(x = scenario, y = value, fill = statistics)) +
    geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
    ylim(0, 1) +
    scale_color_manual(values = c("red", "green", "blue", "purple")) +
    ggtitle(paste0("Classification Result")) +
    theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 25, face = "bold"),
          axis.text.y = element_text(size = 26),
          legend.text = element_text(size = 23, face = "bold"),
          legend.title = element_text(size = 18),
          legend.key.height = unit(31, units = "pt"))

  p.ClassCaseCostRF <- ggplot(data = allClassCaseCostRFDf, aes(x = scenario, y = value, fill = statistics)) +
    geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
    ylim(0, 1) +
    scale_color_manual(values = c("red", "green", "blue", "purple")) +
    ggtitle(paste0("Classification Result")) +
    theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 25, face = "bold"),
          axis.text.y = element_text(size = 26),
          legend.text = element_text(size = 23, face = "bold"),
          legend.title = element_text(size = 18),
          legend.key.height = unit(31, units = "pt"))

  p.ClassUseR2SVM <- ggplot(data = allClassUseR2SVMDf, aes(x = scenario, y = value, fill = statistics)) +
    geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
    ylim(0, 1) +
    scale_color_manual(values = c("red", "green", "blue", "purple")) +
    ggtitle(paste0("Classification Result")) +
    theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 25, face = "bold"),
          axis.text.y = element_text(size = 26),
          legend.text = element_text(size = 23, face = "bold"),
          legend.title = element_text(size = 18),
          legend.key.height = unit(31, units = "pt"))

  p.ClassUseR2RF <- ggplot(data = allClassUseR2RFDf, aes(x = scenario, y = value, fill = statistics)) +
    geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
    ylim(0, 1) +
    scale_color_manual(values = c("red", "green", "blue", "purple")) +
    ggtitle(paste0("Classification Result")) +
    theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 25, face = "bold"),
          axis.text.y = element_text(size = 26),
          legend.text = element_text(size = 23, face = "bold"),
          legend.title = element_text(size = 18),
          legend.key.height = unit(31, units = "pt"))


  fileNameClassCaseCostSVM <- paste0(dirClassificationRes, scriptID,
                                     "_classification_results_case&cost_SVM.png")
  fileNameClassCaseCostRF <- paste0(dirClassificationRes, scriptID,
                                    "_classification_results_case&cost_RF.png")
  fileNameClassUseR2SVM <- paste0(dirClassificationRes, scriptID,
                                  "_classification_results_case&cost&R2_SVM.png")
  fileNameClassUseR2RF <- paste0(dirClassificationRes, scriptID,
                                 "_classification_results_case&cost&R2_RF.png")

  # png(fileNameClassCaseCostSVM, height = 700, width = 800)
  # print(p.ClassCaseCostSVM)
  # dev.off()
  #
  # png(fileNameClassCaseCostRF, height = 700, width = 800)
  # print(p.ClassCaseCostRF)
  # dev.off()
  #
  # png(fileNameClassUseR2SVM, height = 700, width = 800)
  # print(p.ClassUseR2SVM)
  # dev.off()
  #
  # png(fileNameClassUseR2RF, height = 700, width = 800)
  # print(p.ClassUseR2RF)
  # dev.off()
}
