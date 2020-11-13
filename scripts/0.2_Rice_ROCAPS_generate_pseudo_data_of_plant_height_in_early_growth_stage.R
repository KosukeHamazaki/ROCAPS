########################################################################################
###### 0.2_Rice_ROCAPS_generate_pseudo_data_of_plant_height_in_early_growth_stage ######
########################################################################################

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
scriptID <- 0.2



##### 1.2. Import packages #####
require(data.table)
require(MASS)
require(rrBLUP)
require(psych)
require(GGally)
require(EMMREML)

source(paste0("/home/hamazaki/GitHub/research_secret/",
              crop_name, "/Project/", project, "/1.1.0_Rice_ROCAPS_MTM2.R"))


##### 1.3. Setting some parameters #####
trial.no <- 1

PH.0 <- 15
logistic.d <- 5
logistic.a <- 1
logistic.stop.after.heading <- 3
logistic.x.zoom <- 10

eval.start <- 21
eval.end <- 49
eval.space <- 14
eval.day <- seq(eval.start, eval.end, by = eval.space)

base.sd.hand <- 2.5
base.sd.drone <- 10

n.iter <- 20
n.cluster <- 10

n.group <- 3

seed.range <- 1:10000

seeds <- sample(seed.range, 6)
names(seeds) <- c("seed.hand", "seed.drone", "seed.kmeans",
                  "seed.group", paste0("seed.target.", 1:2))


cost.maintain.base <- 10
cost.hand.base <- c(6, 9, 12) 
cost.drone.base <- 1
cost.drone.plus <- 0.5
cost.drone.own <- 30

eval.methods <- c("hand", "drone")
n.eval.method <- length(eval.methods)
n.eval.day <- length(eval.day)


GC.start <- 20
GC.end  <- 100
GC.space <- 20

GC.day <- seq(GC.start, GC.end, by = GC.space)
n.GC <- length(GC.day)

corsG.max <- c(0.3, 0.8, 0.3, 0.8)
corsG.max.point <- c(n.GC, n.GC, n.GC - 2, n.GC - 2)

corsE.max <- 0.4
corsE.max.point <- n.GC

VG.target <- 1
herits.target <- c(0.3, 0.7)

n.scenario <- length(corsG.max)
n.herit <- length(herits.target)
n.sim <- 20
no.T <- 1

dir.target <- paste0("data/phenotype/", scriptID, "_", n.sim,
                     "_simulated_target_values_with", n.herit,
                     "_herits_", n.scenario, "_scenarios/")
dir.create(dir.target)
dir.mid <- paste0("midstream/", project, "_st=", eval.start, "_end=", eval.end,
                  "_by_", eval.space, "days_with_", n.group, "groups_trial=",
                  trial.no, "/")
dir.create(dir.mid, showWarnings = FALSE, recursive = TRUE)

file.params <- paste0(dir.mid, scriptID, "_parameter_informations.RData")
save.image(file = file.params)

file.seeds <- paste0(dir.mid, scriptID, "_seeds.csv")
write.csv(seeds, file = file.seeds, row.names = F)





###### 2. Sigmoid Curve ######
##### 2.1. Read phenotypic data into R #####
load("data/phenotype/0.1_phenotype_selected.RData")
head(pheno.sel)

n.line <- nrow(pheno.sel)
FT <- pheno.sel[, 1]
PH <- pheno.sel[, 3]



##### 2.2. Simulate growth curve #####
file.growth.curve <- paste0(dir.mid, scriptID,"_growth_curve_for_all_accessions.png")

PH.EGS <- matrix(NA, nrow = n.line, ncol = length(eval.day))
png(file.growth.curve, height = 800, width = 900)
curve(logistic(x * logistic.x.zoom / (FT[1] + logistic.stop.after.heading),
               d = logistic.d, a = logistic.a, z = PH[1] - PH.0) + PH.0,
      xlab = "Day After Planting", 
      ylab = "Plant Height (cm)",
      xlim = c(0, max(FT) + 10),
      ylim = c(0, max(PH)))
title(main = "Simulated Growth Curve")
for(i in 1:n.line){
  curve(logistic(x * logistic.x.zoom / (FT[i] + logistic.stop.after.heading),
                 d = logistic.d, a = logistic.a, z = PH[i] - PH.0) + PH.0, add = T)
  PH.EGS.now <- logistic(eval.day * logistic.x.zoom / (FT[i] + logistic.stop.after.heading),
                         d = logistic.d, a = logistic.a, z = PH[i] - PH.0) + PH.0
  PH.EGS[i, ] <- PH.EGS.now
}
dev.off()



##### 2.3. Add error to simulated plant height (evaluated with hand and drone) #####
PH.EGS.hand.sd <- apply(PH.EGS, 2, function(x) mean(x) * base.sd.hand / 100)
set.seed(seeds[1])
PH.EGS.hand.error <- sapply(PH.EGS.hand.sd, function(x) rnorm(n.line, mean = 0, sd = x))
PH.EGS.hand <- PH.EGS + PH.EGS.hand.error

PH.EGS.drone.sd <- apply(PH.EGS, 2, function(x) mean(x) * base.sd.drone / 100)
set.seed(seeds[2])
PH.EGS.drone.error <- sapply(PH.EGS.drone.sd, function(x) rnorm(n.line, mean = 0, sd = x))
PH.EGS.drone <- PH.EGS + PH.EGS.drone.error


PH.all <- data.frame(PH, PH.EGS, PH.EGS.hand, PH.EGS.drone)
colnames(PH.all) <- c("PH", paste0("PH_", eval.day, "DAP"),
                      paste0("PH_hand_", eval.day, "DAP"),
                      paste0("PH_drone_", eval.day, "DAP"))


file.pairs.plot <- paste0(dir.mid, scriptID, "_pairs_plot_for_plant_height.png")
png(file.pairs.plot, height = 1200, width = 1200)
ggpairs(PH.all)
dev.off()




###### 3. Divide accessions into n.group groups ######
load("data/genotype/0.1_genotype_matrix.RData")
pca.res <- prcomp(geno.mat)

kmeans.clusters <- kmeans.evals <- NULL
set.seed(seeds[3])
for(i in 1:n.iter){
  kmeans.res <- kmeans(x = geno.mat, centers = n.cluster, iter.max = 100)
  kmeans.clusters <- c(kmeans.clusters, list(kmeans.res$cluster))
  kmeans.evals <- c(kmeans.evals, list(kmeans.res$betweenss / kmeans.res$totss))
}

which.max.cluster <- which(unlist(kmeans.evals) == max(unlist(kmeans.evals)))
kmeans.cluster <- kmeans.clusters[[which.max.cluster[1]]]

n.line.each.group.mat <- matrix(rep(table(kmeans.cluster) %/% n.group, each = n.group), nrow = n.group)
plus.no <- table(kmeans.cluster) %% n.group

n.line.each.group.plus.mat <- matrix(0, nrow = n.group, ncol = n.cluster)
save.last <- 0
for(i in 1:n.cluster){
  if(plus.no[i] == 0){
    next
  }else{
    next.no <- (save.last + 1:plus.no[i]) %% n.group
    next.no[next.no == 0] <- n.group
    n.line.each.group.plus.mat[next.no, i] <- 1
    save.last <- next.no[length(next.no)]
  }
}


n.line.each.group.mat <- n.line.each.group.mat + n.line.each.group.plus.mat
rownames(n.line.each.group.mat) <- paste0("group_", 1:n.group)
colnames(n.line.each.group.mat) <- paste0("cluster_", 1:n.cluster)

splits <- NULL
groups <- kmeans.cluster
set.seed(seeds[4])
for(i in 1:n.cluster){
  n.cluster.now <- table(kmeans.cluster)[i]
  split.now <- sample(1:n.cluster.now) %% n.group
  split.now[split.now == 0] <- n.group
  
  if(any(n.line.each.group.plus.mat[, i] != 0)){
    bango.max.now <- which(table(split.now) == as.numeric(rownames(table(table(split.now)))[which.min(table(table(split.now)))]))
    bango.max <- which(n.line.each.group.plus.mat[, i] == as.numeric(rownames(table(n.line.each.group.plus.mat[, i]))[which.min(table(n.line.each.group.plus.mat[, i]))]))
    
    split.now[!is.na(match(split.now, bango.max))] <- n.group + bango.max
    
    split.now[!is.na(match(split.now, bango.max.now))] <- bango.max
    split.now[!is.na(match(split.now, bango.max + n.group))] <- bango.max.now
  }
  splits <- c(splits, list(split.now))
  groups[kmeans.cluster == i] <- split.now
}

group.data <- data.frame(NSFTV = names(kmeans.cluster),
                         cluster = kmeans.cluster,
                         group = groups)

file.group <- paste0("data/extra/", scriptID, "_accession_group_information.csv")
write.csv(group.data, file.group, row.names = F)



Y.all <- data.frame(pheno.sel[, c(2, 4)], PH.all)
colnames(Y.all)

file.pheno.csv <- paste0("data/phenotype/", scriptID, "_observed_and_simulated_phenotype_for_ROCAPS.csv")
file.pheno.RData <- paste0("data/phenotype/", scriptID, "_observed_and_simulated_phenotype_for_ROCAPS.RData")

write.csv(Y.all, file = file.pheno.csv)
save(Y.all, file = file.pheno.RData)





###### 4. All possible combinations ######
n.combination <- n.eval.method * n.eval.day
list.0_1.0 <- rep(list(0:n.group), n.combination)
all.case.0 <- as.matrix(expand.grid(list.0_1.0))
colnames(all.case.0) <- paste0(rep(eval.methods, n.eval.day),
                               "_", rep(eval.day, each = n.eval.method))

all.case.cross.0 <- all.case.0
for(i in 1:n.combination){
  for(j in i:n.combination){
    all.case.cross.now.0 <- all.case.0[, i] * all.case.0[, j]
    all.case.cross.0 <- cbind(all.case.cross.0, all.case.cross.now.0)
  }
}

colnames(all.case.cross.0) <- 1:((n.combination + 1) * n.combination / 2 + n.combination)

n.possible <- nrow(all.case.0)
rownames(all.case.0) <- rownames(all.case.cross.0) <- 1:n.possible

file.case <- paste0(dir.mid, scriptID, "_all_case.RData")
file.case.cross <- paste0(dir.mid, scriptID, "_all_case_with_cross.RData")

save(all.case.0, file = file.case)
save(all.case.cross.0, file = file.case.cross)




###### 5. Costs data ######
all.case.hand <- all.case.0[, seq(1, (2 * n.eval.day - 1), by = 2)]
cost.hand.eff <- cost.hand.base
cost.hand <- all.case.hand %*% cost.hand.eff


all.case.drone <- all.case.0[, seq(2, (2 * n.eval.day), by = 2)]
cost.drone.each <- apply(all.case.drone, 2, function(x){
  pmax(x - 1, 0) * cost.drone.plus + ifelse(x >= 1, cost.drone.base, 0)
})

cost.drone <- apply(cost.drone.each, 1, sum)
cost.drone[cost.drone != 0] <- cost.drone[cost.drone != 0] + cost.drone.own

maintain.last.0 <- apply(all.case.0, 1, function(x) max(which(x != 0)))
maintain.last.0[1] <- 0
maintain.last <- ((maintain.last.0 + 1) %/% 2)
maintain.last.week <- maintain.last * 2 + 1
maintain.last.week[1] <- 0

cost.maintain <- maintain.last.week * cost.maintain.base

cost.all <- c(cost.hand) + c(cost.drone) + cost.maintain

all.case.cost <- data.frame(all.case.0, cost = cost.all)
head(all.case.cost)

file.case.cost <- paste0(dir.mid, scriptID, "_all_case_with_cost.csv")
write.csv(all.case.cost, file = file.case.cost)




###### 5. Simulate phenotypic values of target trait ######
##### 5.1. Estimate genetic correlations of PH.GC #####
PH.GC.mat <- matrix(NA, nrow = n.line, ncol = n.GC)
for(i in 1:n.line){
  PH.GC.now <- logistic(GC.day * logistic.x.zoom / (FT[i] + logistic.stop.after.heading),
                        d = logistic.d, a = logistic.a, z = PH[i] - PH.0) + PH.0
  PH.GC.mat[i, ] <- PH.GC.now
}
rownames(PH.GC.mat) <- rownames(pheno.sel)
colnames(PH.GC.mat) <- paste0("GC_", GC.day)


load("data/genotype/0.1_additive_genetic_relationship_matrix.RData")


X <- t(matrix(rep(1, n.line)))
colnames(X) <- rownames(pheno.sel)


Z <- diag(n.line)
rownames(Z) <- colnames(Z) <- rownames(pheno.sel)

# EMM.res <- emmremlMultivariate(Y = t(PH.GC.mat), X = X, Z = Z, K = K.A)
# VG.mat <- EMM.res$Vg
# VE.mat <- EMM.res$Ve
# sigmaGs <- sqrt(diag(VG.mat))
# CG.mat <- VG.mat / tcrossprod(sigmaGs)
# sigmaEs <- sqrt(diag(VE.mat))
# CE.mat <- VE.mat / tcrossprod(sigmaEs)
# U <- t(EMM.res$Gpred)
# E <- PH.GC.mat - U - crossprod(X, t(EMM.res$Bhat))


dir.MTM <- paste0(dir.mid, scriptID, "_MTM_results/")
dir.create(dir.MTM)
saveAt <- dir.MTM
MTM.res <- MTM2(Y = PH.GC.mat,
                K = list(
                  list(
                    K = K.A,
                    COV = list(
                      type = 'UN',
                      df0 = n.GC,
                      S0 = diag(n.GC)
                    )
                  )
                ),
                resCov = list(type = "UN", df0 = n.GC, S0 = diag(n.GC)), 
                nIter = 6000, burnIn = 1200, saveAt = saveAt, thin = 5)

file.MTM.res <- paste0(dir.MTM, "MTM_results_summary.RData")
save(MTM.res, file = file.MTM.res)


VG.mat <- MTM.res$K[[1]]$G
VE.mat <- MTM.res$resCov$R
sigmaGs <- sqrt(diag(VG.mat))
CG.mat <- VG.mat / tcrossprod(sigmaGs)
sigmaEs <- sqrt(diag(VE.mat))
CE.mat <- VE.mat / tcrossprod(sigmaEs)
U <- MTM.res$K[[1]]$U
E <- PH.GC.mat - U - crossprod(X, t(MTM.res$mu))




##### 5.2. Define genetic correlation between target trait and PH.GC #####
scenario.names <- paste0("Scenario_", 1:n.scenario)

# corsG.start <- rep(0, 4)
# corsG.max <- c(0.8, 0.4, 0.9, 0.24)
# corsG.end <- c(0.8, 0.4, 0.6, 0.16)
# corsG.max.point <- c(n.GC, n.GC, 6, 6)
# corsG.params <- rbind(corsG.start, corsG.max, corsG.end, corsG.max.point)
# colnames(corsG.params) <- scenario.names
# 
# corsG.target <- t(apply(corsG.params, 2, function(x) {
#   corsG.start <- x[1]
#   corsG.max <- x[2]
#   corsG.end <- x[3]
#   corsG.max.point <- x[4]
#   res <- c(1, seq(corsG.start, corsG.max, length.out = corsG.max.point + 1),
#                   seq(corsG.max, corsG.end, length.out = n.GC - corsG.max.point + 1)[-1])
#   
#   return(res)
# })[-2, ])
# colnames(corsG.target) <- c("Target", paste0("GC_", GC.day))


corsG.params <- rbind(corsG.max, corsG.max.point)
colnames(corsG.params) <- scenario.names

corsG.target <- t(apply(corsG.params, 2, function(x) {
  c(1, round(CE.mat[x[2],  ], 2) * x[1])
}))


# corsE.start <- 0
# corsE.end <- 0.4
# corsE.target <- c(1, seq(corsE.start, corsE.end, length.out = n.GC + 2)[-c(1:2)])


corsE.target <- c(1, round(CE.mat[corsE.max.point,  ], 2) * corsE.max)
names(corsE.target) <- c("Target", paste0("GC_", GC.day))


VE.target <- (1 - herits.target) / herits.target


cors.target <- rbind(corsG.target, corsE.target)
rownames(cors.target) <- c(paste0("genetic_Scenario_", 1:n.scenario),
                           "residual")
file.cors.target <- paste0(dir.mid, scriptID,
                           "_genetic_and_residual_correlations_with_target_trait.csv")
write.csv(cors.target, file = file.cors.target)




##### 5.3. Generate phenotypic values of target traits #####
# target.sim.list.0 <- rep(list(NULL), 4)
# names(target.sim.list.0) <- scenario.names
# target.sim.list <- rep(list(target.sim.list.0), 2)
# names(target.sim.list) <- paste0("herit_", herits.target)
# target.sim.list.all <- rep(list(target.sim.list), 3)
# names(target.sim.list.all) <- c("geno", "resid", "pheno")

target.sim.array <- array(NA, dim = c(n.line, n.sim, n.herit, n.scenario, 3),
                          dimnames = list(Line = rownames(geno.mat),
                                          Sim = paste0("sim_", 1:n.sim),
                                          Herit = paste0("herit_", herits.target),
                                          Scenario = scenario.names,
                                          Value = c("geno", "resid", "pheno")))


for (herit.no in 1:n.herit){
  for(scenario.no in 1:n.scenario){
    CG.all <- rbind(corsG.target[scenario.no, ],
                    cbind(corsG.target[scenario.no, -1], CG.mat))
    rownames(CG.all)[no.T] <- "Target"
    
    sigmaGs.all <- c(Target = VG.target, sigmaGs)
    VG.all <- tcrossprod(sigmaGs.all) * CG.all
    
    gene.cor.list <- list(correlation = CG.all, 
                          variance = sigmaGs.all)
    file.gene.cor <- paste0(dir.target, scriptID,
                            "_genetic_correlations_and_variance_herit=", 
                            herits.target[herit.no], "scenario_",
                            scenario.no, ".RData")
    save(gene.cor.list, file = file.gene.cor)
    
    VG.all.inv <- solve(VG.all)
    
    SigmaG.target <- K.A / VG.all.inv[no.T, no.T]
    
    muG.target <- -kronecker(t(as.matrix(VG.all.inv[no.T, -no.T])) / VG.all.inv[no.T, no.T],
                             diag(n.line)) %*% c(U)
    set.seed(seeds[5])
    targetG.sim <- t(mvrnorm(n = n.sim, mu = muG.target, Sigma = SigmaG.target))
    
    
    
    CE.all <- rbind(corsE.target,
                    cbind(corsE.target[-1], CE.mat))
    rownames(CE.all)[no.T] <- "Target"
    
    sigmaEs.all <- c(Target = VE.target[herit.no], sigmaEs)
    VE.all <- tcrossprod(sigmaEs.all) * CE.all
    
    resid.cor.list <- list(correlation = CE.all, 
                          variance = sigmaEs.all)
    file.resid.cor <- paste0(dir.target, scriptID,
                            "_residual_correlations_and_variance_herit=", 
                            herits.target[herit.no], "scenario_",
                            scenario.no, ".RData")
    save(resid.cor.list, file = file.resid.cor)
    
    
    VE.all.inv <- solve(VE.all)
    
    SigmaE.target <- diag(n.line) / VE.all.inv[no.T, no.T]
    
    muE.target <- -kronecker(t(as.matrix(VE.all.inv[no.T, -no.T])) / VE.all.inv[no.T, no.T],
                             diag(n.line)) %*% c(E)
    set.seed(seeds[6])
    targetE.sim <- t(mvrnorm(n = n.sim, mu = muE.target, Sigma = SigmaE.target))
    
    
    
    target.sim <- targetG.sim + targetE.sim
    
    
    target.sim.array[, , herit.no, scenario.no, 1] <- targetG.sim
    target.sim.array[, , herit.no, scenario.no, 2] <- targetE.sim
    target.sim.array[, , herit.no, scenario.no, 3] <- target.sim
    
    file.targetG <- paste0(dir.target, scriptID,
                           "_simulated_genotypic_values_of_target_traits_herit=",
                           herits.target[herit.no], "_scenario_", scenario.no, ".csv")
    file.targetE <- paste0(dir.target, scriptID,
                           "_simulated_residuals_of_target_traits_herit=",
                           herits.target[herit.no], "_scenario_", scenario.no, ".csv")
    file.target <- paste0(dir.target, scriptID,
                          "_simulated_phenotypic_values_of_target_traits_herit=",
                          herits.target[herit.no], "_scenario_", scenario.no, ".csv")
    
    write.csv(targetG.sim, file = file.targetG)
    write.csv(targetE.sim, file = file.targetE)
    write.csv(target.sim, file = file.target)
  }
}
file.target.array <- paste0("data/phenotype/", scriptID,
                            "_all_simulated_values_of_target_traits.RData")
save(target.sim.array, file = file.target.array)