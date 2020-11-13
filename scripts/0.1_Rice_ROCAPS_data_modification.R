###############################################
###### 0.1_Rice_ROCAPS_data_modification ######
###############################################

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
scriptID <- 0.1



##### 1.2. Import packages #####
require(data.table)
require(rrBLUP)
require(GGally)



##### 1.3. Setting some parameters #####




###### 2. Modification of raw data ######
##### 2.1. Read original files into R #####
geno.raw <- data.frame(fread(file = "raw_data/genotype/geno.csv", check.names = FALSE), check.names = FALSE)
pheno.raw <- as.matrix(read.csv(file = "raw_data/phenotype/pheno.csv", row.names = 1))

dim(geno.raw)
dim(pheno.raw)

head(geno.raw)[, 1:6]    ### codings with (-1, 0, 1)
head(pheno.raw)[, 1:6]   ### including NA




##### 2.2. Modification of genotype data (including the step which removes MAF <= 0.05) #####
colnames(geno.raw)[-c(1:3)] <- paste0("NSFTV_", colnames(geno.raw)[-c(1:3)])
x.0 <- t(geno.raw[, -c(1:3)])
map.0 <- geno.raw[, 1:3]


freq <- apply(X = x.0, MARGIN = 2,
              FUN = function(x){
                mean(x + 1) / 2
              })
MAF <- pmin(freq, 1 - freq)
MAF.005 <- MAF >= 0.05

geno.mat <- x.0[, MAF.005]
map <- map.0[MAF.005, ]
colnames(map) <- c("marker", "chr", "pos")
colnames(geno.mat) <- map[, 1]




file.map.csv <- paste0("data/genotype/", scriptID,"_genotype_map.csv")
write.csv(x = map, file = file.map.csv, row.names = T)



##### 2.3. Modification of phenotype data and select some traits #####
rownames(pheno.raw) <- paste0("NSFTV_", rownames(pheno.raw))
sum(rownames(geno.mat) != rownames(pheno.raw))  ### Check if phenotype and genotype data match.

(trait.names.raw <- colnames(pheno.raw))  ### Show trait names
trait.names.sel.relax <- trait.names.raw[c(1:3, 8, 12, 13)]
trait.names.sel <- trait.names.raw[c(1, 8, 12, 13)]
pheno.sel.relax <- pheno.raw[, trait.names.sel.relax]
pheno.sel <- pheno.raw[, trait.names.sel]
ggpairs(as.data.frame(pheno.sel.relax), progress = FALSE)   ### Show pairs plot


colnames(pheno.sel) <- c("FT", "FLL", "PH", "PL")


nonNA <- apply(pheno.sel, 1, function(x) !any(is.na(x)))
sum(nonNA)


##### 2.4. Remove NA lines and save genotype and phenotype #####
geno.mat <- geno.mat[nonNA, ]
file.geno.mat.csv <- paste0("data/genotype/", scriptID, "_genotype_matrix.csv")
file.geno.mat.RData <- paste0("data/genotype/", scriptID, "_genotype_matrix.RData")

write.csv(x = geno.mat, file = file.geno.mat.csv, row.names = T)
save(geno.mat, file = file.geno.mat.RData)


pheno.sel <- pheno.sel[nonNA, ]
file.pheno.sel.csv <- paste0("data/phenotype/", scriptID, "_phenotype_selected.csv")
file.pheno.sel.RData <- paste0("data/phenotype/", scriptID, "_phenotype_selected.RData")

write.csv(x = pheno.sel, file = file.pheno.sel.csv, row.names = T)
save(pheno.sel, file = file.pheno.sel.RData)





##### 2.5. Calculate additive genetic relationship matrix #####
K.A <- A.mat(geno.mat)
file.amat.csv <- paste0("data/genotype/", scriptID, "_additive_genetic_relationship_matrix.csv")
file.amat.RData <- paste0("data/genotype/", scriptID, "_additive_genetic_relationship_matrix.RData")

write.csv(x = K.A, file = file.amat.csv, row.names = T)
save(K.A, file = file.amat.RData)


