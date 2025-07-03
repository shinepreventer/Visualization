
#dat$Method <- as.factor(dat$Method)

library(plspm)

comdata <- list(
  Method      = c("Method") ,
  inputC      =  c("input"),
  metabolism1 = c("Com114pos","Com8351pos"),
  metabolism2 = c("Com34743pos", "Com10922pos"), 
  metabolism3 = c("Com24162neg","Com14299pos"),
  Soil = c("soil","Nitrification","DenitrificationIntensity","Mineralization"),
  Microbe = c("Diversity", "CCycle", "DissimilatoryProcesses", "AssimilatoryProcesses","Denitrification", "NitrogenTransformationFixation", "NitrogenAssimilationandMine", "UptakeProcesses") , 
  Plantnutrients = c("RN", "STN", "LEN"),
  Plantgrowth = c("Root",  "Growth")
)

latent_vars <- c("Method", "inputC", "metabolism1", "metabolism2", "metabolism3","Soil","Microbe", "Plantnutrients","Plantgrowth")

model_formulas <- list(
  Soil = ~ metabolism1 + metabolism2 + metabolism3 + inputC + Method,
  Microbe = ~ metabolism1 + metabolism2 + metabolism3+ inputC + Soil +Method,
  Plantgrowth = ~  metabolism1 + metabolism2 + metabolism3 + inputC + Soil + Microbe+ Plantnutrients + Method,
  Plantnutrients = ~ metabolism1 + metabolism2 + metabolism3 + inputC + Soil + Microbe + Method,
  metabolism2 = ~ inputC + Method,
  metabolism3 = ~ inputC + metabolism2 + Method,
  metabolism1 = ~ inputC + Method ,
  inputC = ~ Method
)


dat_path <- matrix(0, nrow = length(latent_vars), ncol = length(latent_vars),                  
                   dimnames = list(latent_vars, latent_vars))
for (end_var in names(model_formulas)) { 
  rhs_vars <- all.vars(model_formulas[[end_var]])  
  dat_path[end_var, rhs_vars] <- 1
}
dat_path <- dat_path * lower.tri(dat_path)
print(dat_path)

dat_modes <- rep('A', 9)

dat <- na.omit(dat)  
dat_modes
dat_pls <- plspm(dat, dat_path, comdata, modes = dat_modes)
dat_pls$gof
dat_pls
summary(dat_pls)


save.image("1.RData")
rm(list = ls())


##################################################################################################


dat <-read.csv("all.csv", header=TRUE) 
dat$Method <- as.factor(dat$Method)
library(plspm)
comdata <- list(
  Method      = c("Method") ,
  inputC      =  c("input"),
  metabolism = c("Com114pos","Com34743pos", "Com10922pos","Com24162neg","Com14299pos","Com8351pos"),
  Soil = c(   "soil","Nitrification","DenitrificationIntensity","Mineralization"),
  Microbe = c("Diversity", "CCycle", "DissimilatoryProcesses", "AssimilatoryProcesses","Denitrification", "NitrogenTransformationFixation", "NitrogenAssimilationandMine", "UptakeProcesses") ,
  Plantnutrients = c("RN", "STN", "LEN"),
  Plantgrowth = c("Root",  "Growth")
)
latent_vars <- c("Method" , "inputC", "metabolism","Soil","Microbe","Plantnutrients","Plantgrowth")
model_formulas <- list(
  Soil = ~ metabolism +  inputC + Method,
  Microbe = ~ metabolism+ inputC+ Soil + Method,
  Plantgrowth = ~  metabolism + inputC + Soil + Microbe + Plantnutrients + Method,
  Plantnutrients = ~ metabolism + inputC + Soil + Microbe + Method,
  inputC = ~ Method ,
  metabolism = ~ Method
)
dat_path <- matrix(0, nrow = length(latent_vars), ncol = length(latent_vars),                  
                   dimnames = list(latent_vars, latent_vars))
for (end_var in names(model_formulas)) { 
  rhs_vars <- all.vars(model_formulas[[end_var]])  
  dat_path[end_var, rhs_vars] <- 1
}
dat_path <- dat_path * lower.tri(dat_path)
print(dat_path)

dat_modes <- rep('A', 7)

dat <- na.omit(dat)  
dat_modes
dat_pls <- plspm(dat, dat_path, comdata, modes = dat_modes)
dat_pls$gof
dat_pls
summary(dat_pls)

save.image("2.RData")
rm(list = ls())


