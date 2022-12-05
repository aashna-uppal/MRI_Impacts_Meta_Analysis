#-------------------------------------
#-------------------------------------
#meta-analysis for MRI impacts study
#modified on 2022-12-05 by A Uppal
#-------------------------------------
#-------------------------------------

#load libraries
library(metafor)
library(meta)
library(readxl)
library(tidyverse)
library(here)
library(cowplot)

#locate code
here::i_am("meta_analysis_MRI.R")

#read in data
dat <- read_xlsx("MRI impact data.xlsx")

#-------------------------------------
#-------------------------------------

#create function with 3 inputs
meta_analysis <- function(category, variable, surveilance){
  
  if(surveilance == "all"){dat <- dat}
  else if(surveilance == "low"){dat <- dat %>% filter(surv_intensity == "low")}
  else if(surveilance == "high"){dat <- dat %>% filter(surv_intensity == "high")}
  
  #set up metaanalysis variables
  ies_clinical <- escalc(xi = dat[, variable] %>% pull(), 
                         ni = dat$LR, 
                         data = dat, 
                         measure = "PR") #indicating no transformation
  #meta analysis
  pes <- rma(yi, vi, data = ies_clinical, method = "REML")
  
  #print pooled estimates
  paste0("Category: ",  category, " - ",
         "Pooled Estimate: ", round(pes$beta, 2), 
         ", Confidence Interval: (", round(pes$ci.lb, 2), 
         ", ", round(pes$ci.ub, 2), ")")
}

#-------------------------------------
#-------------------------------------

#run function for different categories
#3 = LR_clinical, 4 = LR_imaging, 5 = LR_MRI
print <- c(meta_analysis(category = "Clinical, Overall", variable = 3, surveilance = "all"),
meta_analysis(category = "Imaging, Overall", variable = 4, surveilance = "all"),
meta_analysis(category = "MRI, Overall", variable = 5, surveilance = "all"),

meta_analysis(category = "Clinical, High Intesity", variable = 3, surveilance = "high"),
meta_analysis(category = "Imaging, High Intesity", variable = 4, surveilance = "high"),
meta_analysis(category = "MRI, High Intesity", variable = 5, surveilance = "high"),

meta_analysis(category = "Clinical, Low Intensity", variable = 3, surveilance = "low"),
meta_analysis(category = "Imaging, Low Intensity", variable = 4, surveilance = "low"),
meta_analysis(category = "MRI, Low Intensity", variable = 5, surveilance = "low"))

#save into csv file
write.csv(print, "Results.csv")

#-------------------------------------
#forest plots-------------------------
#-------------------------------------

#subset data to get rid of indeterminate surv intensity
dat <- dat %>% 
  filter(surv_intensity != "indeterminate") %>%
  rename(`Surveillance Intensity` = surv_intensity)

#for MRI analysis only, subset data to get rid of studies with no value for prop detected by MRI
dat2 <- dat %>% 
  filter(!(is.na(LR_MRI)))

#-------------------------------------
#-------------------------------------

#Clinical
pes.summary <- metaprop(LR_clinical, LR, study, data = dat, sm = "PRAW", subgroup = `Surveillance Intensity`)
forest(pes.summary, 
       common = F, 
       text.random = "Pooled Proportion Detected Clinically", 
       leftcols=c("studlab", "event", "n"),
       leftlabs=c("Study", "LR Detected Clinically", "Total LR"), 
       rightlabs=c("Proportion", "95% CI", "Weights"),
       pooled.totals=FALSE, 
       weight.study = "random")

#Imaging
pes.summary.2 <- metaprop(LR_imaging, LR, study, data = dat, sm = "PRAW", subgroup = `Surveillance Intensity`)
forest(pes.summary.2, 
       common = F, 
       text.random = "Pooled Proportion Detected by Imaging", 
       leftcols=c("studlab", "event", "n"),
       leftlabs=c("Study", "LR Detected by Imaging", "Total LR"), 
       rightlabs=c("Proportion", "95% CI", "Weights"), 
       pooled.totals=FALSE, 
       weight.study = "random")

#MRI
pes.summary.3 <- metaprop(LR_MRI, LR, study, data = dat2, sm = "PRAW", subgroup = `Surveillance Intensity`)
forest(pes.summary.3, 
       common = F, 
       text.random = "Pooled Proportion Detected by MRI", 
       leftcols=c("studlab", "event", "n"),
       leftlabs=c("Study", "LR Detected by MRI", "Total LR"), 
       rightlabs=c("Proportion", "95% CI", "Weights"), 
       pooled.totals=FALSE, 
       weight.study = "random")


