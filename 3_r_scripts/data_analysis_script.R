###
# Written by: Conor Taff, Jenny Uehling
# Last updated: 4/11/2021
# Run under R Studio XX on XX

# This is the main analysis script for analyzing data after processing them
# in the data_filtering_script.R (located in this repo).

################################################################################
# Load libraries ----

pacman::p_load("tidyverse", "plyr", "dplyr", "here", "ggpubr", "data.table", "car", 
               "sjPlot", "gtools", "lme4", "MuMIn", "gridExtra")

# tidyverse, plyr, and dplyr for data wrangling
# here for file reference by relative paths
# ggpubr for plotting
# data.table for creating data tables
# car for ??
# sjPlot for ??
# gtools for using logit transformation
# lme4 and MuMIn for modeling and AIC calculation
# gridExtra for making AIC tables

################################################################################
# Description of data frames

# These are the major data frames we will use:

# "aquatic" contains a row for every sample (prefixes: L19N or P19N)
aquatic <- read.csv("2_modified_data/aquatic.csv")

# "aquatic_nestlingsadults" contains a row for every nestling, and information
# about each nestling's mother in the same row.
aquatic_nestlingsadults <- read.csv("2_modified_data/aquatic_nestlingsadults.csv")

################################################################################
# Question 2: What predicts aquatic insect content in diet? ----
# Question 2A: Do age, site, or an interaction between the two predict
# insect content in diet?
################################################################################

# We will use the aquatic data frame, but use only adult samples from birds 
# captured during provisioning.
aquatic_prov_adF <- aquatic[aquatic$age == "Adult" & aquatic$sex == "F" ,] # Select females
aquatic_prov_adF <- aquatic_prov_adF[aquatic_prov_adF$cap_num == "3" ,] # Select females from third capture only (during provisioning)
# aquatic_prov_adM <- aquatic[aquatic$age == "Adult" & aquatic$sex == "M" ,] # Select males
# Note that all males were captured during provisioning so no need to filter any further
aquatic_prov_nestlings <- aquatic[aquatic$age == "Nestling" ,] # Select nestlings

# Bind together all of these samples from provisioning
aquatic_prov <- rbind(aquatic_prov_adF, aquatic_prov_nestlings)

#############
# MODEL: Site, age, and percent aquatic
#############

# Use gtools package to perform a logit transformation on percentage data
aquatic_prov$percent_aquatic_ra_logit <- logit(aquatic_prov$percent_aquatic_ra, min = -0.0001, max = 1.0001)

# Check the structure of each variable before modeling
str(aquatic_prov$age)
str(aquatic_prov$site)
str(aquatic_prov$site_box_year)
str(aquatic_prov$percent_aquatic_ra_logit)
str(aquatic_prov$exp_treat)

q2a <- lmer(percent_aquatic_ra_logit ~ age*site + exp_treat + (1|site_box_year), data = aquatic_prov, REML = FALSE)

# Checking conditions
hist(resid(q2a))
plot(q2a)
tab_model(q2a, file = here("3_r_scripts/q2a_mod_tab.html"))

#############
# PLOT: Site, age, and percent aquatic
#############

p_q2a <- ggplot(aquatic_prov) +
  geom_boxplot(width = 0.5, (aes(x=age, y=percent_aquatic_ra, fill=age))) + 
  facet_wrap(~ site, ncol = 2) +
  xlab("Age") +
  ylab("Percent Aquatic Insects in Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

ggsave(here("3_r_scripts/lh_age_days_site.png"), p_q2a, width = 10, height = 8, device = "png")

################################################################################
# Question 2: What predicts aquatic insect content in diet? ----
# Question 2B: Does a mother's CORT and mass predict the percent aquatic in the diet of her nestlings?
################################################################################

# We will use the aquatic_nestlingsadults data frame because we're looking at
# the relationship between mom's phenotypic traits and her nestlings' diets.

# Use gtools package to perform a logit transformation on percentage data
aquatic_nestlingsadults$n_percent_aquatic_ra_logit <- logit(aquatic_nestlingsadults$n_percent_aquatic_ra, min = -0.0001, max = 1.0001)

# Calculate adult stress response
aquatic_nestlingsadults$a_stress_response <- aquatic_nestlingsadults$ad_Bleed2_CORT_corrected_value - aquatic_nestlingsadults$ad_Bleed1_CORT_corrected_value

# Note that, because we are using model selection for this question, we need to
# make sure that there are no "NA" values for any of the variables we will use
# in the models. We also need to visually inspect the data for outliers or abnormalities.

#############
# PLOT: Mom's CORT + percent aquatic
#############

p_base <- ggplot(aquatic_nestlingsadults) +
  geom_point((aes(x=ad_Bleed1_CORT_corrected_value, y=n_percent_aquatic_ra))) + 
  xlab("Female CORT") +
  ylab("Percent Aquatic Insects in Nestling Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

# This plot shows that 4.15 is a big outlier (30 for baseline CORT). We will remove that individual.
# Take out 4.15
aquatic_nestlingsadults_min4.15 <- aquatic_nestlingsadults[aquatic_nestlingsadults$site_box_year != "Unit_4_15_2019" ,]

p_base <- ggplot(aquatic_nestlingsadults_min4.15) +
  geom_point((aes(x=ad_Bleed1_CORT_corrected_value, y=n_percent_aquatic_ra))) + 
  xlab("Female CORT") +
  ylab("Percent Aquatic Insects in Nestling Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))
# This plot looks better; no outliers.

# How many NAs?
count(is.na(aquatic_nestlingsadults_min4.15$ad_Bleed1_CORT_corrected_value)) # 8 NAs, 271 values
# Take out all rows with NA for baseline CORT
aquatic_nestlingsadults_min4.15 <- aquatic_nestlingsadults_min4.15[!is.na(aquatic_nestlingsadults_min4.15$ad_Bleed1_CORT_corrected_value) ,]

p_sr <- ggplot(aquatic_nestlingsadults_min4.15) +
  geom_point((aes(x=ad_Bleed2_CORT_corrected_value-ad_Bleed1_CORT_corrected_value, y=n_percent_aquatic_ra))) + 
  xlab("Female CORT") +
  ylab("Percent Aquatic Insects in Nestling Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))
# How many NAs?
count(is.na(aquatic_nestlingsadults_min4.15$ad_Bleed2_CORT_corrected_value)) # 119 NAs, 160 values
# There are a lot of NAs here because we did not take stress and dex bleeds from Unit 4 and Turkey Hill birds. We will not use stress response in the models.

#############
# PLOT: Mom's mass + percent aquatic
#############

p_mass <- ggplot(aquatic_nestlingsadults_min4.15) +
  geom_point((aes(x=ad_Mass, y=n_percent_aquatic_ra))) + 
  xlab("Female mass") +
  ylab("Percent Aquatic Insects in Nestling Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))
# No major outliers or biologically unrealistic values

# How many NAs?
count(is.na(aquatic_nestlingsadults_min4.15$ad_Mass)) # None!

#############
# PLOT: Mom's wing + percent aquatic
#############

p_wing <- ggplot(aquatic_nestlingsadults_min4.15) +
  geom_point((aes(x=ad_Flat_Wing, y=n_percent_aquatic_ra))) + 
  xlab("Female flatwing") +
  ylab("Percent Aquatic Insects in Nestling Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))
# No major outliers or biologically unrealistic values

# How many NAs?
count(is.na(aquatic_nestlingsadults_min4.15$ad_Flat_Wing)) # None!

#############
# PLOT: Mom's treatment + percent aquatic
#############

p_exp_treat <- ggplot(aquatic_nestlingsadults_min4.15) +
  geom_boxplot((aes(x=n_exp_treat, y=n_percent_aquatic_ra))) + 
  xlab("Experimental Treatment") +
  ylab("Percent Aquatic Insects in Nestling Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

# How many NAs?
count(is.na(aquatic_nestlingsadults_min4.15$n_exp_treat)) # None!

#############
# PLOT: site + percent aquatic
#############

p_site <- ggplot(aquatic_nestlingsadults_min4.15) +
  geom_boxplot((aes(x=site, y=n_percent_aquatic_ra))) + 
  xlab("Site") +
  ylab("Percent Aquatic Insects in Nestling Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

# How many NAs?
count(is.na(aquatic_nestlingsadults_min4.15$site)) # None!

ggsave(here("3_r_scripts/site_n_diet.png"), p_site, width = 7, height = 7, device = "png")

#############
# PLOT: brood size + percent aquatic
#############

p_brood <- ggplot(aquatic_nestlingsadults_min4.15) +
  geom_point((aes(x=n_Brood_Size_Day6, y=n_percent_aquatic_ra))) + 
  xlab("Brood Size at Day 6") +
  ylab("Percent Aquatic Insects in Nestling Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

# How many NAs?
count(is.na(aquatic_nestlingsadults_min4.15$n_Brood_Size_Day6)) # None!

#############
# PLOT: nestling age + percent aquatic
#############

p_age <- ggplot(aquatic_nestlingsadults_min4.15) +
  geom_boxplot((aes(x=as.factor(n_age.1), y=n_percent_aquatic_ra))) + 
  xlab("Nestling Age") +
  ylab("Percent Aquatic Insects in Nestling Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

# How many NAs?
count(is.na(aquatic_nestlingsadults_min4.15$n_age.1)) # None!


#############
# MODEL SELECTION
#############

# First check structure of each variable in the model to make sure it is correct
str(aquatic_nestlingsadults_min4.15$n_percent_aquatic_ra_logit)
str(aquatic_nestlingsadults_min4.15$ad_Bleed1_CORT_corrected_value)
str(aquatic_nestlingsadults_min4.15$site)
str(aquatic_nestlingsadults_min4.15$n_age.1) # Currently classifying this as an integer; should it be a character?
str(aquatic_nestlingsadults_min4.15$n_exp_treat)
str(aquatic_nestlingsadults_min4.15$ad_Mass)
str(aquatic_nestlingsadults_min4.15$ad_Flat_Wing)
str(aquatic_nestlingsadults_min4.15$n_Brood_Size_Day6)

## Read in list of models to be fit
mods <- read.csv(here("3_r_scripts/model_lists/mod_list_q2b.csv"))
mods<-as.character(mods[,2])

## Loop through model list and fit each model, saving fit model as object in a list
mod.output <- list()
for(i in 1:length(mods)){
  model <- paste("fit.mod",i,sep="_")
  mod <- lmer(as.formula(mods[i]), data = aquatic_nestlingsadults_min4.15, REML = FALSE)
  assign(model[1],mod)
  mod.output[[i]] <- get(model)
}

# Create AIC table based on the loop of models above
AICtab <- as.data.frame(mods)
for(i in 1:nrow(AICtab)){
  AICtab$AICc[i] <- AICc(mod.output[[i]])
  AICtab$LL[i] <- logLik(mod.output[[i]])
  AICtab$K[i] <- attr(logLik(mod.output[[i]]),"df")
}
minim <- min(AICtab$AICc)
AICtab$D.AICc <- AICtab$AICc - minim
AICtab$rellik <- exp(-.5* AICtab$D.AICc)
sumAICs <- sum(AICtab$rellik)
for(i in 1:nrow(AICtab)){
  AICtab$Wi[i] <- round(AICtab$rellik[i]/sumAICs,2)
}
AICtab <- AICtab[order(-AICtab$rellik),]

# save the AIC table to a csv
dev.off()
write.csv(AICtab, here("3_r_scripts/model_outputs/AIC_Output_q2b.csv"))

### Model 1
q2b_mod1 <- lmer(n_percent_aquatic_ra_logit ~ ad_Mass + site + + n_age.1 + n_exp_treat + (1|site_box_year),
                 data = aquatic_nestlingsadults_min4.15, REML = FALSE)

# Checking conditions
hist(resid(q2b_mod1))
plot(q2b_mod1)
tab_model(q2b_mod1, file = here("3_r_scripts/model_outputs/q2b_mod1_tab.html"))

### Model 2
q2b_mod2 <- lmer(n_percent_aquatic_ra_logit ~ n_Brood_Size_Day6 + site + n_age.1 + n_exp_treat + (1|site_box_year),
                 data = aquatic_nestlingsadults_min4.15, REML = FALSE)

# Checking conditions
hist(resid(q2b_mod2))
plot(q2b_mod2)
tab_model(q2b_mod2, file = here("3_r_scripts/model_outputs/q2b_mod2_tab.html"))

### Model 3
q2b_mod3 <- lmer(n_percent_aquatic_ra_logit ~ ad_Flat_Wing + site + n_age.1 + n_exp_treat + (1|site_box_year),
                data = aquatic_nestlingsadults, REML = FALSE)

# Checking conditions
hist(resid(q2b_mod3))
plot(q2b_mod3)
tab_model(q2b_mod3, file = here("3_r_scripts/model_outputs/q2b_mod3_tab.html"))

### Model 4
q2b_mod4 <- lmer(n_percent_aquatic_ra_logit ~ ad_Bleed1_CORT_corrected_value + site + n_age.1 + n_exp_treat + (1|site_box_year),
                 data = aquatic_nestlingsadults, REML = FALSE)

# Checking conditions
hist(resid(q2b_mod4))
plot(q2b_mod4)
tab_model(q2b_mod4, file = here("3_r_scripts/model_outputs/q2b_mod4_tab.html"))

################################################################################
# Question 2: What predicts aquatic insect content in diet? ----
# Question 2C: Does a female's CORT and mass predict the percent aquatic in her diet?
################################################################################

aquatic_ad <- aquatic[aquatic$age == "Adult" ,]
aquatic_adF <- aquatic_ad[aquatic_ad$sex == "F" ,]

# Use gtools package to perform a logit transformation on percentage data
aquatic_adF$percent_aquatic_ra_logit <- logit(aquatic_adF$percent_aquatic_ra, min = -0.0001, max = 1.0001)

# Some third capture females do not have wing measurements, so import the wing measurements from their first capture
# However, some of these first wing measurements are recorded during captures when there was not a
# fecal sample taken. We need to go back into the raw data to extract what we need.
captures <- read.csv(here("1_raw_data", "Captures_Hormone_Bleeding_Blood_DNA_11.18.2020.csv"))
captures <- captures[!is.na(captures$Exp_Year), ]
captures <- captures[captures$Exp_Year == "2019" ,]
captures <- captures[captures$Adult_or_Nestling == "Adult" ,]
captures <- captures[captures$Sex == "F" ,]

captures_wing <- captures[captures$Capture_Number == "1" ,]
cols_needed <- as.vector(c("Individual_Band", "Flat_Wing"))
captures_wing <- captures_wing[,cols_needed]

# Rename column so that it can be merged with aquatic_adF
captures_wing <- dplyr::rename(captures_wing, band = Individual_Band)

# Remove duplicated bands
captures_wing <- captures_wing[!duplicated(captures_wing$band) ,]

# Add in new column with wing measurements from first captures
aquatic_adF <- left_join(aquatic_adF, captures_wing)

# Delete unneeded columns
aquatic_adF <- aquatic_adF %>% select (-c(species, buty, bolus_nb, fate, fwing))

#############
# PLOT: CORT + percent aquatic
#############

p_base <- ggplot(aquatic_adF) +
  geom_point((aes(x=cort1_correct, y=percent_aquatic_ra))) + 
  xlab("Female CORT") +
  ylab("Percent Aquatic Insects in Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))
# There are a few outliers here; notably, two individuals have baseline CORT levels above 30. Removed these.
aquatic_adF_nooutliers <- aquatic_adF[aquatic_adF$blood_nb != "190790" & aquatic_adF$blood_nb != "190105" ,]

p_base <- ggplot(aquatic_adF_nooutliers) +
  geom_point((aes(x=cort1_correct, y=percent_aquatic_ra))) + 
  xlab("Female CORT") +
  ylab("Percent Aquatic Insects in Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))
# There are a few outliers here; notably, two individuals have baseline CORT levels above 30. Removed these.

# How many NAs?
count(is.na(aquatic_adF_nooutliers$cort1_correct)) # 3 NAs, 172 values
# Take out all rows with NA for baseline CORT
aquatic_adF_nooutliers <- aquatic_adF_nooutliers[!is.na(aquatic_adF_nooutliers$cort1_correct) ,]
count(is.na(aquatic_adF_nooutliers$cort1_correct)) # No missing values now.

#############
# PLOT: mass + percent aquatic
#############

p_mass <- ggplot(aquatic_adF_nooutliers, (aes(x=mass, y=percent_aquatic_ra_logit))) +
  geom_point() + 
  geom_smooth(method=lm , se=TRUE) +
  xlab("Mass") +
  ylab("Percent Aquatic Insects in Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))
# No major outliers or biologically unrealistic values

# How many NAs?
count(is.na(aquatic_adF_nooutliers$mass)) # None!

ggsave(here("3_r_scripts/mom_mass_mom_diet.png"), p_mass, width = 7, height = 7, device = "png")

#############
# PLOT: wing + percent aquatic
#############

p_wing <- ggplot(aquatic_adF_nooutliers, (aes(x=Flat_Wing, y=percent_aquatic_ra))) +
  geom_point() + 
  geom_smooth(method=lm , se=TRUE) +
  xlab("Female flatwing") +
  ylab("Percent Aquatic Insects in Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))
# No major outliers or biologically unrealistic values

# How many NAs?
count(is.na(aquatic_adF_nooutliers$Flat_Wing)) # None!

ggsave(here("3_r_scripts/mom_wing_mom_diet.png"), p_wing, width = 7, height = 7, device = "png")

#############
# PLOT: experimental treatment + percent aquatic
#############

p_exp <- ggplot(aquatic_adF_nooutliers) +
  geom_boxplot((aes(x=exp_treat, y=percent_aquatic_ra))) + 
  xlab("Treatment") +
  ylab("Percent Aquatic Insects in Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

#############
# PLOT: site + percent aquatic
#############

p_site <- ggplot(aquatic_adF_nooutliers) +
  geom_boxplot((aes(x=site, y=percent_aquatic_ra_logit))) + 
  xlab("Site") +
  ylab("Percent Aquatic Insects in Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

ggsave(here("3_r_scripts/mom_site_mom_diet.png"), p_site, width = 7, height = 7, device = "png")

#############
# PLOT: brood size + percent aquatic
#############

p_brood <- ggplot(aquatic_adF_nooutliers) +
  geom_point((aes(x=Brood_Size_Day6, y=percent_aquatic_ra_logit))) + 
  xlab("Site") +
  ylab("Percent Aquatic Insects in Diet") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

#############
# MODEL SELECTION
#############

## Read in list of models to be fit
mods <- read.csv(here("3_r_scripts/model_lists/mod_list_q2c.csv"))
mods <- as.character(mods[,2])

## Loop through model list and fit each model, saving fit model as object in a list
mod.output <- list()
for(i in 1:length(mods)){
  model <- paste("fit.mod",i,sep="_")
  mod <- lmer(as.formula(mods[i]), data = aquatic_adF_nooutliers, REML = FALSE)
  assign(model[1],mod)
  mod.output[[i]] <- get(model)
}

# Create AIC table based on the loop of models above
AICtab <- as.data.frame(mods)
for(i in 1:nrow(AICtab)){
  AICtab$AICc[i] <- AICc(mod.output[[i]])
  AICtab$LL[i] <- logLik(mod.output[[i]])
  AICtab$K[i] <- attr(logLik(mod.output[[i]]),"df")
}
minim <- min(AICtab$AICc)
AICtab$D.AICc <- AICtab$AICc - minim
AICtab$rellik <- exp(-.5* AICtab$D.AICc)
sumAICs <- sum(AICtab$rellik)
for(i in 1:nrow(AICtab)){
  AICtab$Wi[i] <- round(AICtab$rellik[i]/sumAICs,2)
}
AICtab <- AICtab[order(-AICtab$rellik),]

# save the AIC table to a csv
dev.off()
write.csv(AICtab, here("3_r_scripts/model_outputs/AIC_Output_q2c.csv"))

### Run Model 1
q2c_mod1 <- lmer(percent_aquatic_ra_logit ~ mass + site + exp_treat + (1|site_box_year),
            data = aquatic_adF_nooutliers, REML = FALSE)

# Checking conditions
hist(resid(q2c_mod1))
plot(q2c_mod1)
tab_model(q2c_mod1, file = here("3_r_scripts/model_outputs/q2c_mod1_tab.html"))

### Run Model 2
q2c_mod2 <- lmer(percent_aquatic_ra_logit ~ Flat_Wing + site + exp_treat + (1|site_box_year),
                 data = aquatic_adF_nooutliers, REML = FALSE)

# Checking conditions
hist(resid(q2c_mod2))
plot(q2c_mod2)
tab_model(q2c_mod2, file = here("3_r_scripts/model_outputs/q2c_mod2_tab.html"))

################################################################################
# Question 3A: What are the consequences of variation in diet?
# Does aquatic insect content affect nestling morphology?
################################################################################

aquatic_n <- aquatic[aquatic$age == "Nestling" ,]
aquatic_n <- aquatic_n[aquatic_n$age.1 == "12" | aquatic_n$age.1 == "15" ,]
aquatic_n_day12 <- aquatic_n[aquatic_n$age.1 == "12" ,]

q3 <- lmer(mass ~ percent_aquatic_ra*site + (1|site_box_year), data = aquatic_n_day12, REML = FALSE)

tab_model(q3, file = here("3_r_scripts/model_outputs/q3_mod_tab.html"))

p <- ggplot(aquatic_n_day12) +
  geom_point((aes(x=percent_aquatic_ra, y=mass))) + 
  xlab("Nestling Percent Aquatic") +
  ylab("Mass") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

ggsave(here("3_r_scripts/nd12_aq_mass.png"), p, width = 7, height = 7, device = "png")

p <- ggplot(aquatic_n) +
  geom_boxplot((aes(x=fate, y=percent_aquatic_ra))) + 
  xlab("Fate") +
  ylab("Nestling Percent Aquatic") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

ggsave(here("3_r_scripts/nd12_aq_mass.png"), p, width = 7, height = 7, device = "png")

################################################################################
# Question 3B: What are the consequences of variation in diet?
# Does aquatic insect content of nestling predict fledge success?
################################################################################

p <- ggplot(aquatic_n_day12) +
  geom_boxplot((aes(x=fate, y=percent_aquatic_ra))) + 
  xlab("Nestling Percent Aquatic") +
  ylab("Fate") +
  theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
  theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14))

################################################################################
# Question 3C: What are the consequences of variation in diet?
# Does aquatic insect content of adult predict fledge success?
################################################################################

