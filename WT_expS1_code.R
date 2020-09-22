####### Evaluations of Empathizers Depend on the Target of Empathy ########
#######                       Experiment S1                         #######
#######                    Data Analysis Code                      ########
######################  Written by Y. Andre Wang  #########################



# Prepare Packages --------------------------------------------------------

# Name the packages needed
packages <- c("ez", "apaTables", "tidyverse", "effsize", "psych", "rstatix", 
              "lsr", "lavaan")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == F)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = T))



# Load and prepare data ---------------------------------------------------

dat <- read.csv("WT_expS1_data.csv", stringsAsFactors = T)

# Recode factor levels of target valence condition
dat$TargetValence <- factor(dat$TargetValence, 
                            levels(dat$TargetValence)[c(2, 1)])



# Manipulation checks -----------------------------------------------------

# *- Manipulation check on target valence ---------------------------------

t.test(imp_tar ~ TargetValence, data = dat)
effsize::cohen.d(dat$imp_tar, dat$TargetValence, na.rm = T)

# Mean and SD
dat %>%
  group_by(TargetValence) %>%
  summarise(mean(imp_tar, na.rm = T), round(sd(imp_tar, na.rm = T), 2))



# *- Manipulation check on empathy ----------------------------------------

t.test(empathize ~ ResponseType, data = dat)
effsize::cohen.d(dat$empathize, dat$ResponseType, na.rm = T)

# Mean and SD
dat %>%
  group_by(ResponseType) %>%
  summarise(mean(empathize, na.rm = T), round(sd(empathize, na.rm = T), 2))



# Confirmatory factor analysis --------------------------------------------

# *- Two-factor solution --------------------------------------------------

# Specify model
mod_2f <- "
imp1 =~ liking + respect + trust + friend
imp2 =~ understanding + kind + cold_r + caring"

# Fit model
fit_2f <- cfa(mod_2f, data = dat, std.lv = T)

# Fit indices
round(fitMeasures(fit_2f)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_2f)["pvalue"], 3)

# Full model output
summary(fit_2f, fit.measures = T, standardized = T)


# *- One-factor solution --------------------------------------------------

# Specify model
mod_1f <- "
imp =~ liking + respect + trust + friend +
understanding + kind + cold_r + caring"

# Fit model
fit_1f <- cfa(mod_1f, data = dat, std.lv = T)

# Fit indices
round(fitMeasures(fit_1f)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_1f)["pvalue"], 3)

# Full model output
summary(fit_1f, fit.measures = T, standardized = T)

# Model comparison with the two-factor solution
anova(fit_2f, fit_1f)



# *- Create composite scores ----------------------------------------------

dat$imp1 <- apply(dat[, c("liking", "respect", "trust", "friend")], 
                 1, mean, na.rm = T)
dat$imp2 <- apply(dat[, c("understanding", "kind", "cold_r", "caring")], 
                 1, mean, na.rm = T)

# Reliabilities
psych::alpha(dat[, c("liking", "respect", "trust", "friend")])
psych::alpha(dat[, c("understanding", "kind", "cold_r", "caring")])



# Key analyses ------------------------------------------------------------


# *- ANOVA on respect/liking ----------------------------------------------

aov_imp1 <- ezANOVA(dat, dv = imp1, wid = ID, 
                    between = .(TargetValence, ResponseType), type = 3)
aov_imp1_dat <- as.data.frame(aov_imp1$ANOVA)



# *--- 2 x 2 ANOVA --------------------------------------------------------

apa.ezANOVA.table(aov_imp1)

# 90% CIs for target valence
get.ci.partial.eta.squared(
  aov_imp1_dat$F[1], aov_imp1_dat$DFn[1], aov_imp1_dat$DFd[1])

# 90% CIs for response type
get.ci.partial.eta.squared(
  aov_imp1_dat$F[2], aov_imp1_dat$DFn[2], aov_imp1_dat$DFd[2])

# 90% CIs for interaction
get.ci.partial.eta.squared(
  aov_imp1_dat$F[3], aov_imp1_dat$DFn[3], aov_imp1_dat$DFd[3])



# *--- Simple main effects of response type -------------------------------

sme_imp1_RT <- 
  dat %>% group_by(TargetValence) %>%
  anova_test(imp1 ~ ResponseType, 
             error = lm(imp1 ~ TargetValence * ResponseType, data = dat)) %>%
  as.data.frame()

sme_imp1_RT

# 90% CI for positive target
get.ci.partial.eta.squared(
  sme_imp1_RT$F[1], sme_imp1_RT$DFn[1], sme_imp1_RT$DFd[1])

# 90% CI for negative target
get.ci.partial.eta.squared(
  sme_imp1_RT$F[2], sme_imp1_RT$DFn[2], sme_imp1_RT$DFd[2])



# *--- Descriptives -------------------------------------------------------

# Mean and SD by target valence
dat %>%
  group_by(TargetValence) %>%
  summarise(mean(imp1), sd(imp1))

# Mean and SD by response type
dat %>%
  group_by(ResponseType) %>%
  summarise(mean(imp1), sd(imp1))

# Mean and SD per cell
dat %>%
  group_by(TargetValence, ResponseType) %>%
  summarise(mean(imp1), sd(imp1))



# *- ANOVA on warmth ------------------------------------------------------

aov_imp2 <- ezANOVA(dat, dv = imp2, wid = ID, 
                    between = .(TargetValence, ResponseType), type = 3)
aov_imp2_dat <- as.data.frame(aov_imp2$ANOVA)



# *--- 2 x 2 ANOVA --------------------------------------------------------

apa.ezANOVA.table(aov_imp2)

# 90% CIs for target valence
get.ci.partial.eta.squared(
  aov_imp2_dat$F[1], aov_imp2_dat$DFn[1], aov_imp2_dat$DFd[1])

# 90% CIs for response type
get.ci.partial.eta.squared(
  aov_imp2_dat$F[2], aov_imp2_dat$DFn[2], aov_imp2_dat$DFd[2])

# 90% CIs for interaction
get.ci.partial.eta.squared(
  aov_imp2_dat$F[3], aov_imp2_dat$DFn[3], aov_imp2_dat$DFd[3])



# *--- Simple main effects of response type -------------------------------

sme_imp2_RT <- 
  dat %>% group_by(TargetValence) %>%
  anova_test(imp2 ~ ResponseType, 
             error = lm(imp2 ~ TargetValence * ResponseType, data = dat)) %>%
  as.data.frame()

sme_imp2_RT

# 90% CI for positive target
get.ci.partial.eta.squared(
  sme_imp2_RT$F[1], sme_imp2_RT$DFn[1], sme_imp2_RT$DFd[1])

# 90% CI for negative target
get.ci.partial.eta.squared(
  sme_imp2_RT$F[2], sme_imp2_RT$DFn[2], sme_imp2_RT$DFd[2])



# *--- Descriptives -------------------------------------------------------

# Mean and SD by target valence
dat %>%
  group_by(TargetValence) %>%
  summarise(mean(imp2), sd(imp2))

# Mean and SD by response type
dat %>%
  group_by(ResponseType) %>%
  summarise(mean(imp2), sd(imp2))

# Mean and SD per cell
dat %>%
  group_by(TargetValence, ResponseType) %>%
  summarise(mean(imp2), sd(imp2))



# Exploratory analyses ----------------------------------------------------


# For internal meta-analysis: Covariances between DVs ---------------------

dat %>%
  filter(ResponseType == "Empathic" & TargetValence == "Positive Target") %>%
  select(imp1, imp2) %>% cov(use = "complete.obs")

dat %>%
  filter(ResponseType == "Empathic" & TargetValence == "Negative Target") %>%
  select(imp1, imp2) %>% cov(use = "complete.obs")

dat %>%
  filter(ResponseType == "Nonempathic" & TargetValence == "Positive Target") %>%
  select(imp1, imp2) %>% cov(use = "complete.obs")

dat %>%
  filter(ResponseType == "Nonempathic" & TargetValence == "Negative Target") %>%
  select(imp1, imp2) %>% cov(use = "complete.obs")
