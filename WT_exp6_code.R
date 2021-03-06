####### Evaluations of Empathizers Depend on the Target of Empathy ########
#######                       Experiment 6                         ########
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

dat <- read.csv("WT_exp6_data.csv", stringsAsFactors = T)

# Recode factor levels of response type condition
dat$ResponseType <- factor(dat$ResponseType, levels(dat$ResponseType)[c(2, 1)])



# Manipulation checks -----------------------------------------------------

# *- Manipulation check on target valence ---------------------------------

t.test(dat$imp_tar, mu = 4)
psych::d.ci(lsr::cohensD(dat$imp_tar, mu = 4), n1 = length(dat$imp_tar))
round(sd(dat$imp_tar), 2)



# *- Manipulation check on empathy ----------------------------------------

t.test(empathize ~ ResponseType, data = dat)
effsize::cohen.d(dat$empathize, dat$ResponseType, na.rm = T)

# Mean and SD
dat %>%
  group_by(ResponseType) %>%
  summarise(mean(empathize, na.rm = T), round(sd(empathize, na.rm = T), 2))



# *- Manipulation check on target affect ----------------------------------

t.test(dat$tar_feel, mu = 4)
round(psych::d.ci(lsr::cohensD(dat$tar_feel, mu = 4), 
                  n1 = length(dat$tar_feel)), 2)

# SD
round(sd(dat$tar_feel), 2)



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
                    between = .(CharGender, ResponseType), type = 3)
aov_imp1_dat <- as.data.frame(aov_imp1$ANOVA)



# *--- 2 x 2 ANOVA --------------------------------------------------------

apa.ezANOVA.table(aov_imp1)

# 90% CIs for character gender
get.ci.partial.eta.squared(
  aov_imp1_dat$F[1], aov_imp1_dat$DFn[1], aov_imp1_dat$DFd[1])

# 90% CIs for response type
get.ci.partial.eta.squared(
  aov_imp1_dat$F[2], aov_imp1_dat$DFn[2], aov_imp1_dat$DFd[2])

# 90% CIs for interaction
get.ci.partial.eta.squared(
  aov_imp1_dat$F[3], aov_imp1_dat$DFn[3], aov_imp1_dat$DFd[3])



# *--- Descriptives -------------------------------------------------------

# Mean and SD by character gender
dat %>%
  group_by(CharGender) %>%
  summarise(mean(imp1), sd(imp1))

# Mean and SD by response type
dat %>%
  group_by(ResponseType) %>%
  summarise(mean(imp1), sd(imp1))

# Mean and SD per cell
dat %>%
  group_by(CharGender, ResponseType) %>%
  summarise(mean(imp1), sd(imp1))



# *- ANOVA on warmth ------------------------------------------------------

aov_imp2 <- ezANOVA(dat, dv = imp2, wid = ID, 
                    between = .(CharGender, ResponseType), type = 3)
aov_imp2_dat <- as.data.frame(aov_imp2$ANOVA)



# *--- 2 x 2 ANOVA --------------------------------------------------------

apa.ezANOVA.table(aov_imp2)

# 90% CIs for character gender
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
  dat %>% group_by(CharGender) %>%
  anova_test(imp2 ~ ResponseType, 
             error = lm(imp2 ~ CharGender * ResponseType, data = dat)) %>%
  as.data.frame()

sme_imp2_RT

# 90% CI for female target
get.ci.partial.eta.squared(
  sme_imp2_RT$F[1], sme_imp2_RT$DFn[1], sme_imp2_RT$DFd[1])

# 90% CI for male target
get.ci.partial.eta.squared(
  sme_imp2_RT$F[2], sme_imp2_RT$DFn[2], sme_imp2_RT$DFd[2])



# *--- Descriptives -------------------------------------------------------

# Mean and SD by character gender
dat %>%
  group_by(CharGender) %>%
  summarise(mean(imp2), sd(imp2))

# Mean and SD by response type
dat %>%
  group_by(ResponseType) %>%
  summarise(mean(imp2), sd(imp2))

# Mean and SD per cell
dat %>%
  group_by(CharGender, ResponseType) %>%
  summarise(mean(imp2), sd(imp2))



# For internal meta-analysis: Covariances between DVs ---------------------

dat %>%
  filter(ResponseType == "Empathic") %>%
  select(imp1, imp2) %>% cov(use = "complete.obs")

dat %>%
  filter(ResponseType == "Condemning") %>%
  select(imp1, imp2) %>% cov(use = "complete.obs")
