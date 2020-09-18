####### Evaluations of Empathizers Depend on the Target of Empathy ########
#######                       Experiment 1                         ########
#######                    Data Analysis Code                      ########
######################  Written by Y. Andre Wang  #########################



# Prepare Packages --------------------------------------------------------

# Name the packages needed
packages <- c("ez", "apaTables", "tidyverse", "effsize", "psych", "rstatix", 
              "lsr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == F)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = T))



# Load and prepare data ---------------------------------------------------

dat <- read.csv("WT_exp1_data.csv", stringsAsFactors = T)

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



# Exploratory factor analyses ---------------------------------------------

# Create dataframe for DVs only
s1_DV <- dat[, c(6:11, 13:14)]


# *- Two-factor solution --------------------------------------------------


efa2 <- factanal(~ like + respect + trust + friend + understanding + kind + 
                   cold_r + caring, factors = 2, 
                 covmat = cov(s1_DV, use = "pairwise.complete.obs"), 
                 n.obs = dim(s1_DV)[1], rotation = "promax")
efa2



# *- Solutions with one factor or 3+ factors ------------------------------

efa1 <- factanal(~ like + respect + trust + friend + understanding + kind + 
                   cold_r + caring, factors = 1, 
                 covmat = cov(s1_DV, use = "pairwise.complete.obs"), 
                 n.obs = dim(s1_DV)[1], rotation = "promax")
efa1

efa3 <- factanal(~ like + respect + trust + friend + understanding + kind + 
                   cold_r + caring, factors = 3, 
                 covmat = cov(s1_DV, use = "pairwise.complete.obs"), 
                 n.obs = dim(s1_DV)[1], rotation = "promax")
efa3

efa4 <- factanal(~ like + respect + trust + friend + understanding + kind + 
                   cold_r + caring, factors = 4, 
                 covmat = cov(s1_DV, use = "pairwise.complete.obs"), 
                 n.obs = dim(s1_DV)[1], rotation = "promax")
efa4



# *- Create composite scores ----------------------------------------------

dat$imp1 <- apply(dat[, c("liking", "respect", "trust", "friend")], 
                 1, mean, na.rm = T)
dat$imp2 <- apply(dat[, c("understanding", "kind", "cold_r", "caring")], 
                 1, mean, na.rm = T)

# Reliabilities and correlation
psych::alpha(s1_DV[ , c(1:4)]); psych::alpha(s1_DV[ , c(5:8)])
cor.test(dat$imp1, dat$imp2)



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


# *- ANOVA on similarity --------------------------------------------------

apa.ezANOVA.table(ezANOVA(dat, dv = similar, wid = ID, 
                          between = .(TargetValence, ResponseType), type = 3))



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
