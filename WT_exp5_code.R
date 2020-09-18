####### Evaluations of Empathizers Depend on the Target of Empathy ########
#######                       Experiment 5                         ########
#######                    Data Analysis Code                      ########
######################  Written by Y. Andre Wang  #########################



# Prepare Packages --------------------------------------------------------

# Name the packages needed
packages <- c("ez", "apaTables", "tidyverse", "effsize", "psych", "rstatix", 
              "lsr", "lavaan", "emmeans", "semTools", "JSmediation")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == F)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = T))



# Load and prepare data ---------------------------------------------------

dat <- read.csv("WT_exp5_data.csv", stringsAsFactors = T)

# Recode factor levels of response type condition
dat$ResponseType <- factor(dat$ResponseType, levels(dat$ResponseType)[c(2, 1)])



# Manipulation checks -----------------------------------------------------

# *- Manipulation check on target valence ---------------------------------

t.test(dat$imp_tar, mu = 4)
psych::d.ci(lsr::cohensD(dat$imp_tar, mu = 4), n1 = length(dat$imp_tar))
round(sd(dat$imp_tar), 2)



# *- Manipulation check on empathy ----------------------------------------

t.test(empathize ~ ResponseType, data = dat)
effsize::cohen.d(dat$empathize, dat$ResponseType)

# Mean and SD
dat %>%
  group_by(ResponseType) %>%
  summarise(mean(empathize, na.rm = T), sd(empathize, na.rm = T))



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



# t-test on respect/liking ------------------------------------------------

t.test(imp1 ~ ResponseType, data = dat)
effsize::cohen.d(dat$imp1, dat$ResponseType)

# Means and SDs
dat %>%
  group_by(ResponseType) %>%
  summarise(mean(imp1), sd(imp1))



# t-test on warmth --------------------------------------------------------

t.test(imp2 ~ ResponseType, data = dat)
effsize::cohen.d(dat$imp2, dat$ResponseType)

# Means and SDs
dat %>%
  group_by(ResponseType) %>%
  summarise(mean(imp2), sd(imp2))



# For internal meta-analysis: Covariances between DVs ---------------------

dat %>%
  filter(ResponseType == "Empathic") %>%
  dplyr::select(imp1, imp2) %>% cov(use = "complete.obs")

dat %>%
  filter(ResponseType == "Condemning") %>%
  dplyr::select(imp1, imp2) %>% cov(use = "complete.obs")
