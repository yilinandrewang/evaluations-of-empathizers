####### Evaluations of Empathizers Depend on the Target of Empathy ########
#######                       Experiment 4                         ########
#######                    Data Analysis Code                      ########
######################  Written by Y. Andre Wang  #########################



# Prepare Packages --------------------------------------------------------

# Name the packages needed
packages <- c("ez", "apaTables", "tidyverse", "effsize", "psych", "rstatix", 
              "lsr", "lavaan", "emmeans", "semTools")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == F)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = T))



# Load and prepare data ---------------------------------------------------

dat <- read.csv("WT_exp4_data.csv", stringsAsFactors = T)

# Recode factor levels of target valence condition
dat$TargetValence <- factor(dat$TargetValence, 
                            levels(dat$TargetValence)[c(2, 1)])

# Recode factor levels of response type condition
dat$ResponseType <- factor(dat$ResponseType, 
                           levels(dat$ResponseType)[c(1, 3, 2)])



# Manipulation checks -----------------------------------------------------

# *- Manipulation check on target valence ---------------------------------

t.test(imp_tar ~ TargetValence, data = dat, alternative = "greater")
effsize::cohen.d(dat$imp_tar, dat$TargetValence, na.rm = T)

# Mean and SD
dat %>%
  group_by(TargetValence) %>%
  summarise(mean(imp_tar, na.rm = T), sd(imp_tar, na.rm = T))



# *- Manipulation check on empathy ----------------------------------------

# Empathic vs. Positive nonempathic response
t.test(dat$empathize[dat$ResponseType == "Empathic"], 
       dat$empathize[dat$ResponseType == "Positive Nonempathic"])
effsize::cohen.d(dat$empathize[dat$ResponseType == "Empathic"], 
                 dat$empathize[dat$ResponseType == "Positive Nonempathic"])

# Empathic vs. Neutral nonempathic response
t.test(dat$empathize[dat$ResponseType == "Empathic"], 
       dat$empathize[dat$ResponseType == "Neutral Nonempathic"])
effsize::cohen.d(dat$empathize[dat$ResponseType == "Empathic"], 
                 dat$empathize[dat$ResponseType == "Neutral Nonempathic"])

# Mean and SD
dat %>%
  group_by(ResponseType) %>%
  summarise(mean(empathize, na.rm = T), sd(empathize, na.rm = T))



# *- Manipulation check on target affect ----------------------------------

t.test(dat$tar_feel, mu = 4, alternative = "less")
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



# ANOVA on respect/liking -------------------------------------------------

# ANOVA table as dataframe
aov_imp1 <- ezANOVA(dat, dv = imp1, wid = ID, 
                    between = .(TargetValence, ResponseType), type = 3)
aov_imp1_dat <- as.data.frame(aov_imp1$ANOVA)

# Get aov object
aov_imp1_obj <- ezANOVA(dat, dv = imp1, wid = ID, 
                        between = .(TargetValence, ResponseType), type = 3,
                        return_aov = T)

# Marginal means
emm_imp1 <- emmeans(aov_imp1_obj$aov, ~ TargetValence * ResponseType)



# *- 3 x 2 ANOVA ----------------------------------------------------------

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



# *- Descriptives ---------------------------------------------------------

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

# Cell size
dat %>%
  group_by(TargetValence, ResponseType) %>%
  count()


# *- Planned contrasts: Positive target -----------------------------------

# Tests
contrast(emm_imp1, list(
  PosTar_Empathic_PositiveNonempathic = c(1, 0, -1, 0, 0, 0),
  PosTar_Empathic_NeutralNonempathic = c(1, 0, 0, 0, -1, 0),
  PosTar_PositiveNonempathic_NeutralNonempathic = c(0, 0, 1, 0, -1, 0)))

# Calculate ds and CIs for contrast between empathic and positive nonempathic
round(d.ci(t2d(5.518, n1 = 139, n2 = 129), n1 = 139, n2 = 129), 2)

# Calculate ds and CIs for contrast between empathic and neutral nonempathic
round(d.ci(t2d(8.762, n1 = 139, n2 = 111), n1 = 139, n2 = 111), 2)

# Calculate ds and CIs for contrast between positive and neutral nonempathic
round(d.ci(t2d(3.404, n1 = 129, n2 = 111), n1 = 129, n2 = 111), 2)



# *- Post hoc pairwise comparisons: Negative target -----------------------

# Tests
contrast(emm_imp1, list(
  NegTar_Empathic_PositiveNonempathic = c(0, 1, 0, -1, 0, 0), 
  NegTar_Empathic_NeutralNonempathic = c(0, 1, 0, 0, 0, -1),
  NegTar_PositiveNonempathic_NeutralNonempathic = c(0, 0, 0, 1, 0, -1)), 
  adjust = "bonferroni")

# Calculate ds and CIs for contrast between positive and neutral nonempathic
round(d.ci(t2d(2.551, n1 = 123, n2 = 128), n1 = 123, n2 = 128, 
           alpha = .05/3), 2)



# ANOVA on warmth ---------------------------------------------------------

# ANOVA table as dataframe
aov_imp2 <- ezANOVA(dat, dv = imp2, wid = ID, 
                    between = .(TargetValence, ResponseType), type = 3)
aov_imp2_dat <- as.data.frame(aov_imp2$ANOVA)

# Get aov object
aov_imp2_obj <- ezANOVA(dat, dv = imp2, wid = ID, 
                        between = .(TargetValence, ResponseType), type = 3,
                        return_aov = T)

# Marginal means
emm_imp2 <- emmeans(aov_imp2_obj$aov, ~ TargetValence * ResponseType)



# *- 3 x 2 ANOVA ----------------------------------------------------------

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



# *- Descriptives ---------------------------------------------------------

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



# *- Planned contrasts: Positive target -----------------------------------

# Tests
contrast(emm_imp2, list(
  PosTar_Empathic_PositiveNonempathic = c(1, 0, -1, 0, 0, 0),
  PosTar_Empathic_NeutralNonempathic = c(1, 0, 0, 0, -1, 0),
  PosTar_PositiveNonempathic_NeutralNonempathic = c(0, 0, 1, 0, -1, 0)))

# Calculate ds and CIs for contrast between empathic and positive nonempathic
round(d.ci(t2d(7.089, n1 = 139, n2 = 129), n1 = 139, n2 = 129), 2)

# Calculate ds and CIs for contrast between empathic and neutral nonempathic
round(d.ci(t2d(12.864, n1 = 139, n2 = 111), n1 = 139, n2 = 111), 2)

# Calculate ds and CIs for contrast between positive and neutral nonempathic
round(d.ci(t2d(5.953, n1 = 129, n2 = 111), n1 = 129, n2 = 111), 2)



# *- Post hoc pairwise comparisons: Negative target -----------------------

# Tests
contrast(emm_imp2, list(
  NegTar_Empathic_PositiveNonempathic = c(0, 1, 0, -1, 0, 0), 
  NegTar_Empathic_NeutralNonempathic = c(0, 1, 0, 0, 0, -1),
  NegTar_PositiveNonempathic_NeutralNonempathic = c(0, 0, 0, 1, 0, -1)), 
  adjust = "bonferroni")

# Calculate ds and CIs for contrast between empathic and positive nonempathic
round(d.ci(t2d(0.533, n1 = 110, n2 = 123), n1 = 110, n2 = 123, 
           alpha = .05/3), 2)

# Calculate ds and CIs for contrast between empathic and neutral nonempathic
round(d.ci(t2d(5.520, n1 = 110, n2 = 128), n1 = 110, n2 = 128, 
           alpha = .05/3), 2)

# Calculate ds and CIs for contrast between empathic and neutral nonempathic
round(d.ci(t2d(6.238, n1 = 123, n2 = 128), n1 = 123, n2 = 128, 
           alpha = .05/3), 2)



# Moderated mediation -----------------------------------------------------


# *- Data preparation -----------------------------------------------------

# Effect-code target valence (Z)
dat$Z <- recode(dat$TargetValence, 
                "Negative Target" = -1, "Positive Target" = 1)

# Effect-code response type by empathy (X_e)
dat$X_e <- recode(dat$ResponseType, 
                  "Empathic" = 2/3, "Positive Nonempathic" = -1/3,
                  "Neutral Nonempathic" = -1/3)

# Effect-code response type by positivity/valence (X_v)
dat$X_v <- recode(dat$ResponseType, 
                  "Empathic" = 1/3, "Positive Nonempathic" = 1/3,
                  "Neutral Nonempathic" = -2/3)

# Effect-code response type by empathic vs. positive nonempathic (X_p)
dat$X_p <- recode(dat$ResponseType, 
                  "Empathic" = 1, "Positive Nonempathic" = -1)

# Create interaction term of X and Z
dat$XZ_e <- dat$X_e * dat$Z
dat$XZ_v <- dat$X_v * dat$Z
dat$XZ_p <- dat$X_p * dat$Z

# Reliability of the mediator (M)
psych::alpha(dat[, c("res_like", "res_pos", "res_unfav_r")])

# Create product indicators of MZ
dat <- indProd(dat, var1 = "Z", var2 = c("res_like", "res_pos", "res_unfav_r"), 
               match = F, meanC = T, residualC = T, doubleMC = F,
               namesProd = c("mz1", "mz2", "mz3"))

# Define parameters for first-stage and second-stage mediation to get 95% CIs
med1 <- 'aMod*b'
med2 <- 'a*bMod'
abneg <- '(a + aMod*(-1))*(b + bMod*(-1))'
abpos <- '(a + aMod*(+1))*(b + bMod*(+1))'

# Subset data for Models 5 and 6
dat_p <- subset(dat, is.na(X_p) == F)



# *- Model 1: Respect/liking, Empathy -------------------------------------

mod_f1_e <- "
M =~ p*res_like + p*res_pos + res_unfav_r
Y =~ liking + respect + trust + friend
MZ =~ q*mz1 + q*mz2 + mz3

M ~ a*X_e + Z + aMod*XZ_e
Y ~ cp*X_e + Z + cpMod*XZ_e + b*M + bMod*MZ

c := cp + (a*b)
cMod := cpMod + (aMod * bMod)
aModb := aMod*b
abMod := a*bMod

aneg := a + aMod*(-1)
apos := a + aMod*(+1)
bneg := b + bMod*(-1)
bpos := b + bMod*(+1)
abneg := aneg * bneg
abpos := apos * bpos

res_like ~~ res_pos
mz1 ~~ mz2
"

# Fit model
fit_f1_e <- sem(model = mod_f1_e, data = dat, std.lv = T)

# Summary
summary(fit_f1_e, standardized = T, fit.measures = T, ci = T)

# Fit indices
round(fitMeasures(fit_f1_e)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_f1_e)["pvalue"], 3)

# Get 95% CI for first-stage moderated mediation (aMod*b)
set.seed(42); monteCarloMed(med1, object = fit_f1_e)

# Get 95% CI for second-stage moderated mediation (a*bMod)
set.seed(42); monteCarloMed(med2, object = fit_f1_e)

# Get 95% CI for indirect effect for positive target (apos*bpos)
set.seed(42); monteCarloMed(abpos, object = fit_f1_e)

# Get 95% CI for indirect effect for negative target (aneg*bneg)
set.seed(42); monteCarloMed(abneg, object = fit_f1_e)



# *- Model 2: Warmth, Empathy -------------------------------------

mod_f2_e <- "
M =~ p*res_like + p*res_pos + res_unfav_r
Y =~ understanding + kind + cold_r + caring
MZ =~ q*mz1 + q*mz2 + mz3

M ~ a*X_e + Z + aMod*XZ_e
Y ~ cp*X_e + Z + cpMod*XZ_e + b*M + bMod*MZ

c := cp + (a*b)
cMod := cpMod + (aMod * bMod)
aModb := aMod*b
abMod := a*bMod

aneg := a + aMod*(-1)
apos := a + aMod*(+1)
bneg := b + bMod*(-1)
bpos := b + bMod*(+1)
abneg := aneg * bneg
abpos := apos * bpos

res_like ~~ res_pos
mz1 ~~ mz2
"

# Fit model
fit_f2_e <- sem(model = mod_f2_e, data = dat, std.lv = T)

# Summary
summary(fit_f2_e, standardized = T, fit.measures = T, ci = T)

# Fit indices
round(fitMeasures(fit_f2_e)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_f2_e)["pvalue"], 3)

# Get 95% CI for first-stage moderated mediation (aMod*b)
set.seed(42); monteCarloMed(med1, object = fit_f2_e)

# Get 95% CI for second-stage moderated mediation (a*bMod)
set.seed(42); monteCarloMed(med2, object = fit_f2_e)

# Get 95% CI for indirect effect for positive target (apos*bpos)
set.seed(42); monteCarloMed(abpos, object = fit_f2_e)

# Get 95% CI for indirect effect for negative target (aneg*bneg)
set.seed(42); monteCarloMed(abneg, object = fit_f2_e)



# *- Model 3: Respect/liking, Positivity ----------------------------------

mod_f1_v <- "
M =~ p*res_like + p*res_pos + res_unfav_r
Y =~ liking + respect + trust + friend
MZ =~ q*mz1 + q*mz2 + mz3

M ~ a*X_v + Z + aMod*XZ_v
Y ~ cp*X_v + Z + cpMod*XZ_v + b*M + bMod*MZ

c := cp + (a*b)
cMod := cpMod + (aMod * bMod)
aModb := aMod*b
abMod := a*bMod

aneg := a + aMod*(-1)
apos := a + aMod*(+1)
bneg := b + bMod*(-1)
bpos := b + bMod*(+1)
abneg := aneg * bneg
abpos := apos * bpos

res_like ~~ res_pos
mz1 ~~ mz2
"

# Fit model
fit_f1_v <- sem(model = mod_f1_v, data = dat, std.lv = T)

# Summary
summary(fit_f1_v, standardized = T, fit.measures = T, ci = T)

# Fit indices
round(fitMeasures(fit_f1_v)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_f1_v)["pvalue"], 3)

# Get 95% CI for first-stage moderated mediation (aMod*b)
set.seed(42); monteCarloMed(med1, object = fit_f1_v)

# Get 95% CI for second-stage moderated mediation (a*bMod)
set.seed(42); monteCarloMed(med2, object = fit_f1_v)

# Get 95% CI for indirect effect for positive target (apos*bpos)
set.seed(42); monteCarloMed(abpos, object = fit_f1_v)

# Get 95% CI for indirect effect for negative target (aneg*bneg)
set.seed(42); monteCarloMed(abneg, object = fit_f1_v)



# *- Model 4: Warmth, Positivity ------------------------------------------

mod_f2_v <- "
M =~ p*res_like + p*res_pos + res_unfav_r
Y =~ understanding + kind + cold_r + caring
MZ =~ q*mz1 + q*mz2 + mz3

M ~ a*X_v + Z + aMod*XZ_v
Y ~ cp*X_v + Z + cpMod*XZ_v + b*M + bMod*MZ

c := cp + (a*b)
cMod := cpMod + (aMod * bMod)
aModb := aMod*b
abMod := a*bMod

aneg := a + aMod*(-1)
apos := a + aMod*(+1)
bneg := b + bMod*(-1)
bpos := b + bMod*(+1)
abneg := aneg * bneg
abpos := apos * bpos

res_like ~~ res_pos
mz1 ~~ mz2
"

# Fit model
fit_f2_v <- sem(model = mod_f2_v, data = dat, std.lv = T)

# Summary
summary(fit_f2_v, standardized = T, fit.measures = T, ci = T)

# Fit indices
round(fitMeasures(fit_f2_v)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_f2_v)["pvalue"], 3)

# Get 95% CI for first-stage moderated mediation (aMod*b)
set.seed(42); monteCarloMed(med1, object = fit_f2_v)

# Get 95% CI for second-stage moderated mediation (a*bMod)
set.seed(42); monteCarloMed(med2, object = fit_f2_v)

# Get 95% CI for indirect effect for positive target (apos*bpos)
set.seed(42); monteCarloMed(abpos, object = fit_f2_v)

# Get 95% CI for indirect effect for negative target (aneg*bneg)
set.seed(42); monteCarloMed(abneg, object = fit_f2_v)



# *- Model 5: Respect/liking, Empathic vs. Positive Nonempathic -----------

mod_f1_p <- "
M =~ p*res_like + p*res_pos + res_unfav_r
Y =~ liking + respect + trust + friend
MZ =~ q*mz1 + q*mz2 + mz3

M ~ a*X_p + Z + aMod*XZ_p
Y ~ cp*X_p + Z + cpMod*XZ_p + b*M + bMod*MZ

c := cp + (a*b)
cMod := cpMod + (aMod * bMod)
aModb := aMod*b
abMod := a*bMod

aneg := a + aMod*(-1)
apos := a + aMod*(+1)
bneg := b + bMod*(-1)
bpos := b + bMod*(+1)
abneg := aneg * bneg
abpos := apos * bpos

res_like ~~ res_pos
mz1 ~~ mz2
"

# Fit model
fit_f1_p <- sem(model = mod_f1_p, data = dat, std.lv = T)

# Summary
summary(fit_f1_p, standardized = T, fit.measures = T, ci = T)

# Fit indices
round(fitMeasures(fit_f1_p)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_f1_p)["pvalue"], 3)

# Get 95% CI for first-stage moderated mediation (aMod*b)
set.seed(42); monteCarloMed(med1, object = fit_f1_p)

# Get 95% CI for second-stage moderated mediation (a*bMod)
set.seed(42); monteCarloMed(med2, object = fit_f1_p)

# Get 95% CI for indirect effect for positive target (apos*bpos)
set.seed(42); monteCarloMed(abpos, object = fit_f1_p)

# Get 95% CI for indirect effect for negative target (aneg*bneg)
set.seed(42); monteCarloMed(abneg, object = fit_f1_p)



# *- Model 6: Warmth, Empathic vs. Positive Nonempathic -------------------

mod_f2_p <- "
M =~ p*res_like + p*res_pos + res_unfav_r
Y =~ understanding + kind + cold_r + caring
MZ =~ q*mz1 + q*mz2 + mz3

M ~ a*X_p + Z + aMod*XZ_p
Y ~ cp*X_p + Z + cpMod*XZ_p + b*M + bMod*MZ

c := cp + (a*b)
cMod := cpMod + (aMod * bMod)
aModb := aMod*b
abMod := a*bMod

aneg := a + aMod*(-1)
apos := a + aMod*(+1)
bneg := b + bMod*(-1)
bpos := b + bMod*(+1)
abneg := aneg * bneg
abpos := apos * bpos

res_like ~~ res_pos
mz1 ~~ mz2
"

# Fit model
fit_f2_p <- sem(model = mod_f2_p, data = dat, std.lv = T)

# Summary
summary(fit_f2_p, standardized = T, fit.measures = T, ci = T)

# Fit indices
round(fitMeasures(fit_f2_p)[c("chisq", "df", "cfi", "tli", "rmsea")], 2)
round(fitMeasures(fit_f2_p)["pvalue"], 3)

# Get 95% CI for first-stage moderated mediation (aMod*b)
set.seed(42); monteCarloMed(med1, object = fit_f2_p)

# Get 95% CI for second-stage moderated mediation (a*bMod)
set.seed(42); monteCarloMed(med2, object = fit_f2_p)

# Get 95% CI for indirect effect for positive target (apos*bpos)
set.seed(42); monteCarloMed(abpos, object = fit_f2_p)

# Get 95% CI for indirect effect for negative target (aneg*bneg)
set.seed(42); monteCarloMed(abneg, object = fit_f2_p)



# Exploratory analyses ----------------------------------------------------


# *- ANOVA on similarity --------------------------------------------------

apa.ezANOVA.table(ezANOVA(dat, dv = similar, wid = ID, 
                          between = .(TargetValence, ResponseType), type = 3))



# For internal meta-analysis: Covariances between DVs ---------------------

dat %>%
  filter(ResponseType == "Empathic" & 
           TargetValence == "Positive Target") %>%
  select(imp1, imp2) %>% cov(use = "complete.obs")

dat %>%
  filter(ResponseType == "Empathic" & 
           TargetValence == "Negative Target") %>%
  select(imp1, imp2) %>% cov(use = "complete.obs")

dat %>%
  filter(ResponseType == "Neutral Nonempathic" & 
           TargetValence == "Positive Target") %>%
  select(imp1, imp2) %>% cov(use = "complete.obs")

dat %>%
  filter(ResponseType == "Neutral Nonempathic" & 
           TargetValence == "Negative Target") %>%
  select(imp1, imp2) %>% cov(use = "complete.obs")

dat %>%
  filter(ResponseType == "Positive Nonempathic" & 
           TargetValence == "Positive Target") %>%
  select(imp1, imp2) %>% cov(use = "complete.obs")

dat %>%
  filter(ResponseType == "Positive Nonempathic" & 
           TargetValence == "Negative Target") %>%
  select(imp1, imp2) %>% cov(use = "complete.obs")
