
# ANALYSIS SCRIPT FOR PAPER: Does the institutional environment shape human life-histories? The
# historical case of Augustus Smith and the Isles of Scilly

# LIFESPAN

# Setup -------------------------------------------------------

## Install/load packages --------------------------------------
packages <- c(
  "knitr", "flextable", "ggplot2", "broom", "devtools", "lme4",
  "gridExtra", "merTools", "mgcv", "tidyverse", "gratia", "ggpubr",
  "glmm", "ggeffects", "gamlss", "performance", "cowplot", "patchwork"
)

# Install any packages that aren't already installed
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if (length(to_install)) install.packages(to_install)

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))


## Function for getting relevant info for comparing smooths in GAM analysis -------------
smooth_diff <- function(model, newdata, f1, f2, var, alpha = 0.05,
                        unconditional = FALSE) {
  xp <- predict(model, newdata = newdata, type = 'lpmatrix')
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other lochs
  X[, ! (c1 | c2)] <- 0
  ## zero out the parametric cols
  X[, !grepl('^s\\(', colnames(xp))] <- 0
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  data.frame(pair = paste(f1, f2, sep = '-'),
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr)
}


## Load data -----------------------------------------------------

#LIFESPAN
deaths <- read.csv("data/deaths.csv")
summary(deaths)
deaths$Gender <- as.factor(deaths$Gender)
deaths$DeathYear.factor <- as.factor(as.character(deaths$DeathYear.factor))
deaths$Group <- as.factor(as.character(deaths$Group))

#Get end of relevant period for gam = end of intervention period + two generations, 
#where the length of a generation is the average age at first reproduction
#Load parents data
parents <- read.csv("data/parents.csv")
time.line <- 1872 + (2 * mean(parents$mother_age_first_repro, na.rm = T)) #25.87



# Lifespans overall ITS --------------------------------------------
## Fit model -------------------------------------------------------
ts.ls <- lmer(DeathAge ~ DeathYear + Group + TimeElapsed  + (1|DeathYear), data = deaths)
#saveRDS(ts.ls, file = "output/lifespans/ts.ls.rds")




## Results table ---------------------------------------------------
tab.ls <- cbind(data.frame(Predictors = c("(Intercept)", "Death year", "Group1 (post)", "Time elapsed")), as.data.frame(coef(summary(ts.ls))), confint(ts.ls)[3:6, ], row.names = NULL)
tab.ls <- tab.ls %>% mutate_if(is.numeric, round, digits = 3)
#saveRDS(tab.ls, file = "output/lifespans/tab.ls.rds")




## Plot -------------------------------------------------------------

#Get predicted values for plot
#Counterfactual line
plot.ls.its.c <- ggpredict(ts.ls, c("DeathYear[1798:2002]", "TimeElapsed[0]", "Group[0]"))
plot.ls.its.c <- plot.ls.its.c %>% rename(DeathYear = x, TimeElapsed = group, Group = facet)
plot.ls.its.c$TimeElapsed <- as.numeric(as.character(plot.ls.its.c$TimeElapsed))
plot.ls.its.c$ITS <- "Counterfactual"


#Predicted line
plot.ls.its.p <- ggpredict(ts.ls, c("DeathYear[1798:2002]","TimeElapsed[0:169]", "Group"))
plot.ls.its.p <- plot.ls.its.p %>% rename(DeathYear = x, TimeElapsed = group, Group = facet)
#Select rows where year < 1834 & group = 0 & TimeElapsed = 0
plot.ls.its.p <- plot.ls.its.p[which(plot.ls.its.p$DeathYear < 1834 & plot.ls.its.p$Group == "0" & plot.ls.its.p$TimeElapsed == 0 | plot.ls.its.p$DeathYear >= 1834 & plot.ls.its.p$Group == "1" & plot.ls.its.p$TimeElapsed == plot.ls.its.p$DeathYear - 1833), ]
plot.ls.its.p$ITS <- "Predicted"


# Plot

its.ls.plot <- ggplot(aes(x = DeathYear, y = DeathAge), data = deaths) +
  # Points
  geom_point(data = deaths, aes(x = DeathYear, y = DeathAge), colour = "grey", shape = 16, alpha = 0.4) +
  # CIs
  geom_ribbon(data = plot.ls.its.p[which(plot.ls.its.p$DeathYear >= 1834),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.4) +
  geom_ribbon(data = plot.ls.its.c[which(plot.ls.its.c$DeathYear > 1834.5),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.4) +
  geom_ribbon(data = plot.ls.its.c[which(plot.ls.its.c$DeathYear <= 1834.5),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.4) +
  # Lines
  geom_line(data = plot.ls.its.p[which(plot.ls.its.p$DeathYear <= 1833),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.ls.its.p[which(plot.ls.its.p$DeathYear >= 1834),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.ls.its.c[which(plot.ls.its.p$DeathYear >= 1833),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  #scale_x_continuous(breaks = sort(c(seq(1700, 2000, length.out=4), 1834))) +
  # Intervention line
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  # Theme and stop y at 100
  theme_minimal() +
  xlab("Year of death") +
  ylab("Age at death") +
  coord_cartesian(ylim=c(0, 100)) +
  # Legend
  scale_linetype_manual(labels = c("Counterfactual", "Predicted values"), values = c("dashed", "solid")) +
  scale_alpha_manual(labels = c("Counterfactual", "Predicted values"), values = c(0.4, 0.4), guide = "none") +
  theme(legend.position = "bottom", legend.title=element_text(face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("All data")

its.ls.plot
#saveRDS(its.ls.plot, file = "output/lifespans/its.ls.plot.rds")






# Lifespans overall GAM --------------------------------------------------
## Fit GAM ---------------------------------------------------------------
gam.ls <- gam(DeathAge ~ s(DeathYear, bs = "tp", k = -1) + s(DeathYear.factor, bs = "re"), gamma = 1.4, data = deaths, method = "REML")
#saveRDS(gam.ls, file = "output/lifespans/gam.ls.rds")


## Plot trend ------------------------------------------------------------------
# Predicted data for plots
p.ls <- ggpredict(gam.ls, c("DeathYear[1798:2002]"))
p.ls$group <- NULL
p.ls <- p.ls %>% rename(DeathYear = x, pfit = predicted, pse = std.error)

#Get derivative
data.ls <- data.frame(DeathYear = 1798:2002)
data.ls$DeathYear.factor <- as.factor(data.ls$DeathYear)
deriv.ls <- derivatives(gam.ls, term = "s(DeathYear)", data = data.ls)
names(deriv.ls)[names(deriv.ls) == "data"] <- "DeathYear"
deriv.ls <- transform(deriv.ls, pfit = p.ls$pfit)
 
#Make separate deriv dataframes for inc. vs dec.
deriv.ls.inc <- deriv.ls
deriv.ls.inc[(deriv.ls.inc$.lower_ci < 0 & deriv.ls.inc$.upper_ci < 0), c(".derivative", "pfit")] <- NA
deriv.ls.dec <- deriv.ls
deriv.ls.dec[!(deriv.ls.dec$.lower_ci < 0 & deriv.ls.dec$.upper_ci < 0), c(".derivative", "pfit")] <- NA



#Plotting data and estimated trend line with periods of significant change indicated
gam.ls.plot <- deaths %>%
  ggplot (aes (x = DeathYear, y = DeathAge)) +
  labs (x = "Year of death", y= "Age at death") +
  theme_minimal () +
# Points
  geom_point(data = deaths, aes(x = DeathYear, y = DeathAge), colour = "grey", shape = 16, alpha = 0.6) +
# Intervention lines
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
# CIs
  geom_ribbon (data = p.ls, aes(x = DeathYear, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.3, fill = "black") +
# Lines
  geom_line(data = p.ls, aes(x = DeathYear, y = pfit), color="black") +
 geom_line(data = deriv.ls.inc, aes(x = DeathYear, y = pfit, color="blue"), size = 1) +
 geom_line(data = deriv.ls.dec, aes(x = DeathYear, y = pfit, color="red"), size = 1) +
# Legend and axis titles
  scale_color_manual(labels = c("Sig. incr.", "Sig. decr."), values = c("blue", "red"), name = "GAM") +
  theme(legend.position = "bottom", legend.title=element_text(face = "bold")) +
  theme(legend.position = "none") +
  xlab(NULL)
gam.ls.plot 
#saveRDS(gam.ls.plot, file = "output/lifespans/gam.ls.plot.rds")





## Plot derivative change --------------------------------------------------
gam.ls.deriv <-  ggplot(data = deriv.ls, aes(x = DeathYear, y = .derivative)) +
  theme_minimal()  +
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  geom_hline(yintercept = 0) +
  geom_line(data = deriv.ls, aes(x = DeathYear, y = .derivative), colour = "black") +
  geom_ribbon(data = deriv.ls, aes(x = DeathYear, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.3) +
  geom_line(data = deriv.ls.inc, aes(x = DeathYear, y = .derivative), colour = "blue", size = 1) +
  geom_line(data = deriv.ls.dec, aes(x = DeathYear, y = .derivative), colour = "red", size = 1) +
  xlab("Year of death") +
  ylab("Derivative") +
  ylim(0, 0.35)
gam.ls.deriv

#saveRDS(gam.ls.deriv, file = "output/lifespans/gam.ls.deriv.rds")









# Lifespans gender ITS ----------------------------------------------------------

## fit model --------------------------------------------------------------------
ts.gender <- lmer(DeathAge ~ DeathYear + Group + TimeElapsed + Gender + Gender:DeathYear + Gender:Group + Gender:TimeElapsed + (1|DeathYear), data = deaths)
#saveRDS(ts.gender, file = "output/lifespans/ts.gender.rds")



## Results table ----------------------------------------------------------------
tab.ts.gender <- cbind(data.frame(Predictors = c("(Intercept)", "DeathYear", "Group1 (post)", "TimeElapsed", "Gender", "Group1:Gender", "DeathYear:Gender", "TimeElapsed:Gender")), as.data.frame(coef(summary(ts.gender))), confint(ts.gender)[3:10, ], row.names = NULL)
tab.ts.gender <- tab.ts.gender %>% mutate_if(is.numeric, round, digits = 3)
#saveRDS(tab.ts.gender, file = "output/lifespans/tab.ts.gender.rds")


## Plot -------------------------------------------------------------------------


#Predict from model - ggpredict
#Counter
plot.ls.its.gen.c <- ggpredict(ts.gender, c("DeathYear[1798:2002]", "Gender", "TimeElapsed[0]", "Group"))
plot.ls.its.gen.c <- plot.ls.its.gen.c %>% rename(DeathYear = x, Gender = group, TimeElapsed = facet, Group = panel)
#Remove rows where Group = 1
plot.ls.its.gen.c <- plot.ls.its.gen.c[which(plot.ls.its.gen.c$Group == 0), ]
plot.ls.its.gen.c$TimeElapsed <- as.numeric(as.character(plot.ls.its.gen.c$TimeElapsed))
plot.ls.its.gen.c$ITS <- "Counterfactual"


#Predicted
plot.ls.its.gen.p <- ggpredict(ts.gender, c("DeathYear[1798:2002]", "Gender", "TimeElapsed[0:169]", "Group"))
plot.ls.its.gen.p <- plot.ls.its.gen.p %>% rename(DeathYear = x, Gender = group, TimeElapsed = facet, Group = panel)
#Select rows where year < 1834 & group = 0 & TimeElapsed = 0
plot.ls.its.gen.p <- plot.ls.its.gen.p[which(plot.ls.its.gen.p$DeathYear < 1834 & plot.ls.its.gen.p$Group == "0" & plot.ls.its.gen.p$TimeElapsed == 0 | plot.ls.its.gen.p$DeathYear >= 1834 & plot.ls.its.gen.p$Group == "1" & plot.ls.its.gen.p$TimeElapsed == plot.ls.its.gen.p$DeathYear - 1833), ]
plot.ls.its.gen.p$ITS <- "Predicted"


# Plot

its.ls.gender.plot <- ggplot(aes(x = DeathYear, y = DeathAge), data = deaths[!is.na(deaths$Gender),]) +
  # Points
  geom_point(data = deaths[!is.na(deaths$Gender),], aes(x = DeathYear, y = DeathAge), colour = "grey", shape = 16, alpha = 0.4) +
  # CIs
  geom_ribbon(data = plot.ls.its.gen.p[which(plot.ls.its.gen.p$DeathYear >= 1834),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.6) +
  geom_ribbon(data = plot.ls.its.gen.c[which(plot.ls.its.gen.c$DeathYear >= 1834),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.4) +
  geom_ribbon(data = plot.ls.its.gen.c[which(plot.ls.its.gen.c$DeathYear <= 1834.5),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.6) +
  # Lines
  geom_line(data = plot.ls.its.gen.p[which(plot.ls.its.gen.p$DeathYear <= 1833),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.ls.its.gen.p[which(plot.ls.its.gen.p$DeathYear >= 1834),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.ls.its.gen.c[which(plot.ls.its.gen.c$DeathYear >= 1833),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  # Intervention line
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  # Theme and stop y at 100
  theme_minimal() +
  coord_cartesian(ylim=c(0, 100)) +
  # Legend
  scale_linetype_manual(labels = c("Counterfactual", "Predicted values"), values = c("dashed", "solid")) +
  scale_alpha_manual(labels = c("Counterfactual", "Predicted values"), values = c(0.4, 0.4), guide = "none") +
  theme(legend.position = "bottom", legend.title=element_text(face = "bold")) +
  # Remove axis title
  theme(axis.title = element_blank()) +
  #By gender
  facet_wrap(~Gender, labeller = as_labeller(c("F" = "Female", "M" = "Male")))
its.ls.gender.plot


#saveRDS(its.ls.gender.plot, file = "output/lifespans/its.ls.gender.plot.rds")




# Lifespans gender GAM -----------------------------------------------------------------
## Fit model ---------------------------------------------------------------------------
gam.gender <- gam(DeathAge ~ s(DeathYear, bs = "tp", k = -1, by = Gender) + Gender + s(DeathYear.factor, bs = "re"), gamma = 1.4, data = deaths, method = "REML")

#saveRDS(gam.gender, file = "output/lifespans/gam.gender.rds")

## Plot trend -------------------------------------------------------------------------
#Predict from model
p.gen <- as.data.frame(ggpredict(gam.gender, c("DeathYear[1798:2002]", "Gender")))
p.gen <- p.gen %>% rename(DeathYear = x, Gender = group, pfit = predicted, pse = std.error)

#Get derivative
data.gen <- data.frame(DeathYear = p.gen$DeathYear, pfit = p.gen$pfit, pse = p.gen$pse, Gender = p.gen$Gender)
data.gen$Gender <- as.factor(data.gen$Gender)
data.gen$DeathYear.factor <- as.factor(as.character(data.gen$DeathYear))
deriv.gen <- derivatives(gam.gender, term = c("s(DeathYear)"), partial_match = TRUE, data = data.gen)
#Create gender column
deriv.gen$Gender <- ""
deriv.gen[str_detect(deriv.gen$.smooth, "GenderF"), "Gender"] <- "F"
deriv.gen[str_detect(deriv.gen$.smooth, "GenderM"), "Gender"] <- "M"
#Rename data column
#deriv.gen <- deriv.gen %>% rename(DeathYear = data)
#Add in pfit from p.gen
add <- p.gen[, c("DeathYear", "Gender", "pfit")]
deriv.gen <- merge(deriv.gen, add)
deriv.gen$Gender <- as.factor(deriv.gen$Gender)

#Make separate deriv dataframes for inc. vs dec.
deriv.gen.inc <- deriv.gen
deriv.gen.inc[!(deriv.gen.inc$.lower_ci > 0 & deriv.gen.inc$.upper_ci > 0), c(".derivative", "pfit")] <- NA
deriv.gen.dec <- deriv.gen
deriv.gen.dec[!(deriv.gen.dec$.lower_ci < 0 & deriv.gen.dec$.upper_ci < 0), c(".derivative", "pfit")] <- NA


# Overall gam plot with periods of significant change highlighted
gam.ls.gen.plot <- 
  deaths[!is.na(deaths$Gender),] %>%
  ggplot (aes (x = DeathYear, y = DeathAge)) +
  labs (x = "Year of death", y= "Age at death") +
  theme_minimal () +
  # Points
  geom_point(data = deaths[!is.na(deaths$Gender),], aes(x = DeathYear, y = DeathAge), colour = "grey", shape = 16, alpha = 0.6) +
  # Intervention lines
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  # CIs and lines
  geom_ribbon(data = p.gen, aes(x = DeathYear, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.5) +
  geom_line(data = p.gen, aes(x = DeathYear, y = pfit)) +
  #Periods of significant change
  geom_line(data = deriv.gen.inc, aes(x = DeathYear, y = pfit, colour = "blue"), size = 1) +
  geom_line(data = deriv.gen.dec, aes(x = DeathYear, y = pfit, colour = "red"), size = 1) +
  # Legend and axis titles
  coord_cartesian(ylim=c(0, 100)) +
  scale_color_manual(labels = c("Sig. incr.", "Sig. decr."), values = c("blue", "red"), name = "GAM") +
  theme(legend.position = "bottom", legend.title=element_text(face = "bold")) +
  #By gender
  facet_wrap(~Gender, labeller = as_labeller(c("F" = "Female", "M" = "Male")))

gam.ls.gen.plot
#saveRDS(gam.ls.gen.plot, file = "output/lifespans/gam.ls.gen.plot.rds")




## Plot derivative change ----------------------------------------------------
gam.ls.gen.deriv <- ggplot(data = deaths[!is.na(deaths$Gender),], aes(x = DeathYear, y = .derivative)) +
  theme_minimal()  +
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  geom_hline(yintercept = 0) +
  geom_line(data = deriv.gen, aes(x = DeathYear, y = .derivative), colour = "black") +
  geom_ribbon(data = deriv.gen, aes(x = DeathYear, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2) +
  geom_line(data = deriv.gen.inc, aes(x = DeathYear, y = .derivative), colour = "blue", size = 1) +
  geom_line(data = deriv.gen.dec, aes(x = DeathYear, y = .derivative), colour = "red", size = 1) +
  facet_wrap(~Gender, labeller = as_labeller(c("F" = "Female", "M" = "Male")))

gam.ls.gen.deriv
#saveRDS(gam.ls.gen.deriv, file = "output/lifespans/gam.ls.gen.deriv.rds")

  
# Difference between smooths ---------------------------------------------------
pdat <- expand.grid(DeathYear = seq(min(deaths$DeathYear), max(deaths$DeathYear), length = 400), Gender = c("M", "F")) 
pdat$DeathYear.factor <- as.factor(as.character(pdat$DeathYear))
#Using smooth_diff function in setup
comp1 <- smooth_diff(gam.gender, pdat, 'F', 'M', 'Gender')
comp <- cbind(DeathYear = seq(min(deaths$DeathYear), max(deaths$DeathYear), length = 400), comp1)
gam.death.gender.diff <- ggplot(comp, aes(x = DeathYear, y = diff, group = pair)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line() +
  labs(x = 'Year of death', y = 'Difference in gender trend') +
  theme_minimal ()

gam.death.gender.diff
