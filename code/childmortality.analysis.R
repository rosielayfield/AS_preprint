
# ANALYSIS SCRIPT FOR PAPER: Does the institutional environment shape human life-histories? The
# historical case of Augustus Smith and the Isles of Scilly

# CHILD MORTALITY


# Setup -------------------------------------------------------

## Install/load packages --------------------------------------

packages <- c(
  "knitr", "flextable", "ggplot2", "broom", "devtools", "lme4",
  "gridExtra", "merTools", "mgcv", "tidyverse", "gratia", "ggpubr",
  "glmm", "ggeffects", "gamlss", "performance", "cowplot", "patchwork",
  "glmmTMB", "betareg"
)

# Install missing packages
to_install <- setdiff(packages, rownames(installed.packages()))
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
time.line #1923.746


# Create child mortality data -------------------------------------------------------------------

# Child mortality – number of deaths of a child under five divided by number of children under five, calculated per year

childDeaths <- data.frame(DeathYear = unique(deaths$DeathYear))

# Count number of children under 5 who were alive in each year
i <- 1 
childDeaths$DeathYear[i] #1903

for (i in 1:nrow(childDeaths)) {
  sub <- deaths[which(deaths$DeathYear >= childDeaths$DeathYear[i] & deaths$BirthYear <= childDeaths$DeathYear[i]), ]
  sub$currentAge <- childDeaths$DeathYear[i] - sub$BirthYear
  childDeaths$aliveUnder5[i] <- nrow(sub[which(sub$currentAge < 5), ])
}

# Count number of children under 5 who died that year
for (i in 1:nrow(childDeaths)) {
  sub <- deaths[which(deaths$DeathYear == childDeaths$DeathYear[i]), ]
  childDeaths$childDeaths[i] <- nrow(sub[which(sub$DeathAge < 5), ])
}

childDeaths$childMortProp <- childDeaths$childDeaths / childDeaths$aliveUnder5



# Add group and timeelapsed
intervention_year <- 1834
childDeaths <- childDeaths %>%
  mutate(
    group = ifelse(DeathYear < intervention_year, 0, 1),
    timeelapsed = DeathYear - intervention_year
  )

childDeaths[which(childDeaths$DeathYear < 1834), "timeelapsed"] <- 0

#Save
saveRDS(childDeaths, "data/childMortalityData.csv")


## Data for gender interaction models -------------------------------------------
# Create a data frame with unique DeathYear and Gender
childDeaths.g <- deaths %>%
  distinct(DeathYear, Gender) %>%
  arrange(DeathYear, Gender)  # Ensure sorting

# Count number of children under 5 who were alive in each year, grouped by gender
childDeaths.g <- childDeaths.g %>%
  rowwise() %>%
  mutate(
    aliveUnder5 = sum(
      deaths$BirthYear <= DeathYear & 
        deaths$DeathYear >= DeathYear & 
        deaths$Gender == Gender & 
        (DeathYear - deaths$BirthYear) < 5, na.rm = TRUE
    )
  ) %>%
  ungroup()


# Count number of children under 5 who died that year, grouped by gender
childDeaths.g <- childDeaths.g %>%
  rowwise() %>%
  mutate(
    childDeaths = sum(
      deaths$DeathYear == DeathYear & 
        deaths$Gender == Gender & 
        deaths$DeathAge < 5, na.rm = TRUE
    )
  ) %>%
  ungroup()

# Compute child mortality proportion
childDeaths.g <- childDeaths.g %>%
  mutate(childMortProp = ifelse(aliveUnder5 > 0, childDeaths / aliveUnder5, NA))  # Avoid division by zero

# Define intervention year
intervention_year <- 1834

# Add group and timeelapsed
childDeaths.g <- childDeaths.g %>%
  mutate(
    group = ifelse(DeathYear < intervention_year, 0, 1),
    timeelapsed = pmax(DeathYear - intervention_year, 0)  # Ensures no negative values
  )

childDeaths.g <- childDeaths.g %>%
  filter(!is.na(Gender))

childDeaths.g$Gender <- as.factor(childDeaths.g$Gender)

saveRDS(childDeaths.g, "data/childMortalityGenderData.csv")



# ITS binomial ------------------------------------------------------------------------
ts.cm.binomial <- glm(cbind(childDeaths, aliveUnder5 - childDeaths) ~ DeathYear + group + timeelapsed, 
                      family = binomial(), 
                      data = childDeaths)

summary(ts.cm.binomial)
confint(ts.cm.binomial)




## Plot ------------------------------------------------------------------------

#Get predicted values for plot
#Counterfactual line
plot.cm.its.c <- as.data.frame(ggpredict(ts.cm.binomial, c("DeathYear[1798:1981]", "timeelapsed[0]", "group[0]")))
plot.cm.its.c <- plot.cm.its.c %>% rename(DeathYear = x, timeelapsed = group, group = facet)
plot.cm.its.c$timeelapsed <- as.numeric(as.character(plot.cm.its.c$timeelapsed))
plot.cm.its.c$ITS <- "Counterfactual"


#Predicted line
plot.cm.its.p <- as.data.frame(ggpredict(ts.cm.binomial, c("DeathYear[1798:1981]","timeelapsed[0:168]", "group")))
plot.cm.its.p <- plot.cm.its.p %>% rename(DeathYear = x, timeelapsed = group, group = facet)
#Select rows where year < 1834 & group = 0 & TimeElapsed = 0
plot.cm.its.p <- plot.cm.its.p[which(plot.cm.its.p$DeathYear < 1834 & plot.cm.its.p$group == "0" & plot.cm.its.p$timeelapsed == 0 | plot.cm.its.p$DeathYear >= 1834 & plot.cm.its.p$group == "1" & plot.cm.its.p$timeelapsed == plot.cm.its.p$DeathYear - 1833), ]
plot.cm.its.p$ITS <- "Predicted"


# Plot
ts.cm.plot <- ggplot(aes(x = DeathYear, y = childDeaths), data = childDeaths) +
  # Points
  #geom_point(data = childDeaths, aes(x = DeathYear, y = childMortProp), colour = "grey", shape = 16, alpha = 0.4) +
  # CIs
  geom_ribbon(data = plot.cm.its.p[which(plot.cm.its.p$DeathYear >= 1834),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_ribbon(data = plot.cm.its.c[which(plot.cm.its.c$DeathYear > 1834.5),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_ribbon(data = plot.cm.its.c[which(plot.cm.its.c$DeathYear <= 1834.5),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  # Lines
  geom_line(data = plot.cm.its.p[which(plot.cm.its.p$DeathYear <= 1833),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.cm.its.p[which(plot.cm.its.p$DeathYear >= 1834),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.cm.its.c[which(plot.cm.its.p$DeathYear >= 1833),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  #scale_x_continuous(breaks = sort(c(seq(1700, 2000, length.out=4), 1834))) +
  # Intervention line
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  # Theme and stop y at 100
  theme_minimal() +
  xlab("Year of Death") +
  ylab("Child mortality\nprobability") +
  coord_cartesian(ylim=c(0, 0.12), xlim=c(1798,1981)) +
  # Legend
  scale_linetype_manual(labels = c("Counterfactual", "Predicted values"), values = c("dashed", "solid")) +
  scale_alpha_manual(labels = c("Counterfactual", "Predicted values"), values = c(0.4, 0.4), guide = "none") +
  theme(legend.position = "bottom", legend.title=element_text(face = "bold")) +
  theme(legend.position = "none") +
  ggtitle("All data")
ts.cm.plot 

#saveRDS(ts.cm.plot, file = "output/lifespans/ts.cm.plot.rds")





# GAM binomial ------------------------------------------------------------------------
# Fit a binomial GAM
cm.gam <- gam(cbind(childDeaths, aliveUnder5 - childDeaths) ~ s(DeathYear, bs = "tp", k = -1), 
              family = binomial(), gamma = 1.4, method = "REML",
              data = childDeaths)


# View model summary
summary(cm.gam)

#saveRDS(cm.gam, file = "output/lifespans/cm.gam.rds")



## Plot binomial GAM --------------------------------------------------

# Predicted data for plots
p.ls <- as.data.frame(ggpredict(cm.gam, c("DeathYear[1798:1981]")))
p.ls$group <- NULL
p.ls <- p.ls %>% rename(DeathYear = x, pfit = predicted, pse = std.error)


# Get derivative
data.ls <- data.frame(DeathYear = 1798:1981, aliveUnder5_log = log(10))
deriv.ls <- derivatives(cm.gam, data = data.ls)
deriv.ls <- derivatives(cm.gam, select = "s(DeathYear)", data = data.ls)
names(deriv.ls)[names(deriv.ls) == "data"] <- "DeathYear"
deriv.ls <- transform(deriv.ls, pfit = p.ls$pfit)


#Make separate deriv dataframes for inc. vs dec.
deriv.ls.inc <- deriv.ls
deriv.ls.inc[!(deriv.ls.inc$.lower_ci > 0 & deriv.ls.inc$.upper_ci > 0), c(".derivative", "pfit")] <- NA
deriv.ls.dec <- deriv.ls
deriv.ls.dec[!(deriv.ls.dec$.lower_ci < 0 & deriv.ls.dec$.upper_ci < 0), c(".derivative", "pfit")] <- NA



# Plot overall GAM

#Plotting data and estimated trend line with periods of significant change indicated
gam.cm.plot <- childDeaths %>%
  ggplot (aes (x = DeathYear)) +
  labs (x = "Year of Death", y= "Child mortality\nprobability") +
  theme_minimal () +
  # Points
  #geom_point(data = childDeaths, aes(x = DeathYear, y = childMortProp), colour = "grey", shape = 16, alpha = 0.6) +
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
  coord_cartesian(xlim = c(1798, 1981), ylim = c(0, 0.8)) +
  xlab(NULL)

gam.cm.plot 

#saveRDS(gam.cm.plot, file = "output/lifespans/gam.cm.plot.rds")


## Plot derivative change ------------------------------------------------------------------

#Plot derivative change
gam.cm.deriv <-  ggplot(data = deriv.ls, aes(x = DeathYear, y = .derivative)) +
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
  coord_cartesian(xlim = c(1798, 1981), ylim = c(-0.25, 0.15))
gam.cm.deriv

#saveRDS(gam.cm.deriv, file = "output/lifespans/gam.cm.deriv.rds")



# ITS interaction ------------------------------------------------------------------

ts.cm.gen.binomial <- glm(cbind(childDeaths, aliveUnder5 - childDeaths) ~ 
                            DeathYear +
                            group +
                            timeelapsed +
                            Gender +
                            DeathYear:Gender + 
                            group:Gender + 
                            timeelapsed:Gender, 
                          family = binomial(), 
                          data = childDeaths.g)

summary(ts.cm.gen.binomial)



## Plot ------------------------------------------------------------------------
#Predict from model - ggpredict
#Counter
plot.cm.its.gen.c <- as.data.frame(ggpredict(ts.cm.gen.binomial, c("DeathYear[1798:1981]", "Gender", "timeelapsed[0]", "group")))
plot.cm.its.gen.c <- plot.cm.its.gen.c %>% rename(DeathYear = x, Gender = group, TimeElapsed = facet, Group = panel)
#Remove rows where Group = 1
plot.cm.its.gen.c <- plot.cm.its.gen.c[which(plot.cm.its.gen.c$Group == 0), ]
plot.cm.its.gen.c$TimeElapsed <- as.numeric(as.character(plot.cm.its.gen.c$TimeElapsed))
plot.cm.its.gen.c$ITS <- "Counterfactual"


#Predicted
plot.cm.its.gen.p <- as.data.frame(ggpredict(ts.cm.gen.binomial, c("DeathYear[1798:1981]", "Gender", "timeelapsed[0:168]", "group")))
plot.cm.its.gen.p <- plot.cm.its.gen.p %>% rename(DeathYear = x, Gender = group, TimeElapsed = facet, Group = panel)
#Select rows where year < 1834 & group = 0 & TimeElapsed = 0
plot.cm.its.gen.p <- plot.cm.its.gen.p[which(plot.cm.its.gen.p$DeathYear < 1834 & plot.cm.its.gen.p$Group == "0" & plot.cm.its.gen.p$TimeElapsed == 0 | plot.cm.its.gen.p$DeathYear >= 1834 & plot.cm.its.gen.p$Group == "1" & plot.cm.its.gen.p$TimeElapsed == plot.cm.its.gen.p$DeathYear - 1833), ]
plot.cm.its.gen.p$ITS <- "Predicted"


# Plot

its.cm.girls.plot <- ggplot(aes(x = DeathYear, y = childMortProp), data = childDeaths.g[which(childDeaths.g$Gender == "F"), ]) +
  # Points
  #geom_point(data = childDeaths.g, aes(x = DeathYear, y = childMortProp), colour = "grey", shape = 16, alpha = 0.4) +
  # CIs
  geom_ribbon(data = plot.cm.its.gen.p[which(plot.cm.its.gen.p$DeathYear >= 1834 & plot.cm.its.gen.p$Gender == "F"),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_ribbon(data = plot.cm.its.gen.c[which(plot.cm.its.gen.c$DeathYear >= 1834 & plot.cm.its.gen.c$Gender == "F"),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_ribbon(data = plot.cm.its.gen.c[which(plot.cm.its.gen.c$DeathYear <= 1834.5 & plot.cm.its.gen.c$Gender == "F"),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  # Lines
  geom_line(data = plot.cm.its.gen.p[which(plot.cm.its.gen.p$DeathYear <= 1833 & plot.cm.its.gen.p$Gender == "F"),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.cm.its.gen.p[which(plot.cm.its.gen.p$DeathYear >= 1834 & plot.cm.its.gen.p$Gender == "F"),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.cm.its.gen.c[which(plot.cm.its.gen.c$DeathYear >= 1833 & plot.cm.its.gen.c$Gender == "F"),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  # Intervention line
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  xlab("Year of Death") +
  ylab("Child mortality\nprobability") +
  # Theme and stop y at 100
  theme_minimal() +
  coord_cartesian(ylim=c(0, 0.12), xlim=c(1800,1981)) +
  # Legend
  scale_linetype_manual(labels = c("Counterfactual", "Predicted values"), values = c("dashed", "solid")) +
  scale_alpha_manual(labels = c("Counterfactual", "Predicted values"), values = c(0.4, 0.4), guide = "none") +
  theme(legend.position = "bottom", legend.title=element_text(face = "bold"), strip.text = element_blank()) +
  # Remove axis title
  #theme(axis.title = element_blank()) +
  #By gender
  facet_wrap(~Gender, labeller = as_labeller(c("F" = "Girls", "M" = "Boys"))) +
  ggtitle("Girls")
its.cm.girls.plot

#saveRDS(its.cm.girls.plot, file = "output/lifespans/its.cm.girls.plot.rds")



its.cm.boys.plot <- ggplot(aes(x = DeathYear, y = childMortProp), data = childDeaths.g[which(childDeaths.g$Gender == "M"), ]) +
  # Points
  #geom_point(data = childDeaths.g, aes(x = DeathYear, y = childMortProp), colour = "grey", shape = 16, alpha = 0.4) +
  # CIs
  geom_ribbon(data = plot.cm.its.gen.p[which(plot.cm.its.gen.p$DeathYear >= 1834 & plot.cm.its.gen.p$Gender == "M"),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_ribbon(data = plot.cm.its.gen.c[which(plot.cm.its.gen.c$DeathYear >= 1834 & plot.cm.its.gen.c$Gender == "M"),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_ribbon(data = plot.cm.its.gen.c[which(plot.cm.its.gen.c$DeathYear <= 1834.5 & plot.cm.its.gen.c$Gender == "M"),], aes(x = DeathYear, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  # Lines
  geom_line(data = plot.cm.its.gen.p[which(plot.cm.its.gen.p$DeathYear <= 1833 & plot.cm.its.gen.p$Gender == "M"),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.cm.its.gen.p[which(plot.cm.its.gen.p$DeathYear >= 1834 & plot.cm.its.gen.p$Gender == "M"),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.cm.its.gen.c[which(plot.cm.its.gen.c$DeathYear >= 1833 & plot.cm.its.gen.c$Gender == "M"),], aes(x=DeathYear, y = predicted, linetype = ITS), colour = "black") +
  # Intervention line
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  xlab("Year of Death") +
  ylab("Child mortality\nprobability") +
  # Theme and stop y at 100
  theme_minimal() +
  coord_cartesian(ylim=c(0, 0.12), xlim=c(1800,1981)) +
  # Legend
  scale_linetype_manual(labels = c("Counterfactual", "Predicted values"), values = c("dashed", "solid")) +
  scale_alpha_manual(labels = c("Counterfactual", "Predicted values"), values = c(0.4, 0.4), guide = "none") +
  theme(legend.position = "none", strip.text = element_blank()) +
  # Remove axis title
  #theme(axis.title = element_blank()) +
  #By gender
  facet_wrap(~Gender, labeller = as_labeller(c("F" = "Girls", "M" = "Boys"))) +
    ggtitle("Boys")
its.cm.boys.plot

#saveRDS(its.cm.boys.plot, file = "output/lifespans/its.cm.boys.plot.rds")




# GAM interaction -----------------------------------------------------------

# Fit a binomial GAM
cm.gam.gen <- gam(cbind(childDeaths, aliveUnder5 - childDeaths) ~ s(DeathYear, bs = "tp", k = -1, by = Gender) + Gender, 
                  family = binomial(), gamma = 1.4, method = "REML",
                  data = childDeaths.g)

#saveRDS(cm.gam.gen, file = "output/lifespans/cm.gam.gen.rds")



## Plot ----------------------------------------------------------------------

#Predict from model
p.gen <- as.data.frame(ggpredict(cm.gam.gen, c("DeathYear[1798:1981]", "Gender")))
p.gen <- p.gen %>% rename(DeathYear = x, Gender = group, pfit = predicted, pse = std.error)

#Get derivative
data.gen <- data.frame(DeathYear = p.gen$DeathYear, pfit = p.gen$pfit, pse = p.gen$pse, Gender = p.gen$Gender)
data.gen$Gender <- as.factor(data.gen$Gender)
data.gen$DeathYear.factor <- as.factor(as.character(data.gen$DeathYear))
deriv.gen <- derivatives(cm.gam.gen, term = c("s(DeathYear)"), partial_match = TRUE, data = data.gen)
#Create gender column
deriv.gen$Gender <- ""
deriv.gen[str_detect(deriv.gen$.smooth, "GenderF"), "Gender"] <- "F"
deriv.gen[str_detect(deriv.gen$.smooth, "GenderM"), "Gender"] <- "M"

#Add in pfit from p.gen
add <- p.gen[, c("DeathYear", "Gender", "pfit")]
deriv.gen <- merge(deriv.gen, add)
deriv.gen$Gender <- as.factor(deriv.gen$Gender)

#Make separate deriv dataframes for inc. vs dec.
deriv.gen.inc <- deriv.gen
deriv.gen.inc[!(deriv.gen.inc$.lower_ci > 0 & deriv.gen.inc$.upper_ci > 0), c(".derivative", "pfit")] <- NA
deriv.gen.dec <- deriv.gen
deriv.gen.dec[!(deriv.gen.dec$.lower_ci < 0 & deriv.gen.dec$.upper_ci < 0), c(".derivative", "pfit")] <- NA


# Plot 
gam.cm.gen.female.plot <- 
childDeaths.g %>%
  ggplot (aes (x = DeathYear, y = childMortProp)) +
  labs (x = "Year of death", y = "Child mortality\nprobability") +
  theme_minimal () +
  # Points
  #geom_point(data = deaths[!is.na(deaths$Gender) & deaths$Gender == "F",], aes(x = DeathYear, y = DeathAge), colour = "grey", shape = 16, alpha = 0.6) +
  # Intervention lines
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  # CIs and lines
  geom_ribbon(data = p.gen[p.gen$Gender == "F", ], aes(x = DeathYear, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.3) +
  geom_line(data = p.gen[p.gen$Gender == "F", ], aes(x = DeathYear, y = pfit)) +
  #Periods of significant change
  geom_line(data = deriv.gen.inc[deriv.gen.inc$Gender == "F", ], aes(x = DeathYear, y = pfit, colour = "blue"), size = 1) +
  geom_line(data = deriv.gen.dec[deriv.gen.dec$Gender == "F", ], aes(x = DeathYear, y = pfit, colour = "red"), size = 1) +
  # Legend and axis titles
  scale_color_manual(labels = c("Sig. incr.", "Sig. decr."), values = c("blue", "red"), name = "GAM") +
  theme(legend.position = "bottom", legend.title=element_text(face = "bold")) +
  theme(legend.position = "none") +
  xlab(NULL) +
  ylab(NULL) +
  coord_cartesian(xlim = c(1798, 1981), ylim = c(0, 0.8))
gam.cm.gen.female.plot

#saveRDS(gam.cm.gen.female.plot, file = "output/lifespans/gam.cm.gen.female.plot.rds")


gam.cm.gen.male.plot <- 
childDeaths.g %>%
  ggplot (aes (x = DeathYear, y = childMortProp)) +
  labs (x = "Year of death", y = "Child deaths prop") +
  theme_minimal () +
  # Points
  #geom_point(data = deaths[!is.na(deaths$Gender) & deaths$Gender == "M",], aes(x = DeathYear, y = DeathAge), colour = "grey", shape = 16, alpha = 0.6) +
  # Intervention lines
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  # CIs and lines
  geom_ribbon(data = p.gen[p.gen$Gender == "M", ], aes(x = DeathYear, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.3) +
  geom_line(data = p.gen[p.gen$Gender == "M", ], aes(x = DeathYear, y = pfit)) +
  #Periods of significant change
  geom_line(data = deriv.gen.inc[deriv.gen.inc$Gender == "M", ], aes(x = DeathYear, y = pfit, colour = "blue"), size = 1) +
  geom_line(data = deriv.gen.dec[deriv.gen.dec$Gender == "M", ], aes(x = DeathYear, y = pfit, colour = "red"), size = 1) +
  # Legend and axis titles
  scale_color_manual(labels = c("Sig. incr.", "Sig. decr."), values = c("blue", "red"), name = "GAM") +
  theme(legend.position = "bottom", legend.title=element_text(face = "bold"))+
  theme(legend.position = "none") +
  xlab(NULL) +
  ylab(NULL) +
  coord_cartesian(xlim = c(1795, 1981), ylim = c(0, 0.8))

gam.cm.gen.male.plot 
#saveRDS(gam.cm.gen.male.plot, file = "output/lifespans/gam.cm.gen.male.plot.rds")




## Plot derivative change ------------------------------------------------------------

gam.cm.gen.female.deriv <- ggplot(data = deaths[!is.na(deaths$Gender) & deaths$Gender == "F",], aes(x = DeathYear, y = .derivative)) +
  theme_minimal()  +
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  geom_hline(yintercept = 0) +
  geom_line(data = deriv.gen[deriv.gen$Gender == "F", ], aes(x = DeathYear, y = .derivative), colour = "black") +
  geom_ribbon(data = deriv.gen[deriv.gen$Gender == "F", ], aes(x = DeathYear, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.3) +
  geom_line(data = deriv.gen.inc[deriv.gen.inc$Gender == "F", ], aes(x = DeathYear, y = .derivative, colour = "blue"), size = 1) +
  geom_line(data = deriv.gen.dec[deriv.gen.dec$Gender == "F", ], aes(x = DeathYear, y = .derivative, colour = "red"), size = 1) +
  xlab("Year of death") +
  ylab(NULL) +
  coord_cartesian(xlim = c(1798, 1981), ylim = c(-0.25, 0.15)) +
  scale_color_manual(labels = c("Significant increase", "Significant decrease"), values = c("blue", "red"), name = "GAM") +
  theme(legend.position = "bottom", legend.title=element_text(face = "bold"))
gam.cm.gen.female.deriv


#saveRDS(gam.cm.gen.female.deriv, file = "output/lifespans/gam.cm.gen.female.deriv.rds")


gam.cm.gen.male.deriv <- ggplot(data = deaths[!is.na(deaths$Gender) & deaths$Gender == "M",], aes(x = DeathYear, y = .derivative)) +
  theme_minimal()  +
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  geom_hline(yintercept = 0) +
  geom_line(data = deriv.gen[deriv.gen$Gender == "M", ], aes(x = DeathYear, y = .derivative), colour = "black") +
  geom_ribbon(data = deriv.gen[deriv.gen$Gender == "M", ], aes(x = DeathYear, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.3) +
  geom_line(data = deriv.gen.inc[deriv.gen.inc$Gender == "M", ], aes(x = DeathYear, y = .derivative), colour = "blue", size = 1) +
  geom_line(data = deriv.gen.dec[deriv.gen.dec$Gender == "M", ], aes(x = DeathYear, y = .derivative), colour = "red", size = 1)+
  xlab("Year of death") +
  coord_cartesian(xlim = c(1798, 1981), ylim = c(-0.25, 0.15)) +
  ylab(NULL)

gam.cm.gen.male.deriv
#saveRDS(gam.cm.gen.male.deriv, file = "output/lifespans/gam.cm.gen.male.deriv.rds")

