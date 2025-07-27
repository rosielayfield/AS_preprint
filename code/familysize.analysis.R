
# Setup -------------------------------------------------------

## Install/load packages --------------------------------------
packages <- c(
  "knitr", "flextable", "ggplot2", "broom", "devtools", "lme4",
  "gridExtra", "merTools", "mgcv", "tidyverse", "gratia", "ggpubr",
  "glmm", "ggeffects", "boot", "patchwork"
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

## Function to check for overdispersion ----------------------------------------------------
#Function to calculate a point estimate of overdispersion from a mixed model object
od.point<-function(modelobject){
  x<-sum(resid(modelobject,type="pearson")^2)
  rdf<-summary(modelobject)$AICtab[5]
  return(x/rdf)
}

#Function to pass to parametric bootstrap function 'bootMer' that calculates the sum of squared Pearson residuals (required for 'od' function)
FUN <- function(fit) {
  #return(fixef(fit))
  x<-resid(fit,type="pearson")
  return(sum(x^2))
}	
#Function To Calculate Ratio of Model SS to Mean Parametric Bootstrap SS ('bias')
od<-function(bootobject){
  biasvals<-bootobject $t0/bootobject[2]$t
  bias<-mean(biasvals,na.rm=T)
  intervals<-quantile(biasvals,c(0.025,0.975),na.rm=T)
  dat<-c(bias,intervals)
  return(dat)
}


## Load data -----------------------------------------------------

#Data - parents
parents <- read.csv("data/parents.csv")
parents$Group <- as.factor(as.character(parents$Group))
parents$first.birth.factor <- as.factor(as.character(parents$first.birth))
parents[is.na(parents$HCLASS.cat), "HCLASS.cat"] <- "Unknown"
parents$HCLASS.cat <- as.factor(parents$HCLASS.cat)
parents$Island <- as.factor(parents$Island)
parents$OLRE <-seq(nrow(parents))
summary(parents)

#Get end of relevant period for gam = end of intervention period + two generations, 
#where the length of a generation is the average age at first reproduction
time.line <- 1872 + (2 * mean(parents$mother_age_first_repro, na.rm = T)) #25.87



# Family size ITS ----------------------------------------------------------------------
## Check for overdispersion ------------------------------------------------------------
#ts.main.fam.old <- glmer(children ~ first.birth + Group + TimeElapsed  + (1|first.birth), family = "poisson", 
#                         control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000000)),data = parents)

#Check for overdispersion - point estimate of overdispersion from a mixed model object
#od.point(ts.main.fam.old)
#1.494


# Refit model using observation level random effects to deal with overdispersion
#ts.main.fam.olre <- glmer(children ~ first.birth + Group + TimeElapsed  + (1|first.birth) + (1|OLRE),
#                          family = "poisson", data = parents)
#od.point(ts.main.fam.olre)
#0.501

## Fit model ---------------------------------------------------------------------------
# Fit model using observation level random effects to deal with overdispersion
ts.main.fam <- glmer(children ~ first.birth + Group + TimeElapsed  + (1|first.birth) + (1|OLRE), family = "poisson",
                     data = parents, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000000)))

#od.point(ts.main.fam)
#AIC(ts.main.fam.old, ts.main.fam)

#saveRDS(ts.main.fam, file = "output/lifespans/ts.main.fam.rds")
summary(ts.main.fam)
confint(ts.main)



## Results table ------------------------------------------------------
#Get confidence intervals
coeftbl <- as.data.frame(coef(summary(ts.main.fam)))
conf.ints <- with(coeftbl,
                  +      Estimate + outer(`Std. Error`, c(lower=-1, upper=1)) * sqrt(qchisq(0.95, 1)))


tab.ts.main.fam <- cbind(data.frame(Predictors = c("(Intercept)", "Year first birth", "Group", "Time elapsed")), 
                         as.data.frame(coef(summary(ts.main.fam))), row.names = NULL)

tab.ts.main.fam <- cbind(tab.ts.main.fam, conf.ints)
names(tab.ts.main.fam)[names(tab.ts.main.fam) == "lower"] <- "2.5%"
names(tab.ts.main.fam)[names(tab.ts.main.fam) == "upper"] <- "97.5%"

#saveRDS(tab.ts.main.fam, file = "output/lifespans/tab.ts.main.fam.rds")


## Plot ------------------------------------------------------------------
#Predict from model
#Counter
range(parents$first.birth)
plot.fam.its.c <- as.data.frame(ggpredict(ts.main.fam, c("first.birth[1754:1915]", "TimeElapsed[0]")))
plot.fam.its.c <- plot.fam.its.c %>% rename(first.birth = x, TimeElapsed = group)
plot.fam.its.c$TimeElapsed <- as.numeric(as.character(plot.fam.its.c$TimeElapsed))
plot.fam.its.c$ITS <- "Counterfactual"

#Predicted
plot.fam.its.p <- as.data.frame(ggpredict(ts.main.fam, c("first.birth[1754:1915]", "TimeElapsed[0:82]", "Group")))
plot.fam.its.p <- plot.fam.its.p %>% rename(first.birth = x, TimeElapsed = group, Group = facet)
#Select rows where year < 1834 & group = 0 & TimeElapsed = 0
plot.fam.its.p <- plot.fam.its.p[which(plot.fam.its.p$first.birth < 1834 & plot.fam.its.p$Group == "0" & plot.fam.its.p$TimeElapsed == 0 | plot.fam.its.p$first.birth >= 1834 & plot.fam.its.p$Group == "1" & plot.fam.its.p$TimeElapsed == plot.fam.its.p$first.birth - 1833), ]
plot.fam.its.p$ITS <- "Predicted"


# Plot

its.fam.plot <- ggplot(aes(x = first.birth, y = children), data = parents) +
  # Points
  geom_point(data = parents, aes(x = first.birth, y = children), colour = "grey", shape = 16, alpha = 0.6) +
  # CIs
  geom_ribbon(data = plot.fam.its.p[which(plot.fam.its.p$first.birth >= 1834), ], aes(x = first.birth, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.4) +
  geom_ribbon(data = plot.fam.its.c[which(plot.fam.its.c$first.birth >= 1834), ], aes(x = first.birth, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.4) +
  geom_ribbon(data = plot.fam.its.c[which(plot.fam.its.c$first.birth <= 1834.5), ], aes(x = first.birth, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.4) +
  # Lines
  geom_line(data = plot.fam.its.p[which(plot.fam.its.p$first.birth < 1834), ], aes(x=first.birth, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.fam.its.p[which(plot.fam.its.p$first.birth >=1834), ], aes(x=first.birth, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.fam.its.c, aes(x=first.birth, y = predicted, linetype = ITS), colour = "black") + 
  #scale_x_continuous(breaks = sort(c(seq(1700, 2000, length.out=4), 1834))) +
  # Intervention line
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  # Theme and stop y at 15
  theme_minimal() +
  coord_cartesian(ylim=c(0, 15)) +
  # Legend
  scale_linetype_manual(labels = c("Counterfactual", "Predicted values"), values = c("dashed", "solid")) +
  scale_alpha_manual(labels = c("Counterfactual", "Predicted values"), values = c(0.4, 0.4), guide = "none") +
  theme(legend.position = "bottom", legend.title=element_text(face = "bold")) +
  xlab("Year of first birth in family") +
  ylab("Number of children") +
  ggtitle("All data")

its.fam.plot


#saveRDS(its.fam.plot, file = "output/lifespans/its.fam.plot.rds")








# Family size GAM ----------------------------------------------------------
## Check for overdispersion ------------------------------------------------
#gam.fam.old <- gam(children ~ s(first.birth, bs = "tp", k = -1) + s(first.birth.factor, bs = "re"), family = "poisson", gamma = 1.4, data = parents, method = "REML")
#gam.check(gam.fam)

#Check for overdispersion
#sum(residuals(gam.fam.old, type = "pearson")^2) / df.residual(gam.fam.old)

# Fix overdispersion
#parents$OLRE <-seq(nrow(parents))
#gam.fam.p.olre <- gam(children ~ s(first.birth, bs = "tp", k = -1) + s(first.birth.factor, bs = "re") + s(OLRE, bs = "re"#), family = "poisson", gamma = 1.4, data = parents, method = "REML")

#sum(residuals(gam.fam.p.olre, type = "pearson")^2) / df.residual(gam.fam.p.olre)

## Fit model ----------------------------------------------------------------
gam.fam <- gam(children ~ s(first.birth, bs = "tp", k = -1) + s(first.birth.factor, bs = "re") + s(OLRE, bs = "re"), 
               family = "poisson", gamma = 1.4, data = parents, method = "REML")

summary(gam.fam)
#saveRDS(gam.fam, file = "output/lifespans/gam.fam.rds")


## Plot --------------------------------------------------------------------

#Predict from model
pred <- ggpredict(gam.fam, c("first.birth[1754:1915]"))
pred$group <- NULL
pred <- pred %>% rename(first.birth = x, pfit = predicted, pse = std.error)

#Get derivative
data.fam <- data.frame(first.birth = 1754:1915)
data.fam$first.birth.factor <- as.factor(as.character(data.fam$first.birth))
data.fam$OLRE <- 1
deriv.fam <- derivatives(gam.fam, term = "s(first.birth)", data = data.fam)
names(deriv.fam)[names(deriv.fam) == "data"] <- "first.birth"
deriv.fam <- transform(deriv.fam, pfit = pred$pfit)

#Make separate deriv dataframes for inc. vs dec.
deriv.fam.inc <- deriv.fam
deriv.fam.inc[!(deriv.fam.inc$.lower_ci > 0 & deriv.fam.inc$.upper_ci > 0), c(".derivative", "pfit")] <- NA
deriv.fam.dec <- deriv.fam
deriv.fam.dec[!(deriv.fam.dec$.lower_ci < 0 & deriv.fam.dec$.upper_ci < 0), c(".derivative", "pfit")] <- NA



# Plot

# Plot overall GAM
fam.gam.plot <-  ggplot(data = pred) +
# Points
  geom_point(data = parents, aes(x = first.birth, y = children), colour = "grey", shape = 16, aplha = 0.6) +
  geom_ribbon (aes(x = first.birth, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.4) +
  geom_line(aes(x = first.birth, y = pfit)) +
  geom_line(data = deriv.fam.inc, aes(x = first.birth, y = pfit), colour = "blue") +
  geom_line(data = deriv.fam.dec, aes(x = first.birth, y = pfit), colour = "red") +
  labs (x = "Year of first birth", y= "Children") +
  # Intervention lines
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  coord_cartesian(ylim=c(0, 15)) +
  xlab(NULL) +
  ylab("Number of children")+
  theme_minimal()
fam.gam.plot 
#saveRDS(fam.gam.plot, file = "output/lifespans/fam.gam.plot.rds")


## Plot derivative change ----------------------------------------------
deriv.fam.plot <-  ggplot(data = deriv.fam, aes(x = first.birth, y = .derivative)) +
  theme_minimal()  +
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim=c(-0.08, 0.05)) +
  geom_line(data = deriv.fam, aes(x = first.birth, y = .derivative), colour = "black") +
  geom_ribbon(data = deriv.fam, aes(x = first.birth, ymin = .lower_ci, ymax = .upper_ci), alpha = 0.3) +
  geom_line(data = deriv.fam.inc, aes(x = first.birth, y = .derivative, colour = "blue"), size = 1) +
  geom_line(data = deriv.fam.dec, aes(x = first.birth, y = .derivative, colour = "red"), size = 1) +
    scale_color_manual(labels = c("Sig. incr.", "Sig. decr."), values = c("blue", "red"), name = "GAM") +
    theme(legend.position = "bottom", legend.title=element_text(face = "bold")) +
    xlab("Year of first birth in family") +
    ylab("Derivative")
    
deriv.fam.plot 
#saveRDS(deriv.fam.plot, file = "output/lifespans/deriv.fam.plot.rds")



# Family size occupation ITS ---------------------------------------------------------------
## Fit model -------------------------------------------------------------------------------
ts.fam.ses.cat <- glmer(children ~ first.birth + Group + TimeElapsed  + 
                          HCLASS.cat + HCLASS.cat:first.birth + HCLASS.cat:Group + HCLASS.cat:TimeElapsed +
                          (1|first.birth) + (1|OLRE), family = "poisson", data = parents, 
                        control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=20000000)))

summary(ts.fam.ses.cat)

#saveRDS(ts.fam.ses.cat, file = "output/lifespans/ts.fam.ses.cat.rds")


## Relevel to get SE of the slope for Non-manual and unknown ------------------------------
### Non-Manual ----------------------------------------------------------------------------
parents$HCLASS.cat <- relevel(parents$HCLASS.cat, ref = "Non-Manual")

ts.nm <- glmer(children ~ first.birth + Group + TimeElapsed  + 
                 HCLASS.cat + HCLASS.cat:first.birth + HCLASS.cat:Group + HCLASS.cat:TimeElapsed +
                 (1|first.birth) + (1|OLRE), data = parents, family = poisson(link = "log"),
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e7)),
               nAGQ = 0)
summary(ts.nm)

#saveRDS(ts.nm, file = "output/lifespans/ts.fam.nm.rds")


### Unknown ----------------------------------------------------------------------
parents$HCLASS.cat <- relevel(parents$HCLASS.cat, ref = "Unknown")

ts.un <- glmer(children ~ first.birth + Group + TimeElapsed  + 
                 HCLASS.cat + HCLASS.cat:first.birth + HCLASS.cat:Group + HCLASS.cat:TimeElapsed +
                 (1|first.birth) + (1|OLRE), data = parents, family = poisson(link = "log"),
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e7)),
               nAGQ = 0)
summary(ts.un)

#saveRDS(ts.un, file = "output/lifespans/ts.fam.un.rds")



## Plot -------------------------------------------------------------------------
#Counterfactual
plot.fam.its.ses.c <- as.data.frame(ggpredict(ts.fam.ses.cat, c("first.birth[1754:1915]", "HCLASS.cat", "TimeElapsed[0]")))
plot.fam.its.ses.c <- plot.fam.its.ses.c %>% rename(first.birth = x, HCLASS.cat = group, TimeElapsed = facet)
plot.fam.its.ses.c$TimeElapsed <- as.numeric(as.character(plot.fam.its.ses.c$TimeElapsed))
plot.fam.its.ses.c$ITS <- "Counterfactual"


#Predicted
plot.fam.its.ses.p <- as.data.frame(ggpredict(ts.fam.ses.cat, c("first.birth[1754:1915]", "HCLASS.cat", "TimeElapsed[0:82]", "Group")))
plot.fam.its.ses.p  <- plot.fam.its.ses.p  %>% rename(first.birth = x, HCLASS.cat = group, TimeElapsed = facet, Group = panel)
#Select rows where year < 1834 & group = 0 & TimeElapsed = 0
plot.fam.its.ses.p <- plot.fam.its.ses.p[which(plot.fam.its.ses.p$first.birth < 1834 & plot.fam.its.ses.p$Group == "0" & plot.fam.its.ses.p$TimeElapsed == 0 | 
                                                 plot.fam.its.ses.p$first.birth >= 1834 & plot.fam.its.ses.p$Group == "1" & plot.fam.its.ses.p$TimeElapsed == plot.fam.its.ses.p$first.birth - 1833), ]
plot.fam.its.ses.p$ITS <- "Predicted"


# PLot

its.fam.ses.plot <-
  ggplot(aes(x = first.birth, y = children), data = parents) +
  # Points
  geom_point(data = parents, aes(x = first.birth, y = children), colour = "grey", shape = 16, alpha = 0.4) +
  # CIs
  geom_ribbon(data = plot.fam.its.ses.p[which(plot.fam.its.ses.p$first.birth >= 1834),], aes(x = first.birth, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.6) +
  geom_ribbon(data = plot.fam.its.ses.c[which(plot.fam.its.ses.c$first.birth >= 1834),], aes(x = first.birth, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.6) +
  geom_ribbon(data = plot.fam.its.ses.c[which(plot.fam.its.ses.c$first.birth <= 1834.5),], aes(x = first.birth, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.6) +
  # Lines
  geom_line(data = plot.fam.its.ses.p[which(plot.fam.its.ses.p$first.birth <= 1833), ], aes(x=first.birth, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.fam.its.ses.p[which(plot.fam.its.ses.p$first.birth >= 1834), ], aes(x=first.birth, y = predicted, linetype = ITS), colour = "black") +
  geom_line(data = plot.fam.its.ses.c[which(plot.fam.its.ses.c$first.birth > 1833), ], aes(x=first.birth, y = predicted, linetype = ITS), colour = "black") +
  # Intervention line
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  # Theme and stop y at 100
  theme_minimal() +
  coord_cartesian(ylim=c(0, 15)) +
  # Legend
  scale_linetype_manual(labels = c("Counterfactual", "Predicted values"), values = c("dashed", "solid")) +
  scale_alpha_manual(labels = c("Counterfactual", "Predicted values"), values = c(0.4, 0.4), guide = "none") +
  theme(legend.position = "bottom", legend.title=element_text(face = "bold")) +
  # Remove axis title
  theme(axis.title = element_blank()) +
  #By gender
  facet_wrap(~HCLASS.cat)

its.fam.ses.plot

#saveRDS(its.fam.ses.plot, file = "output/lifespans/its.fam.ses.plot.rds")



# Family size occupation cat GAM ----------------------------------------------------

## Check for overdisperation ----------------------------------------------------------------------
#poisson
#gam.fam.ses.cat1 <- gam(children ~ s(first.birth, bs = "tp", k = -1, by = HCLASS.cat) + 
#                         HCLASS.cat + s(first.birth, bs = "re"), gamma = 1.4, 
#                       family = "poisson", data = parents, method = "REML")

#sum(residuals(gam.fam.ses.cat1, type = "pearson")^2) / df.residual(gam.fam.ses.cat1)
#Overdispersed 1.51

#OLRE
#gam.fam.ses.cat2 <- gam(children ~ s(first.birth, bs = "tp", k = -1, by = HCLASS.cat) + 
#                          HCLASS.cat + s(first.birth, bs = "re")+ s(OLRE, bs = "re"), 
#                        gamma = 1.4, family = "poisson", data = parents, method = "REML")

#sum(residuals(gam.fam.ses.cat2, type = "pearson")^2) / df.residual(gam.fam.ses.cat2)
#Overdispersed 1.39

## Fit model --------------------------------------------------------------------------------
gam.fam.ses.cat <- gam(children ~ s(first.birth, bs = "tp", k = -1, by = HCLASS.cat) + 
                          HCLASS.cat + s(first.birth.factor, bs = "re") + s(OLRE, bs = "re"), 
                        gamma = 1.4, family = "poisson", data = parents, method = "REML")

summary(gam.fam.ses.cat)

#saveRDS(gam.fam.ses.cat, file = "output/lifespans/gam.fam.ses.cat.rds")

## Plot --------------------------------------------------------------------------------
#Predict from model
df.fam.ses <- as.data.frame(ggpredict(gam.fam.ses.cat, c("first.birth[1754:1915]", "HCLASS.cat")))
df.fam.ses <- df.fam.ses %>% rename(first.birth = x, HCLASS.cat = group, pfit = predicted, pse = std.error)

#Get derivative
data.fam.ses <- as.data.frame(df.fam.ses[, c("first.birth", "pfit", "pse", "HCLASS.cat")])
data.fam.ses$HCLASS.cat <- as.factor(data.fam.ses$HCLASS.cat)
data.fam.ses$OLRE <- 1
data.fam.ses$first.birth.factor <- as.factor(as.character(data.fam.ses$first.birth))
colnames(data.fam.ses)
deriv.fam.ses <- derivatives(gam.fam.ses.cat, term = c("s(first.birth)"), partial_match = TRUE, data = data.fam.ses)

#Create HCLASS.cat column column
deriv.fam.ses$HCLASS.cat <- ""
deriv.fam.ses[str_detect(deriv.fam.ses$.smooth, "Manual"), "HCLASS.cat"] <- "Manual"
deriv.fam.ses[str_detect(deriv.fam.ses$.smooth, "Non-Manual"), "HCLASS.cat"] <- "Non-Manual"
deriv.fam.ses[str_detect(deriv.fam.ses$.smooth, "Unknown"), "HCLASS.cat"] <- "Unknown"
#Add in pfit from p.gen
add <- df.fam.ses[, c("first.birth", "HCLASS.cat", "pfit")]
deriv.fam.ses <- merge(deriv.fam.ses, add)
deriv.fam.ses$HCLASS.cat <- as.factor(deriv.fam.ses$HCLASS.cat)


#Make separate deriv dataframes for inc. vs dec.
deriv.fam.ses.inc <- deriv.fam.ses
deriv.fam.ses.inc[!(deriv.fam.ses.inc$.lower_ci > 0 & deriv.fam.ses.inc$.upper_ci > 0), c(".derivative", "pfit")] <- NA
deriv.fam.ses.dec <- deriv.fam.ses
deriv.fam.ses.dec[!(deriv.fam.ses.dec$.lower_ci < 0 & deriv.fam.ses.dec$.upper_ci < 0), c(".derivative", "pfit")] <- NA


#### Overall gam plot with periods of significant change highlighted
gam.fam.ses.plot <- 
  parents %>%
  ggplot (aes (x = first.birth, y =children)) +
  labs (x = "Year of first birth", y= "Children") +
  theme_minimal () +
  # Points
  geom_point(data = parents, aes(x = first.birth, y = children), colour = "grey", shape = 16, alpha = 0.6) +
  # Intervention lines
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  # CIs and lines
  geom_ribbon(data = df.fam.ses, aes(x = first.birth, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.5) +
  geom_line(data = df.fam.ses, aes(x = first.birth, y = pfit)) +
  #Periods of significant change
  geom_line(data = deriv.fam.ses.inc, aes(x = first.birth, y = pfit, colour = "blue"), size = 1) +
  geom_line(data = deriv.fam.ses.dec, aes(x =first.birth, y = pfit, colour = "red"), size = 1) +
  # Legend and axis titles
  scale_color_manual(labels = c("Sig. incr.", "Sig. decr."), values = c("blue", "red"), name = "GAM") +
  theme(legend.position = "bottom", legend.title=element_text(face = "bold")) +
  coord_cartesian(ylim=c(0, 15)) +
  #By gender
  facet_wrap(~HCLASS.cat)
gam.fam.ses.plot


### Plot each group separately ------------------------------------------

# Subset data for each HCLASS.cat level
parents_manual <- subset(parents, HCLASS.cat == "Manual")
parents_non_manual <- subset(parents, HCLASS.cat == "Non-Manual")
parents_unknown <- subset(parents, HCLASS.cat == "Unknown")
df_fam_ses_manual <- subset(df.fam.ses, HCLASS.cat == "Manual")
df_fam_ses_non_manual <- subset(df.fam.ses, HCLASS.cat == "Non-Manual")
df_fam_ses_unknown <- subset(df.fam.ses, HCLASS.cat == "Unknown")
deriv_fam_ses_inc_manual <- subset(deriv.fam.ses.inc, HCLASS.cat == "Manual")
deriv_fam_ses_inc_non_manual <- subset(deriv.fam.ses.inc, HCLASS.cat == "Non-Manual")
deriv_fam_ses_inc_unknown <- subset(deriv.fam.ses.inc, HCLASS.cat == "Unknown")
deriv_fam_ses_dec_manual <- subset(deriv.fam.ses.dec, HCLASS.cat == "Manual")
deriv_fam_ses_dec_non_manual <- subset(deriv.fam.ses.dec, HCLASS.cat == "Non-Manual")
deriv_fam_ses_dec_unknown <- subset(deriv.fam.ses.dec, HCLASS.cat == "Unknown")
deriv_fam_ses_manual <- subset(deriv.fam.ses, HCLASS.cat == "Manual")
deriv_fam_ses_non_manual <- subset(deriv.fam.ses, HCLASS.cat == "Non-Manual")
deriv_fam_ses_unknown <- subset(deriv.fam.ses, HCLASS.cat == "Unknown")


# Plot for "Manual"
gam.fam.ses.man.plot <- ggplot(aes(x = first.birth, y = children), data = parents_manual) +
  labs(x = "Year of first birth", y = "Children") +
  theme_minimal() +
  geom_point(data = parents_manual, aes(x = first.birth, y = children), colour = "grey", shape = 16, alpha = 0.6) +
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  geom_ribbon(data = df_fam_ses_manual, aes(x = first.birth, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.4) +
  geom_line(data = df_fam_ses_manual, aes(x = first.birth, y = pfit)) +
  geom_line(data = deriv_fam_ses_inc_manual, aes(x = first.birth, y = pfit, colour = "blue"), size = 1) +
  geom_line(data = deriv_fam_ses_dec_manual, aes(x = first.birth, y = pfit, colour = "red"), size = 1) +
  scale_color_manual(labels = c("Sig. incr.", "Sig. decr."), values = c("blue", "red"), name = "GAM") +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(0, 15)) +
  xlab(NULL) +
  ylab(NULL)
gam.fam.ses.man.plot
#saveRDS(gam.fam.ses.man.plot, file = "output/lifespans/gam.fam.ses.man.plot.rds")



# Plot for "Non-Manual"
gam.fam.ses.nonMan.plot <- ggplot(aes(x = first.birth, y = children), data = parents_non_manual) +
  labs(x = "Year of first birth", y = "Children") +
  theme_minimal() +
  geom_point(data = parents_non_manual, aes(x = first.birth, y = children), colour = "grey", shape = 16, alpha = 0.6) +
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  geom_ribbon(data = df_fam_ses_non_manual, aes(x = first.birth, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.4) +
  geom_line(data = df_fam_ses_non_manual, aes(x = first.birth, y = pfit)) +
  geom_line(data = deriv_fam_ses_inc_non_manual, aes(x = first.birth, y = pfit, colour = "blue"), size = 1) +
  geom_line(data = deriv_fam_ses_dec_non_manual, aes(x = first.birth, y = pfit, colour = "red"), size = 1) +
  scale_color_manual(labels = c("Sig. incr.", "Sig. decr."), values = c("blue", "red"), name = "GAM") +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(0, 15)) +
  xlab(NULL) +
  ylab(NULL)
gam.fam.ses.nonMan.plot
#saveRDS(gam.fam.ses.nonMan.plot, file = "output/lifespans/gam.fam.ses.nonMan.plot.rds")



# Plot for "Unknown"
gam.fam.ses.un.plot <- ggplot(aes(x = first.birth, y = children), data = parents_unknown) +
  labs(x = "Year of first birth", y = "Children") +
  theme_minimal() +
  geom_point(data = parents_unknown, aes(x = first.birth, y = children), colour = "grey", shape = 16, alpha = 0.6) +
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  geom_ribbon(data = df_fam_ses_unknown, aes(x = first.birth, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.4) +
  geom_line(data = df_fam_ses_unknown, aes(x = first.birth, y = pfit)) +
  geom_line(data = deriv_fam_ses_inc_unknown, aes(x = first.birth, y = pfit, colour = "blue"), size = 1) +
  geom_line(data = deriv_fam_ses_dec_unknown, aes(x = first.birth, y = pfit, colour = "red"), size = 1) +
  scale_color_manual(labels = c("Sig. incr.", "Sig. decr."), values = c("blue", "red"), name = "GAM") +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  coord_cartesian(ylim = c(0, 15)) +
  xlab(NULL) +
  ylab(NULL)

gam.fam.ses.un.plot
#saveRDS(gam.fam.ses.un.plot, file = "output/lifespans/gam.fam.ses.un.plot.rds")







## Plot derivative change --------------------------------------------------------

# Plot for "Manual"
gam.fam.ses.man.deriv.plot <- ggplot(aes(x = first.birth, y = children), data = parents_manual) +
  labs(x = "Year of first birth", y = "Children") +
  theme_minimal() +
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  geom_hline(yintercept = 0) +
  geom_ribbon(data = deriv_fam_ses_manual, aes(x = first.birth, ymin = .lower_ci, ymax = .upper_ci), inherit.aes = FALSE, alpha = 0.3) +
  geom_line(data = deriv_fam_ses_manual, aes(x = first.birth, y = .derivative)) +
  geom_line(data = deriv_fam_ses_inc_manual, aes(x = first.birth, y = .derivative, colour = "blue"), size = 1) +
  geom_line(data = deriv_fam_ses_dec_manual, aes(x = first.birth, y = .derivative, colour = "red"), size = 1) +
  scale_color_manual(labels = c("Sig. incr.", "Sig. decr."), values = c("blue", "red"), name = "GAM") +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  xlab("Year of first birth in family") +
  ylab(NULL) +
  ylim(c(-0.08, 0.05))
gam.fam.ses.man.deriv.plot
#saveRDS(gam.fam.ses.man.deriv.plot, file = "output/lifespans/gam.fam.ses.man.deriv.plot.rds")


# Plot for "Non-Manual"
gam.fam.ses.nonMan.deriv.plot <- ggplot(aes(x = first.birth, y = children), data = parents_non_manual) +
  labs(x = "Year of first birth", y = "Children") +
  theme_minimal() +
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  geom_hline(yintercept = 0) +
  geom_ribbon(data = deriv_fam_ses_non_manual, aes(x = first.birth, ymin = .lower_ci, ymax = .upper_ci), inherit.aes = FALSE, alpha = 0.3) +
  geom_line(data = deriv_fam_ses_non_manual, aes(x = first.birth, y = .derivative)) +
  geom_line(data = deriv_fam_ses_inc_non_manual, aes(x = first.birth, y = .derivative, colour = "blue"), size = 1) +
  geom_line(data = deriv_fam_ses_dec_non_manual, aes(x = first.birth, y = .derivative, colour = "red"), size = 1) +
  scale_color_manual(labels = c("Sig. incr.", "Sig. decr."), values = c("blue", "red"), name = "GAM") +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  xlab("Year of first birth in family") +
  ylab(NULL) +
  ylim(c(-0.08, 0.05))
gam.fam.ses.nonMan.deriv.plot
#saveRDS(gam.fam.ses.nonMan.deriv.plot, file = "output/lifespans/gam.fam.ses.nonMan.deriv.plot.rds")



# Plot for "Unknown
gam.fam.ses.unknown.deriv.plot <- ggplot(aes(x = first.birth, y = children), data = parents_unknown) +
  labs(x = "Year of first birth", y = "Children") +
  theme_minimal() +
  geom_vline(xintercept = 1834, linetype = "dashed", colour = "darkred") +
  geom_vline(xintercept = time.line, linetype = "dashed", colour = "darkred") +
  geom_hline(yintercept = 0) +
  geom_ribbon(data = deriv_fam_ses_unknown, aes(x = first.birth, ymin = .lower_ci, ymax = .upper_ci), inherit.aes = FALSE, alpha = 0.3) +
  geom_line(data = deriv_fam_ses_unknown, aes(x = first.birth, y = .derivative)) +
  geom_line(data = deriv_fam_ses_inc_unknown, aes(x = first.birth, y = .derivative, colour = "blue"), size = 1) +
  geom_line(data = deriv_fam_ses_dec_unknown, aes(x = first.birth, y = .derivative, colour = "red"), size = 1) +
  scale_color_manual(labels = c("Sig. incr.", "Sig. decr."), values = c("blue", "red"), name = "GAM") +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  xlab("Year of first birth in family") +
  ylab(NULL) +
  ylim(c(-0.08, 0.05))
gam.fam.ses.unknown.deriv.plot
#saveRDS(gam.fam.ses.unknown.deriv.plot, file = "output/lifespans/gam.fam.ses.unknown.deriv.plot.rds")
