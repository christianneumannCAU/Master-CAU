## Set libraries ##
library(ggplot2)
library(ggpubr)
library(moments)
library(nlme)

## read data ##
rg_tab <- read.csv("../02_data/04_final/regression_table.csv")
tt_tab <- read.csv("../02_data/04_final/ttest_table.csv")

## Add new variable for distance to target (got determined by MRT beforehand - the Depth 0 is the target)##
rg_tab$distance <- abs(rg_tab$DEPTH)


## Assign data types ##
rg_tab$ID <- as.factor(rg_tab$ID)
rg_tab$SIDE <- as.factor(rg_tab$SIDE)
rg_tab$CHANNEL <- as.factor(rg_tab$CHANNEL)


## Test for normal distribution is bad: gets significant automatically with big samplesize ##
shapiro.test(rg_tab$AP_EXPONENT)


## visualizing data ##

# AP_EXPONENT
hist(rg_tab$AP_EXPONENT)
qqnorm(rg_tab$AP_EXPONENT)
qqline(rg_tab$AP_EXPONENT)

ggdensity(rg_tab, x = "AP_EXPONENT", fill = "lightgray", title = "AP Exponent") +
  scale_x_continuous(limits = c(-2, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

# THETA_POWER
hist(rg_tab$THETA_POWER)
qqnorm(rg_tab$THETA_POWER)
qqline(rg_tab$THETA_POWER)

ggdensity(rg_tab, x = "THETA_POWER", fill = "lightgray", title = "Theta") +
  scale_x_continuous(limits = c(-1, 2)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

# ALPHA_POWER
hist(rg_tab$ALPHA_POWER)
qqnorm(rg_tab$ALPHA_POWER)
qqline(rg_tab$ALPHA_POWER)

ggdensity(rg_tab, x = "ALPHA_POWER", fill = "lightgray", title = "Alpha") +
  scale_x_continuous(limits = c(-1, 2)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

# BETA_POWER
hist(rg_tab$BETA_POWER)
qqnorm(rg_tab$BETA_POWER)
qqline(rg_tab$BETA_POWER)

ggdensity(rg_tab, x = "BETA_POWER", fill = "lightgray", title = "Beta") +
  scale_x_continuous(limits = c(-1, 2)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

# root mean square
hist(rg_tab$root_mean_square)
qqnorm(rg_tab$root_mean_square)
qqline(rg_tab$root_mean_square)

ggdensity(rg_tab, x = "root_mean_square", fill = "lightgray", title = "root_mean_square") +
  scale_x_continuous(limits = c(-2, 30)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")


## visualizing z data ##

# AP_EXPONENT
hist(rg_tab$z_exp)
qqnorm(rg_tab$z_exp)
qqline(rg_tab$z_exp)

ggdensity(rg_tab, x = "z_exp", fill = "lightgray", title = "AP Exponent") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

# THETA_POWER
hist(rg_tab$z_theta)
qqnorm(rg_tab$z_theta)
qqline(rg_tab$z_theta)

ggdensity(rg_tab, x = "z_theta", fill = "lightgray", title = "Theta") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

# ALPHA_POWER
hist(rg_tab$z_alpha)
qqnorm(rg_tab$z_alpha)
qqline(rg_tab$z_alpha)

ggdensity(rg_tab, x = "z_alpha", fill = "lightgray", title = "Alpha") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

# BETA_POWER
hist(rg_tab$z_beta)
qqnorm(rg_tab$z_beta)
qqline(rg_tab$z_beta)

ggdensity(rg_tab, x = "z_beta", fill = "lightgray", title = "Beta") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

# root mean square
hist(rg_tab$z_rms)
qqnorm(rg_tab$z_rms)
qqline(rg_tab$z_rms)

ggdensity(rg_tab, x = "z_rms", fill = "lightgray", title = "root mean square") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")


### That made it better and good, exept for rms


## try log-transformation ##
rg_tab$l_theta <- log10(rg_tab$THETA_POWER)
rg_tab$l_alpha <- log10(rg_tab$ALPHA_POWER)
rg_tab$l_beta <- log10(rg_tab$BETA_POWER)


## visualizing log data ##

# THETA_POWER
ggdensity(rg_tab, x = "l_theta", fill = "lightgray", title = "Theta") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

# ALPHA_POWER
ggdensity(rg_tab, x = "l_alpha", fill = "lightgray", title = "Alpha") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

# BETA_POWER
ggdensity(rg_tab, x = "l_beta", fill = "lightgray", title = "Beta") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

### that made it worse


######## H1 #########
## compare rms and beta near target with far from target

# check normal distribution of difference for paired t-test
dif_beta <- tt_tab$far_beta - tt_tab$near_beta
dif_rms <- tt_tab$far_rms - tt_tab$near_rms

shapiro.test(dif_beta)
shapiro.test(dif_rms)

# we can assume normal distribution
ttest_beta <- t.test(tt_tab$near_beta, tt_tab$far_beta, paired = T, "greater")
ttest_beta
ttest_rms <- t.test(tt_tab$near_rms, tt_tab$far_rms, paired = T, "greater")
ttest_rms

######## H2 #########
## compare aperiodic exponent near target with far from target

# check normal distribution of difference for paired t-test
dif_exp <- tt_tab$far_exp - tt_tab$near_exp

shapiro.test(dif_exp)

# we can assume normal distribution
ttest_exp <- t.test(tt_tab$near_exp, tt_tab$far_exp, paired = T, "two.sided")
ttest_exp

## correlation between Depth and aperiodic exponent

# is depth normal distributed?
ggdensity(rg_tab, x = "DEPTH", fill = "lightgray", title = "Depth of electrode") +
  scale_x_continuous(limits = c(-20, 20)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")
# we can't assume normal distribution for depth -> kendall correlation which is non-parametric

# scatter plot
ggscatter(rg_tab, x = "DEPTH", y = "z_exp", add = "reg.line", 
          conf.int = T, cor.coef = T, cor.method = "kendall", xlab = "Depth of electrode", ylab = "aperiodic exponent")

# correlation test
r <- cor.test(rg_tab$DEPTH,rg_tab$z_exp, "two.sided", "kendall")
r

######## Exploration ######### 
## correlation matrix

# kendall correlation which is non-parametric (spearman has issues with ties)
short <- data.frame(rg_tab$distance, rg_tab$DEPTH, rg_tab$z_exp, rg_tab$z_theta, rg_tab$z_alpha, rg_tab$z_beta, rg_tab$z_rms)
cortab <- cor(short, method = "kendall")
cortab

#visualize found correlations

ggscatter(rg_tab, x = "DEPTH", y = "z_rms", add = "reg.line", 
          conf.int = T, cor.coef = T, cor.method = "kendall", xlab = "Depth of electrode", ylab = "root mean square")

ggscatter(rg_tab, x = "DEPTH", y = "z_theta", add = "reg.line", 
          conf.int = T, cor.coef = T, cor.method = "kendall", xlab = "Depth of electrode", ylab = "theta power")

ggscatter(rg_tab, x = "DEPTH", y = "z_alpha", add = "reg.line", 
          conf.int = T, cor.coef = T, cor.method = "kendall", xlab = "Depth of electrode", ylab = "alpha power")

# regression is robust (source), data is close to normal distribution
full_model_depth <- lme(fixed=DEPTH ~ z_exp + z_rms + z_alpha + z_theta, random=~1|ID, data=rg_tab)
summary(full_model_depth)
anova(full_model_depth)

####### discussion #########
## t-test for other variables
# normal distribution?

dif_theta <- tt_tab$far_theta - tt_tab$near_theta
dif_alpha <- tt_tab$far_alpha - tt_tab$near_alpha

shapiro.test(dif_theta)
shapiro.test(dif_alpha)

# we can assume normal distribution
ttest_theta <- t.test(tt_tab$near_theta, tt_tab$far_theta, paired = T, "less")
ttest_alpha <- t.test(tt_tab$near_alpha, tt_tab$far_alpha, paired = T, "greater")

# what if we distinct between low and high beta?
# normal distribution?
dif_lbeta <- tt_tab$far_lbeta - tt_tab$near_lbeta
dif_hbeta <- tt_tab$far_hbeta - tt_tab$near_hbeta

shapiro.test(dif_lbeta)
shapiro.test(dif_hbeta)

ttest_lbeta <- t.test(tt_tab$near_lbeta, tt_tab$far_lbeta, paired = T, "greater")
ttest_hbeta <- t.test(tt_tab$near_hbeta, tt_tab$far_hbeta, paired = T, "greater")