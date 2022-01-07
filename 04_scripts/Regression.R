## Set libraries ##
library(ggplot2)
library(ggpubr)
library(moments)
library(nlme)

## read data ##
rg_tab <- read.csv("C:/Users/acer/Desktop/Master CAU/02_data/04_final/regression_table.csv")


## Add new variable for distance to STN ##
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


## try z-transformation ##
for(i in 1:30) {
  idx <-  rg_tab$ID == i
  rg_tab$z_exp[idx] <- scale(rg_tab$AP_EXPONENT[idx])
  rg_tab$z_theta[idx] <- scale(rg_tab$THETA_POWER[idx])
  rg_tab$z_alpha[idx] <- scale(rg_tab$ALPHA_POWER[idx])
  rg_tab$z_beta[idx] <- scale(rg_tab$BETA_POWER[idx])
}


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

### That made it better and good


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


## correlation between distance and beta ##

# we can't assume normal distribution for distance
hist(rg_tab$distance)
qqnorm(rg_tab$distance)
qqline(rg_tab$distance)

ggdensity(rg_tab, x = "distance", fill = "lightgray", title = "Distance to STN") +
  scale_x_continuous(limits = c(-20, 20)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

#visualize data
ggscatter(rg_tab, x = "distance", y = "z_beta", add = "reg.line", 
          conf.int = T, cor.coef = T, cor.method = "spearman", xlab = "distance to STN", ylab = "beta power")

# kendall correlation which is non-parametric (spearman has issues with ties)
H <- cor.test(rg_tab$distance,rg_tab$z_beta, method = "kendall")
H
### no significant correlation between beta and distance to STN 

## exploring ##
# regression is robust (source), data is close to normal distribution
full_model <- lme(fixed=DEPTH ~ z_exp + z_theta + z_alpha + z_beta, random=~1|ID, data=rg_tab)
summary(full_model)
anova(full_model)
coef(full_model)