## Set libraries ##
library(grid)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(moments)
library(nlme)

## read data ##
rg_tab <- read.csv("../02_data/04_final/regression_table.csv")
tt_tab <- read.csv("../02_data/04_final/ttest_table.csv")
rgd <- read.csv("../02_data/04_final/regression_table_discussion.csv")
ttd <- read.csv("../02_data/04_final/ttest_table_discussion.csv")
beta_depth_nf <- read.csv("../02_data/04_final/beta_depth_nf.csv")
or_beta_depth_nf <- read.csv("../02_data/04_final/or_beta_depth_nf.csv")
rg_fd <- read.csv("../02_data/04_final/regression_table_fd.csv")
tt_fd <- read.csv("../02_data/04_final/ttest_table_fd.csv")
fbeta_id <- read.csv("../02_data/04_final/beta_ID_fd.csv")

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

# THETA_POWER
hist(rg_tab$THETA_POWER)
qqnorm(rg_tab$THETA_POWER)
qqline(rg_tab$THETA_POWER)

# ALPHA_POWER
hist(rg_tab$ALPHA_POWER)
qqnorm(rg_tab$ALPHA_POWER)
qqline(rg_tab$ALPHA_POWER)

# BETA_POWER
hist(rg_tab$BETA_POWER)
qqnorm(rg_tab$BETA_POWER)
qqline(rg_tab$BETA_POWER)

# root mean square
hist(rg_tab$root_mean_square)
qqnorm(rg_tab$root_mean_square)
qqline(rg_tab$root_mean_square)

#density plots for all 5
A <- ggdensity(rg_tab, x = "AP_EXPONENT", fill = "lightgray", ylab = "Dichte", xlab = "aperoidischer Exponent") +
  scale_x_continuous(limits = c(-2, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

B <- ggdensity(rg_tab, x = "THETA_POWER", fill = "lightgray", ylab = "Dichte", xlab = "Thetapower (µV)") +
  scale_x_continuous(limits = c(-0.5, 1.5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

C <- ggdensity(rg_tab, x = "ALPHA_POWER", fill = "lightgray", ylab = "Dichte", xlab = "Alphapower (µV)") +
  scale_x_continuous(limits = c(-1, 1.5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

D <- ggdensity(rg_tab, x = "BETA_POWER", fill = "lightgray", ylab = "Dichte", xlab = "Betapower (µV)") +
  scale_x_continuous(limits = c(-0.5, 1)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

E <- ggdensity(rg_tab, x = "root_mean_square", fill = "lightgray", ylab = "Dichte", xlab = "Quadratisches Mittel") +
  scale_x_continuous(limits = c(-2, 30)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

ggarrange(A,B,C,D,E,ncol = 2, nrow = 3)

## visualizing z data ##

# AP_EXPONENT
hist(rg_tab$z_exp)
qqnorm(rg_tab$z_exp)
qqline(rg_tab$z_exp)

# THETA_POWER
hist(rg_tab$z_theta)
qqnorm(rg_tab$z_theta)
qqline(rg_tab$z_theta)

# ALPHA_POWER
hist(rg_tab$z_alpha)
qqnorm(rg_tab$z_alpha)
qqline(rg_tab$z_alpha)

# BETA_POWER
hist(rg_tab$z_beta)
qqnorm(rg_tab$z_beta)
qqline(rg_tab$z_beta)

# root mean square
hist(rg_tab$z_rms)
qqnorm(rg_tab$z_rms)
qqline(rg_tab$z_rms)

# density plots for z-transformed variables
Az <- ggdensity(rg_tab, x = "z_exp", fill = "lightgray", ylab = "Dichte", xlab = "aperoidischer Exponent") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

Bz <- ggdensity(rg_tab, x = "z_theta", fill = "lightgray", ylab = "Dichte", xlab = "Thetapower (µV)") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

Cz <- ggdensity(rg_tab, x = "z_alpha", fill = "lightgray", ylab = "Dichte", xlab = "Alphapower (µV)") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

Dz <- ggdensity(rg_tab, x = "z_beta", fill = "lightgray", ylab = "Dichte", xlab = "Betapower (µV)") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

Ez <- ggdensity(rg_tab, x = "z_rms", fill = "lightgray", ylab = "Dichte", xlab = "Quadratisches Mittel") +
  scale_x_continuous(limits = c(-5, 5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

ggarrange(Az,Bz,Cz,Dz,Ez,ncol = 2, nrow = 3)

### That made it better and good, exept for rms


## try log-transformation ##
rg_tab$l_theta <- log10(rg_tab$THETA_POWER)
rg_tab$l_alpha <- log10(rg_tab$ALPHA_POWER)
rg_tab$l_beta <- log10(rg_tab$BETA_POWER)


## visualizing log data ##

# THETA_POWER
Bl <- ggdensity(rg_tab, x = "l_theta", fill = "lightgray", ylab = "Dichte", xlab = "Thetapower (µV)") +
  scale_x_continuous(limits = c(-3, 1.5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

# ALPHA_POWER
Cl <- ggdensity(rg_tab, x = "l_alpha", fill = "lightgray", ylab = "Dichte", xlab = "Alphapower (µV)") +
  scale_x_continuous(limits = c(-3, 1.5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

# BETA_POWER
Dl <- ggdensity(rg_tab, x = "l_beta", fill = "lightgray", ylab = "Dichte", xlab = "Betapower (µV)") +
  scale_x_continuous(limits = c(-3, 1.5)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")

ggarrange(Bl,Cl,Dl,ncol = 2, nrow = 2)

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
#boxplot 
beta_col1 <- data.frame(tt_tab$near_beta,c(rep('"nah"',30)))
colnames(beta_col1) <- c("Betapower", "Bedingung")
beta_col2 <- data.frame(tt_tab$far_beta,c(rep('"fern"',30)))
colnames(beta_col2) <- c("Betapower", "Bedingung")
t_beta_plot <- rbind(beta_col1,beta_col2)
ggplot(t_beta_plot, aes(x=Bedingung, y=Betapower)) + 
  geom_boxplot() + labs(y = "Betapower (µV)")

ttest_rms <- t.test(tt_tab$near_rms, tt_tab$far_rms, paired = T, "greater")
ttest_rms
#boxplot 
rms_col1 <- data.frame(tt_tab$near_rms,c(rep('"nah"',30)))
colnames(rms_col1) <- c("Quadratisches_Mittel", "Bedingung")
rms_col2 <- data.frame(tt_tab$far_rms,c(rep('"fern"',30)))
colnames(rms_col2) <- c("Quadratisches_Mittel", "Bedingung")
t_rms_plot <- rbind(rms_col1,rms_col2)
ggplot(t_rms_plot, aes(x=Bedingung, y=Quadratisches_Mittel)) + 
  geom_boxplot() + labs(y = "Quadratisches Mittel")

######## H2 #########
## compare aperiodic exponent near target with far from target

# check normal distribution of difference for paired t-test
dif_exp <- tt_tab$far_exp - tt_tab$near_exp

shapiro.test(dif_exp)

# we can assume normal distribution
ttest_exp <- t.test(tt_tab$near_exp, tt_tab$far_exp, paired = T, "two.sided")
ttest_exp
#boxplot 
exp_col1 <- data.frame(tt_tab$near_exp,c(rep('"nah"',30)))
colnames(exp_col1) <- c("Aperiodischer_Exponent", "Bedingung")
exp_col2 <- data.frame(tt_tab$far_exp,c(rep('"fern"',30)))
colnames(exp_col2) <- c("Aperiodischer_Exponent", "Bedingung")
t_exp_plot <- rbind(exp_col1,exp_col2)
ggplot(t_exp_plot, aes(x=Bedingung, y=Aperiodischer_Exponent)) + 
  geom_boxplot() + labs(y = "Aperiodischer Exponent")

## correlation between Depth and aperiodic exponent

# is depth normal distributed?
ggdensity(rg_tab, x = "DEPTH", fill = "lightgray", xlab = "Tiefe (mm)", ylab = "Dichte") +
  scale_x_continuous(limits = c(-20, 20)) +
  stat_overlay_normal_density(color = "red", linetype = "dashed")
# we can't assume normal distribution for depth -> kendall correlation which is non-parametric

# scatter plot
ggscatter(rg_tab, x = "DEPTH", y = "z_exp", add = "reg.line", 
          conf.int = T, xlab = "Tiefe (mm)", ylab = "Aperiodischer Exponent")

# correlation test
r <- cor.test(rg_tab$DEPTH,rg_tab$z_exp, "two.sided", "kendall")
r

######## Exploration ######### 
## correlation matrix

# kendall correlation which is non-parametric (spearman has issues with ties)
short <- data.frame(rg_tab$DEPTH, rg_tab$z_exp, rg_tab$z_theta, rg_tab$z_alpha, rg_tab$z_beta, rg_tab$z_rms)
cortab <- round(cor(short, method = "kendall"),3)
rownames(cortab) <- c("Tiefe", "Aperiodischer Exponent", "Thetapower", "Alphapower", "Betapower", "Quadratisches Mittel")
colnames(cortab) <- c("Tiefe", "Aperiodischer Exponent", "Thetapower", "Alphapower", "Betapower", "Quadratisches Mittel")
cortab

#visualize found correlations

ggscatter(rg_tab, x = "DEPTH", y = "z_rms", add = "reg.line", 
          conf.int = T, xlab = "Tiefe (mm)", ylab = "Quadratisches Mittel")

ggscatter(rg_tab, x = "DEPTH", y = "z_theta", add = "reg.line", 
          conf.int = T, xlab = "Tiefe (mm)", ylab = "Thetapower (µV)")

ggscatter(rg_tab, x = "DEPTH", y = "z_alpha", add = "reg.line", 
          conf.int = T, xlab = "Tiefe (mm)", ylab = "Alphapower (µV)")


# regression is robust (source), data is close to normal distribution
full_model_depth <- lme(fixed=DEPTH ~ z_exp + z_rms + z_alpha + z_theta, random=~1|ID, data=rg_tab)
summary(full_model_depth)
anova(full_model_depth)

## t-test for other variables
# normal distribution?

dif_theta <- tt_tab$far_theta - tt_tab$near_theta
dif_alpha <- tt_tab$far_alpha - tt_tab$near_alpha

shapiro.test(dif_theta)
shapiro.test(dif_alpha)

# we can assume normal distribution
ttest_theta <- t.test(tt_tab$near_theta, tt_tab$far_theta, paired = T, "greater")
ttest_theta

ttest_alpha <- t.test(tt_tab$near_alpha, tt_tab$far_alpha, paired = T, "less")
ttest_alpha


####### discussion #########

## correlation between beta and depth but for depth < 4
cor(rg_tab$DEPTH[rg_tab$DEPTH < 4],rg_tab$z_beta[rg_tab$DEPTH < 4],method =  "kendall")

## t-tests near target vs far target for low-beta and high beta
# low-beta
dif_lbeta <- tt_tab$far_lbeta - tt_tab$near_lbeta
shapiro.test(dif_lbeta)
ttest_lbeta <- t.test(tt_tab$near_lbeta, tt_tab$far_lbeta, paired = T, "greater")
ttest_lbeta #not significant

# high-beta
dif_hbeta <- tt_tab$far_hbeta - tt_tab$near_hbeta
shapiro.test(dif_hbeta)
ttest_hbeta <- t.test(tt_tab$near_hbeta, tt_tab$far_hbeta, paired = T, "greater")
ttest_hbeta # not significant

## correlation between depth and low-beta/ hight-beta
# low-beta
cor(rg_tab$DEPTH, rg_tab$z_lbeta, method = "kendall")
# high-beta
cor(rg_tab$DEPTH, rg_tab$z_hbeta, method = "kendall")

## same for depth < 4
# low-beta
cor(rg_tab$DEPTH[rg_tab$DEPTH < 4], rg_tab$z_lbeta[rg_tab$DEPTH < 4], method = "kendall")
cor.test(rg_tab$DEPTH[rg_tab$DEPTH < 4], rg_tab$z_lbeta[rg_tab$DEPTH < 4],  "less" ,"kendall")
# high-beta
cor(rg_tab$DEPTH[rg_tab$DEPTH < 4], rg_tab$z_hbeta[rg_tab$DEPTH < 4], method = "kendall")

## compare beta near target with far from target again but with original powerspectrum
# check normal distribution of difference for paired t-test
dif_beta_d <- ttd$far_beta - ttd$near_beta
shapiro.test(dif_beta_d)
# we can assume normal distribution
ttest_beta_d <- t.test(ttd$near_beta, ttd$far_beta, paired = T, "greater")
ttest_beta_d #significant

#boxplot 
d_beta_col1 <- data.frame(ttd$near_beta,c(rep('"nah"',30)))
colnames(d_beta_col1) <- c("Betapower", "Bedingung")
d_beta_col2 <- data.frame(ttd$far_beta,c(rep('"fern"',30)))
colnames(d_beta_col2) <- c("Betapower", "Bedingung")
t_d_beta_plot <- rbind(d_beta_col1,d_beta_col2)
ggplot(t_d_beta_plot, aes(x=Bedingung, y=Betapower)) + 
  geom_boxplot() + labs(y = "Betapower (µV)")

# compare depth of channels near target for original powerspectrum vs powerspectrum after fooof
par(mfrow=c(2,2))
plot(or_beta_depth_nf$near_beta, ttd$near_beta, xlab = 'Tiefe "nah" (mm)', ylab = "Betapower (µV)", main = 'Originales Powerspektrum')
plot(beta_depth_nf$near_beta, tt_tab$near_beta, xlab = 'Tiefe "nah" (mm)', ylab = "Betapower (µV)", main = 'Powerspektrum nach FOOOF-Algorithmus')
plot(or_beta_depth_nf$far_beta, ttd$far_beta, xlab = 'Tiefe "fern" (mm)', ylab = "Betapower (µV)", main = 'Originales Powerspektrum')
plot(beta_depth_nf$far_beta, tt_tab$far_beta, xlab = 'Tiefe "fern" (mm)', ylab = "Betapower (µV)", main = 'Powerspektrum nach FOOOF-Algorithmus')
# compare power of channels near target for original powerspectrum vs powerspectrum after fooof
mean(ttd$near_beta)
mean(tt_tab$near_beta)
mean(ttd$far_beta)
mean(tt_tab$far_beta)

## t-tests near target vs far target for beta without aperiodic component, but without cleaning of data
dif_fbeta <- tt_fd$far_beta - tt_fd$near_beta
shapiro.test(dif_fbeta)
ttest_fbeta <- t.test(tt_fd$near_beta, tt_fd$far_beta, paired = T, "greater")
ttest_fbeta 

mean(tt_fd$near_beta)
mean(tt_fd$far_beta)
