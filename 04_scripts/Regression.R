## read data
rg_tab <- read.csv("C:/Users/acer/Desktop/Master CAU/02_data/04_final/regression_table.csv")

## Add new variable for distance to STN
rg_tab$distance <- abs(rg_tab$DEPTH)

## Test for normal distribution is bad: gets significant automatically with big samplesize
shapiro.test(rg_tab$AP_EXPONENT)

## visualizing data
hist(rg_tab$AP_EXPONENT)
qqnorm(rg_tab$AP_EXPONENT)
qqline(rg_tab$AP_EXPONENT)

hist(rg_tab$THETA_POWER)
qqnorm(rg_tab$THETA_POWER)
qqline(rg_tab$THETA_POWER)

hist(rg_tab$ALPHA_POWER)
qqnorm(rg_tab$ALPHA_POWER)
qqline(rg_tab$ALPHA_POWER)

hist(rg_tab$BETA_POWER)
qqnorm(rg_tab$BETA_POWER)
qqline(rg_tab$BETA_POWER)

## try z-transformation
for(i in 1:30) {
  idx <-  rg_tab$ID == i
  rg_tab$z_exp[idx] <- scale(rg_tab$AP_EXPONENT[idx])
  rg_tab$z_theta[idx] <- scale(rg_tab$THETA_POWER[idx])
  rg_tab$z_alpha[idx] <- scale(rg_tab$ALPHA_POWER[idx])
  rg_tab$z_beta[idx] <- scale(rg_tab$BETA_POWER[idx])
}

## visualizing z data
hist(rg_tab$z_exp)
qqnorm(rg_tab$z_exp)
qqline(rg_tab$z_exp)
shapiro.test(rg_tab$z_exp)

hist(rg_tab$z_theta)
qqnorm(rg_tab$z_theta)
qqline(rg_tab$z_theta)

hist(rg_tab$z_alpha)
qqnorm(rg_tab$z_alpha)
qqline(rg_tab$z_alpha)

hist(rg_tab$z_beta)
qqnorm(rg_tab$z_beta)
qqline(rg_tab$z_beta)

# try transforming exept ap_exp
rg_tab$t_exp <- sqrt(max(rg_tab$z_exp+1) - rg_tab$z_exp) # square-root/ negatively skewed data
rg_tab$t_theta <- 1/(rg_tab$z_theta+4) # inverse/ positively skewed data
rg_tab$t_alpha <- 1/(rg_tab$z_alpha+6) # inverse/ positively skewed data
rg_tab$t_beta <- log10(rg_tab$z_beta+4) # log/ positively skewed data

## visualizing log data
hist(rg_tab$t_exp)
qqnorm(rg_tab$t_exp)
qqline(rg_tab$t_exp)
shapiro.test(rg_tab$t_exp)

hist(rg_tab$t_theta)
qqnorm(rg_tab$t_theta)
qqline(rg_tab$t_theta)
shapiro.test(rg_tab$t_theta)

hist(rg_tab$t_alpha)
qqnorm(rg_tab$t_alpha)
qqline(rg_tab$t_alpha)
shapiro.test(rg_tab$t_alpha)

hist(rg_tab$t_beta)
qqnorm(rg_tab$t_beta)
qqline(rg_tab$t_beta)
shapiro.test(rg_tab$t_beta)

# correlation between distance and beta

# Regression
testmodel <- lm(DEPTH~AP_EXPONENT, data = rg_tab)
summary(testmodel)