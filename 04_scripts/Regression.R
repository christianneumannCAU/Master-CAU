rg_tab <- read.csv("C:/Users/acer/Desktop/Master CAU/02_data/04_final/regression_table.csv")

testmodel <- lm(DEPTH~AP_EXPONENT, data = rg_tab)

summary(testmodel)
