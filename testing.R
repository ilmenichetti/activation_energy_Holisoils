
dim(post_bytreat_indA_Moy$amplitude)
dim(post_bytreat_indA_Moy$peak_day)
Ea = mean(post_bytreat_indA_Moy$Ea)
A = mean(post_bytreat_indA_Moy$A)
a = mean(post_bytreat_indA_Moy$a)
b = mean(post_bytreat_indA_Moy$b)
T_0= 227.13

moist = validation_data$Soil_moist
temp =  validation_data$T1_soil

xi_moist =   a * M - b * moist^2;
xi_temp= A*(exp(-Ea / ((temp + 273.15) - T_0)));

model_resp = (xi_temp * xi_moist);


plot(model_resp)
lines(validation_data$CO2_flux_norm, col="red")

plot(model_resp, validation_data$CO2_flux_norm)
abline(a=0, b=1, col="black", lty=2)




plot(xi_moist)

lines(validation_data$CO2_flux_hour, col="red")

