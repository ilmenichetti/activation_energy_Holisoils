
beta1 = 0.08;
beta2 = 0.3;
beta3 = 0.3;
beta4 = 0.01;
beta5 = 0.04;
beta6 = 0.02;
beta7 = 0.06;
beta8 = 0.09;

M=processed_data_filtered$Soil_moist[1]
BD=1.2
Clay=0.2
SOC=0.03

xi_moist  =  beta1*M + beta2*M^2 + beta3*M^3 + beta4*BD + beta5*M*BD + beta6*Clay + beta7*M*Clay + beta8*SOC

Ea = 400
logA=1

model_resp_temp = exp(logA - Ea / (R * temp))

temp=15
Ea = 300
T_0 =seq(from=210.13, to=220, by=0.1)
A = 1
model_resp_temp = A*exp(-Ea/((temp+273.15)-T_0))
plot(T_0, model_resp_temp)



k_b = 1.380649*10^-23 #kJ K^-1
h=6.62607015*10^-34 #kJ S^-1
R=8.314462618 # mol^-1 K^-1
Dg_cat=51*10^3 #kJ
DH_eq=92*10^3 #kJ
T_eq=25+273 #K
temp=seq(-20:100)+273 #Kß
E_0=1

K_cat = ((k_b*temp)/h)*exp(-(Dg_cat/(R*temp)));
plot(temp-273, K_cat)

K_eq = exp((DH_eq/R)*((1/T_eq)-(1/temp)));

model_resp_temp = (K_cat * E_0)/(1+K_eq);

plot(temp-273, model_resp_temp)
plot(temp-273, K_cat)
plot(temp-273, K_eq)




Dg_cat=51*10^3 #kJ
DH_eq=92*10^3 #kJ
T_eq=25+273 #K
temp=seq(-20:100)+273 #Kß
E_0=1
A=0.0001

K_cat = ((k_b*temp)/h)*exp(-(Dg_cat/(R*temp)));
K_eq = exp((DH_eq/R)*((1/T_eq)-(1/temp)));

model_resp_temp = A* (K_cat * E_0)/(1+K_eq);

plot(temp-273, model_resp_temp)


# Constants
h <- 6.62607015e-34  # Planck's constant, Joule*seconds
kB <- 1.380649e-23   # Boltzmann's constant, Joule/Kelvin
R <- 8.314           # Universal gas constant, J/(mol*K)

# Activation parameters (example values)
deltaH_dagger <- 200000  # Activation enthalpy, J/mol (example value)
deltaS_dagger <- 100     # Activation entropy, J/(mol*K) (example value)

# Temperature range: 250K to 350K
temperature <- seq(250, 350, 10)

# Calculate Gibbs free energy of activation
deltaG_dagger <- deltaH_dagger - temperature * deltaS_dagger

# Eyring-Polanyi equation using Gibbs free energy to calculate rate constant k
k <- (kB * temperature / h) * exp(-deltaG_dagger / (R * temperature))

# Output the rate constants
plot(temperature, k)

x=seq(from=0, to=1, by=0.01)
logistic <- function(L, k, x0, x){
  L/(1+exp(-k*(x-x0)))
}
plot(x, logistic(L=1, k=10, x0=0.5, x=x), type="l", col="red")
lines(x, logistic(L=1, k=-10, x0=0.5, x=x), type="l", col="blue")
legend("topleft", c("L=1, k=10, x0=0.5", "L=1, k=-10, x0=0.5"), col=c("red", "blue"), lty=1)
