

data {
  int<lower=1> N;        // Number of data points
  int<lower=1> Tr;       // Number of treatments
  int<lower=1> Pl;       // Number of plot ids
  int<lower=1> treatment[N];  // Treatment group for each observation
  int<lower=1> plot_id[N];    // Plot ID for each observation
  int<lower=1> day_year[N];   // Day of the year for each observation
  vector[N] temp;        // Temperature data
  vector[N] resp;        // Respiration rates
  vector[N] M;           // Soil moisture (volumetric, cm3 cm^-3)
}

transformed data {
  real T_0 = 227.13;  // Lloyd-Taylor parameter
  real sigma =  0.2712706;  // Total SD of the dataset
}

parameters {
  // real<lower=0> T_0;

  vector<lower=0>[Pl] A;           // Log of the pre-exponential factor
  vector<lower=0>[Pl] Ea;          // Ea/R for each treatment

  vector<lower=3.11*0.25, upper=3.11*1.75>[Tr] a;  // Scaling of Moyano function
  vector<lower=2.42*0.25, upper=2.42*1.75>[Tr] b;  // Scaling of Moyano function

  vector<lower= -0.5, upper=0.5>[Tr] amplitude;   // Amplitude for each treatment
  vector<lower=0, upper=365>[Tr] peak_day;  // Peak day for each treatment

}

model {
  vector[N] model_resp;
  vector[N] xi_moist;
  vector[N] xi_temp;
  vector[N] sine_wave;


  // Model prediction
  for (i in 1:N) {
    xi_moist[i] =   a[treatment[i]] * M[i] - b[treatment[i]] * M[i]^2;
    xi_temp[i] = A[plot_id[i]]*(exp(-Ea[plot_id[i]] / ((temp[i] + 273.15) - T_0)));
    sine_wave[i] = amplitude[treatment[i]] * cos((2 * pi() / 365) * day_year[i] + (2 * pi() / 365) * (peak_day[treatment[i]] - 1) - pi() / 2);
    model_resp[i] = sine_wave[i] +  (xi_temp[i] * xi_moist[i]);
    // model_resp[i] = (xi_temp[i] * xi_moist[i]);
    }

  // Priors
  // T_0 ~ normal(227.13, 10); //now this is also calibrated but general

  A ~ normal(400,100); #this prior is hard to find
  Ea ~ normal(398.5, 50);
  a ~ normal(3.11, 0.25);
  b ~ normal(2.42, 0.25);
  amplitude ~ normal(0, 0.25);
  peak_day ~ normal(196, 50); #centered on the 15th of July

  // Likelihood
  resp ~ normal(model_resp, sigma);
}



generated quantities {
  vector[N] model_resp;
  vector[N] xi_moist;
  vector[N] xi_temp;
  vector[N] sine_wave;
  vector[N] res;

  // Model prediction
  for (i in 1:N) {
    xi_moist[i] =   a[treatment[i]] * M[i] - b[treatment[i]] * M[i]^2;
    xi_temp[i] = A[plot_id[i]]*(exp(-Ea[plot_id[i]] / ((temp[i] + 273.15) - T_0)));
    sine_wave[i] = amplitude[treatment[i]] * cos((2 * pi() / 365) * day_year[i] + (2 * pi() / 365) * (peak_day[treatment[i]] - 1) - pi() / 2);
    model_resp[i] = sine_wave[i] +  (xi_temp[i] * xi_moist[i]);
    // model_resp[i] =   (xi_temp[i] * xi_moist[i]);
    res[i] = model_resp[i] - resp[i];
  }
}
