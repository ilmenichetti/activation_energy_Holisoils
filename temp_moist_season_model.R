
temp_moist_season <- function(temp, M, day_year, a, b, Ea, A, peak_day, amplitude, treatment, plot_id){

  N = length(temp)
  xi_moist <- c()
  xi_temp <- c()
  sine_wave <- c()
  model_resp <- c()

  T_0 = 227.13

  for (i in 1:N) {
    xi_moist[i] =   a[treatment[i]] * M[i] - b[treatment[i]] * M[i]^2;
    xi_temp[i] = A[plot_id[i]] * (exp(-Ea[treatment[i]] / ((temp[i] + 273.15) - T_0)));
    sine_wave[i] = amplitude[treatment[i]] * cos((2 * pi / 365) * day_year[i] + (2 * pi / 365) * (peak_day[treatment[i]] - 1) - pi / 2);
    model_resp[i] = sine_wave[i] +  (xi_temp[i] * xi_moist[i]);
  }

return(model_resp)

}

