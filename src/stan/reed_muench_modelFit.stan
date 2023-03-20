functions{
  real logistic(real a, real b, real concentration){
    real prob = 1 / (1 + exp(-(a + b * concentration)));
    return prob;
  }
}

data {
  int N;
  int nreplicates[N];
  int survival[N];
  int dilution[N];
  int is_log;
  
  // posterior predictive
  int ndilutionssim;
  vector[ndilutionssim] dilutionssim;
}

parameters {
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> phi;
}

model {
  for(i in 1:N) {
    real concentration;
    if(is_log == 1)
      concentration = log(phi / dilution[i]);
    else
      concentration = phi / dilution[i];
      
    survival[i] ~ binomial(nreplicates[i], logistic(a, b, concentration));
  }
  
  // priors
  a ~ cauchy(0, 10);
  b ~ cauchy(0, 10);
  phi ~ lognormal(2, 0.5);
}

generated quantities {
  vector[ndilutionssim] prob_survive;
  real ed50 = phi * exp(-a)^(-1.0 / b);
  {
    real concentration;
    for(i in 1:ndilutionssim) {
      if(is_log == 1)
        concentration = log(phi / dilutionssim[i]);
      else
        concentration = phi / dilutionssim[i];
        prob_survive[i] = logistic(a, b, concentration);
    }
  }
}
