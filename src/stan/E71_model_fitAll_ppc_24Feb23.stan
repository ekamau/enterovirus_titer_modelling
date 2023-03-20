// each individual phi estimate

functions{
  vector logistic(real a, real b, real c, int N_len, vector concentration){
    vector[N_len] log_denom = c * log(1 + exp(-(a + b * concentration)));
    return(1 ./ exp(log_denom));
  }
}

data {
  int N_panels;
  int phi_panel_index[N_panels];
  int N;
  int N_sim;
  int nreplicates[N];
  int survival[N];
  vector[N] dilution;
  vector[N_sim] dilution_sim;
  int nsample;
  int sample[N];
  int sample_sim[N_sim];
}

parameters {
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> c;
  vector<lower=0>[nsample] phi;
}

model {
  vector[N] concentration = log(phi[sample] ./ dilution);
  survival ~ binomial(nreplicates, logistic(a, b, c, N, concentration));
  
  // priors
  a ~ cauchy(0, 10);
  b ~ cauchy(0, 10);
  c ~ cauchy(0, 10);
  phi ~ cauchy(0, 100);
}

generated quantities{
  vector[N_sim] prob;
  {
    vector[N_panels] phi_sim = phi[phi_panel_index];
    vector[N_sim] concentration = log(phi_sim[sample_sim] ./ dilution_sim);
    prob = logistic(a, b, c, N_sim, concentration);
  }
}
