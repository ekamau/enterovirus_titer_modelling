// each individual phi estimate

functions{
  vector logistic(real a, real b, int N_len, vector log_concentration){
    vector[N_len] denom = 1 + exp(-(a + b * log_concentration));
    return(1 ./ denom);
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
  real a;
  real b;
}

parameters {
  vector<lower=0>[nsample] phi;
}

model {
  vector[N] concentration = log(phi[sample] ./ dilution);
  survival ~ binomial_logit(nreplicates, a + b * concentration);
  
  // priors
  phi ~ normal(0, 2000);
}

generated quantities{
  vector[N_sim] prob;
  {
    vector[N_panels] phi_sim = phi[phi_panel_index];
    vector[N_sim] concentration1 = log(phi_sim[sample_sim] ./ dilution_sim);
    prob = logistic(a, b, N_sim, concentration1);
  }
}
