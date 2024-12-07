data {
    int<lower=1> N;          // Number of observations
    array[N] real y;         // Observed data
}
parameters {
    real<lower=0> theta;     // Parameter of interest
}
transformed parameters {
    vector[N] L;             // Likelihood
    vector[N] pro;           // Normalized probabilities
    real prior_num;          // Custom prior numerator
    real C = 1e6;            // Normalization constant

    // Custom prior numerator
    prior_num = (1 / (theta^4 + theta^3 + 2*theta^2 + 6*theta)) *
                sqrt(theta^6 + 4*theta^5 + 18*theta^4 + 
                     96*theta^3 + 72*theta^2 + 96*theta + 144);

    // Compute likelihood and probabilities
    for (i in 1:N) {
        L[i] = ((theta^4) / (theta^3 + theta^2 + 2*theta + 6)) *
               (1 + y[i] + y[i]^2 + y[i]^3) * exp(-theta * y[i]);
        pro[i] = L[i] / C;
    }
}
model {
    // Prior
    target += log(prior_num);  // Log-prior for numerical stability

    // Likelihood
    for (i in 1:N) {
        target += bernoulli_lpmf(1 | pro[i]); // Log-likelihood
    }
}

