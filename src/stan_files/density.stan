data{
  int M;
  vector[M] x;
  vector[M] y;
  real rho;
  real sigma;
}

transformed data{
  vector[M] Mu = rep_vector(0,M);
}

parameters{
  vector[M] f;
}

transformed parameters{
  simplex[M] p;
  cov_matrix[M] Cov;
  real log_lik;
  log_lik = 0;

  //p
  p = softmax(f);

  //kernel function
  for(i in 1:M){
    for(j in 1:M){
      Cov[i,j] = sigma*exp(-(x[i]-x[j])^2/rho);
      if(i==j) Cov[i,j] = Cov[i,j] + 0.01;
    }
  }

  //likelihood function
  {
    for(i in 1:M){
      log_lik += y[i]*log(p[i]);
    }
  }
}

model{
  f ~ multi_normal(Mu,Cov);
  //likelihood
  target += log_lik;
}
