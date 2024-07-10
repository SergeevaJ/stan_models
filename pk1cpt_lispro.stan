// One compartment model with zero- and first-order absorption for insulin lispro pharmacokinetic data
functions{
  vector ode_rhs(real t, vector x, array[] real parms, array[] real
                 x_r, array[] int x_i){
    real Ka = parms[1];
    real CL = parms[2];
    real V = parms[3];
    real tlag = parms[4];
    real tinf = parms[5];
    real frac = parms[6];
    real dose = parms[7];
    //dose must be multiplied by frac
    real CL_V = CL / V;
    real k0 = (1-frac)*dose/(tinf*V);

    if (t>= tinf) {
       k0 = 0;
    } 
       
    vector[2] y;
    y[1] = -Ka*x[1];
    y[2] = Ka*x[1] + k0 - CL_V*x[2];
    return y;
  }
}
data{
  int<lower = 1> nt; //n - number of all events
  int<lower = 1> nObs; //n-k - number of measurements
  int<lower = 1> nSubjects; // k - number of studies 
  array[nObs] int<lower = 1> iObs; //n-k, index of observations - array of measurement indexes  
  array[nSubjects] int<lower = 1> start; //k, array of indexes where start every new subject 
  array[nSubjects] int<lower = 1> end; //k, array of indexes where end every new subject 
  array[nt] int<lower = 1> cmt; //CMT
  array[nt] int evid; //EVID
  array[nt] int addl; // n, zeros
  array[nt] int ss; // n, zeros
  array[nt] real amt; //AMT
  array[nt] real time; //TIME
  array[nt] real rate; // n, zeros
  array[nt] real ii; // n, zeros
  vector<lower = 0>[nObs] cObs; //DV
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int nTheta = 6; // number of estimated parameters
  int nCmt = 1; // number of compartments
  array[nSubjects] int nti; // number of obs for every subject
  for (i in 1:nSubjects) nti[i] = end[i] - start[i] + 1;
  real biovar = 1;
  }
parameters{
  real<lower = 0> Ka_lis_hat;
  real<lower = 0> CL_lis_hat;
  real<lower = 0> V_lis_hat;
  real<lower = 0> tlag_lis_hat;
  real<lower = 0> tinf_lis_hat;
  real<lower = 0, upper=1>  frac_lis_hat;
  real <lower=0> sigma;
}

transformed parameters{
  array[nTheta-1] real<lower = 0> thetaHat;
  array[nTheta] real<lower = 0> theta_d;
  real<lower=0> tlag;
  vector[nt] amt_mod; 
  matrix<lower = 0>[nCmt, nt] x;
  row_vector<lower = 0>[nt] cHat; // estimation of DV
  row_vector<lower = 0>[nObs] cHatObs; //

  thetaHat[1] = Ka_lis_hat;
  thetaHat[2] = CL_lis_hat;
  thetaHat[3] = V_lis_hat;
  tlag = tlag_lis_hat;
  thetaHat[4] = tinf_lis_hat;
  thetaHat[5] = frac_lis_hat;
  amt_mod = amt
  amt_mod = amt_mod*thetaHat[5];
  for(j in 1:nSubjects)
  { theta_d[1:5] = thetaHat[1:5]; 
    theta_d[6] = amt[start[j]];

    x[, start[j]:end[j]] = pmx_solve_rk45(ode_rhs, 1, time[start[j]:end[j]], 
                                            amt_mod[start[j]:end[j]],
                                            rate[start[j]:end[j]],
                                            ii[start[j]:end[j]],
                                            evid[start[j]:end[j]],
                                            cmt[start[j]:end[j]],
                                            addl[start[j]:end[j]],
                                            ss[start[j]:end[j]],
theta_d, biovar, tlag, 1e-5, 1e-8, 1e5);
                                       
    cHat[start[j]:end[j]] = x[2, start[j]:end[j]] ./ theta_d[3]; // divide by V_lis
  }

  cHatObs  = cHat[iObs];
}

model{
  // Prior
  Ka_lis_hat ~ lognormal(log(0.989), 0.03);
  CL_lis_hat ~ lognormal(log(30.5), 0.5);
  V_lis_hat ~ lognormal(log(43), 0.5);
  tlag_lis_hat ~ lognormal(log(0.265), 0.01);
  tinf_lis_hat ~ lognormal(log(1.06), 0.01);
  frac_lis_hat ~ uniform(0, 1); 
  sigma ~ cauchy(0, 1);

  logCObs ~ normal(log(cHatObs), sigma);
}



