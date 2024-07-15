// One compartment model with zero- and first-order absorption for insulin lispro pharmacokinetic data
functions{
  vector ode_rhs(real t, vector x, array[] real parms, array[] real
                 x_r, array[] int x_i){
    real Ka = parms[1];
    real CL = parms[2];
    real V = parms[3];
  
    real tinf = parms[4];
    real frac = parms[5];
    real dose = parms[6];
    real CL_V = CL / V;
    //(1-frac)*dose goes to zero-order absorption
    real k0 = (1-frac)*dose/tinf;

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
  array[nSubjects] int<lower = 1> start; //k, array of indexes where starts every new subject 
  array[nSubjects] int<lower = 1> end; //k, array of indexes where ends every new subject 
  array[nt] int<lower = 1> cmt; //CMT, 1 - dose "compartment", 2 - main compartment
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
  vector[nObs] logCObs = log(cObs+10^-10); //added small value to cObs, because at t=0 cObs is zero and log(cObs) is -inf
  int nTheta = 6; // number of estimated parameters
  int nCmt = 2; // number of compartments +1 (dose "compartment")
  array[nSubjects] int nti; // number of obs for every subject
  for (i in 1:nSubjects) nti[i] = end[i] - start[i] + 1;
  array[nCmt] real biovar; 
  biovar[1] = 1;
  biovar[2] = 1;
  }
parameters{
  real<lower = 0> Ka_lis_hat;
  real<lower = 0> CL_lis_hat;
  real<lower = 0> V_lis_hat;
  real<lower = 0> tlag_lis_hat;
  real<lower = 0> tinf_lis_hat;
  real<lower = 0, upper=1>  frac_lis_hat;
  real <lower=0, upper=10000> sigma;
}

transformed parameters{
  array[nTheta-1] real<lower = 0> thetaHat;
  array[nTheta] real<lower = 0> theta_d;
  array[nCmt] real<lower=0> tlag;
  array[nt] real amt_mod;
  matrix<lower = 0>[nCmt, nt] x;
  row_vector<lower = 0>[nt] cHat; // estimation of DV
  row_vector<lower = 0>[nObs] cHatObs; 

  thetaHat[1] = Ka_lis_hat;
  thetaHat[2] = CL_lis_hat;
  thetaHat[3] = V_lis_hat;
  tlag[1] = tlag_lis_hat;
  tlag[2] = tlag_lis_hat;
  thetaHat[4] = tinf_lis_hat;
  thetaHat[5] = frac_lis_hat;
  // multiply amt (dose) by frac, because this part goes through 1-order absorption
  amt_mod=amt;
  for (ind in start)
  {amt_mod[ind]=amt[ind]*frac_lis_hat;}


  for(j in 1:nSubjects)
  { theta_d[1:5] = thetaHat[1:5]; 
    theta_d[6] = amt[start[j]];

    x[, start[j]:end[j]] = pmx_solve_bdf(ode_rhs, 2, time[start[j]:end[j]], 
                                            amt_mod[start[j]:end[j]],
                                            rate[start[j]:end[j]],
                                            ii[start[j]:end[j]],
                                            evid[start[j]:end[j]],
                                            cmt[start[j]:end[j]],
                                            addl[start[j]:end[j]],
                                            ss[start[j]:end[j]],
theta_d, biovar, tlag, 1e-5, 1e-8, 1e5);
                                       
    cHat[start[j]:end[j]] = x[2, start[j]:end[j]] ./ theta_d[3]; // divide by V
  }

  cHatObs  = cHat[iObs]+10^-10;
}

model{
  // Prior
  Ka_lis_hat ~ lognormal(log(0.989), 0.03);
  CL_lis_hat ~ lognormal(log(30.5), 0.5);
  V_lis_hat ~ lognormal(log(43), 0.5);
  tlag_lis_hat ~ lognormal(log(0.265), 0.01);
  tinf_lis_hat ~ lognormal(log(1.06), 0.01);
  frac_lis_hat ~ uniform(0, 1); 
  sigma ~ uniform(0, 10000);

  logCObs ~ normal(log(cHatObs), sigma);
}



