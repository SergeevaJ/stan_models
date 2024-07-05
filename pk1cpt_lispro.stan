// One compartment model with zero- and first-order absorption for insulin lispro pharmokinetic data
functions{
//dose must be included into params as 7th param!!!!!
  vector ode_rhs(real t, vector x, array[] real parms, array[] real
                 x_r, array[] int x_i){
    real Ka = parms[1];
    real CL = parms[2];
    real V = parms[3];
    real tlag = parms[4];
    real tinf = parms[5];
    real frac = parms[6];
    real dose =parms[7];

    real CL_V = CL / V;
    real abs_ = Ka*dose*frac;
    real k0 = (1-frac)*dose/(tinf);
    real abs_cond=abs_;
    real k0_cond=k0;

    if (t>tlag) {
      abs_cond = abs_ / V;
    }
    else { 
      abs_cond = 0;
    }

    if (t<= tinf) {
       k0_cond = k0 / V; 
    } 
    else { 
      k0_cond = 0;
    }
       
    vector[1] y;

    y[1] = abs_cond + k0_cond - CL_V*x[1];
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
  int nCmt = 1; // number of compartments +1
  array[nSubjects] int nti; // number of obs for every subject
  for (i in 1:nSubjects) nti[i] = end[i] - start[i] + 1;
  }
parameters{
  real<lower = 0, upper = 200> Ka_lis_hat;
  real<lower = 0, upper = 200> CL_lis_hat;
  real<lower = 0, upper = 200> V_lis_hat;
  real<lower = 0, upper = 200> tlag_lis_hat;
  real<lower = 0, upper = 200> tinf_lis_hat;
real<lower = 0, upper = 1> frac_lis_hat;
}

transformed parameters{
 // array[nTheta] real<lower = 0> theta;
array[nTheta+1] real<lower = 0> theta_d;
  matrix<lower = 0>[nCmt, nt] x;
  row_vector<lower = 0>[nt] cHat; // estimation of DV
  row_vector<lower = 0>[nObs] cHatObs; //

  thetaHat[1] = Ka_lis_hat;
  thetaHat[2] = CL_lis_hat;
  thetaHat[3] = V_lis_hat;
  thetaHat[4] = tlag_lis_hat;
  thetaHat[5] = tinf_lis_hat;
  thetaHat[6] = frac_lis_hat;
// how to define theta for all 
  thetaM = rep_matrix(thetaHat, nSubjects);
  for(j in 1:nSubjects)
  {
    theta[1] = thetaM[j, 1]; 
    theta[2] = thetaM[j, 2];
    theta[3] = thetaM[j, 3]; 

    theta[4] = thetaM[j, 4]; 
    theta[5] = thetaM[j, 5]; 
    theta[6] = thetaM[j, 6]; 
    theta_d[1:6] = theta[1:6];
    theta_d[7] = amt[start[j]]; 

    x[, start[j]:end[j]] = pmx_solve_rk45(ode_rhs, 1, time[start[j]:end[j]], 
                                            amt[start[j]:end[j]],
                                            rate[start[j]:end[j]],
                                            ii[start[j]:end[j]],
                                            evid[start[j]:end[j]],
                                            cmt[start[j]:end[j]],
                                            addl[start[j]:end[j]],
                                            ss[start[j]:end[j]],
theta_d, 1e-5, 1e-8, 1e5);
                                       
    cHat[start[j]:end[j]] = x[2, start[j]:end[j]] ./ theta[3]; // divide by V_lis
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

// how to transform cHat, how to define ThetaM.
//how to deal with population
//define differential equations

}



