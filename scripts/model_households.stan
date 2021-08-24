/* Household final size model, as used in DFM Reukers, M van Boven, A Meijer,*/
/* N Rots, C Reusken, I Roof, AB van Gageldonk-Lafeber, W van der Hoek,      */
/* S van den Hof, High infection secondary attack rates of SARS-CoV-2 in     */
/* Dutch households revealed by dense sampling, Clinical Infectious Diseases */
/* (2021), https://doi.org/10.1093/cid/ciab237.                              */
/* Notation is mostly taken from F Ball (1986) Advances in Applied           */
/* Probability 18, 289-310, and later paers.                                 */
/* Let a, n and j be the number of initial infections, the number of         */
/* uninfected individuals in the household, and the number of new            */
/* infections after the outbreak. Throughout we assume there are 4 types.    */
/* Let q(j|a,n,beta) denote the probability of final size j+a, given         */
/* a and n and parameters beta etcetera, then                                */
/* q(j|a,n,beta) = \binom{n}{j} p(j|a,n,beta) where p is defined by          */
/* 1 = \sum_{\omega \leq j} \binom{j}{\omega} p(\omega | a, n) \times        */
/* \prod_{k=1}^4 \phi(\sum_{l=1}^4 \beta_{lk}(n_l-j_l)/|a+n|)^{-omega_i-a_i} */
/* Here, phi is the Laplace transform of the pdf of the infectious period.   */
/* For instance, in case of a fixed infections period of 1 time unit, we get */
/* \phi(x) = \int_0^{\infty} e^{-tx} \delta_{t-1} dt = e^{-x}                */
/* Also included is a probability of external infection (i.e. 1-escape).     */
/* Code is drafted by Chris van Dorp and Michiel van Boven (10/2020), and    */
/* licensed with BSD3 clause (use but mention).                              */

functions {
  /* Laplace transform of the scaled infectious period (E(T_I)=1) with realistic  */
  /* gamma distribution (alpha=beta=50; alt par k=50, theta=1/50). ~95% coverage  */
  /* is 0.75-1.25, ie 6-10 days when mean is 8 days                               */
  /* Infectious period of fixed duration yields virtually identical results.      */
  real phi(real x) {
    return (1+x/50)^-50;                                                                                
  }
  
  /* Product of powers of Laplace transforms and external infection escape probabilities */
  /* take san=1 in function below for a model with density-dependent transmission        */
  /* (instead of frequency-dependent transmission). Also note that a possibility for     */
  /* external infection is included (or rather external escape): b[i]^(n[i] - j[i])      */
  real phi_prod(int[] omega, int[] j, int[] a, int[] n, matrix beta, vector b) {
    int san = sum(a) + sum(n);  // adapt for density-dependent transmission (san=1)                      
    real F = 1.0;                                     
    for ( i in 1:4 ) {
      real Bx = 0;
      for ( k in 1:4 ) {
        Bx += beta[k,i]*(n[k] - j[k]);
      }
      F *= phi(Bx/san)^(omega[i] + a[i]) * b[i]^(n[i] - j[i]); 
    }
    return F;
  }
  
  /* binomial coefficient for integer vectors as product of components */
  int binom_prod(int[] n, int[] k) {
    int B = 1;
    for ( i in 1:num_elements(n) ) {
      B *= choose(n[i], k[i]);
    }
    return B;
  }
  
  /* iterative calculation of the houshold infection probabilities, */
  real prob_infect_pattern(int[] j, int[] a, int[] n, matrix beta, vector b) {
    int san = sum(a) + sum(n);
    // cache array: we need space for 0,...,j elements (i.e. j+1)
    real P[j[1]+1, j[2]+1, j[3]+1, j[4]+1];
    for ( i1 in 1:j[1]+1 ) {
      for ( i2 in 1:j[2]+1 ) {
        for ( i3 in 1:j[3]+1 ) {
          for ( i4 in 1:j[4]+1 ) {
            // Stan starts indexing at 1: correct this.
            int omega[4] = {i1-1, i2-1, i3-1, i4-1};
            // compute value P[i1, i2, i3, i4] using previous values
            real S = 0;
            for ( k1 in 1:i1 ) {
              for ( k2 in 1:i2 ) {
                for ( k3 in 1:i3 ) {
                  for ( k4 in 1:i4 ) {
                    // exclude the point we want to compute...
                    if ( k1 < i1 || k2 < i2 || k3 < i3 || k4 < i4 ) {
                      int wou[4] = {k1-1, k2-1, k3-1, k4-1};
                      S += P[k1, k2, k3, k4] * binom_prod(omega, wou) / phi_prod(wou, omega, a, n, beta, b);
                    }
                  } // k4
                } // k3
              } // k2
            } // k1
            /* compute the missing element of P using the sum S. */
            /* notice that the binomial coefficients equal 1     */
            P[i1, i2, i3, i4] = (1-S) * phi_prod(omega, omega, a, n, beta, b);
            // inner loops
          } // i4
        } // i3
      } // i2
    } // i1
    // outer loops
    return P[j[1]+1, j[2]+1, j[3]+1, j[4]+1] * binom_prod(n, j); 
  }
  
  /* include conditioning for households with a non-primary index case */
  real prob_pos_infect(int c, int[] a, int[] n, matrix beta, vector b) {
    int sel[4] = rep_array(1, 4); // sel selects what for loops to use below
    real S = 0.0;
    if ( c > 0 ) {
      sel[c] = 0; // if select[c] is zero, the c-th for loop only has one term.
      for ( w1 in 0:n[1]*sel[1] ) {
        for ( w2 in 0:n[2]*sel[2] ) {
          for ( w3 in 0:n[3]*sel[3] ) {
            for ( w4 in 0:n[4]*sel[4] ) {
              int wau[4] = {w1, w2, w3, w4};
              S += prob_infect_pattern(wau, a, n, beta, b);
            } // w4
          } // w3
        } // w2
      } // w1
    } // if c > 0 (else function returns 1.0)
    return 1.0 - S; 
  }
}

data {
  int<lower=0> H;                                     // number of households
  int<lower=0> N[H, 4];                               // initial uninfected individuals
  int<lower=0> A[H, 4];                               // number of initial cases
  int<lower=0> J[H, 4];                               // number of household infections
  int<lower=0, upper=4> conditioning[H];              // conditioning in case of single index
  int <lower = 0, upper = 1> mode;                    // 0 = estimation, 1 = WBIC calculation
 
  /* fix selected parameters for various scenarios */
  vector<lower = 0>[2] susceptibility;                // susceptibility as data; streamline!
  //real<lower = 0> susceptibility_children;
  vector<lower = 0>[2] infectivity;                   // infectivity as data
  vector<lower = 0, upper = 1>[4] external_escape;    // probabilities of escape from external infection
}

transformed data {	
  /* sampling temperature */  
  real<lower = 0, upper = 1> watanabe_beta;           // 0 = parameter estimation; 1 = WBIC
  
  /* sampling mode. 0 = normal; 1 = WBIC */
  if ( mode == 0 ) {
    watanabe_beta = 1.0;
  }
  else { // mode == 1
    watanabe_beta = 1.0/log(H);                       // determines the sampling temperature
  }
}

parameters {
  //vector<lower = 0>[2] susceptibility;              // susceptibility as parameter
  //real<lower = 0> susceptibility_children;          // susceptibility of children as parameter
  //vector<lower = 0>[2] infectivity;                 // infectivity as parameter
  //real<lower = 0, upper = 1> ext_escape;            // probability of escape from external infection
  real<lower = 0> beta;                               // transmission rate in reference class
  //vector<lower = 0, upper = 1>[4] external_escape;  // probabilities of escape from external infection (\approx 1)
}

transformed parameters {
  matrix<lower = 0>[4,4] transmission_rate;           // transmission rates between different types
  //vector<lower = 0>[2] susceptibility;
  //vector<lower = 0, upper = 1>[4] external_escape;  // probabilities of escape from external infection (\approx 1)
  vector<lower = 0>[4] sus;                           // susceptibility with 1 padded for reference class (last class)
  row_vector<lower = 0>[4] inf;                       // susceptibility with 1 padded for reference class (last class)
  vector[H] log_lik;                                  // household log-likelihood contributions
  real log_like;                                      // sum of log-likelihoods
  
  /* scenario with only variable susceptbility of children */
  //susceptibility[1] = susceptibility_children;
  //susceptibility[2] = 1;

  /* inclusion of external infection */
  /*
  for ( i in 1:4 ) {
	  external_escape[i] = ext_escape;
  }
  */

  /* calculate transmission rates from underlying assumptions */
  /* notice that the third class is now the reference group   */
  sus = append_row(append_row(susceptibility, 1.0), 0.0); 
  inf = append_row(append_row(infectivity, 1.0), 0.0)';  
  transmission_rate = beta * sus * inf;               // type-to-type transmission rates per infectious period
  
  /* log-likelihood contributions */
  for ( i in 1:H ) {
	  if ( conditioning[i] == 0 ) { // no conditioning, index case is also primary case
	    log_lik[i] = log(prob_infect_pattern(J[i,:], A[i,:], N[i,:], transmission_rate, external_escape));
    } else { // conditioning for non-primary index case 
	    log_lik[i] = log(prob_infect_pattern(J[i,:], A[i,:], N[i,:], transmission_rate, external_escape) / 
	        prob_pos_infect(conditioning[i], A[i,:], N[i,:], transmission_rate, external_escape));
	  }
  }
  
  /* WBIC calculation (mode = 1) */
  log_like = sum(log_lik);
}

model {
  /* prior distributions */
  //ext_escape ~ beta(1,1);                           // flat prior for external infection
  //ext_escape ~ beta(29,1);                          // weekly informative; 95% prior coverage >0.9
  
  /* log-likelihood */
  target += watanabe_beta * sum(log_lik);
}

generated quantities {
  // generate random household infection data for scenario analyses -> extensions possible
  // include sample size analyses; how many households are needed?
}

