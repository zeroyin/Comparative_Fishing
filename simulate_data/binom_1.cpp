// various binomial models
// various smoothing assumptions

#include <TMB.hpp>

template<class Type>
// Template function code rom R package glmmTMB
// verify re-parameterization into mu and phi
Type dbetabinom(Type y, Type s1, Type s2, Type n, int give_log=0)
{
    /*
      Wikipedia:
      f(k|n,\alpha,\beta) =
      \frac{\Gamma(n+1)}{\Gamma(k+1)\Gamma(n-k+1)}
      \frac{\Gamma(k+\alpha)\Gamma(n-k+\beta)}{\Gamma(n+\alpha+\beta)}
      \frac{\Gamma(\alpha+\beta)}{\Gamma(\alpha)\Gamma(\beta)}
    */
	Type logres =
	lgamma(n + 1) - lgamma(y + 1)     - lgamma(n - y + 1) +
	lgamma(y + s1) + lgamma(n - y + s2) - lgamma(n + s1 + s2) +
	lgamma(s1 + s2) - lgamma(s1)         - lgamma(s2) ;
	if(!give_log) return exp(logres);
	else return logres;
}

template<class Type>
// prob of at least one vessel has pos catch: p 
// should include stations where both catch zero?
Type dzerobinom(Type y, Type n, Type mu, Type p, int give_log=0)
{
	Type logres = 0;
	if(y == 0){
		logres = log((1 - p)/2);
	}else if(y == n){
		logres = log((1 - p)/2);
	}else{
		logres = log(p) + dbinom(y, n, mu, true);
	}
	if(!give_log) return exp(logres);
	else return logres;
}

template<class Type>
Type objective_function<Type>::operator() () {
	
	using namespace density;
	
	// data:
	// filter out zeros?
	DATA_MATRIX(A); // catch number by A, by station by length
	DATA_MATRIX(B); // catch number by B, by station by length
	DATA_MATRIX(offset); // offset for each tow
	DATA_MATRIX(Xf); // design matrix for fixed effect of mu, phi
	DATA_MATRIX(Xr); // design matrix for random effect of mu, phi
	DATA_VECTOR(d); // positive eigenvalues of penalty matrix S
	DATA_INTEGER(idist); // distribution type: BI 0; BB 1.

	// parameters:
	PARAMETER_VECTOR(beta); // coeff for fixed effect for mu
	PARAMETER_VECTOR(b); // coeff for random effect for mu
	PARAMETER(log_s_b); // smooth term for b
	PARAMETER_VECTOR(gamma); // coeff for fixed effect for phi
	PARAMETER_VECTOR(g); // coeff for random effect for phi
	PARAMETER(log_s_g); // smooth term for g
	PARAMETER_MATRIX(delta); // overdispersion for beta
	PARAMETER_VECTOR(chol_delta); // chol decomposition and vectorized cov matrix for delta
	PARAMETER_MATRIX(epsilon); // overdispersion for b
	PARAMETER(log_s_epsilon); // smooth term for epsilon
	PARAMETER(beta_0); 
    PARAMETER(gamma_0);
    PARAMETER_VECTOR(delta_0);
    PARAMETER(log_sigma_delta_0);
    PARAMETER_MATRIX(p);
    PARAMETER(log_p_s1);
    PARAMETER(log_p_s2);


	// transformation
	const int n_len = A.cols(); // number of length bins
	const int n_s = A.rows(); // number of stations
	const int n_f = beta.size(); // default 2
	const int n_r = b.size(); // default 8

	matrix<Type> N = A + B; // total catch

	Type s_b = exp(log_s_b); // smooth term for b
	Type s_g = exp(log_s_g); // smooth term for g
	Type s_epsilon = exp(log_s_epsilon); // smooth term for epsilon

	// Covariance matrix for delta
	matrix<Type> L_delta(n_f,n_f); L_delta.setZero(); 
	int k_f = 0;
	for(int i_f = 0; i_f < n_f; i_f++){
		for(int j_f = 0; j_f <= i_f; j_f++){
			L_delta(i_f, j_f) = chol_delta(k_f);
			k_f++;
		}
	}
	matrix<Type> C_delta(n_f,n_f); 
	C_delta = L_delta * L_delta.transpose();


	// set up objective fn components: negative log likelihood
	vector<Type> nll(14); nll.setZero(); // initialize


	// random effects with smooth penalty
	for(int i_r = 0; i_r < n_r; i_r++){
		nll(0) -= dnorm(b(i_r), Type(0), d(i_r)/s_b, true);
		nll(1) -= dnorm(g(i_r), Type(0), d(i_r)/s_g, true);
	}

	// station-level over-dispersion
	for (int i_s = 0; i_s < n_s; i_s++){
		vector<Type> tmp_delta = delta.row(i_s);
		nll(2) += MVNORM(C_delta)(tmp_delta);
		vector<Type> tmp_epsilon = epsilon.row(i_s);
		nll(3) -= sum(dnorm(tmp_epsilon, Type(0), d/s_epsilon, true));
	}
	// option for delta_0
	nll(4) -= sum(dnorm(delta_0, Type(0), exp(log_sigma_delta_0), true));
	// option for ZB
	for (int i_s = 0; i_s < n_s; i_s++){
		for (int i_len = 0; i_len < n_len; i_len++){
			nll(5) -= dbeta(p(i_s,i_len), exp(log_p_s1), exp(log_p_s2), true);
		}
	}


	// Observation likelihood
	matrix<Type> eta_mu(n_s, n_len); eta_mu.setZero();
	matrix<Type> eta_phi(n_s, n_len); eta_phi.setZero();
	matrix<Type> mu(n_s, n_len);
	matrix<Type> phi(n_s, n_len);
	matrix<Type> rho(n_s, n_len);

	for(int i_s = 0; i_s < n_s; i_s++){
		for(int i_len = 0; i_len < n_len; i_len++){
			// linear predictors
			eta_mu(i_s, i_len) += beta_0 + delta_0(i_s);
			eta_phi(i_s, i_len) += gamma_0;
			for(int i_f = 0; i_f < n_f; i_f++){
				eta_mu(i_s, i_len) += Xf(i_len, i_f) * (beta(i_f) + delta(i_s, i_f));
				eta_phi(i_s, i_len) += Xf(i_len, i_f) * gamma(i_f);
			}
			for(int i_r = 0; i_r < n_r; i_r++){
				eta_mu(i_s, i_len) += Xr(i_len, i_r) * (b(i_r) + epsilon(i_s, i_r));
				eta_phi(i_s, i_len) += Xr(i_len, i_r) * g(i_r);
			}
			// link function
			mu(i_s, i_len) = invlogit(eta_mu(i_s, i_len) + offset(i_s, i_len));
			phi(i_s, i_len) = exp(eta_phi(i_s, i_len));
			// conversion rate calculation alternative: 
			rho(i_s, i_len) = exp(eta_mu(i_s, i_len));

			// transformation to BB shape parameter
				Type s1 = mu(i_s, i_len)*phi(i_s, i_len);
				Type s2 = (Type(1)-mu(i_s, i_len))*phi(i_s, i_len);

			// observation likelihood: unsuccessful using switch statement
			// if(A(i_s, i_len) > 0 | B(i_s, i_len) > 0){
			if(idist == 0){
				nll(10) -= dbinom(A(i_s, i_len), N(i_s, i_len), mu(i_s, i_len), true);
			}else if(idist == 1){
				nll(11) -= dbetabinom(A(i_s, i_len), s1, s2, N(i_s, i_len), true);
			}else if(idist == 2){
				nll(12) -= dzerobinom(A(i_s, i_len), N(i_s, i_len), mu(i_s, i_len), p(i_s, i_len), true);
			}
			// }
		}
	}


	// derived quantities  
	vector<Type> mean_mu = mu.colwise().sum()/mu.rows(); // column means: mu
	ADREPORT(mean_mu);
	vector<Type> mean_phi = phi.colwise().sum()/phi.rows(); // column means: phi
	ADREPORT(mean_phi);
	vector<Type> mean_log_rho = eta_mu.colwise().sum()/eta_mu.rows(); // use log rho
	ADREPORT(mean_log_rho);

	// // sdreport
	// ADREPORT(mu);
	// ADREPORT(phi);
	// ADREPORT(rho);

	// report
	REPORT(mu);
	REPORT(phi);
	REPORT(rho);
	REPORT(eta_mu);
	REPORT(eta_phi);
	REPORT(beta);
	REPORT(b);
	REPORT(gamma);
	REPORT(g);
	REPORT(delta);
	REPORT(epsilon);
	REPORT(C_delta);
	REPORT(p);
	REPORT(mean_mu);
	REPORT(mean_phi);
	REPORT(mean_log_rho);

	// sum up nll as joint nll: 
	// computation is sacrificed for code readability  
	Type jnll = nll.sum();

	// report for diagnostics
	REPORT(nll);
	REPORT(jnll);

	return jnll;
}






