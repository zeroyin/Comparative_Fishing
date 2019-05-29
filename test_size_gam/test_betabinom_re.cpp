// test nagative multinomial model

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
Type objective_function<Type>::operator() () {
	
	using namespace density;
	
	// data:
	DATA_MATRIX(A); // catch number by A, by station by length
	DATA_MATRIX(B); // catch number by B, by station by length
	DATA_MATRIX(Xf); // design matrix for fixed effect
	DATA_MATRIX(Xr); // design matrix for random effect
	DATA_VECTOR(d); // positive eigenvalues of penalty matrix S

	// parameters:
	PARAMETER_VECTOR(beta); // coeff for fixed effect for mu
	PARAMETER_VECTOR(b); // coeff for random effect for mu
	PARAMETER_VECTOR(gamma); // coeff for fixed effect for phi
	PARAMETER_VECTOR(g); // coeff for random effect for phi
	PARAMETER_MATRIX(delta); // overdispersion for beta
	PARAMETER_MATRIX(epsilon); // overdispersion for b
	PARAMETER(log_s_b); // smooth term for b
	PARAMETER(log_s_g); // smooth term for g
	PARAMETER(log_s_epsilon); // smooth term for epsilon
	PARAMETER_MATRIX(C_delta); // Covariance matrix for delta

	// transformation
	matrix<Type> N = A + B; // total catch
	Type s_b = exp(log_s_b); // smooth term for b
	Type s_g = exp(log_s_g); // smooth term for g
	Type s_epsilon = exp(log_s_epsilon); // smooth term for epsilon


	// set up objective fn 
	const int n_len = A.cols(); // number of length bins
	const int n_s = A.rows(); // number of stations
	const int n_f = beta.size();
	const int n_r = b.size();
	vector<Type> nll(5); nll.setZero(); // initialize negative log likelihood

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

	// Observation likelihood
	matrix<Type> mu(n_s, n_len);
	matrix<Type> phi(n_s, n_len);
	for(int i_s = 0; i_s < n_s; i_s++){
		for(int i_len = 0; i_len < n_len; i_len++){
			// linear predictors
			Type eta_mu = Type(0);
			Type eta_phi = Type(0);
			for(int i_f = 0; i_f < n_f; i_f++){
				eta_mu += Xf(i_len, i_f) * (beta(i_f) + delta(i_s, i_f));
				eta_phi += Xf(i_len, i_f) * gamma(i_f);
			}
			for(int i_r = 0; i_r < n_r; i_r++){
				eta_mu += Xr(i_len, i_r) * (b(i_r) + epsilon(i_s, i_r));
				eta_phi += Xr(i_len, i_r) * g(i_r);
			}
			// link function
			mu(i_s, i_len) = invlogit(eta_mu);
			phi(i_s, i_len) = exp(eta_phi);

	        Type s1 = mu(i_s, i_len)*phi(i_s, i_len); // s1 = mu(i) * mu(i) / phi(i);
	        Type s2 = (Type(1)-mu(i_s, i_len))*phi(i_s, i_len); // phi(i) / mu(i);

			nll(4) -= dbetabinom(A(i_s, i_len), s1, s2, N(i_s, i_len), true);
		}
	}



	// report
	REPORT(mu);
	REPORT(phi);
	REPORT(beta);
	REPORT(b);
	REPORT(gamma);
	REPORT(g);
	REPORT(delta);
	REPORT(epsilon);


	// derived quantities  
	// vector<Type> rho = mu/(1-mu);
	// REPORT(rho);

	// sdreport
	// ADREPORT(mu);
	// ADREPORT(rho);
	
	Type jnll = nll.sum();
	
	REPORT(nll);
	REPORT(jnll);

	return jnll;
}






