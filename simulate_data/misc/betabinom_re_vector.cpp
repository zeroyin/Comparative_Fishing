// test betabinomial model, length GAM, with random effects

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
	DATA_VECTOR(A);
	DATA_VECTOR(B);
	DATA_IVECTOR(indstn);
	DATA_IVECTOR(indlen);
	DATA_VECTOR(offset);
	DATA_MATRIX(Xf); // design matrix for fixed effect of mu, phi
	DATA_MATRIX(Xr); // design matrix for random effect of mu, phi
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
	PARAMETER_VECTOR(chol_delta); // chol decomposition and vectorized cov matrix for delta

	// transformation
	vector<Type> N = A + B; // total catch
	const int n_obs = N.size();
	const int n_f = beta.size();
	const int n_r = b.size();
	const int n_s = delta.rows();


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
	vector<Type> nll(5); nll.setZero(); // initialize

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
	vector<Type> eta_mu(n_obs); eta_mu.setZero();
	vector<Type> eta_phi(n_obs); eta_phi.setZero();
	vector<Type> mu(n_obs);
	vector<Type> phi(n_obs);
	vector<Type> rho(n_obs);
	for(int i = 0; i < n_obs; i++){
			// linear predictors
		for(int i_f = 0; i_f < n_f; i_f++){
			eta_mu(i) += Xf(indlen(i), i_f) * (beta(i_f) + delta(indstn(i), i_f));
			eta_phi(i) += Xf(indlen(i), i_f) * gamma(i_f);
		}
		for(int i_r = 0; i_r < n_r; i_r++){
			eta_mu(i) += Xr(indlen(i), i_r) * (b(i_r) + epsilon(indstn(i), i_r));
			eta_phi(i) += Xr(indlen(i), i_r) * g(i_r);
		}

			// link function
		mu(i) = invlogit(eta_mu(i) + offset(i));
		phi(i) = exp(eta_phi(i));

			// conversion rate calculation alternative: 
		rho(i) = exp(eta_mu(i));
			// conversion rate derived from proportion:
			// rho(i) = mu(i)/(Type(1)-mu(i))*exp(-offset(i));

			// transformation to shape parameter
			Type s1 = mu(i)*phi(i); // s1 = mu(i) * mu(i) / phi(i);
			Type s2 = (Type(1)-mu(i))*phi(i); // phi(i) / mu(i);

			// observation likelihood
			nll(4) -= dbetabinom(A(i), s1, s2, N(i), true);

		}


	// sdreport
		ADREPORT(mu);
		ADREPORT(phi);
		ADREPORT(rho);

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

	// sum up nll as joint nll
		Type jnll = nll.sum();

	// report for diagnostics
		REPORT(nll);
		REPORT(jnll);

		return jnll;
	}






