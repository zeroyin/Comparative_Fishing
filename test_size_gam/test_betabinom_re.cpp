// test nagative multinomial model

#include <TMB.hpp>

template<class Type>
// Template function code rom R package glmmTMB
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
	
	// data:
	DATA_VECTOR(n_A);
	DATA_VECTOR(n_B);
	DATA_MATRIX(Xf);
	DATA_MATRIX(Xr);

	vector<Type> n_N = n_A + n_B; 

	// parameters:
	PARAMETER_VECTOR(beta);
	PARAMETER_VECTOR(b);
	PARAMETER_VECTOR(gamma);
	PARAMETER_VECTOR(g);
	PARAMETER_VECTOR(delta);
	PARAMETER_VECTOR(epsilon);
	PARAMETER(lambda, tau)

	// set up objective fn 
	const int nlen = n_N.size();
	Type nll = Type(0.0); // initialize negative log likelihood


    
	// linear predictors
	vector<Type> eta = Xf * (beta + delta) + Xr * (b + epsilon);
	vector<Type> etad = Xf * gamma + Xr * g;

	// link function
	vector<Type> mu = invlogit(eta);
	vector<Type> phi = exp(etad);


	// random effects
	for(int i=0;i<nlen;i++){
		nll -= dnorm(b(i), 0, [pos eigen val of S]/lambda);
	}
	for(int i=0;i<nlen;i++){
		nll -= dnorm(g(i), 0, [pos eigen val of S]/tau);
	}
	

	// over-dispersion
	nll -= dmvnorm(delta, 0, C1);
	nll -= dmvnorm(epsilon, 0, C2)


	// Observation likelihood
	for(int i=0;i<nlen;i++){
        Type s1 = mu(i)*phi(i); // s1 = mu(i) * mu(i) / phi(i);
        Type s2 = (Type(1)-mu(i))*phi(i); // phi(i) / mu(i);
		nll -= dbetabinom(n_A(i), s1, s2, n_N(i), TRUE);
	}



	// report
	REPORT(mu);
	REPORT(phi);

	vector<Type> rho = mu/(1-mu);
	REPORT(rho);


	return nll;
}






