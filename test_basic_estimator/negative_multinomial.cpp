// test nagative multinomial model

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
	
	// data:
	DATA_VECTOR(n_A);
	DATA_VECTOR(n_B);
    
	// parameters:
	PARAMETER(log_d);
	PARAMETER(log_r);
	PARAMETER_VECTOR(log_dens);
	PARAMETER(log_q_A);
	PARAMETER(log_q_B);

	Type d = exp(log_d);
	Type r = exp(log_r);
	vector<Type> dens = exp(log_dens);
	Type q_A = exp(log_q_A);
	Type q_B = exp(log_q_B);
	
	vector<Type> mu_A = q_A*dens;
	vector<Type> mu_B = q_B*dens;
	Type rho = q_B/q_A;

	// obj fn 
	const int n_i = dens.size();
	Type nll = Type(0.0); // initialize negative log likelihood

	for(int i=0;i<n_i;i++){
		nll -= dgamma(dens(i), r, d/r , TRUE); // shape and scale from assumed mean and var
	}

	for(int i=0;i<n_i;i++){
		nll -= dpois(n_A(i), mu_A(i), TRUE);
		nll -= dpois(n_B(i), mu_B(i), TRUE);
	}


	// report
	// REPORT(dens);
	// REPORT(rho);
	ADREPORT(rho);


	return nll;
}






