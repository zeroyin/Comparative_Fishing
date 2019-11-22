// Spatial model, basic version +
// Options for temporal dependence of ST effect:
// Options(0)==2: AR(1) for ST effect
// Depth included as covariate
#include <TMB.hpp>


template<class Type>
Type dbetabinom(Type y, Type n, Type mu, Type phi, int give_log=0)
{
	Type s1 = mu*phi;
	Type s2 = (Type(1)-mu)*phi;
	Type logres =
	lgamma(n + 1) - lgamma(y + 1)     - lgamma(n - y + 1) +
	lgamma(y + s1) + lgamma(n - y + s2) - lgamma(n + s1 + s2) +
	lgamma(s1 + s2) - lgamma(s1)         - lgamma(s2) ;
	if(!give_log) return exp(logres);
	else return logres;
}


template<class Type>
Type objective_function<Type>::operator() (){
    
    using namespace density;
    using namespace Eigen;
    using namespace R_inla;

    // data:
    DATA_INTEGER(nobs);
    DATA_INTEGER(nsite);
    DATA_INTEGER(nyear);
    DATA_INTEGER(nvessel);
    DATA_INTEGER(nsurvey);
    DATA_IVECTOR(knot); // spatial knot number
    DATA_IVECTOR(site); // site number
    DATA_IVECTOR(year); // year number
    DATA_IVECTOR(vessel); // vessel number
    DATA_IVECTOR(survey); // survey number
    DATA_VECTOR(offset); // offset for tow standardization
    DATA_VECTOR(C); // catch number
    DATA_STRUCT(spde, spde_t); // INLA SPDE 
    DATA_IVECTOR(Options); // User options


    // parameters:
	PARAMETER_VECTOR(log_rho); // relative catchability
	PARAMETER_VECTOR(log_N); // total abundance by year
 	PARAMETER_VECTOR(log_S); // spatial effect
 	PARAMETER_ARRAY(log_ST); // spatiotemporal effect 
    PARAMETER(log_kappa_s);
    PARAMETER(log_tau_s);
    PARAMETER(log_kappa_st);
    PARAMETER(log_tau_st);
    PARAMETER(beta);
    PARAMETER(log_nbk);
    PARAMETER(log_phi);



    // ---------------------------
    // Joint negative log-likelihood
    vector<Type> nll(5); nll.setZero();
    
	// process equation
    vector<Type> mu(nobs); mu.setZero();
    for(int i = 0; i < nobs; i++){
    	Type q = 1; // catchability: vessel == 0, q = 1
    	if(vessel(i) > 0) // relative to vessel == 0
    		q = exp(log_rho(vessel(i)-1));
		// expected catch for each tow
    	mu(i) = q*exp(log_N(year(i))+log_S(knot(i))+log_ST(knot(i),year(i))+offset(i));
    }

    // among knot: spatial and spatiotemporal effect
    SparseMatrix<Type> Q_s = Q_spde(spde, exp(log_kappa_s));
    nll(0) += SCALE(GMRF(Q_s), 1/exp(log_tau_s))(log_S);

	SparseMatrix<Type> Q_st = Q_spde(spde, exp(log_kappa_st));
    if(Options(0) == 0){
    	for(int i_year = 0; i_year < nyear; i_year++){
    		nll(1) += SCALE(GMRF(Q_st), 1/exp(log_tau_st))(vector<Type>(log_ST.col(i_year)));
    	}
    }else if(Options(0) == 1){
    	nll(1) += SCALE(GMRF(Q_st), 1/exp(log_tau_st))(vector<Type>(log_ST.col(0)));
    	for(int i_year = 1; i_year < nyear; i_year++){
    		nll(1) += SCALE(GMRF(Q_st), 1/exp(log_tau_st))(vector<Type>(log_ST.col(i_year)-log_ST.col(i_year-1)));
    	}
    }else if(Options(0) == 2){
    	nll(1) += SEPARABLE(AR1(beta), SCALE(GMRF(Q_st), 1/exp(log_tau_st)))(log_ST);
    }

	// among site: negative binomial for catch
	for(int i_site = 0; i_site < nsite; i_site++){
		Type mu_site = 0;
		Type catch_site = 0;
		// Type zip_site = 0;
		for(int i = 0; i < nobs; i++){
			if(i_site == site(i)){
				mu_site += mu(i);
				catch_site += C(i);
			}
		}
		// combined mu and combined catch for each site
		nll(2) -= dnbinom2(catch_site, mu_site, mu_site+mu_site*mu_site/exp(log_nbk), true);
	}


	// within-site: beta-binomial for pair
	for(int j = 0; j < nobs; j++){
		for(int i = 0; i < nobs; i++){
			if((site(j) == site(i)) & (vessel(j) < vessel(i))){
				// expected catch proportion for each pair
				nll(3) -= dbetabinom(C(j), C(i)+C(j), mu(j)/(mu(i)+mu(j)), exp(log_phi), true);
			}
		}
	}



    // ---------------------------
    // reports:
	REPORT(mu);
	REPORT(log_rho);
	REPORT(log_N);
	REPORT(log_S);
	REPORT(log_ST);
	REPORT(nll);




	
    Type jnll = nll.sum();
    return(jnll);

}
