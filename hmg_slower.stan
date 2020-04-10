data{

	// num_trials: number of trials
	int<lower=1> num_trials ;

	// outcome: outcome on each trial
	vector[num_trials] outcome ;

	// num_subjects: number of subjects
	int<lower=1> num_subjects ;

	// subj_label: integer label assigning trials to subjects
	int<lower=1> subj_label[num_trials] ;

	// num_cols_W: num cols in within predictor matrix W
	int<lower=1> num_cols_W ;

	// W: within predictor matrix
	matrix[num_trials,num_cols_W] W ;

	// num_cols_B: num cols in between predictor matrix B
	int<lower=1> num_cols_B ;

	// B: between predictor matrix
	matrix[num_subjects,num_cols_B] B ;


}
transformed data{

	// outcome_mean: mean outcome value
	real outcome_mean = mean(outcome) ;

	// outcome_sd: sd of outcomes
	real outcome_sd = sd(outcome) ;

	// scaled_outcome: outcomes scaled to have zero mean and unit variance
	vector[num_trials] scaled_outcome = (outcome-outcome_mean)/outcome_sd ;

}
parameters{

	// multi_normal_helper: a helper variable for implementing non-centered parameterization
	matrix[num_cols_W,num_subjects] multi_normal_helper ;

	// cor_helper: population-level correlations (on cholesky factor scale) amongst within-subject predictors
	cholesky_factor_corr[num_cols_W] cor_helper ;

	// scaled_sds: population-level sds for each within-subject predictor
	vector<lower=0>[num_cols_W] scaled_sds ;

	// scaled_coef: coefficients for between and within subject predictors
	matrix[num_cols_B,num_cols_W] scaled_coef ;

	// scaled_noise: measurement noise
	real<lower=0> scaled_noise ;

}
model{

	////
	// Priors
	////

	// multi_normal_helper must have normal(0,1) prior for non-centered parameterization
	to_vector(multi_normal_helper) ~ std_normal() ;

	// relatively flat prior on correlations
	cor_helper ~ lkj_corr_cholesky(2) ;

	// normal(0,1) priors on all scaled_sds
	scaled_sds ~ std_normal() ;

	// normal(0,1) priors on all coefficients
	to_vector(scaled_coef) ~ std_normal() ;

	// low-near-zero prior on measurement noise
	scaled_noise ~ weibull(2,1) ; // weibull(2,1) is peaked around .8


	{// local environment to permit putting variable declarations closer to their computation

		// subj_devs: subject-by-subject deviations from mean for each within-subject predictor
		matrix[num_subjects,num_cols_W] subj_devs = transpose(diag_pre_multiply(scaled_sds,cor_helper) * multi_normal_helper) ;

		// compute coefficients for each subject/condition
		matrix[num_subjects,num_cols_W] subj_coefs = B * scaled_coef + subj_devs ;

		// Likelihood
		scaled_outcome ~ normal(
			rows_dot_product(
				subj_coefs[subj_label]
				, W
			)
			, scaled_noise
		) ;

	}

}
generated quantities{

	// cor: correlation matrix for the full set of within-subject predictors
	corr_matrix[num_cols_W] cor = multiply_lower_tri_self_transpose(cor_helper) ;

	// sds: population-level sds for each within-subject predictor
	vector[num_cols_W] sds = scaled_sds*outcome_sd ;

	// coef: coefficients for between and within subject predictors
	matrix[num_cols_B,num_cols_W] coef = scaled_coef*outcome_sd ;

	// noise: measurement noise
	real noise = scaled_noise*outcome_sd ;

	// adding back the mean to the intercept only
	coef[1,1] += outcome_mean ;

	// tweak cor
	for(i in 1:num_cols_W){
		cor[i,i] += uniform_rng(1e-16, 1e-15) ;
	}

}
