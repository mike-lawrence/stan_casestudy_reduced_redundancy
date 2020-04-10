library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)

library(ezStan) #for some useful functions; install via remotes::install_github('mike-lawrence/ezStan')

# Read in the data and look at it ----
dat = readr::read_csv('dat.csv')
print(dat)
summary(dat)


# Get within-subjects contrast matrix ----
dat %>%
	ezStan::get_contrast_matrix(
		formula = ~ a
		, contrast_kind = ezStan::halfsum_contrasts
	) ->
	W_full
head(W_full)
nrow(W_full)



# Get between-subjects contrast matrix ----

#collapse to one row per id, keeping any between-subjects variables
dat %>%
	dplyr::group_by(
		id
		, group
	) %>%
	dplyr::summarize() ->
	dat_for_B

#get the full contrast matrix
dat_for_B %>%
	ezStan::get_contrast_matrix(
		formula = ~ group
		, contrast_kind = halfsum_contrasts
	) ->
	B_full
head(B_full)
nrow(B_full)


# package for stan ----
data_for_stan = list(

	# num_trials: number of trials
	  num_trials = nrow(dat)

	# outcome: outcome on each trial
	, outcome = dat$value

	# num_subjects: number of subjects
	, num_subjects = length(unique(dat$id))

	# subj_label: integer label assigning trials to subjects
	, subj_label = as.numeric(factor(dat$id)) #trick to ensure a list of consecutive integers

	# num_cols_W: num cols in within predictor matrix W
	, num_cols_W = ncol(W_full)

	# W: within predictor matrix
	, W = W_full

	# num_cols_B: num cols in between predictor matrix B
	, num_cols_B = ncol(B_full)

	# B: between predictor matrix
	, B = B_full

)

#compile and sample ----
slower_mod = ezStan::build_stan('hmg_slower.stan')
ezStan::start_stan(
	data = data_for_stan
	, mod = slower_mod
	, iter = 4e3
	, seed_start = 1
	, pars = c(
		#list of variables to include in output
		'cor'
		, 'sds'
		, 'coef'
		, 'noise'
	)
	, control = list(
		max_treedepth = 12
	)
)
ezStan::watch_stan()

#collect posterior and summarize ----
post = ezStan::collect_stan()

#how long did it take?
summary(rowSums(get_elapsed_time(post)))

#check diagnostics:
check_hmc_diagnostics(post)

#look at the coefficients
ezStan::stan_summary(
	from_stan = post
	, par = 'coef'
	, W = W_full
	, B = B
)

#look at the sds
ezStan::stan_summary(
	from_stan = post
	, par = 'sds'
	, X = W_full
)

#look at the correlations
ezStan::stan_summary(
	from_stan = post
	, par = 'cor'
	, X = W_full
	, is_cor = T
)
