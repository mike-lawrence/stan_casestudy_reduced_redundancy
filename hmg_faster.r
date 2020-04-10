library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)

library(ezStan) #for some useful functions; install via remotes::install_github('mike-lawrence/ezStan')

# Read in the data and look at it ----
dat = readr::read_csv('dat.csv')
print(dat)
summary(dat)


# Get within-subjects contrast matrix and helpers ----
dat %>%
	ezStan::get_contrast_matrix(
		formula = ~ a
		, contrast_kind = ezStan::halfsum_contrasts
	) ->
	W_full
head(W_full)
nrow(W_full)

# get the unique entries in W
W_full %>%
	tibble::as_tibble() %>%
	dplyr::distinct() %>%
	as.matrix() ->
	W_unique
print(W_unique)

#for each unique condition specified by W_unique, the stan model will
# work out values for that condition for each subject, and we'll need to index
# into the resulting subject-by-condition matrix. So we need to create our own
# subject-by-condition matrix and get the indices of the observed data into a
# the array produced when that matrix is flattened.

#first copy W_unique for each id, yielding the same order as the
# the flattened subject-by-condition matrix that we'll eventually index into
W_unique %>%
	tibble::as_tibble() %>%
	#first repeat the matrix so there's a copy for each subject
	dplyr::slice(
		rep(
			dplyr::row_number()
			, length(unique(dat$id))
		)
	) %>%
	#now add the subject labels
	dplyr::mutate(
		id = rep(sort(unique(dat$id)),each=nrow(W_unique))
	) ->
	W_unique_per_id

#Now find the row index in W_unique_per_id corresponding to each trial
W_unique_per_id %>%
	dplyr::mutate(
		row = 1:n()
	) %>%
	dplyr::right_join(
		tibble::as_tibble(W_full) %>%
			dplyr::mutate(id=dat$id)
	) %>%
	dplyr::pull(row) ->
	trial_index_in_flattened_value_for_subj_cond

# Get between-subjects contrast matrix and helpers ----

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

# get the unique entries in B_full
B_full %>%
	tibble::as_tibble() %>%
	dplyr::distinct() %>%
	as.matrix() ->
	B_unique
print(B_unique)

#Now find the row index in B_unique_per_id corresponding to each subject
B_unique %>%
	tibble::as_tibble() %>%
	dplyr::mutate(
		row = 1:n()
	) %>%
	dplyr::right_join(
		tibble::as_tibble(B_full)
	) %>%
	dplyr::pull(row) ->
	subject_row_in_B_unique


# package for stan ----
data_for_stan = list(

	# num_subjects: number of subjects
	num_subjects = length(unique(dat$id))

	# num_trials: number of trials
	, num_trials = nrow(dat)

	# outcome: outcome on each trial
	, outcome = dat$value

	# num_rows_W: num rows in within predictor matrix W
	, num_rows_W = nrow(W_unique)

	# num_cols_W: num cols in within predictor matrix W
	, num_cols_W = ncol(W_unique)

	# W: within predictor matrix
	, W = W_unique

	# trial_index_in_flattened_value_for_subj_cond: index of each trial in flattened version of subject-by-condition value matrix
	, trial_index_in_flattened_value_for_subj_cond = trial_index_in_flattened_value_for_subj_cond

	# num_rows_B: num rows in between predictor matrix B
	, num_rows_B = nrow(B_unique)

	# num_cols_B: num cols in between predictor matrix B
	, num_cols_B = ncol(B_unique)

	# B: between predictor matrix
	, B = B_unique

	# subject_row_in_B_unique: row for each subject in between-predictor matrix B
	, subject_row_in_B_unique = subject_row_in_B_unique

)

#compile and sample ----
faster_mod = ezStan::build_stan('hmg_faster.stan')
ezStan::start_stan(
	data = data_for_stan
	, mod = faster_mod
	, iter = 4e3
	, seed_start = 1
	, pars = c(
		#list of variables to include in output
		'cor'
		, 'sds'
		, 'coef'
		, 'noise'
	)
	# , iter = 1e
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
	, W = W_unique
	, B = B_unique
)

#look at the sds
ezStan::stan_summary(
	from_stan = post
	, par = 'sds'
	, X = W_unique
)

#look at the correlations
ezStan::stan_summary(
	from_stan = post
	, par = 'cor'
	, X = W_unique
	, is_cor = T
)
