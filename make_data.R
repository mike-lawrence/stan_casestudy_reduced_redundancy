library(rstan)
rstan_options(auto_write = TRUE)
nSubj = 100
nY = 100
data_parameters = list()
data_parameters$nX = 2
data_parameters$nY = nY
data_parameters$noise = 100
data_parameters$coef_means = c(400,100)
data_parameters$coef_sds = c(100,30)
data_parameters$cor_mat = matrix(c(1,0,0,1),nrow=2,ncol=2)
# 	1,0,0,0,
# 	0,1,0,0,
# 	0,0,1,0,
# 	0,0,0,1
# ),nrow=4,ncol=4)
data_parameters$contrasts = expand.grid(
	trial = 1:nY
	, intercept = 1
	, a = c(-.5,.5)
	#, b = c(-.5,.5)
)
data_parameters$contrasts = data_parameters$contrasts[,-1]
#data_parameters$contrasts$ab = data_parameters$contrasts$a*data_parameters$contrasts$b
data_parameters$contrasts = as.matrix(data_parameters$contrasts)

from_stan = rstan::stan(
	file = 'make_data.stan'
	, data = data_parameters
	, iter = nSubj
	, seed = 1
	, chains = 1
)
dat = rstan::extract(from_stan,inc_warmup=T,permute=F,par='dat')
a = data.frame(matrix(dat,ncol=3)[,2:3])
names(a) = c('a','value')
a$id = 1:nSubj
a$a = factor(a$a)
levels(a$a) = c('a1','a2')


#add a group column & main effect for fun
a %>%
	dplyr::mutate(
		group = case_when(
			id%%2==1 ~ 'group1'
			, T ~ 'group2'
		)
		, value = case_when(
			group=='group1' ~ value+100
			, T ~ value-100
		)
	) %>%
	dplyr::arrange(
		id
		, a
	) %>%
	dplyr::select(
		id
		, group
		, a
		, value
	) ->
	a

readr::write_csv(a,path='dat.csv')

