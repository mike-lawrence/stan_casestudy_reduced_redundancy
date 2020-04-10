data{
	int nX ;
	int nY ;
	real noise ;
	vector[nX] coef_means ;
	vector[nX] coef_sds ;
	corr_matrix[nX] cor_mat ;
	matrix[nX*nY,nX] contrasts ;
}
parameters{
  vector[nX] coefs ; #observed intercept & effect
}
model{
  #observed intercept & effect as multivariate normal
  coefs ~ multi_normal(
      coef_means
    , quad_form_diag( cor_mat , coef_sds )
  ) ;
}
generated quantities{
  matrix[nX*nY,nX+1] dat ;
  dat[,1:nX] = contrasts ; #save the contrasts used
  for(i in 1:(nX*nY)){
	  dat[i,nX+1] = normal_rng(contrasts[i,] * coefs,noise) ; #compute observed vals
  }
}
