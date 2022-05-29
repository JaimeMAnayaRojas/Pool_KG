data{
  
  // Survival
  int<lower=1> N_survG;
  int<lower=1> N_growG;
  int<lower=1> N_repG;
  int<lower=1> N_recrG;
  int<lower=1> N_survK;
  int<lower=1> N_growK;
  int<lower=1> N_repK;
  int<lower=1> N_recrK;
  int<lower=1> N_stream;
  
  // Response variables
  int  Surv_G[N_survG];
  real  z_survG[N_survG];
  real  NK_survG[N_survG];
  int  stream_survG[N_survG];
  real  canopy_survG[N_survG];
  real FB_survG[N_survG];
  
  
  real  z1_G[N_growG];
  real  z_growG[N_growG];
  real  NK_growG[N_growG];
  int  stream_growG[N_growG];
  real  canopy_growG[N_growG];
  real FB_growG[N_growG];
  
  
  
  int  Recr_G[N_recrG];
  real  z_recrG[N_recrG];
  real  NK_recrG[N_recrG];
  int  stream_recrG[N_recrG];
  real  canopy_recrG[N_recrG];
  real FB_recrG[N_recrG];
  
  
  int  Surv_K[N_survK];
  real  z_survK[N_survK];
  real  NG_survK[N_survK];
  int  stream_survK[N_survK];
  real  canopy_survK[N_survK];
  real FB_survK[N_survK];
  
  
  
  real  z1_K[N_growK];
  real  z_growK[N_growK];
  real  NG_growK[N_growK];
  int  stream_growK[N_growK];
  real  canopy_growK[N_growK];
  real  FB_growK[N_growK];
  
  
  int  Recr_K[N_recrK];
  real  z_recrK[N_recrK];
  real  NG_recrK[N_recrK];
  int  stream_recrK[N_recrK];
  real  canopy_recrK[N_recrK];
  real  FB_recrK[N_recrK];
  
}

transformed data{
  
  real  grow_G[N_growG];
  real  grow_K[N_growK];
  real  z2_growG[N_growG];
  real  z2_growK[N_growK];
  
  for(i in 1:N_growG){
    grow_G[i] = log(z1_G[i]/z_growG[i]);
    z2_growG[i] = (18 + z_growG[i])^2 - 18^2;
  }
  
  for(i in 1:N_growK){
    grow_K[i] = log(z1_K[i]/z_growG[i]);
    z2_growK[i] = (18 + z_growK[i])^2 - 18^2;
  }
  
  
  
  
}

parameters{
  
  // surv guppies 
  real Intercept_survG;
  // real<lower=1.33, upper=4.86>  Intercept_survG;
  real b_NK_survG;
  //  real<lower=0.04, upper=0.28> b_z_survG;
  real b_z_survG;
  real b_zNK_survG;
  real b_FB_survG;
  real b_canopy_survG;
  vector[N_stream] v_Intercept_survG;
  real<lower=0> sigma_stream_G;
  
  
  // growth guppy
  //  real<lower= (18+2.83), upper= (18+3.62)> Intercept_growG;
  real Intercept_growG;
  real b_NK_growG;
  real b_z_growG;
  real b_zNK_growG;
  vector[N_stream] v_Intercept_growG;
  real<lower=0> sigma_growG;
  real b_canopy_growG;
  real b_FB_growG;
  
  // recruitment
  
  real<lower=0.61, upper=1.8> Intercept_recrG;
  real b_NK_recrG;
  real b_z_recrG;
  real b_zNK_recrG;
  vector[N_stream] v_Intercept_recrG;
  real b_canopy_recrG;
  real b_FB_recrG;  
  
  // Killifish
  
  // surv  
  real<lower=1.39, upper=3.43>  Intercept_survK; // restricted to the survival observed by Bassar 2017
  real b_NG_survK;
  real b_z_survK;
  real b_zNG_survK;
  vector[N_stream] v_Intercept_survK;
  real<lower=0> sigma_stream_K;
  real b_canopy_survK;
  real b_FB_survK;
  
  // growth 
  real Intercept_growK;
  real b_NG_growK;
  real b_z_growK;
  real b_z2_growK;
  real b_zNG_growK;
  vector[N_stream] v_Intercept_growK;
  real<lower=0> sigma_growK;
  real b_canopy_growK;
  real b_FB_growK;
  
  
  real<lower=-1.76, upper= 0.71> Intercept_recrK;
  real b_NG_recrK;
  real b_z_recrK;
  real b_zNG_recrK;
  vector[N_stream] v_Intercept_recrK;
  real b_canopy_recrK;
  real b_FB_recrK;  
  
}


model{
  
  vector[N_survG] p_survG;
  vector[N_growG] mu_growG;
  vector[N_recrG] lambda_G;
  vector[N_survK] p_survK;
  vector[N_growK] mu_growK;
  vector[N_recrK] lambda_K;
  
  // Survival
  
  sigma_stream_G ~ cauchy( 0 , 1 );
  v_Intercept_survG ~ normal( 0 , sigma_stream_G );
  b_zNK_survG ~ normal( 0 , 1 );
  b_z_survG ~ normal( 0 , 1 );
  b_NK_survG ~ normal( 0 , 1 );
  Intercept_survG ~ normal( 2.6 , 3 ); //  
  b_FB_survG ~ normal( 0 , 1 );
  
  // I use priors from Bassar 2017 evolution
  
  
  for ( i in 1:N_survG ) {
    p_survG[i] = Intercept_survG + b_NK_survG * NK_survG[i] + b_z_survG * z_survG[i] + 
    b_zNK_survG* NK_survG[i] * z_survG[i] + v_Intercept_survG[stream_survG[i]] + 
    b_canopy_survG * canopy_survG[i] + 
    b_FB_survG * FB_survG[i];
    
    p_survG[i] = inv_logit(p_survG[i]);
  }
  Surv_G ~ binomial( 1 , p_survG);
  
  // Growth
  
  sigma_growG ~ cauchy( 0.8 , 2 );
  
  v_Intercept_growG ~ normal( 0 , sigma_stream_G );
  b_zNK_growG ~ normal( 0 , 1 );
  b_z_growG ~ normal( 0 , 1 );
  b_NK_growG ~ normal( 0 , 1 );
  Intercept_growG ~ normal( (18) , 10); // Bassar et al., 2017
  
  
  for ( i in 1:N_growG ) {
    mu_growG[i] = Intercept_growG + b_NK_growG * NK_growG[i] + 
    b_z_growG * z_growG[i] + b_zNK_growG * NK_growG[i] * z_growG[i] + 
    v_Intercept_growG[stream_growG[i]]+ 
    b_canopy_growG * canopy_growG[i] + 
    b_FB_growG*FB_growG[i];
  }
  
  z1_G ~ normal( mu_growG , sigma_growG );
  
  
  // Reproduction
  
  v_Intercept_recrG ~ normal( 0 , sigma_stream_G );
  b_zNK_recrG ~ normal( 0 , 1 );
  b_z_recrG ~ normal( 0 , 1 );
  b_NK_recrG ~ normal( 0 , 1 );
  Intercept_recrG ~ normal( 0 , 1 );
  b_FB_recrG ~ normal( 0 , 1 );

  for ( i in 1:N_recrG ) {
    lambda_G[i] = Intercept_recrG + b_NK_recrG * NK_recrG[i] + b_z_recrG * z_recrG[i] + 
    b_zNK_recrG * NK_recrG[i] * z_recrG[i] + 
    v_Intercept_recrG[stream_recrG[i]] + 
    b_canopy_recrG * canopy_recrG[i] + 
    b_FB_recrG * FB_recrG[i];
  }
  
  Recr_G ~ poisson_log(lambda_G);  
  
  
  
  // Killifish
  
  // Survival
  
  sigma_stream_K ~ cauchy( 0 , 2 );
  v_Intercept_survK ~ normal( 0 , sigma_stream_K );
  b_zNG_survK ~ normal( 0 , 1 );
  b_z_survK ~ normal( 0 , 1 );
  b_NG_survK ~ normal( 0 , 1 );
  Intercept_survK ~ normal( 2, 2 ); // killifish survival is atleast that of guppies
  b_FB_survK ~ normal( 0 , 1 );
  
  for ( i in 1:N_survK ) {
    p_survK[i] = Intercept_survK + b_NG_survK * NG_survK[i] + b_z_survK * z_survK[i] + 
    b_zNG_survK* NG_survK[i] * z_survK[i] + v_Intercept_survK[stream_survK[i]] + 
    b_canopy_survK * canopy_survK[i] + 
    b_FB_survK * FB_survK[i];
    
    p_survK[i] = inv_logit(p_survK[i]);
  }
  
  Surv_K ~ binomial( 1 , p_survK);
  
  // K growth
  
  sigma_growK ~ cauchy( 0 , 1 );
  v_Intercept_growK ~ normal( 0 , sigma_stream_K );
  b_zNG_growK ~ normal( 0 , 1 );
  b_z_growK ~ normal( 0 , 1 );
  b_z2_growK ~ normal( 0 , 1 );
  b_NG_growK ~ normal( 0 , 1 );
  Intercept_growK ~ normal( 18 , 10 );
  b_FB_growK ~ normal( 0 , 1 );

  for ( i in 1:N_growK ) {
    mu_growK[i] = Intercept_growK + b_NG_growK * NG_growK[i] + 
    b_z_growK * z_growK[i] + b_z2_growK * z2_growK[i]  + b_zNG_growK * NG_growK[i] * z_growK[i] + 
    v_Intercept_growK[stream_growK[i]] + 
    b_canopy_growK * canopy_growK[i] + 
    b_FB_growK*FB_growK[i];
  }
  
  z1_K ~ normal( mu_growK , sigma_growK );
  
  
  // Reproduction
  
  v_Intercept_recrK ~ normal( 0 , sigma_stream_K );
  b_zNG_recrK ~ normal( 0 , 1 );
  b_z_recrK ~ normal( 0 , 1 );
  b_NG_recrK ~ normal( 0 , 1 );
  Intercept_recrK ~ normal( 0 , 2 );
  b_FB_recrK ~ normal( 0 , 1 );
  
  for ( i in 1:N_recrK ) {
    lambda_K[i] = Intercept_recrK + b_NG_recrK * NG_recrK[i] + b_z_recrK * z_recrK[i] + 
    b_zNG_recrK * NG_recrK[i] * z_recrK[i] + 
    v_Intercept_recrK[stream_recrK[i]] + 
    b_canopy_recrK * canopy_recrK[i] + 
    b_FB_recrK * FB_recrK[i];
  }
  
  Recr_K ~ poisson_log(lambda_K);
  
}
