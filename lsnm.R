source("microb_skew.R")







lsnm=function(W_count,G){
  
  #algorithm to get the result
  const=NULL
  for (i in 1:dim(W_count)[1]) {
    const[i]=constant_fun(W_count[i,],dim(W_count)[2]-1)
  }

    initial_value=initial_skew_lnm(W_count,G,1)
    pi_g=initial_value$new_pi_g
    z_ig=initial_value$new_z_ig
    m=initial_value$new_m
    V=initial_value$new_V
    ti_m=initial_value$new_ti_m
    ti_V=initial_value$new_ti_V
    mu_g=initial_value$new_mu_g
    sig_g=initial_value$new_sig_g
    lam_g=initial_value$new_lam_g
    res=Mico_skew_lnm(W=W_count,G=g,z_ig=z_ig,pi_g = pi_g,m=m,V=V,ti_m = ti_m,ti_V=ti_V,
                      mu_g = mu_g,sig_g=sig_g,lam_g = lam_g,const=const)
    return(res)
}
