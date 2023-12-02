library(mvtnorm)
library(condMVNorm)
library(gtools)
library(mclust)
library(stringr)
library(MASS)
library(foreach)
library(doParallel)
library(pgmm)
library(sn)






initial_skew_lnm=function(df,G,q){
  n=dim(df)[1]
  K=dim(df)[2]-1
  # if(para_lamb=="D"){
  #   q=K
  # }else{q=5}
    
  
  V=vector("list",G)
  m=vector("list",G)
  ti_V=vector("list",G)
  ti_m=vector("list",G)
  z_ig=matrix(0,nrow=n,ncol=G)
  pi_g=rep(0,G)
  mu_g=list()
  sig_g=list()
  lam_g=list()
  


    Z=df/rowSums(df)
    Z[Z==0]=0.001
    Y=as.matrix(log(Z[,1:K]/c(Z[,K+1])),ncol=K)
    
    #pgmm start
    # res_class=pgmmEM(Y,rG=G,rq=1:round(K+0.5-sqrt(2*K+0.25)),zstart=2,
    #                  modelSubset=c("UUU"),relax=TRUE)
    
    res_class=NULL
    
    if(is.null(res_class)|is.function(res_class)){
      res_class=kmeans(Y,G)
      lab=res_class$cluster
    }else{lab=res_class$map}
    
    
    #initial for z_ig
    for(i in 1:n){
      for(g in 1:G){
        z_ig[i,lab[i]]=1
      }
    }
    
    
    #initial for pi_g
    pi_g=table(lab)/n
    
    
    #initial for mu_g
    for (g in 1:G) {
      #mu_g[[g]]=rep(0,Q_g[g])
      mu_g[[g]]=colMeans(Y[lab==g,])
    }
    
    #initial for sig_g
    for (g in 1:G) {
      sig_g[[g]]=cov(Y[lab==g,])
    }
    
    
    #initial for lam_g
    set.seed(1029)
    for(g in 1:G){
      lam_g[[g]]=matrix(rnorm(q*K,0,0.1),nrow=K)
    }
    set.seed(NULL)
    
    
    
  
    #initial for m
    for (g in 1:G) {
      m[[g]]=t(Y)
    }
    
    
    
    #inital for V
    for (g in 1:G) {
      #vg=sig_g[[g]]
      for (i in 1:n) {
        #V[[g]][[i]]=diag(diag(vg))
        V[[g]][[i]]=diag(0.1,K)
      }
    }
    
    
    
    #initial for ti_V
    for(g in 1:G){
      bbb=solve(diag(1,q)+t(lam_g[[g]])%*%ginv(sig_g[[g]])%*%lam_g[[g]])
      for(i in 1:n){
      ti_V[[g]][[i]]=bbb
      }
    }
    
    
    
    
    #initial for ti_m
    for(g in 1:G){
      ti_m[[g]]=t(lam_g[[g]])%*%ginv(sig_g[[g]]+lam_g[[g]]%*%t(lam_g[[g]]))%*%
        (m[[g]]-c(mu_g[[g]]))
    }
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    return(list(new_pi_g=pi_g,
                new_z_ig=z_ig,
                new_m=m,
                new_V=V,
                new_ti_m=ti_m,
                new_ti_V=ti_V,
                new_mu_g=mu_g,
                new_sig_g=sig_g,
                new_lam_g=lam_g))
  
  
}





#main algorithm
Mico_skew_lnm=function(W,G,pi_g,z_ig,mu_g,sig_g,lam_g,V,m,ti_V,ti_m,const){
  n=dim(W)[1]
  K=dim(W)[2]-1
  q=dim(ti_V[[1]][[1]])[1]
  ti_sig=vector("list",G)

  # const=NULL
  # for (i in 1:dim(W)[1]) {
  #   const[i]=constant_fun(W[i,],dim(W)[2]-1)
  # }
  
  overall_loglik=rep(-Inf,2)
  
  h=2
  
  repeat{
    #print(h)
    #update m and V
    for (g in 1:G) {
      for (i in 1:n) {
        mv=GD_mvxi(V[[g]][[i]],ti_V[[g]][[i]],m[[g]][,i],ti_m[[g]][,i],sig_g[[g]],
                   lam_g[[g]],W[i,],mu_g[[g]],K)
        m[[g]][,i]=c(mv$m)
        V[[g]][[i]]=mv$V
      }
    }

    
    #print("GD_mvxi")
    # print(p_comp_logw(W,m,ti_m,V,ti_V,mu_g,pi_g,sig_g,lam_g,G,n,K,const,z_ig))
    for(g in 1:G){
      #ti_sig[[g]]=ginv(diag(1,q)+t(lam_g[[g]])%*%ginv(sig_g[[g]])%*%lam_g[[g]])
      ti_sig[[g]]=diag(1,q)-t(lam_g[[g]])%*%ginv(sig_g[[g]]+lam_g[[g]]%*%t(lam_g[[g]]))%*%lam_g[[g]]
    }
    
    
    
    #update ti_m and ti_V
    for (g in 1:G) {
      for (i in 1:n) {
        ti_mv=Eti_mv(V[[g]][[i]],ti_V[[g]][[i]],m[[g]][,i],ti_m[[g]][,i],sig_g[[g]],
                     lam_g[[g]],W[i,],mu_g[[g]],K,ti_sig[[g]])
        ti_m[[g]][,i]=c(ti_mv$ti_m)
        ti_V[[g]][[i]]=ti_mv$ti_V

      }
      #print(t(ti_m[[g]]))
    }
    
    
    
    
    #calculate tilde_F
    log_w_approx=matrix(0,nrow = n,ncol=G)
    for (g in 1:G) {
      for (i in 1:n) {
        log_w_approx[i,g]=
          lb_loglik(V[[g]][[i]],ti_V[[g]][[i]],K,sig_g[[g]],lam_g[[g]],m[[g]][,i],
                    ti_m[[g]][,i],mu_g[[g]],W[i,],const[i],ti_sig[[g]])
      }
    }
    #update z_ig
    z_ig=z_ig_update(pi_g,n,G,log_w_approx)
    
    # print("z_ig")
    # print(p_comp_logw(W,m,ti_m,V,ti_V,mu_g,pi_g,sig_g,lam_g,G,n,K,const,z_ig))
    #update m and V and ti_m and ti_V
    
    
    
    

    
    
    
    
    
    #update pi_g
    pi_g=colMeans(z_ig)
    label=apply(z_ig, 1, which.max)
    if(length(pi_g[pi_g==0])>0|any(table(label)<3)|length(unique(label))<G){break}

    # if(overall_loglik[h-1]<
    #    p_comp_logw(W,m,ti_m,V,ti_V,mu_g,new_pi_g,sig_g,lam_g,G,n,K,const,z_ig))
    # {
    #   pi_g=new_pi_g
    # }
    
    # print("pi_g")
    #print(p_comp_logw(W,m,ti_m,V,ti_V,mu_g,pi_g,sig_g,lam_g,G,n,K,const,z_ig))
    
    
    
    
    
    
    #update mu_g
    for (g in 1:G) {
      mu_g[[g]]=(m[[g]]-lam_g[[g]]%*%ti_m[[g]])%*%z_ig[,g]/sum(z_ig[,g])
    }
    
    # #use true mu
    # set.seed(123)
    # mu1=rnorm(K,-0.3,0.01)
    # set.seed(223)
    # mu2=rnorm(K,0.3,0.01)
    # mu_g[[1]]=mu1;mu_g[[2]]=mu2
    
    # if(overall_loglik[h-1]<
    #    p_comp_logw(W,m,ti_m,V,ti_V,new_mu_g,pi_g,sig_g,lam_g,G,n,K,const,z_ig))
    # {
    #   mu_g=new_mu_g
    # }

    # print("mu_g")
    # print(p_incomp_logw(W,m,ti_m,V,ti_V,mu_g,pi_g,sig_g,lam_g,G,n,K,const))
    
    
    
    
    
    
    #update lam_g
    for(g in 1:G){
    fenzi=0
    fenmu=0
    for(i in 1:n){
      fenzi=fenzi+z_ig[i,g]*(m[[g]][,i]-mu_g[[g]])%*%t(ti_m[[g]][,i])
      fenmu=fenmu+z_ig[i,g]*(ti_m[[g]][,i]%*%t(ti_m[[g]][,i])+ti_V[[g]][[i]])
    }
    # if(dim(lam_g[[g]])[1]==dim(lam_g[[g]])[2]){
    #   lam_g[[g]]=diag(c(ginv(ginv(sig_g[[g]])*fenmu)%*%rowSums(ginv(sig_g[[g]])*t(fenzi))))
    #   p_lamb=G*K
    # }else{
      lam_g[[g]]=fenzi%*%ginv(fenmu)
      
    #}
    }
    p_lamb=G*K*q
    
    
    # #try true lambda
    # set.seed(123)
    # Lamb1=matrix(rnorm(K*q,0,0.3),nrow=K)
    # set.seed(223)
    # Lamb2=matrix(rnorm(K*q,0.8,0.1),nrow=K)
    # lam_g[[1]]=Lamb1;lam_g[[2]]=Lamb2
    
    
    # if(overall_loglik[h-1]<
    #    p_comp_logw(W,m,ti_m,V,ti_V,mu_g,pi_g,sig_g,new_lam_g,G,n,K,const,z_ig))
    # {
    #   lam_g=new_lam_g
    # }
    
    # print("lam_g")
    # print(p_incomp_logw(W,m,ti_m,V,ti_V,mu_g,pi_g,sig_g,lam_g,G,n,K,const))
    

    
    
    
    
    
    #update sig_g
    for(g in 1:G){
    u_s=0
    for(i in 1:n){
      u_s=u_s+z_ig[i,g]*(V[[g]][[i]]+
                           (m[[g]][,i]-mu_g[[g]]-lam_g[[g]]%*%ti_m[[g]][,i])%*%
                           t(m[[g]][,i]-mu_g[[g]]-lam_g[[g]]%*%ti_m[[g]][,i])+
                           lam_g[[g]]%*%ti_V[[g]][[i]]%*%t(lam_g[[g]]))
    }
      sig_g[[g]]=u_s/sum(z_ig[,g])
    }
    p_sig=G*(K+1)*K/2
    
    
    # #try true sig
    # set.seed(123)
    # sig1=diag(runif(K,0,0.03))+0.03
    # set.seed(223)
    # sig2=diag(runif(K,0,0.05))+0.01
    # sig_g[[1]]=sig1;sig_g[[2]]=sig2
    
    
    #print(sig_g)
    
    # if(overall_loglik[h-1]<
    #    p_comp_logw(W,m,ti_m,V,ti_V,mu_g,pi_g,new_sig_g,lam_g,G,n,K,const,z_ig))
    # {
    #   sig_g=new_sig_g
    # }
    
    # print("sig_g")
    # print(p_incomp_logw(W,m,ti_m,V,ti_V,mu_g,pi_g,sig_g,lam_g,G,n,K,const))
    

    
    
    
    o_g=order(pi_g,decreasing = T)
    mu_g=mu_g[o_g]
    pi_g=pi_g[o_g]
    z_ig=z_ig[,o_g]
    m=m[o_g]
    V=V[o_g]
    ti_m=ti_m[o_g]
    ti_V=ti_V[o_g]
    sig_g=sig_g[o_g]
    lam_g=lam_g[o_g]
    ti_sig=ti_sig[o_g]

    if(is.vector(z_ig)){lab_rec=rep(1,n)}
    else{lab_rec=apply(z_ig,1,which.max)}
    
    overall_loglik[h]=p_incomp_logw(W,m,ti_m,V,ti_V,mu_g,pi_g,sig_g,lam_g,G,n,K,const,ti_sig)

    EN2=0
    for(j in 1:G){
    for(i in 1:n){
      EN2_new=2*sum(W[i,])*
        log(sum(exp(mu_g[[j]]+lam_g[[j]]*ti_m[[j]][,i]+
              as.vector(lam_g[[j]])^2*ti_V[[j]][[i]]/2))/
              sum(exp(mu_g[[j]]+lam_g[[j]]*ti_m[[j]][,i]))
            )
      EN2=EN2+EN2_new
    }
    }
    
    
    #print(overall_loglik[h])
    
    if(is.na(overall_loglik[h])){overall_loglik[h]=overall_loglik[h-1]}
    
    if(h>3){
      a=(overall_loglik[h]-overall_loglik[h-1])/(overall_loglik[h-1]-overall_loglik[h-2])
      L_inf=overall_loglik[h-1]+(overall_loglik[h]-overall_loglik[h-1])/(1-a)}else{L_inf=Inf}
    
    diff=L_inf-overall_loglik[h]
    if(is.na(diff)){diff=Inf}
    if((abs(diff)<10^(-2))&h>100|h>200|length(pi_g[pi_g==0])>0){break}
    h=h+1
    
    
  }
  
  if(length(pi_g[pi_g==0])>0|any(table(label)<3)|length(unique(label))<G){return("fewer class than define")}
  else{
    ICL=overall_loglik[h]*2-log(n)*(G*K+G-1+(p_lamb+p_sig))+2*sum(as.vector(z_ig)*log(as.vector(z_ig)))
    BIC=overall_loglik[h]*2-log(n)*(G*K+G-1+(p_lamb+p_sig))
    AIC=overall_loglik[h]*2-2*(G*K+G-1+(p_lamb+p_sig))
    HIC=overall_loglik[h]*2-log(n)*(G*K+G-1+(p_lamb+p_sig)+EN2)
    
    return(list(z_ig=z_ig,
                cluster=lab_rec,
                mu_g=mu_g,
                pi_g=pi_g,
                sig_g=sig_g,
                lam_g=lam_g,
                m=m,
                V=V,
                overall_loglik=overall_loglik,
                ICL=ICL,
                BIC=BIC,
                AIC=AIC,
                HIC=HIC))
  }
  
  
}




#calculate approx log liklyhood for i and g observation
#constant for each i
constant_fun=function(W,K){
  funforlb <- function(k){if(W[k]==0){return(0)}else{return(sum(log(1:W[k])))}}
  const <- sum(log(1:sum(W)))-sum(sapply(1:(K+1),funforlb))
  const
}


#tilda_F
lb_loglik=function(V,ti_V,K,sig,lam,m,ti_m,mu,W,const,ti_sig){
  
  #### q(y) is gaussian, but q(tau) is truncated normal ####
  #parmeter for truncated normal
  ti_mu=t(lam)%*%ginv(sig+lam%*%t(lam))%*%(m-mu)
  q=length(ti_mu)
  alpha=pmvnorm(lower=rep(0,q),upper=Inf, mean=c(ti_mu),sigma=ti_sig)

 wml=
  log(det(ti_sig))/2+log(2*pi)*q/2+log(alpha)+
   tr((diag(1,q)+t(lam)%*%ginv(sig)%*%lam)%*%
        (ti_V+ti_m%*%t(ti_m)+ti_mu%*%t(ti_mu)-2*ti_mu%*%t(ti_m)))/2+
    log(det(V))/2+K/2+K*log(2*pi)/2-
    0.5*(log(det(sig))+K*log(2*pi)+tr(ginv(sig)%*%V)+
           t(m-mu-lam%*%ti_m)%*%ginv(sig)%*%(m-mu-lam%*%ti_m)+
           tr(t(lam)%*%ginv(sig)%*%lam%*%ti_V))+
    K*log(2/pi)/2-0.5*(t(ti_m)%*%ti_m+tr(ti_V))+
    W[1:K]%*%m-sum(W)*(log(sum(c(exp(m+diag(V)/2),1))))+const
  #
  
  
  if(is.na(wml)){wml=-Inf}
  return(wml)
  
}


#incomplete log likelihood
elbo_fun = function(pi_g,G,n,log_w_approx){
  forelbo = matrix(nrow=n,ncol=G)
  for(i in 1:n){
    for(g in 1:G){
      if(pi_g[g]*exp(log_w_approx[i,g])==0){
        forelbo[i,g]<-pi_g[g]*.Machine$double.xmin
      }else if(pi_g[g]*exp(log_w_approx[i,g])==Inf){
        forelbo[i,g]<-pi_g[g]*.Machine$double.xmax
        }else{forelbo[i,g]=pi_g[g]*exp(log_w_approx[i,g])}
      
    }
    if(sum(forelbo[i,])==0){
      forelbo[i,]<-pi_g*.Machine$double.xmin
    }else if(sum(forelbo[i,])==Inf){
      forelbo[i,]<-pi_g*.Machine$double.xmax
    }
    
  }
  
  mix_den=rowSums(forelbo)
  return(sum(log(mix_den)))
}



#calculate incomplete likelihood at the end of each iteration
p_incomp_logw=function(W,m,ti_m,V,ti_V,mu_g,pi_g,sig_g,lam_g,G,n,K,const,ti_sig){
  #calculate overall log likelyhood.
  log_w_approx=matrix(0,nrow = n,ncol=G)
  for (g in 1:G) {
    for (i in 1:n) {
      log_w_approx[i,g]=
        lb_loglik(V[[g]][[i]],ti_V[[g]][[i]],K,sig_g[[g]],lam_g[[g]],m[[g]][,i],ti_m[[g]][,i],
                  mu_g[[g]],W[i,],const[i],ti_sig[[g]])
    }
  }
  return(elbo_fun(pi_g,G,n,log_w_approx))
}




#calculate complete likelihood at the end of each parameter updates
p_comp_logw=function(W,m,ti_m,V,ti_V,mu_g,pi_g,sig_g,lam_g,G,n,K,const,z_ig,ti_sig){
  #calculate overall log likelyhood.
  log_w_approx=matrix(0,nrow = n,ncol=G)
  for (g in 1:G) {
    for (i in 1:n) {
      log_w_approx[i,g]=
        lb_loglik(V[[g]][[i]],ti_V[[g]][[i]],K,sig_g[[g]],lam_g[[g]],m[[g]][,i],ti_m[[g]][,i],
                  mu_g[[g]],W[i,],const[i],ti_sig[[g]])
    }
  }
  forelbo = matrix(nrow=n,ncol=G)
  for(i in 1:n){
    for(g in 1:G){
      if(pi_g[g]*exp(log_w_approx[i,g])==0){forelbo[i,g]<-pi_g[g]*.Machine$double.xmin}
      else{forelbo[i,g]=pi_g[g]*exp(log_w_approx[i,g])}
      
    }
  }
  
  
  return(sum(log(forelbo)*z_ig))
}





#one step newtown raphson to update m and V, and close form for ti_m and ti_V
GD_mvxi=function(V,ti_V,m,ti_m,sig_g,lam_g,W,mu_g,K){
  V_old=sqrt(V)
  m_old=m
  
  
  f_V=diag(V_old)^(-1)-diag(ginv(sig_g))*diag(V_old)-
    c(sum(W)*diag(V_old)*c(exp(m_old+diag(V_old^2)/2)/
                                   sum(c(exp(m_old+diag(V_old^2)/2),1))))
  s_V=-1/diag(V_old)^2-diag(ginv(sig_g))-(diag(V_old)^2+1)*
    sum(W)*exp(m_old+diag(V_old)^2/2)/sum(c(exp(m_old+diag(V_old^2)/2),1))
  
  
  if(any(is.na(f_V/s_V))){
    V_new=V_old
    }else{
    V_new=V_old-diag(f_V/s_V)
  }
  
  
  f_m=c(W[1:K])-ginv(sig_g)%*%(m_old-mu_g)-
    c(sum(W)*c(exp(m_old+diag(V_old^2)/2)/
                                       sum(c(exp(m_old+diag(V_old^2)/2),1))))+
    ginv(sig_g)%*%lam_g%*%ti_m
  
  s_m=-ginv(sig_g)-
    sum(W)*(diag(exp(m_old+diag(V_old^2)/2))/
                    sum(c(exp(m_old+diag(V_old^2)/2),1)))
  
  if(det(s_m)==0|is.na(det(s_m))|det(s_m)==Inf){
    m_new=m_old
  }else{
    m_new=m_old-ginv(s_m)%*%f_m
  }
  
  list(m=m_new,V=V_new^2)
}



Eti_mv=function(V,ti_V,m,ti_m,sig_g,lam_g,W,mu_g,K,ti_sig){
  ti_V_old=ti_V
  ti_m_old=ti_m
  
  
  #### q(y) is gaussian, but q(tau) is truncated normal ####
  #where ti_sig and ti_mu are parameters for it.
  ti_mu=t(lam_g)%*%ginv(sig_g+lam_g%*%t(lam_g))%*%(m-mu_g)
  q=length(ti_mu)
  alpha=pmvnorm(lower=rep(0,q),upper=Inf,mean=c(ti_mu),sigma=ti_sig)
  

  if(q==1){
    qq=dnorm(0,mean=ti_mu,sd=sqrt(ti_sig))
  }else{
    qq=sapply(c(1:q),
                function(x){
                  dnorm(0,mean=ti_mu[x],sd=sqrt(ti_sig[x,x]))*
                  pcmvnorm(lower=0,upper=Inf,
                          mean=ti_mu,
                          sigma=round(ti_sig,7),
                          dependent.ind = c(1:q)[-x],
                          given.ind = x,X.given = 0)
                })
  }
  
  ti_m_new=ti_mu+ti_sig%*%qq/alpha
  
  
  if(any(is.na(ti_m_new)|ti_m_new==Inf)){
    ti_m_new=ti_m_old
  }

  

  if(q==1){
    value=0
    }else if(q==2){
    value=mvtnorm::dmvnorm(rep(0,q),mean = ti_mu,sigma=ti_sig)
  }else{
  value=apply(combinations(q,2,c(1:q)),1,function(x){
    i=x[1];j=x[2]
    mvtnorm::dmvnorm(c(0,0),mean = c(ti_mu[c(i,j)]),sigma=ti_sig[c(i,j),c(i,j)])*
      pcmvnorm(lower=0,upper=Inf,
              mean=ti_mu,
              sigma=round(ti_sig,7),
              dependent.ind = c(1:q)[-c(i,j)],
              given.ind = c(i,j),X.given = rep(0,2))
  })
  }
  
  

  HH=matrix(0,q,q)
  if(q==1){
    HH=value
  }else{
  HH[lower.tri(HH)]=value
  HH=HH+t(HH)
  }

  DD=diag((-ti_mu*qq-diag(ti_sig%*%HH))/diag(ti_sig))
  
  
  ti_V_new=ti_sig+ti_sig%*%(HH+DD)%*%ti_sig/alpha-(ti_sig%*%qq/alpha)%*%t(ti_sig%*%qq/alpha)
  
  
  #ti_V_new=ti_sig+ti_sig*(-ti_mu)*qq/alpha-(ti_sig%*%qq/alpha)%*%t(ti_sig%*%qq/alpha)
  
  
  if(any(is.na(ti_V_new)|ti_V_new==Inf)){
    ti_V_new=ti_V_old
  }
  
  list(ti_m=ti_m_new,ti_V=ti_V_new)
}



#update z_ig
z_ig_update=function(pi_g,n,G,log_w_approx){
  z_ig = matrix(0,nrow=n,ncol=G)
  for(i in 1:n){
    for (g in 1:G) {
      if(pi_g[g]*exp(log_w_approx[i,g])==0){
        z_ig[i,g]<-pi_g[g]*.Machine$double.xmin
      }else if(pi_g[g]*exp(log_w_approx[i,g])==Inf){
        z_ig[i,g]<-pi_g[g]*.Machine$double.xmax
         }else{z_ig[i,g]=pi_g[g]*exp(log_w_approx[i,g])}
    }
    if(sum(z_ig[i,])==0){
      z_ig[i,]<-pi_g*.Machine$double.xmin
    }else if(sum(z_ig[i,])==Inf){
      z_ig[i,]<-pi_g*.Machine$double.xmax
    }
  }
  
  z_ig=z_ig/rowSums(z_ig)
  return(z_ig)
}




####trace function####
tr=function(x){sum(diag(x))}

#function to store parameter to table
store_table=function(para,name,procid){
  output.filename = paste("simulationOutput_",name,".csv",sep='')
  write.table(t(para),file=output.filename,row.names=paste("simulation",procid),
              col.names=F,sep=",",append = TRUE)
}















