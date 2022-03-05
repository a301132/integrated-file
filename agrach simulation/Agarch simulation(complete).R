## 최종완성본


library(fGarch)
library(ggplot2)


########################################################################################################################

sample_agarch<-function(omega0,alpha0,beta0,gamma0,omega1,alpha1,beta1,gamma1,n,eta_dis,changepoint){
  
  ran<-matrix(0,n+500,2)
  x<-rep(0,n+500);x[1]<-1
  sigma2<-rep(0,n+500);sigma2[1]<-1
  cp<-500+changepoint
  
  if (eta_dis=="normal") {eta<-rnorm(n+500,0,1)}
  if (eta_dis=="t5") {eta<-sqrt(3/5)*rt(n+500,5)}
  if (eta_dis=="t10") {eta<-sqrt(4/5)*rt(n+500,10)}
  
  
  for (j in 2:cp){
    sigma2[j]<-omega0+alpha0*(abs(x[j-1])-gamma0*x[j-1])^2+beta0*sigma2[j-1]
    x[j]<-sqrt(sigma2[j])*eta[j]
  }
  for (j in (cp+1):(n+500)){
    sigma2[j]<-omega1+alpha1*(abs(x[j-1])-gamma1*x[j-1])^2+beta1*sigma2[j-1]
    x[j]<-sqrt(sigma2[j])*eta[j]
  }
  ran[,1]<-x
  ran[,2]<-sigma2
  return(ran[501:(n+500),])
}

#######################################################################################################################################################


omega0<-0.2;    alpha0<-0.2;    beta0<-0.2;   gamma0<-(-0.2)
omega1<-omega0; alpha1<-alpha0; beta1<-beta0; gamma1<-gamma0



n_seq<-seq(200,3000,200)
Kn<-10
cv1<-2.241
cv2<-2.728
stack_num=1
size_DnkR <- numeric(5*15)
size_DnkS <- numeric(5*15)

for (dodo in (1:5)){
  for (n_num in n_seq){
    n=n_num
    cpn<-0.5;K<-floor(n*Kn);cp<-floor(n*cpn)
    
    
    count0<-0;count_DnkR<-0;count_DnkS<-0
    taun_DnkR<-numeric(RR)
    taun_DnkS<-numeric(RR)
    eta_dis<-"normal"
    
  
    RR<-300
  
    for (R in 1:RR){
      
      repeat {
        agarch_sample<-sample_agarch(omega0,alpha0,beta0,gamma0,omega1,alpha1,beta1,gamma1,n+K,eta_dis,cp)
        y<-agarch_sample[,1]
        sigma2<-agarch_sample[,2]
        
        h_theta<-garchFit(formula =~aparch(1,1),data=y,include.delta=F,include.mean = F,trace = F)@fit$par  
        count0<-count0+1
        if (h_theta[1]>0.00005 & h_theta[2]>0.00005 & h_theta[2]<0.99999 & h_theta[4]>0.00005 & h_theta[4]<0.99999 & h_theta[3]>-0.99999 & h_theta[3]<0.99999) break
      }
      
      w<-h_theta[1];a<-h_theta[2];g<-h_theta[3];b<-h_theta[4]
      h_h<-numeric(n+K);h_h[1]<-1
      
      for (i in 2:(n+K)){
        h_h[i]<-w+a*(abs(y[i-1])-g*y[i-1])^2+b*h_h[i-1]
      }
      
      
      h_eta<-y/sqrt(h_h)
      h_eta2<-h_eta^2
      
      h_nu2<-mean(h_eta[1:n]^4)-1
      h_nu<-sqrt(h_nu2)
      
      DnkR<-numeric(K)
      temp_DnkR<-0
      
      
      for (j in 1:K){
        DnkR[j]<-sqrt(n)/h_nu*abs(mean(h_eta2[1:(n+j)])- mean(h_eta2[1:n]))
        ifelse(DnkR[j]>cv1,break,temp_DnkR<-temp_DnkR+1)
      }
      
      taun_DnkR[R]<-temp_DnkR
      ifelse(temp_DnkR==K,count_DnkR<-count_DnkR+1,NA)
      
      if(R%%10==0){cat("n=",n,"Kn=",Kn,"iterating: ",R,"of",RR,"\n")}
      
      h_eta4<-h_eta^4
      
      row_omega <- 0; row_alpha <- 0 ; row_beta <- 0 ; row_gamma <- 0
      
      
      for (i in (2:(n+K))){
        row_omega[i] <- 1 + b * row_omega[i-1]
        row_alpha[i] <- (abs(y[i-1]) - g * y[i-1])^2 + b * row_alpha[i-1]
        row_beta[i] <- h_h[i-1] + b * row_beta[i-1]
        row_gamma[i] <- 2 * a * (abs(y[i-1] - g * y[i-1]) * (-y[i-1]) + b * row_gamma[i-1])
      }
      
      h_h_dot <- t(rbind(row_omega, row_alpha, row_beta, row_gamma))
      h_h_dot2 <- h_h_dot/h_h
      In <- ((mean(h_eta4[1:n])-1)/4)*((t(h_h_dot2[1:n,1:4]) %*% h_h_dot2[1:n,1:4])/n)
      
      h_lt_dot <- (h_h_dot/h_h)*(1-(y^2 /h_h))/(-2)
      
      In2 <- eigen(In)$vectors%*% diag((eigen(In)$values)^(-0.5)) %*% t(eigen(In)$vectors)
      DnkS<-numeric(K)
      temp_DnkS<-0
      
      
      for (j in (1:K)){
        DnkS[j]<-max(abs(In2 %*% (colMeans(h_lt_dot[1:(n+j),1:4])-colMeans(h_lt_dot[1:n,1:4])))) *sqrt(n)
        ifelse(DnkS[j]>cv2,break,temp_DnkS<-temp_DnkS+1)
      }
      
      taun_DnkS[R]<-temp
      ifelse(temp_DnkS==K,count_DnkS<-count_DnkS+1,NA)
      
      
      if(R%%10==0){cat("n=",n,"Kn=",Kn,"iterating: ",R,"of",RR,"\n")}
           
    
    }
    size_DnkR[stack_num]<-1-count_DnkR/RR
    size_DnkS[stack_num]<-1-count_DnkS/RR
    stack_num <- stack_num+1
  }
}
