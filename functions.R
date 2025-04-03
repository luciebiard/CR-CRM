library(survival)
library(cmprsk)
library(xtable)


######################################################
# getprior_exp
######################################################


getprior_exp <- function (halfwidth, target, nu, nlevel, model = "empiric", intcpt = 3, tstar=NULL) 
{
  dosescaled <- prior <- rep(NA, nlevel)
  b <- rep(NA, nlevel + 1)
  b[1] <- -Inf
  b[(nlevel + 1)] <- Inf
  
  if (model == "empiric") {
    dosescaled[nu] <- target
    for (k in nu:2) {
      b[k] <- log(log(target + halfwidth)/log(dosescaled[k]))
      if (nu > 1) {
        dosescaled[k - 1] <- exp(log(target - halfwidth)/exp(b[k]))
      }
    }
    if (nu < nlevel) {
      for (k in nu:(nlevel - 1)) {
        b[k + 1] <- log(log(target - halfwidth)/log(dosescaled[k]))
        dosescaled[k + 1] <- exp(log(target + halfwidth)/exp(b[k + 1]))
      }
    }
    val <- dosescaled
  }
  
  else if (model == "logistic") {
    dosescaled[nu] <- log(target/(1 - target)) - intcpt
    for (k in nu:2) {
      b[k] <- log((log((target + halfwidth)/(1 - target - 
                                               halfwidth)) - intcpt)/dosescaled[k])
      if (nu > 1) {
        dosescaled[k - 1] <- (log((target - halfwidth)/(1 - 
                                                          target + halfwidth)) - intcpt)/exp(b[k])
      }
    }
    if (nu < nlevel) {
      for (k in nu:(nlevel - 1)) {
        b[k + 1] <- log((log((target - halfwidth)/(1 - 
                                                     target + halfwidth)) - intcpt)/dosescaled[k])
        dosescaled[k + 1] <- (log((target + halfwidth)/(1 - 
                                                          target - halfwidth)) - intcpt)/exp(b[k + 1])
      }
    }
    val <- {
      1 + exp(-intcpt - dosescaled)
    }^{
      -1
    }
  }
  else if (model == "exponential"){
    
    dosescaled[nu] <- log(-log(1-target)/tstar)
    for (k in nu:2) {
      b[k] <- log(log(-log(1-(target + halfwidth))/tstar)/ dosescaled[k])      
      if (nu > 1) {
        dosescaled[k - 1] <- log(-log(1-(target - halfwidth))/tstar)/ exp(b[k])
      }
    }
    if (nu < nlevel) {
      for (k in nu:(nlevel - 1)) {
        b[k + 1] <- log(log(-log(1-(target - halfwidth))/tstar)/ dosescaled[k])  
        dosescaled[k + 1] <- log(-log(1-(target + halfwidth))/tstar)/ exp(b[k+1])
      }
    }
    val <- 1-exp(-exp(dosescaled)*tstar)
    
    
  }
  val
}

######################################################
# simulations exponential hazard rate exponential
######################################################

get_alpha <- function(pi, obswin, typ="exponential", typ_wb=NULL){
  
  if(typ %in% c("exponential", "clayton")){
    alpha1 <- (- log(1-pi)/obswin)
  } 
  
  
  if(typ=="weibull"){
    
    if(typ_wb=="decreasing") {a=0.3}
    if(typ_wb=="increasing") {a=3}
    if(typ_wb=="constant") {a=1}
    b = exp(1/a * log(-obswin^a/log(1-pi)))
    alpha1 <- cbind(rep(a, length(pi)), b)
  }
  
  
  return(alpha1)
}


######################################################
# simulations time Clayton
######################################################
get_clayton <- function(n,K, lambdaT, lambdaP, phi=1.5, u1, u2){


  a <- (u1)^(-1/phi)-1
  b <- (u1)^(-(phi+1)/phi)
  se <- u2

  t1 <- qexp(1-u1, rate=lambdaT)
  t2 <-  phi/lambdaP * log((se/b)^(1/(-phi-1))-a)
  
  return(cbind(u1=u1, u2=u2, tt=pmin(t1, t2), dd=apply(cbind(t1, t2), 1, which.min)))
}

######################################################
# simulations time exponential
######################################################
get_expo <- function(n, K, lambdaT, lambdaP, u1, u2){
  
  tt <- qexp(u1,  lambdaT + lambdaP)
  dd <- qbinom(u2, 1, lambdaP / (lambdaT + lambdaP))+1
  
  
  
  
  return(cbind(u1=u1, u2=u2, tt=tt, dd=dd))
}

######################################################
# simulations time Weibull
######################################################
temps <- function (x, a1, b1, a2, b2, u){
  
  udiff <- (x/b1)^a1 + (x/b2)^a2 + log(u)
  
  return(udiff)
}


hx <- function(a,b,x) {
  
  hx1 = a/b * (x/b)^(a-1)  
  return(hx1)
  
}

one_weibull_comp <- function(a1, b1, a2, b2, u1, u2){
  
  u_efs <- u1  
  t_efs <- uniroot(temps, a1=a1, b1=b1, a2=a2, b2=b2, u=u_efs, interval= c(1.e-14, 1e04),
                   extendInt="yes")$root
  
  hx1 <- hx(a1, b1, t_efs)
  hx2 <- hx(a2, b2, t_efs)
  
  evt <- qbinom(u2, 1, hx2/(hx1+hx2))+1
  
  return(c(u1=u1, u2=u2, tt=t_efs, dd=evt))
}

get_weibull <- function(n, a1, b1, a2, b2, u1, u2){
  
  if(n !=length(a1)){
    print("unequal length")
  }else {
    
    baz <- cbind( a1, b1, a2, b2, u1, u2)
    res <- t(apply(baz, 1, function(x){one_weibull_comp(x[1], x[2], x[3], x[4], x[5], x[6])}))
    colnames(res) <- c("u1", "u2", "tt", "dd")
  } 
  return(res)
}


######################################################
# simulations exponential EFS
######################################################
get_dataset0 <- function(n=60, alpha1, alpha2, beta1=NULL, beta2=NULL, tstar, K=5, graine=1234, type='exponential') {
  
  set.seed(graine)

  if(type=='exponential'){
    alpha1 <- rep(alpha1, each=n)
    alpha2 <- rep(alpha2, each=n)
    
    u1 <- rep(runif(n), K)
    u2 <- rep(runif(n), K)
    
    res <- as.data.frame(get_expo(n, K, alpha1, alpha2, u1, u2))
  }
  
  if(type=='clayton'){
    alpha1 <- rep(alpha1, each=n)
    alpha2 <- rep(alpha2, each=n)
    
    u1 <- rep(runif(n), K)
    u2 <- rep(runif(n), K)

    res <- as.data.frame(get_clayton(n, K, alpha1, alpha2, phi=1.5, u1, u2))
  }  
  
  if(type=="weibull"){
    u1 <- rep(runif(n), K)
    u2 <- rep(runif(n), K)
    
    res <- as.data.frame(get_weibull(n*K, a1=rep(alpha1[,1], each=n), b1=rep(alpha1[,2], each=n), a2=rep(alpha2[,1], each=n), b2=rep(alpha2[,2],each=n), u1, u2))
    
  }
  
  
  dd <- ifelse(res$tt > tstar, 0, res$dd)
  tt <- ifelse(res$tt > tstar, tstar, res$tt)
  
  baz <- data.frame(id=1:(n*K), u1=res$u1, u2=res$u2, tt=ceiling(tt*7)/7, dd=dd, status=dd, tstar=tstar)
  baz$id <- rep(1:n, K)
  baz$dose <- rep(1:K, each=n)
  baz <- baz[order(baz$id), c('id','u1' ,'u2', 'tt',"status",'dose')]
  colnames(baz) <- c('patient','u1' ,'u2', 'time',"status",'dose')
  baz <- as.matrix(baz)
  
  T_entrance_basic = c(0,rep(tstar,n-1))
  T_entrance_basic = cumsum(T_entrance_basic)
  
  T_entrance_fixrapid = c(0, rep(tstar/4 ,n-1))
  T_entrance_fixrapid = cumsum(T_entrance_fixrapid)
  
  T_entrance_fixslow= c(0, rep(tstar/2 ,n-1))
  T_entrance_fixslow = cumsum(T_entrance_fixslow)
  
  
  T_entrance <- data.frame(T_entrance_basic=T_entrance_basic,
                           T_entrance_fixslow=T_entrance_fixslow,
                           T_entrance_fixrapid=T_entrance_fixrapid
                           
  )
  
  return(list(data_complete=baz, T_entrance=T_entrance) )
  
}

################################################################################
# LIKELIHOODS
################################################################################

logvraisemblance_prog <- function(para, dv, deltav, tv, xref){
  
  lv <- sum( I(deltav==2)*1 * ((para[1]) + (para[2])*xref[dv] + (para[3])*xref[dv]^2) - tv * exp((para[1]) + (para[2])*xref[dv]  + (para)[3]*xref[dv]^2))
  return (lv)
}

logvraisemblance_tox_expo <- function(para, deltav,  dv, tv, xref) {
  
  lv <-  sum(I(deltav==1)*1 * exp(para)*xref[dv]  - tv*exp(exp(para)*xref[dv]))
  
  return(lv)
  
}


################################################################################
# NEXT DOSE
################################################################################

dosesuiv <- function(doseref.ds, objds=0.25, dds, deltads, tds, tstar, stage=2, rando="adapt", fourchette=0.1, ni=0){
  
    if (sum(deltads==1)==0) {
    
      print("No DLT yet")
      print(ni)
    
    } else {
      
      nbds <- length(dds)
  
      bn1 <- optim(par=c(0), logvraisemblance_tox_expo, dv=dds, deltav=deltads, tv=tds, xref=doseref.ds , method="L-BFGS-B", control=list(fnscale=-1),lower=-10, upper= 10)$par
      
      l1 <- exp(exp(bn1)* doseref.ds)
      
      cif_dlt_lat <- 1 - exp(- l1 * tstar)
      
        
     if(stage==1){
       if (sum(cif_dlt_lat <= (objds))==0) {

        nledose <- doseref.ds[1]
      } else {
        nledose <- doseref.ds[which.min(abs(cif_dlt_lat - objds))[1]]
      }
       return(list(nextdoselevel= which(doseref.ds==nledose), nextdoseclear=nledose, 
                   lat_dlt=cif_dlt_lat, lat_prog=NA, 
                   bhat= c(bn1, NA, NA, NA)))
       
       }
    
      if (stage==2) {
        
        bn2 <- optim(par=c(0,0,0), logvraisemblance_prog, dv=dds, deltav=deltads, tv=tds, xref=doseref.ds , method="L-BFGS-B", control=list(fnscale=-1),lower=rep(-10,3), upper=rep(10,3))$par
        l2 <- exp(bn2[1] + bn2[2]* doseref.ds + bn2[3]*doseref.ds^2)
        cif_prog_lat <- 1 - exp(- l2 * tstar) 

        
        if (sum(cif_dlt_lat  <= objds)==0) {
          nledose <- doseref.ds[1]
        } else {
        set_dlt <- unique(c(which(cif_dlt_lat <= objds), which.min(abs(cif_dlt_lat - objds))))   
        #set_dlt <- unique(c(which(cif_dlt_lat <= objds)))   

        
        if(length(set_dlt) ==1) {
          nledose <- doseref.ds[1]
        }else {
          if(sum(cif_prog_lat[set_dlt]==1) == length(set_dlt)) {
            ix <- sample(set_dlt, 1)
          }else{
            fi <- which(cif_prog_lat[set_dlt] - min(cif_prog_lat[set_dlt])<= fourchette)
            if(length(fi)==1){
              ix <- set_dlt[fi]
            } else{
              if(rando=="adapt"){
                
                ix <- sample(set_dlt[fi], 1, prob=(1-cif_prog_lat[set_dlt[fi]])/sum(1-cif_prog_lat[set_dlt[fi]]))
              }
              
              if(rando=="range"){
              
                ix <- sample(set_dlt[fi], 1)
              }
            }
          }
          nledose <- doseref.ds[ix]
        }
        
        
      }
        return(list(nextdoselevel= which(doseref.ds==nledose), nextdoseclear=nledose, 
                    lat_dlt=cif_dlt_lat,lat_prog=cif_prog_lat, 
                    bhat= c(bn1, bn2)))
        
      }
      
      
      if (stage==3) {
        bn2 <- optim(par=c(0,0,0), logvraisemblance_prog, dv=dds, deltav=deltads, tv=tds, xref=doseref.ds , method="L-BFGS-B", control=list(fnscale=-1),lower=rep(-10,3), upper=rep(10,3))$par
        l2 <- exp(bn2[1] + bn2[2]* doseref.ds + bn2[3]*doseref.ds^2)
        cif_prog_lat <- 1 - exp(- l2 * tstar) 
 
        if (sum(cif_dlt_lat <= (objds))==0) {
          nledose <- doseref.ds[1]
        } else {
          
          set_dlt <- unique(c(which(cif_dlt_lat <= objds), which.min(abs(cif_dlt_lat - objds))))   
          
          if(length(set_dlt) ==1 || sum(1-cif_prog_lat[set_dlt])==0) {
            nledose <- doseref.ds[1]
          } else {
            
            ix <- which(cif_prog_lat == min(cif_prog_lat[set_dlt])[1])[1]
            
            nledose <- doseref.ds[ix]
            
           
          }
          
        }
        return(list(nextdoselevel= which(doseref.ds==nledose), nextdoseclear=nledose, 
                    lat_dlt=cif_dlt_lat, lat_prog=cif_prog_lat, 
                    bhat= c(bn1, bn2)))
        
        } 
 
    }
}




initiation <- function(baz, nch, n, doseinit, tstar, type, K){
  
  baztmp <- matrix(baz[baz[,'patient'] %in% 1:nch & baz[, 'dose']==doseinit,  c("patient",'entrance', "time", "status", "dose" )], ncol=5)
  colnames(baztmp) <-  c("patient",'entrance', "time", "status", "dose")
  
  baztmp <- cbind(baztmp, "date_evt"= (baztmp[,'time']) + (baztmp[,'entrance']))
  
  date_assign <- max(baztmp[,'entrance']) + tstar
  
  av_fup <- pmin(baztmp[,'date_evt'], date_assign) - baztmp[,'entrance']
  av_evt <- ifelse(baztmp[,'time'] - av_fup > 10^(-9) , 0, baztmp[,'status'])
  drun <- doseinit
  
  if (type=='both'){
    
    while((sum(av_evt==1)==0 | sum(av_evt==2)==0) & (nrow(baztmp) + nch) < n ) {
      
      if (av_evt[nrow(baztmp)]==2) {drun <- min(drun + 1, K)} # increasing dose until Kth
      if (av_evt[nrow(baztmp)]==1) {drun <- max(drun - 1, 1)}  # decrease dose until 1st
      
      
      tmp <- matrix(baz[baz[,'patient'] %in% c(nrow(baztmp)+1:nch)& baz[,'dose']==drun,  c("patient", 'entrance',"time", "status", "dose"  )], ncol=5)
      colnames(tmp) <-  c("patient",'entrance', "time", "status", "dose")
      
      tmp <- cbind(tmp, 'date_evt'= tmp[,'time'] + tmp[,'entrance'])
      tmp[,'entrance'] <-  date_assign + tmp[,'entrance'] - tmp[1,'entrance']
      tmp[,'date_evt'] <- tmp[,'time'] + tmp[,'entrance']
      
      baztmp <- rbind(baztmp, tmp)
      
      date_assign <- max(baztmp[,'entrance']) + tstar
      av_fup <- pmin(baztmp[,'date_evt'], date_assign) - baztmp[,'entrance']
      av_evt <- ifelse(baztmp[,'time'] - av_fup > 10^(-9) , 0, baztmp[,'status'])
      
    }
  }
  
  if (type=='dlt'){
    nb_prog <- 0
    while(sum(av_evt==1)==0 & (nrow(baztmp) + nch) <= n ) {
      
      if (sum(av_evt==2)>nb_prog) {drun <- min(drun + 1, K);nb_prog <- nb_prog+1} # increasing dose until Kth
      
      
      tmp <- matrix(baz[baz[,'patient'] %in% c(nrow(baztmp)+1:nch)& baz[,'dose']==drun,  c("patient", 'entrance',"time", "status", "dose"  )], ncol=5)
      colnames(tmp) <-  c("patient",'entrance', "time", "status", "dose")
      
      tmp <- cbind(tmp, 'date_evt'= tmp[,'time'] + tmp[,'entrance'])
      tmp[,'entrance'] <-  date_assign + tmp[,'entrance'] - tmp[1,'entrance']
      tmp[,'date_evt'] <- tmp[,'time'] + tmp[,'entrance']
      
      baztmp <- rbind(baztmp, tmp)

      date_assign <- max(baztmp[,'entrance']) + tstar 
      av_fup <- pmin(baztmp[,'date_evt'], date_assign) - baztmp[,'entrance']
      av_evt <- ifelse(baztmp[,'time'] - av_fup > 10^(-9) , 0, baztmp[,'status'])
      
    }
  }
  return(list(baztmp, date_assign=date_assign, av_fup=av_fup, av_evt=av_evt, drun=drun))
}

execute_scenario_surv <- function(baz, nch=3, n=50, tstar, xref, target, cif_prior, doseinit=1, K, nstage1=20, delta, alpha=0.05, rando="adapt", ni)  {
  
        res <- initiation(baz, nch, n=n, doseinit, tstar, type='dlt', K=K)  
        baztmp <- res[[1]]
        
        baz[baz[,'patient']>nrow(baztmp),'entrance'] <- res[[2]] + baz[baz[,'patient']>nrow(baztmp),'entrance'] - min(baz[baz[,'patient']>nrow(baztmp),'entrance'] )

    
    n1 <- nrow(baztmp)
    baztmp <- cbind(baztmp, 'bhat1'=NA, 'bhat20'=NA, 'bhat21'=NA, 'bhat22'=NA)
    baztmp <- cbind(baztmp,'newdose'= NA)
    
    
     
    
    if((sum(baztmp[,'status']==1)==0) & n1 >= n) { 
      
      baztmp <- baztmp[ , c("patient",'entrance', "time", "status" , "dose","date_evt" ,"bhat1","bhat20","bhat21","bhat22",'newdose')]
      print("Heterogeneity issue")
      
      vari <- c(selected_dose=NA, 
                n1=n1, 
                duration=max(baztmp[,'date_evt']),
                fup=sum(baztmp[,"time"]), 
                table(factor(baztmp[,'dose'], 1:K)),
                table(factor(baztmp[,'dose'][baztmp[,'status']==1], 1:K)),
                table(factor(baztmp[,'dose'][baztmp[,'status']==2], 1:K)), 
                rep(NA, 5),
                rep(NA, 5), 
                rep(NA, 4)
      )
      
      
    } else{
      
      if((sum(baztmp[,'status']==1)>0)& sum(res$av_evt==1)==0) { 
        
        date_assign2 <- baztmp[,'date_evt'][baztmp[,'status']==1] 
        av_fup <- pmin(baztmp[, 'date_evt'], date_assign2) - baztmp[,'entrance']
        av_evt <- ifelse(baztmp[,'time'] - av_fup > 10^(-9) , 0, baztmp[,'status'])
        baz[,'entrance'][baz[,'patient']%in% (nrow(baztmp)+1):nrow(baz)] <- baz[,'entrance'][baz[,'patient'] %in% (nrow(baztmp)+1):nrow(baz)] + (date_assign2 - res$date_assign)
        
        } else {av_fup <- res$av_fup; av_evt <- res$av_evt}
      
      
      xxnew <- dosesuiv(doseref.ds=xref, objds=target, dds=baztmp[,'dose'], deltads=av_evt, tds=av_fup, tstar=tstar, stage=1, ni=ni)
      
      
      dose_reco <- min(xxnew$nextdoselevel, min(K, res$drun+1))
      
      baztmp <- cbind(baztmp, matrix(c(rep(c(cif_prior, rep(NA, K)), nrow(baztmp)-1), xxnew$lat_dlt, rep(NA, K)), byrow=T, ncol=K*2))
      colnames(baztmp)[colnames(baztmp)==""] <- as.character(1:(2*K))
    
      baztmp[nrow(baztmp), c('bhat1', 'bhat20','bhat21' , 'bhat22')] <- xxnew$bhat
      baztmp[nrow(baztmp),'newdose'] <- dose_reco
      
      while (nrow(baztmp)< nstage1){ # Stage 1

        tmp <-matrix(baz[baz[,'patient']==nrow(baztmp)+1 & baz[,'dose']==dose_reco, c("patient",'entrance', "time", "status", "dose")], ncol=5)
        colnames(tmp) <-  c("patient",'entrance', "time", "status", "dose")
        tmp <- cbind(tmp, 'date_evt'=tmp[,'time'] + tmp[,'entrance'])
        
        baztmp <- rbind(baztmp, NA) 
        baztmp[nrow(baztmp), c( "patient", 'entrance',"time", "status" ,"dose",  "date_evt") ]<- tmp[, c( "patient", 'entrance',"time", "status" ,"dose",  "date_evt") ]
        
        date_assign <- baz[,'entrance'][baz[,'patient']==nrow(baztmp)+1][1] 
        av_fup <- pmin(baztmp[,'date_evt'], date_assign) - baztmp[,'entrance']
        av_evt <- ifelse( baztmp[,'time'] - av_fup > 10^(-9), 0, baztmp[,'status'])
        
        xxnew <- dosesuiv(doseref.ds=xref, objds=target, dds=baztmp[,'dose'], deltads=av_evt, tds=av_fup, tstar=tstar, stage=1, ni)
        dose_reco <- min(xxnew$nextdoselevel, min(K,dose_reco+1))
        
         
        baztmp[ nrow(baztmp), as.character(1:(K*2))] <- c(xxnew$lat_dlt, rep(NA, K))
          

        baztmp[nrow(baztmp), c('bhat1', 'bhat20','bhat21' , 'bhat22')] <- xxnew$bhat
        baztmp[,'newdose'][nrow(baztmp)] <- dose_reco
          
          
        }

      while (nrow(baztmp) < n){ # STAGE 2

        tmp <-matrix(baz[baz[,'patient']==nrow(baztmp)+1 & baz[,'dose']==dose_reco, c("patient",'entrance', "time", "status", "dose")], ncol=5)
        colnames(tmp) <-  c("patient",'entrance', "time", "status", "dose")
        tmp <- cbind(tmp, 'date_evt'=tmp[,'time'] + tmp[,'entrance'])
        
        baztmp <- rbind(baztmp, NA) 
        baztmp[nrow(baztmp), c( "patient", 'entrance',"time", "status" ,"dose",  "date_evt") ]<- tmp[, c( "patient", 'entrance',"time", "status" ,"dose",  "date_evt") ]
        
        date_assign <- baz[,'entrance'][baz[,'patient']==nrow(baztmp)+1][1] 
        

        if (nrow(baztmp)==n) {
          xxnew2 <-  dosesuiv(doseref.ds=xref, objds=target, dds=baztmp[,'dose'], deltads=baztmp[,'status'], tds=baztmp[,'time'], tstar=tstar, stage=3)
          selected_dose <- min(xxnew2$nextdoselevel, min(K, dose_reco+1))
          fup <- sum(baztmp[,'time'])
          duration <- max(baztmp[,'date_evt'])

          baztmp[ nrow(baztmp), as.character(1:(K*2))] <-c(xxnew2$lat_dlt,xxnew2$lat_prog)
          baztmp[nrow(baztmp), c('bhat1', 'bhat20','bhat21' , 'bhat22')] <- xxnew2$bhat
          baztmp[nrow(baztmp),'newdose'] <- selected_dose
          
            
        } else {
          
          av_fup <- pmin(baztmp[,'date_evt'], date_assign) - baztmp[,'entrance']
          av_evt <- ifelse( baztmp[,'time'] - av_fup > 10^(-9), 0, baztmp[,'status'])
          
         
          xxnew <- dosesuiv(doseref.ds=xref, objds=target, dds=baztmp[,'dose'], deltads=av_evt, tds=av_fup, tstar=tstar, stage=2, rando=rando, fourchette =delta)
          dose_reco <- min(xxnew$nextdoselevel, min(K,dose_reco+1))
          
          baztmp[ nrow(baztmp), as.character(1:(K*2))] <-c(xxnew$lat_dlt,  xxnew$lat_prog)
          baztmp[nrow(baztmp), c('bhat1', 'bhat20','bhat21' , 'bhat22')] <- xxnew$bhat
          baztmp[nrow(baztmp),'newdose'] <- dose_reco
          
        }
        
      }   

      
      
     colnames(baztmp) <-  c("patient", 'entrance',"time", "status" ,"dose", "date_evt" , 'bhat1','bhat20','bhat21','bhat22','newdose',
                            paste0('CIF_dlt_',1:K),paste0('CIF_prog_',1:K)) 
     
     vari <- c(selected_dose=selected_dose, 
               #selected_dose2=paste(selected_dose2, collapse=''),
               n1=n1, 
               duration=duration,
               fup=fup, 
               nperlevel=table(factor(baztmp[,'dose'], 1:K)),
               dltperlevel=table(factor(baztmp[,'dose'][baztmp[,'status']==1], 1:K)),
               progperlevel=table(factor(baztmp[,'dose'][baztmp[,'status']==2], 1:K)), 
               dlt_est=round(xxnew2$lat_dlt,4),
               prog_est=round(xxnew2$lat_prog,4), 
               betas=round(xxnew2$bhat,5)
               )
  }
  
    return(list(baztmp=baztmp[, 1:11], res=vari))
    #return(vari)
  
}




get_nsim_surv <- function(Nsim, n, correct=3, n1=20, K, scenar,  tstar, acc_typ='fixslow',graine=32, delta, rando="adapt", typ="exponential", typ_wb=NULL, half){

  p1 <- list_scenarios6$p1o[[1]]
  p1 <- p1[(K+1-correct): (2*K-correct)]
  

  p2 <- list_scenarios6$p2o[[scenar]]
  
  pi <- getprior_exp(halfwidth=half, target=0.25, nu=3, nlevel=K, model="exponential", tstar=tstar)
  
  xref <- log(-log(1-pi)/tstar)
 

  tox_prior <- pi 
  prog_prior <-rep(NA, K)
  efs_prior <- rep(NA, K)

  trials <- endpoints <- g1hat_list <- g2hat_list <- g1hat_list_bis <- g2hat_list_bis <- list()

  

  for (ni in 1:Nsim) {
    
    
    baz0 <- get_dataset0(n=n,  tstar=tstar, alpha1=get_alpha(p1, obswin=tstar, typ=typ, typ_wb=typ_wb[1]), alpha2=get_alpha(p2, obswin=tstar, typ=typ, typ_wb=typ_wb[2]), K=K, graine=graine+ni, type=typ)
    
    baz <- (baz0$data_complete)
    baz <- cbind(baz, "entrance"=rep(baz0$T_entrance[, paste0('T_entrance_', acc_typ)], each=K))
    

    fit <- execute_scenario_surv(baz, nch=1, n=n, tstar=tstar, nstage1 = n1,    target=0.25, xref=xref,  cif_prior= tox_prior,  doseinit=1,  K=K, delta=delta, rando=rando, ni=ni)
    
    trials[[ni]] <- fit$baztmp
    
    colnames(trials[[ni]]) <-  c("patient", 'entrance',"time", "status" ,"dose", "date_evt" , 'bhat1','bhat20','bhat21','bhat22','newdose')
    
    
    endpoints[[ni]] <- as.numeric(fit$res)
    
 # Benchmark   
    if(sum(baz[,'status']==1)>0){
      g1hat_list[[ni]] <- summary(survfit(Surv(time, I(status==1)*1)~factor(dose), data=as.data.frame(baz), etype=I(status==1)*1), times=tstar)$pstate[,1]
      g1hat_list_bis[[ni]] <- timepoints(cuminc(baz[,'time'], baz[,"status"], factor(baz[,'dose'], 1:5), cencode = 3), 8)$est[6:10]
      }
    
    if(sum(baz[,'status']==2)>0){
      g2hat_list[[ni]] <- summary(survfit(Surv(time, I(status==2)*1)~factor(dose), data=as.data.frame(baz), etype=I(status==2)*1), times=tstar)$pstate[,1]
      g2hat_list_bis[[ni]] <- timepoints(cuminc(baz[,'time'], baz[,"status"], factor(baz[,'dose'], 1:5), cencode = 3), 8)$est[11:15]
      }
    
    
  }
  
  return(list(trials, endpoints, g1hat_list,g1hat_list_bis, g2hat_list,g2hat_list_bis))  
  
}




