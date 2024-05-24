
library("invgamma")
library("SciViews")
library('plyr')
library("dplyr")
library('data.table')
library('matlib')
library('MASS')
library(parallel)
options(scipen=999)
source("/nesi/nobackup/uoa03821/REMAP-CAP/SNAP_function.R")


SimulateAMonthOfAccrualTimes <- function( dPatsPerMonth , dStartMonth )
{
  nQtyPats    <- 1.2 *qpois(0.9999,dPatsPerMonth)
  vTimes      <- cumsum( rexp( nQtyPats, dPatsPerMonth ) )
  vTimes      <- vTimes[ vTimes < 1 ]
  vTimes      <- vTimes + dStartMonth
  return( vTimes )
}

SimulateArrivalTimes <- function( vPatsPerMonth, nMaxQtyPats )
{
  vTimes <- c()
  if( length( vPatsPerMonth ) == 1 )
  {
    vTimes <- cumsum(rexp(nMaxQtyPats ,vPatsPerMonth))
  }
  else
  {
    dStartMonth <- 0
    nMonth     <- 1
    while( length( vTimes ) < nMaxQtyPats  )
    {
      vTimes      <- c( vTimes, SimulateAMonthOfAccrualTimes( vPatsPerMonth[ nMonth ], dStartMonth ))
      dStartMonth <- dStartMonth + 1

      if( nMonth < length( vPatsPerMonth ) )
        nMonth <- nMonth +  1
    }
    vTimes <- vTimes[ 1:nMaxQtyPats ]
  }
  return( vTimes )
}

SimulateOutcomeObservedTime <- function( vStartTime )
{
  vTimeToOutcome <-3 #rnorm( length( vStartTime ),2, 0.3)
  vObsTime <- vStartTime  + vTimeToOutcome
  return( vObsTime )
}

blockrand = function(blocksize,N,armn,armlabel){
  block = rep(1:ceiling(N/blocksize), each = blocksize)
  a1 = data.frame(block, rand=runif(length(block)), envelope= 1: length(block))
  a2 = a1[order(a1$block,a1$rand),]
  a2$arm = rep(armlabel,times = length(block)/armn)
  assign = a2[order(a2$envelope),]
  return(assign[,c("block","arm")])
}



pop<-function(vPatsPerMonth,nMaxQtyPats,enrollrate1){
  populationtotal<-SimulateArrivalTimes (vPatsPerMonth, nMaxQtyPats)
  vStartTime1<-rbinom(nMaxQtyPats,size=1,enrollrate1)
  vStartTime2<- cbind(vStartTime1,populationtotal)
  vStartTime3<-vStartTime2[vStartTime2[,1]==1,]
  return(list(populationtotal,length(populationtotal),as.vector(vStartTime3[,2])))
}


parameters = c('aS','bA','bS','bM','aSBS' )

levellab<- c('S1-A1-S1-M1','S1-A1-S1-M2','S1-A1-S2-M1','S1-A1-S2-M2','S1-A1-S3-M1','S1-A1-S3-M2',
             'S1-A2-S1-M1','S1-A2-S1-M2','S1-A2-S2-M1','S1-A2-S2-M2','S1-A2-S3-M1','S1-A2-S3-M2',
             'S1-A3-S1-M1','S1-A3-S1-M2','S1-A3-S2-M1','S1-A3-S2-M2','S1-A3-S3-M1','S1-A3-S3-M2',
             'S1-A4-S1-M1','S1-A4-S1-M2','S1-A4-S2-M1','S1-A4-S2-M2','S1-A4-S3-M1','S1-A4-S3-M2',
             'S1-A5-S1-M1','S1-A5-S1-M2','S1-A5-S2-M1','S1-A5-S2-M2','S1-A5-S3-M1','S1-A5-S3-M2',

             'S2-A1-S1-M1','S2-A1-S1-M2','S2-A1-S2-M1','S2-A1-S2-M2','S2-A1-S3-M1','S2-A1-S3-M2',
             'S2-A2-S1-M1','S2-A2-S1-M2','S2-A2-S2-M1','S2-A2-S2-M2','S2-A2-S3-M1','S2-A2-S3-M2',
             'S2-A3-S1-M1','S2-A3-S1-M2','S2-A3-S2-M1','S2-A3-S2-M2','S2-A3-S3-M1','S2-A3-S3-M2',
             'S2-A4-S1-M1','S2-A4-S1-M2','S2-A4-S2-M1','S2-A4-S2-M2','S2-A4-S3-M1','S2-A4-S3-M2',
             'S2-A5-S1-M1','S2-A5-S1-M2','S2-A5-S2-M1','S2-A5-S2-M2','S2-A5-S3-M1','S2-A5-S3-M2')
levelss<- c('non-shock','shock')


snap_fun <- function(datapop=sim11,jj=1,x,N2=7500,aaa,sss,mmm,intercept1,interaction1  ) { #N1,

  vStartTime<-sort(datapop[[jj]][[x]][[3]][1:N2], decreasing = FALSE)
  vOutcomeTime<-SimulateOutcomeObservedTime(vStartTime)

  data1<-matrix(NA_real_,nrow=N2,ncol=10)
  data1[,1]<-1:N2
  data1[,2]<-vStartTime
  data1[,3]<-vOutcomeTime
  assign1<-blockrand(blocksize=30*2,N2,armn=30*2,armlabel=levellab)
  armleft<-levellab

  decisiona1<-character(5)
  decisiona2<- character(5)
  decisions1<-character(3)
  decisions2<- character(3)
  decisionm1<-character(2)
  decisionm2<- character(2)

  stoppa1<-rep(NA_real_,5)
  stoppa2<-rep(NA_real_,5)
  stopps1<-rep(NA_real_,3)
  stopps2<-rep(NA_real_,3)
  stoppm1<-rep(NA_real_,2)
  stoppm2<-rep(NA_real_,2)

  stopA1<-c(1,2,3,4,5)
  stopA2<-c(1,2,3,4,5)
  stopS1<-c(1,2,3)
  stopS2<-c(1,2,3)
  stopM1<-c(1,2)
  stopM2<-c(1,2)

  for (jjj in 1:15){ 
    print('ha')
    print(jjj)
    if (jjj==1){
      stopnum=500
    }else{
      stopnum=500*jjj
    }


    for (j in 1:N2){
      total<-sum (max(as.numeric(data1[1:j,3]))>=as.numeric(data1[1:j,2]))

      if (total==stopnum ){
        data11<-as.matrix(data1[which(max(as.numeric(data1[1:j,3]))>=as.numeric(data1[1:j,2])),],ncol=10)
        break
      }
    }


    if (jjj==1){
      data11[,4]<-assign1$arm[1:nrow(data11)]
      data11[,5]<-as.numeric(factor(assign1$arm[1:nrow(data11)],levels =levellab))#as.numeric(factor(assign1$arm[1:total],levels =levellab))
      data11[,6]<-as.numeric( substr(data11[,4], start = 2, stop =2))
    }else{

      lengtha<-nrow(data11[!complete.cases(data11),])
      if (lengtha>0){
        data11[!complete.cases(data11),][,6]<-as.numeric(factor(blockrand(blocksize=4,N=lengtha,armn=2,armlabel=c('non-shock','shock'))[1:lengtha,2],levels =levelss))
        non_shock_num<-sum(data11[!complete.cases(data11),][,6]==1)
        shock_num<-sum(data11[!complete.cases(data11),][,6]==2)


        if (length(stopA1)!=1){
          A_block1<-sample(stopA1,size=non_shock_num,prob =A_alo1,replace = TRUE)
        }else if (length(stopA1)==1){
          A_block1<-rep(stopA1,non_shock_num)
        }

        if (length(stopA2)!=1){
          A_block2<-sample(stopA2,size=shock_num,prob =A_alo2,replace = TRUE)
        }else if (length(stopA2)==1){
          A_block2<-rep(stopA2,shock_num)
        }

        if (length(stopS1)!=1){
           S_block1<-sample(stopS1,size=non_shock_num,prob =S_alo1,replace = TRUE)
        }else if (length(stopS1)==1){
          S_block1<-rep(stopS1,non_shock_num)
        }

        if (length(stopS2)!=1){
          S_block2<-sample(stopS2,size=shock_num,prob =S_alo2,replace = TRUE)
        }else if (length(stopS2)==1){
          S_block2<-rep(stopS2,shock_num)
        }



        if (length(stopM1)!=1){
          M_block1<-sample(stopM1,size=non_shock_num,prob =M_alo1,replace = TRUE)
        }else if (length(stopM1)==1){
          M_block1<-rep(stopM1,non_shock_num)
        }

        if (length(stopM2)!=1){
            M_block2<- sample(stopM2,size=shock_num,prob =M_alo2,replace = TRUE)
        }else if (length(stopM2)==1){
            M_block2<-rep(stopM2,shock_num)
        }

      non_shock_dat<-data11[!complete.cases(data11),][data11[!complete.cases(data11),][,6]==1,]
        shock_dat<-data11[!complete.cases(data11),][data11[!complete.cases(data11),][,6]==2,]

        non_shock_dat[,7]<- A_block1
        shock_dat[,7]<- A_block2
        non_shock_dat[,8]<- S_block1
        shock_dat[,8]<- S_block2
        non_shock_dat[,9]<- M_block1
        shock_dat[,9]<- M_block2

        data11a<-rbind(non_shock_dat,shock_dat)
        data11a[,4]<-paste('S',data11a[,6],'-A',data11a[,7],'-S',data11a[,8],'-M',data11a[,9])
        data11a[,4] <- gsub(" ", "", data11a[,4])
        data11a[,5]<-as.numeric(factor(data11a[,4],levels =levellab))


      }
    }



    if (jjj==1){
      data11[,7]<- as.numeric( substr(data11[,4], start = 5, stop =5))
      data11[,8]<- as.numeric( substr(data11[,4], start = 8, stop =8))
      data11[,9]<- as.numeric( substr(data11[,4], start = 11, stop =11))


      for (i in 1:total){
        if (as.numeric(data11[i,6])==2 & as.numeric(data11[i,8])==2 ){
          p<-exp(intercept1[as.numeric(data11[i,6])]+aaa[as.numeric(data11[i,7])]+sss[as.numeric(data11[i,8])]+mmm[as.numeric(data11[i,9])]+interaction1[1])/
            (1+exp(intercept1[as.numeric(data11[i,6])]+aaa[as.numeric(data11[i,7])]+sss[as.numeric(data11[i,8])]+mmm[as.numeric(data11[i,9])]+interaction1[1]))
          data11[i,10]=rbinom(1,1,p)
        }else if (as.numeric(data11[i,6])==2 & as.numeric(data11[i,8])==3 ){
          p<-exp(intercept1[as.numeric(data11[i,6])]+aaa[as.numeric(data11[i,7])]+sss[as.numeric(data11[i,8])]+mmm[as.numeric(data11[i,9])]+interaction1[2])/
            (1+exp(intercept1[as.numeric(data11[i,6])]+aaa[as.numeric(data11[i,7])]+sss[as.numeric(data11[i,8])]+mmm[as.numeric(data11[i,9])]+interaction1[2]))
          data11[i,10]=rbinom(1,1,p)
        }  else {
          p<-exp(intercept1[as.numeric(data11[i,6])]+aaa[as.numeric(data11[i,7])]+sss[as.numeric(data11[i,8])]+mmm[as.numeric(data11[i,9])])/
            (1+exp(intercept1[as.numeric(data11[i,6])]+aaa[as.numeric(data11[i,7])]+sss[as.numeric(data11[i,8])]+mmm[as.numeric(data11[i,9])]))
          data11[i,10]=rbinom(1,1,p)
        }
      }


      data11[,1]<-as.numeric(data11[,1])

      data1<- merge.data.frame(x=data1, y=data11,    # Merge data
                               by = c('V1'),
                               all=TRUE)

      data1<-subset(data1, select = -c(V2.y, V3.y,V4.x, V5.x, V6.x, V7.x, V8.x, V9.x ,V10.x))
      names(data1) <- gsub(".y", "", names(data1), fixed = TRUE)
      names(data1) <- gsub(".x", "", names(data1), fixed = TRUE)

    }else if (jjj!=1 & lengtha>0){
    
        if (nrow(data11a)>=1){
 
        for (i in 1:nrow(data11a)){
          if (as.numeric(data11a[i,6])==2 & as.numeric(data11a[i,8])==2 ){
            p<-exp(intercept1[as.numeric(data11a[i,6])]+aaa[as.numeric(data11a[i,7])]+sss[as.numeric(data11a[i,8])]+mmm[as.numeric(data11a[i,9])]+interaction1[1])/
              (1+exp(intercept1[as.numeric(data11a[i,6])]+aaa[as.numeric(data11a[i,7])]+sss[as.numeric(data11a[i,8])]+mmm[as.numeric(data11a[i,9])]+interaction1[1]))
            data11a[i,10]=rbinom(1,1,p)
          }else if (as.numeric(data11a[i,6])==2 & as.numeric(data11a[i,8])==3 ){
            p<-exp(intercept1[as.numeric(data11a[i,6])]+aaa[as.numeric(data11a[i,7])]+sss[as.numeric(data11a[i,8])]+mmm[as.numeric(data11a[i,9])]+interaction1[2])/
              (1+exp(intercept1[as.numeric(data11a[i,6])]+aaa[as.numeric(data11a[i,7])]+sss[as.numeric(data11a[i,8])]+mmm[as.numeric(data11a[i,9])]+interaction1[2]))
            data11a[i,10]=rbinom(1,1,p)
          }  else {
            p<-exp(intercept1[as.numeric(data11a[i,6])]+aaa[as.numeric(data11a[i,7])]+sss[as.numeric(data11a[i,8])]+mmm[as.numeric(data11a[i,9])])/
              (1+exp(intercept1[as.numeric(data11a[i,6])]+aaa[as.numeric(data11a[i,7])]+sss[as.numeric(data11a[i,8])]+mmm[as.numeric(data11a[i,9])]))
            data11a[i,10]=rbinom(1,1,p)
          }
        }


        data1$V1<-as.numeric(data1$V1)
        data11a<-as.data.frame(data11a)
        data11a$V1<-as.numeric(data11a$V1)
        data1<- merge.data.frame(x=data1, y=data11a,    # Merge data
                                 by = c('V1'),
                                 all=TRUE)
        data1$V4.x<-ifelse (is.na(data1$V4.x),data1$V4.y,data1$V4.x)
        data1$V5.x<-ifelse (is.na(data1$V5.x),data1$V5.y,data1$V5.x)
        data1$V6.x<-ifelse (is.na(data1$V6.x),data1$V6.y,data1$V6.x)
        data1$V7.x<-ifelse (is.na(data1$V7.x),data1$V7.y,data1$V7.x)
        data1$V8.x<-ifelse (is.na(data1$V8.x),data1$V8.y,data1$V8.x)
        data1$V9.x<-ifelse (is.na(data1$V9.x),data1$V9.y,data1$V9.x)
        data1$V10.x<-ifelse (is.na(data1$V10.x),data1$V10.y,data1$V10.x)

        data1<-subset(data1, select = -c(V2.y, V3.y,V4.y, V5.y, V6.y, V7.y, V8.y, V9.y ,V10.y))
        names(data1) <- gsub(".y", "", names(data1), fixed = TRUE)
        names(data1) <- gsub(".x", "", names(data1), fixed = TRUE)
      }
    }


    data1 <- data1[order(data1$V1),]
    data11a<-data1[complete.cases(data1),  ]

    A_summary<-rename(count(data11a, V6, V7), Freqt = n, shock=V6, A=V7)
    S_summary<-rename(count(data11a, V6, V8), Freqt = n, shock=V6, S=V8)
    M_summary<-rename(count(data11a, V6, V9), Freqt= n, shock=V6, M=V9)

    data11b<-as.data.frame(data11a)
    data11b$V6<-factor(as.numeric(data11b$V6))
    data11b$V7<-factor(as.numeric(data11b$V7))
    data11b$V8<-factor(as.numeric(data11b$V8))
    data11b$V9<-factor(as.numeric(data11b$V9))
    data11b$V10<-as.numeric(data11b$V10)
    data11b$V6 <- relevel(factor(data11b$V6), ref = 1)
    data11b$V7 <- relevel(factor(data11b$V7), ref = 5)
    data11b$V8 <- relevel(factor(data11b$V8), ref = 1)
    data11b$V9 <- relevel(factor(data11b$V9), ref = 1)

    #likelihood
    mod <- glm(V10 ~ 0+V6 +V7+V8+V9+V6*V8, family = "binomial", data =data11b)
    betahat<-summary(mod)$coefficients[,1]
    V<-vcov(mod)


    #prior
    meanprior<-matrix(c(0,0,0,0,0,0,0,0,0,0,0),ncol=1)
    tau2<-diag(x = 0, nrow=11, ncol=11, names = TRUE)
    diag(tau2)<-c(100,100,100,100,100,100,100,100,100,0.15^2,0.15^2)

    mupo<-inv(inv(tau2) + inv(V)) %*% (inv(tau2) %*% meanprior +  inv(V)%*% betahat)

    rownames(mupo)<-c('aS[1]','aS[2]','bA[1]','bA[2]','bA[3]','bA[4]','bS[2]','bS[3]',
                      'bM[2]', 'aSBS[2,2]','aSBS[2,3]')
    sigma2po<-inv(inv(tau2) + inv(V))


    mcmc_samples1 <- as.data.frame(mvrnorm(n = 20000, mu=mupo, Sigma=sigma2po))
    mcmc_samples_var<- mcmc_samples1 %>% mutate(shock_s2=V7+V10, shock_s3=V8+V11) %>%summarise_if(is.numeric, var)%>%
      as.data.frame()

    colnames(mcmc_samples1) <- c('aS[1]','aS[2]','bA[1]','bA[2]','bA[3]','bA[4]','bS[2]','bS[3]',
                                 'bM[2]', 'aSBS[2,2]','aSBS[2,3]')

    colnames(mcmc_samples_var) <- c('aS[1]','aS[2]','bA[1]','bA[2]','bA[3]','bA[4]','bS[2]','bS[3]',
                                 'bM[2]', 'aSBS[2,2]','aSBS[2,3]','shock_s2','shock_s3')


    posterior1<-list()
    posterior2<-list()
    for (i in 1:20000){

      allocation1<-matrix(NA,nrow=4,ncol=2)
      allocation1[,1]<-c('A','S-shock','S-non-shock','M')

      stopA1a<-stopA1[stopA1!=5]
      stopS1a<-stopS1[stopS1!=1]#non-shock
      stopS2a<-stopS2[stopS2!=1]#shock
      stopM1a<-stopM1[stopM1!=1]

      stopA11<- sort(stopA1a+2)
      stopS11<- sort(stopS1a+5)
      stopS22<- sort(stopS2a-1)
      stopM11<- sort(stopM1a)

      if (length(stopA11)!=0 ){
        allocation1[1,2]<-as.numeric(substr(names(which.min(mcmc_samples1[i,stopA11,drop=F])),4,4))
      }else{
        allocation1[1,2]<-0
      }

      if (length(stopS22)==2){#shock
        allocation1[2,2]<-which.min(c(mcmc_samples1[i,7]+mcmc_samples1[i,10],mcmc_samples1[i,8]+mcmc_samples1[i,11])[stopS22])+1
      }else if (length(stopS22)==1){
        allocation1[2,2]<-stopS22+1
      }else{
        allocation1[2,2]<-0
      }

      if (length(stopS11)!=0){#non-shock
        allocation1[3,2]<-as.numeric(substr(names(which.min(mcmc_samples1[i,stopS11,drop=F])),4,4))
      }else{
        allocation1[3,2]<-0
      }

      if (length(stopM11)!=0){
        allocation1[4,2]<-as.numeric(substr(names(which.min(mcmc_samples1[i,stopM11,drop=F])),4,4))
      }else{
        allocation1[4,2]<-0
      }

      posterior1[[i]]<-matrix(allocation1[c(1,2,4),2],nrow=1,ncol=3)#shock
      posterior2[[i]]<-matrix(allocation1[c(1,3,4),2],nrow=1,ncol=3)#non-shock

    }


    rar1<-do.call(rbind,posterior1)#shock
    rar2<-do.call(rbind,posterior2)#non-shock

    dummyA<-data.frame(cbind(A=c(1,2,3,4,5)))
    dummyS<-data.frame(cbind(S=c(1,2,3)))
    dummyM<-data.frame(cbind(M=c(1,2)))

    A_count1<-rename(data.frame(table(rar1[,1])/20000),A=Var1)
    A_count1$A<-as.numeric(as.character(A_count1$A))
    A_count2<-rename(data.frame(table(rar2[,1])/20000),A=Var1)
    A_count2$A<-as.numeric(as.character(A_count2$A))

    S_count1<-rename(data.frame(table(rar1[,2])/20000),S=Var1)
    S_count1$S<-as.numeric(as.character(S_count1$S))
    S_count2<-rename(data.frame(table(rar2[,2])/20000),S=Var1)
    S_count2$S<-as.numeric(as.character(S_count2$S))

    M_count1<-rename(data.frame(table(rar1[,3])/20000),M=Var1)
    M_count1$M<-as.numeric(as.character(M_count1$M))
    M_count2<-rename(data.frame(table(rar2[,3])/20000),M=Var1)
    M_count2$M<-as.numeric(as.character(M_count2$M))


    a_posterior<-rbind(cbind(left_join(dummyA,A_count1,by='A'),shock=2,vars=mcmc_samples_var[,3:6]),#shock  ,useNA = "always"
                       cbind(left_join(dummyA,A_count2,by='A'),shock=1,vars=mcmc_samples_var[,3:6]))#non-shock

    s_posterior<-rbind(cbind(left_join(dummyS,S_count1,by='S'),shock=2,vars.S2=mcmc_samples_var[,12],vars.S3=mcmc_samples_var[,13]),
                       cbind(left_join(dummyS,S_count2,by='S'),shock=1,vars.S2=mcmc_samples_var[,7],vars.S3=mcmc_samples_var[,8]))

    m_posterior<-rbind(cbind(left_join(dummyM,M_count1,by='M'),shock=2,vars=mcmc_samples_var[,9]),
                       cbind(left_join(dummyM,M_count2,by='M'),shock=1,vars=mcmc_samples_var[,9]))
    a_posterior[is.na(a_posterior)] <- 0
    s_posterior[is.na(s_posterior)] <- 0
    m_posterior[is.na(m_posterior)] <- 0

    A_summary<-rename(count(data11a, V6, V7), Freqt = n, shock=V6, A=V7)
    S_summary<-rename(count(data11a, V6, V8), Freqt = n, shock=V6, S=V8)
    M_summary<-rename(count(data11a, V6, V9), Freqt= n, shock=V6, M=V9)

    A_alo <- merge(a_posterior, A_summary, by=c("shock","A"))
    S_alo <- merge(s_posterior, S_summary, by=c("shock","S"))
    M_alo <- merge(m_posterior, M_summary, by=c("shock","M"))



    for (i in 1:nrow(A_alo)){
      if (A_alo$A[i]==1){
        A_alo$alo[i]=sqrt(A_alo$`vars.bA[1]`[i]*A_alo$Freq[i]/A_alo$Freqt[i])
      }else if (A_alo$A[i]==2){
        A_alo$alo[i]=sqrt(A_alo$`vars.bA[2]`[i]*A_alo$Freq[i]/A_alo$Freqt[i])
      }else if (A_alo$A[i]==3){
        A_alo$alo[i]=sqrt(A_alo$`vars.bA[3]`[i]*A_alo$Freq[i]/A_alo$Freqt[i])
      }else if (A_alo$A[i]==4){
        A_alo$alo[i]=sqrt(A_alo$`vars.bA[4]`[i]*A_alo$Freq[i]/A_alo$Freqt[i])
      }else{
        A_alo$alo[i]=NA
      }

    }

    for (i in 1:nrow(S_alo)){
      if (S_alo$S[i]==2){
        S_alo$alo[i]=sqrt(S_alo$`vars.S2`[i]*S_alo$Freq[i]/S_alo$Freqt[i])
      }else if (S_alo$S[i]==3){
        S_alo$alo[i]=sqrt(S_alo$`vars.S3`[i]*S_alo$Freq[i]/S_alo$Freqt[i])
      }else{
        S_alo$alo[i]=NA
      }
    }

    for (i in 1:nrow(M_alo)){
      if (M_alo$M[i]==2){
        M_alo$alo[i]=sqrt(M_alo$vars[i]*M_alo$Freq[i]/M_alo$Freqt[i])
      }else{
        M_alo$alo[i]=NA
      }
    }
    A_alo$alo1<-A_alo$alo
    S_alo$alo1<-S_alo$alo
    M_alo$alo1<-M_alo$alo

    if (!(any(decisiona1=='S')) & !(sum(decisiona1=='I')==4)){
      afun1<-superfun(dat=mcmc_samples1,category1=c('A'),decisionf1=decisiona1,decisionf2=decisiona2,type='superiority',stopnum=stopnum,stopp1=stoppa1,stopp2=stoppa2,jjj=jjj,
                      stopA1=stopA1,stopA2=stopA2,stopS1=stopS1,stopS2=stopS2,stopM1=stopM1,stopM2=stopM2)
      stopA1<-afun1[[1]]
      decisiona1<-afun1[[3]]
      stoppa1<-afun1[[5]]
    }
    if (!(any(decisiona2=='S')) & !(sum(decisiona2=='I')==4)){
      afun2<-afun1
      stopA2<-afun2[[2]]
      decisiona2<-afun2[[4]]
      stoppa2<-afun2[[6]]
    }


    if (!(any(decisions1=='S')) & !(sum(decisions1=='I')==2)){
      sfun1<-superfun(dat=mcmc_samples1,category1=c('S'),decisionf1=decisions1,decisionf2=decisions2,type='superiority',stopnum=stopnum,stopp1=stopps1,stopp2=stopps2,jjj=jjj,
                      stopA1=stopA1,stopA2=stopA2,stopS1=stopS1,stopS2=stopS2,stopM1=stopM1,stopM2=stopM2)
      stopS1<-sfun1[[1]]
      decisions1<-sfun1[[3]]
      stopps1<-sfun1[[5]]
    }
    if (!(any(decisions2=='S')) & !(sum(decisions2=='I')==2)){
      sfun2<- sfun1
      stopS2<-sfun2[[2]]
      decisions2<-sfun2[[4]]
      stopps2<-sfun2[[6]]
    }


    if (!(any(decisionm1=='S')) & !(sum(decisionm1=='I')==1)){
      mfun1<-superfun(dat=mcmc_samples1,category1=c('M'),decisionf1=decisionm1,decisionf2=decisionm2,type='non-inferiority',stopnum=stopnum,stopp1=stoppm1,stopp2=stoppm2,jjj=jjj,
                      stopA1=stopA1,stopA2=stopA2,stopS1=stopS1,stopS2=stopS2,stopM1=stopM1,stopM2=stopM2)
      stopM1<-mfun1[[1]]
      decisionm1<-mfun1[[3]]
      stoppm1<-mfun1[[5]]
    }

    if (!(any(decisionm2=='S')) & !(sum(decisionm2=='I')==1)){
      mfun2<-mfun1
      stopM2<-mfun2[[2]]
      decisionm2<-mfun2[[4]]
      stoppm2<-mfun2[[6]]
    }



    if (any(decisiona1=='S')){
      A_alo1<-c(rep(0.75/length(which(decisiona1=='S')),length(which(decisiona1=='S'))),0.25)
      stopA1<-c(which(decisiona1=='S'),5)
    }else if (sum(decisiona1=='I')==4){
      A_alo1<-c(0.25/4,0.25/4,0.25/4,0.25/4,0.75)
      stopA1<-c(1,2,3,4,5)
    }else{
        A_alo1<-A_alo[A_alo$shock==1 & A_alo$A %in% stopA1,'alo1']
        A_alo1<-c((A_alo1/sum(A_alo1))*(1-1/(length(A_alo1)+1)),1/(length(A_alo1)+1))
        stopA1<-c(stopA1,5)
     }


    if (any(decisiona2=='S')){
      A_alo2<-c(rep(0.75/length(which(decisiona2=='S')),length(which(decisiona2=='S'))),0.25)
      stopA2<-c(which(decisiona2=='S'),5)
    }else if (sum(decisiona2=='I')==4){
      A_alo2<-c(0.25/4,0.25/4,0.25/4,0.25/4,0.75)
      stopA2<-c(1,2,3,4,5)
    }else{
      A_alo2<-A_alo[A_alo$shock==2 & A_alo$A %in% stopA2,'alo1']
      A_alo2<-c((A_alo2/sum(A_alo2))*(1-1/(length(A_alo2)+1)),1/(length(A_alo2)+1))
      stopA2<-c(stopA2,5)
    }



    if (any(decisions1=='S')){
      S_alo1<-c(rep(0.75/length(which(decisions1=='S')),length(which(decisions1=='S'))),0.25)
      stopS1<-c(which(decisions1=='S'),1)
    }else if (sum(decisions1=='I')==2){
      S_alo1<-c(0.25/2,0.25/2,0.75)
      stopS1<-c(2,3,1)
    }else{
      S_alo1<-S_alo[S_alo$shock==1 & S_alo$S %in% stopS1,'alo1']
      S_alo1<-c((S_alo1/sum(S_alo1))*(1-1/(length(S_alo1)+1)),1/(length(S_alo1)+1))
      stopS1<-c(stopS1,1)
    }


    if (any(decisions2=='S')){
      S_alo2<-c(rep(0.75/length(which(decisions2=='S')),length(which(decisions2=='S'))),0.25)
      stopS2<-c(which(decisions2=='S'),1)
    }else if (sum(decisions2=='I')==2){
      S_alo2<-c(0.25/2,0.25/2,0.75)
      stopS2<-c(2,3,1)
    }else{
      S_alo2<-S_alo[S_alo$shock==2 & S_alo$S %in% stopS2,'alo1']
      S_alo2<-c((S_alo2/sum(S_alo2))*(1-1/(length(S_alo2)+1)),1/(length(S_alo2)+1))
      stopS2<-c(stopS2,1)
    }

    if (any(decisionm1=='S')){
      M_alo1<-c(rep(0.75/length(which(decisionm1=='S')),length(which(decisionm1=='S'))),0.25)
      stopM1<-c(which(decisionm1=='S'),1)
    }else if (sum(decisionm1=='I')==1){
       M_alo1<-c(0.75,0.25)
       stopM1<-c(1,2)
    }else{
      M_alo1<-M_alo[M_alo$shock==1 & M_alo$M %in% stopM1,'alo1']
      M_alo1<-c((M_alo1/sum(M_alo1))*(1-1/(length(M_alo1)+1)),1/(length(M_alo1)+1))
      stopM1<-c(stopM1,1)
    }

    if (any(decisionm2=='S')){
      M_alo2<-c(rep(0.75/length(which(decisionm2=='S')),length(which(decisionm2=='S'))),0.25)
      stopM2<-c(which(decisionm2=='S'),1)
    }else if (sum(decisionm2=='I')==1){
      M_alo2<-c(0.75,0.25)
      stopM2<-c(1,2)
    }else{
      M_alo2<-M_alo[M_alo$shock==2 & M_alo$M %in% stopM2,'alo1']
      M_alo2<-c((M_alo2/sum(M_alo2))*(1-1/(length(M_alo2)+1)),1/(length(M_alo2)+1))
      stopM2<-c(stopM2,1)
    }
    
    
    if ((sum(decisiona1=='I')==4 & sum(decisiona2=='I')==4)  &
         (sum(decisions1=='I')==2 & sum(decisions2=='I')==2)  &
         (sum(decisionm1=='I')==1 & sum(decisionm2=='I')==1)){
           break
    }

  }

  data1<-data1[complete.cases(data1),]
  decisiona1[1:4]<-ifelse(decisiona1[1:4]=='','ND',decisiona1[1:4])
  decisiona2[1:4]<-ifelse(decisiona2[1:4]=='','ND',decisiona2[1:4])
  decisions1[2:3]<-ifelse(decisions1[2:3]=='','ND',decisions1[2:3])
  decisions2[2:3]<-ifelse(decisions2[2:3]=='','ND',decisions2[2:3])
  decisionm1[2]<-ifelse(decisionm1[2]=='','ND',decisionm1[2])
  decisionm2[2]<-ifelse(decisionm2[2]=='','ND',decisionm2[2])

  n<-max(length(decisiona1),length(decisiona2),length(decisions1),length(decisions2),
         length(decisionm1),length(decisionm2))
  length(decisiona1)<-n
  length(decisiona2)<-n
  length(decisions1)<-n
  length(decisions2)<-n
  length(decisionm1)<-n
  length(decisionm2)<-n
  decision<-rbind(decisiona1,decisiona2,decisions1,decisions2,decisionm1,decisionm2)
  data1<-data1[complete.cases(data1),]
  Aalof<-rename(count(data1, V6, V7), Freqt = n)#data.frame(table(data1[,7])/nrow(data1))
  Salof<-rename(count(data1, V6, V8), Freqt = n)#data.frame(table(data1[,8])/nrow(data1))
  Malof<-rename(count(data1, V6, V9), Freqt = n)#data.frame(table(data1[,9])/nrow(data1))
  Aalof$alo<-Aalof$Freqt/nrow(data1)
  Salof$alo<-Salof$Freqt/nrow(data1)
  Malof$alo<-Malof$Freqt/nrow(data1)


  stoppa1<-ifelse(is.na(stoppa1),min(7500,nrow(data1)),stoppa1)
  stoppa2<-ifelse(is.na(stoppa2),min(7500,nrow(data1)),stoppa2)
  stopps1<-ifelse(is.na(stopps1),min(7500,nrow(data1)),stopps1)
  stopps2<-ifelse(is.na(stopps2),min(7500,nrow(data1)),stopps2)
  stoppm1<-ifelse(is.na(stoppm1),min(7500,nrow(data1)),stoppm1)
  stoppm2<-ifelse(is.na(stoppm2),min(7500,nrow(data1)),stoppm2)

  nn<-max(length(stoppa1),length(stoppa2),length(stopps1),length(stopps2),
         length(stoppm1),length(stoppm2))
  length(stoppa1)<-nn
  length(stoppa1)<-nn
  length(stopps1)<-nn
  length(stopps2)<-nn
  length(stoppm1)<-nn
  length(stoppm2)<-nn
  stopp<-rbind(stoppa1,stoppa2,stopps1,stopps2,stoppm1,stoppm2)

  return(list(data1[,4],data1[,10],decision,N2,data1,stopp, Aalof,Salof,Malof))

}

ptsinpop1<-function(x,coutcome,datain1,domainn,linen,perct,position,datapop=sim11){

  ps<- datain1[[x]][[6]][linen,1:(domainn-1)]
  stopseq<-rank(ps[1:(domainn-1)],ties.method = 'min')
  decision<-datain1[[x]][[3]][linen,1:(domainn-1)]
  population<-datapop[[jj]][[x]][[1]]
  enrollinpopulation<-as.numeric(datain1[[x]][[5]][,2])
  enrollinpopulation1<-as.numeric(datain1[[x]][[5]][,3])
  trt<-matrix(1:(domainn-1))
  stopdat<-as.data.frame(cbind(trt,ps,stopseq,decision))
  colnames(stopdat)<-c("realarm","number","stopseq","decision")
  stopdat<-stopdat[with(stopdat, order(stopseq, decision)),]
  stoptrial<-as.numeric(unique(stopdat$number))
  stopindex<-vector('list',length(stoptrial))

  if (length( which (decision %in%  'I'))==(domainn-1)){
    stopindex<- lapply(1:length(stoptrial),function(x){
      for (v in 1:(length(population))) {
        if (all.equal(population[v],max(c(enrollinpopulation[stoptrial[x]],enrollinpopulation1[1:stoptrial[[x]]]),na.rm=TRUE))==TRUE){
          return(v)
          break
        }else if (population[v]>max(c(enrollinpopulation[stoptrial[[x]]],enrollinpopulation1[1:stoptrial[[x]]]),na.rm=TRUE)){
          return(v-1)
          break
        }
      }
    })
  }else{

    stopindex11<- function(x){
      for (v in 1:(length(population))) {
        if (all.equal(population[v],max(c(enrollinpopulation[stoptrial[[x]]],enrollinpopulation1[1:stoptrial[[x]]]),na.rm=TRUE))==TRUE){
          return(v)
          break
        } else if (population[v]>max(c(enrollinpopulation[stoptrial[[x]]],enrollinpopulation1[1:stoptrial[[x]]]),na.rm=TRUE)){
          return(v-1)
          break
        }
      }
    }
    stopindex12<-  function(x){
      for (v in 1:(length(population))) {
        if (population[v]>enrollinpopulation1[stoptrial[x]]){
          return(v-1)
          break
        }
      }
    }
    for(l in 1:length(stoptrial)){
      if (l<length(stoptrial)) {
        stopindex[[l]]<-stopindex11(l)
      }else{
        stopindex[[l]]<-stopindex12(l)
      }
    }
  }

  stoppop<- do.call(rbind,lapply(1:length(stoptrial),
                                 function(x){(stopindex[[x]]-stoptrial[x])/2}))

  if (stopindex[[length(stopindex)]]!=length(population)){
    stoppop<-c(stoppop,(length(population)-stopindex[[length(stopindex)]])/2)
  }

  percttot2<-rep(0,length(unique(stopdat$stopseq)))

  if (length(unique(stopdat$stopseq))>1){
    for (i in 1:(length(unique(stopdat$stopseq))-1)) {
      dat1<-subset(stopdat,stopdat$stopseq<=unique(stopdat$stopseq)[i])
      dum1<-  as.numeric(subset(dat1,dat1$decision== "I")$realarm)#+1
      dum2<- subset(stopdat,stopdat$stopseq>unique(stopdat$stopseq)[i])$realarm
      percttot2[which(unique(stopdat$stopseq)==unique(stopdat$stopseq)[as.numeric(i)+1])]= (sum(perct[dum1])/(length(dum2)+1))*
        (stoppop[which(unique(stopdat$stopseq)==unique(stopdat$stopseq)[as.numeric(i)+1])]-
           stoppop[which(unique(stopdat$stopseq)==unique(stopdat$stopseq)[as.numeric(i)])])
    }
  }
  percttot21<-cumsum(percttot2)
  patientintrialt1<-vector('list',domainn)
  patientinpopulation<-vector('list',domainn)

  for (j in 1:domainn){
    patientintrialt1[[j]]<- coutcome[j]#sum(substr(coutcome,position,position)==j)#length(coutcome[coutcome==j])
  }
  patientintrialt<-matrix(do.call(rbind,patientintrialt1),nrow=1,ncol=domainn)


  if (sum(decision %in% "S")>0) {
    for (j in 1:(domainn-1)){
      seqq<-subset(stopdat,stopdat$realarm==(j))$stopseq
      seqq1<-which(unique(stopdat$stopseq) %in% seqq)
      if(decision[j] %in% "S"){
        patientinpopulation[[j]]<-ceiling(perct[j]*stoppop[seqq1]+
                                            percttot21[seqq1]+stoppop[length(stoppop)]/(sum(decision %in% "S")))
      }else if (decision[j] %in% "ND" ){
        patientinpopulation[[j]]<-ceiling(perct[j]*stoppop[seqq1]+
                                            percttot21[seqq1])
      }else if (decision[j] %in% "I"){
        patientinpopulation[[j]]<-ceiling(perct[j]*stoppop[seqq1]+
                                            percttot21[seqq1])

      }
    }
  }else {
    for (j in 1:(domainn-1)){
      seqq<-subset(stopdat,stopdat$realarm==(j))$stopseq
      seqq1<-which(unique(stopdat$stopseq) %in% seqq)
      patientinpopulation[[j]]<-ceiling(perct[j]*stoppop[seqq1]+
                                          percttot21[seqq1])
    }
  }


  if(sum(decision %in% "S")>0){
    patientinpopulation[[domainn]]<-ceiling(perct[domainn]*max(stoppop[1:(length(stoppop)-1)])+max(percttot21))
  } else {
    patientinpopulation[[domainn]]<-ceiling(perct[domainn]*max(stoppop[1:(length(stoppop)-1)])+max(percttot21)+(stoppop[length(stoppop)]-stoppop[length(stoppop)-1]))
  }

  patientinpopulationt1<-do.call(rbind, patientinpopulation)
  patientinpopulationt<-patientinpopulationt1/max(2*sum(stoppop[(length(stoppop)-1):(length(stoppop))]),2*sum( patientinpopulationt1))
  patientinpopulationt<-matrix(patientinpopulationt,nrow=1,ncol=domainn)
  return(list(patientintrialt,patientinpopulationt))
}


ptsinpop2<-function(x,coutcome,datain1,domainn,linen,perct,position,datapop=sim11){
  ps<- datain1[[x]][[6]][linen,2:(domainn)]
  stopseq<-rank(ps[1:(domainn-1)],ties.method = 'min')
  decision<-datain1[[x]][[3]][linen,2:(domainn)]
  population<-datapop[[jj]][[x]][[1]]
  enrollinpopulation<-as.numeric(datain1[[x]][[5]][,2])
  enrollinpopulation1<-as.numeric(datain1[[x]][[5]][,3])
  
  trt<-matrix(1:(domainn-1))
  stopdat<-as.data.frame(cbind(trt,ps,stopseq,decision))
  colnames(stopdat)<-c("realarm","number","stopseq","decision")
  stopdat<-stopdat[with(stopdat, order(stopseq, decision)),]
  stoptrial<-as.numeric(unique(stopdat$number))
  stopindex<-vector('list',length(stoptrial))

  if (length( which (decision %in%  'I'))==(domainn-1)){
    stopindex<- lapply(1:length(stoptrial),function(x){
      for (v in 1:(length(population))) {
        if (all.equal(population[v],max(c(enrollinpopulation[stoptrial[x]],enrollinpopulation1[1:stoptrial[[x]]]),na.rm=TRUE))==TRUE){
          return(v)
          break
        }else if (population[v]>max(c(enrollinpopulation[stoptrial[[x]]],enrollinpopulation1[1:stoptrial[[x]]]),na.rm=TRUE)){
          return(v-1)
          break
        }
      }
    })
  }else{

    stopindex11<- function(x){
      for (v in 1:(length(population))) {
        if (all.equal(population[v],max(c(enrollinpopulation[stoptrial[[x]]],enrollinpopulation1[1:stoptrial[[x]]]),na.rm=TRUE))==TRUE){
          return(v)
          break
        } else if (population[v]>max(c(enrollinpopulation[stoptrial[[x]]],enrollinpopulation1[1:stoptrial[[x]]]),na.rm=TRUE)){
          return(v-1)
          break
        }
      }
    }
    stopindex12<-  function(x){
      for (v in 1:(length(population))) {
        if (population[v]>enrollinpopulation1[stoptrial[x]]){
          return(v-1)
          break
        }
      }
    }
    for(l in 1:length(stoptrial)){
      if (l<length(stoptrial)) {
        stopindex[[l]]<-stopindex11(l)
      }else{
        stopindex[[l]]<-stopindex12(l)
      }
    }
  }

  stoppop<- do.call(rbind,lapply(1:length(stoptrial),
                                 function(x){(stopindex[[x]]-stoptrial[x])/2}))

  if (stopindex[[length(stopindex)]]!=length(population)){
    stoppop<-c(stoppop,(length(population)-stopindex[[length(stopindex)]])/2)
  }

  percttot2<-rep(0,length(unique(stopdat$stopseq)))

  if (length(unique(stopdat$stopseq))>1){
    for (i in 1:(length(unique(stopdat$stopseq))-1)) {
      dat1<-subset(stopdat,stopdat$stopseq<=unique(stopdat$stopseq)[i])
      dum1<-  as.numeric(subset(dat1,dat1$decision== "I")$realarm)+1
      dum2<- subset(stopdat,stopdat$stopseq>unique(stopdat$stopseq)[i])$realarm
      percttot2[which(unique(stopdat$stopseq)==unique(stopdat$stopseq)[as.numeric(i)+1])]= (sum(perct[dum1])/(length(dum2)+1))*
        (stoppop[which(unique(stopdat$stopseq)==unique(stopdat$stopseq)[as.numeric(i)+1])]-
           stoppop[which(unique(stopdat$stopseq)==unique(stopdat$stopseq)[as.numeric(i)])])
    }
  }
  percttot21<-cumsum(percttot2)
  patientintrialt1<-vector('list',domainn)
  patientinpopulation<-vector('list',domainn)

  for (j in 1:domainn){
    patientintrialt1[[j]]<-coutcome[j] #sum(substr(coutcome,position,position)==j)#length(coutcome[coutcome==j])
  }
  patientintrialt<-matrix(do.call(rbind,patientintrialt1),nrow=1,ncol=domainn)


  if (sum(decision %in% "S")>0) {
    for (j in 2:(domainn)){
      seqq<-subset(stopdat,stopdat$realarm==(j-1))$stopseq
      seqq1<-which(unique(stopdat$stopseq) %in% seqq)
      if(decision[j-1] %in% "S"){
        patientinpopulation[[j]]<-ceiling(perct[j]*stoppop[seqq1]+
                                            percttot21[seqq1]+stoppop[length(stoppop)]/(sum(decision %in% "S")))
      }else if (decision[j-1] %in% "ND" ){
        patientinpopulation[[j]]<-ceiling(perct[j]*stoppop[seqq1]+
                                            percttot21[seqq1])
      }else if (decision[j-1] %in% "I"){
        patientinpopulation[[j]]<-ceiling(perct[j]*stoppop[seqq1]+
                                            percttot21[seqq1])

      }
    }
  }else {
    for (j in 2:(domainn)){
      seqq<-subset(stopdat,stopdat$realarm==(j-1))$stopseq
      seqq1<-which(unique(stopdat$stopseq) %in% seqq)
      patientinpopulation[[j]]<-ceiling(perct[j]*stoppop[seqq1]+
                                          percttot21[seqq1])
    }
  }


  if(sum(decision %in% "S")>0){
    patientinpopulation[[1]]<-ceiling(perct[1]*max(stoppop[1:(length(stoppop)-1)])+max(percttot21))
  } else {
    patientinpopulation[[1]]<-ceiling(perct[1]*max(stoppop[1:(length(stoppop)-1)])+max(percttot21)+(stoppop[length(stoppop)]-stoppop[length(stoppop)-1]))
  }

  patientinpopulationt1<-do.call(rbind, patientinpopulation)
  patientinpopulationt<-patientinpopulationt1/max(2*sum(stoppop[(length(stoppop)-1):(length(stoppop))]),2*sum( patientinpopulationt1))
  patientinpopulationt<-matrix(patientinpopulationt,nrow=1,ncol=domainn)

  return(list(patientintrialt,patientinpopulationt))
}


simbrarp<-function(datain1,datapop,jj,x){

  decision<-datain1[[x]][[3]]
  population<-datapop[[jj]][[x]][[1]]
  enrollinpopulation<-as.numeric(datain1[[x]][[5]][,2])
  enrollinpopulation1<-as.numeric(datain1[[x]][[5]][,3])
  triallen<-max(enrollinpopulation1)

  aa1<-ptsinpop1(x=x,coutcome=datain1[[x]][[7]][1:5,'Freqt'],datain1=datain1,domainn=5,linen=1,perct=c(0,0,0,0,0.5),position=5)
  aa2<-ptsinpop1(x=x,coutcome=datain1[[x]][[7]][6:10,'Freqt'],datain1=datain1,domainn=5,linen=2,perct=c(0,0,0,0,0.5),position=5)
  ss1<-ptsinpop2(x=x,coutcome=datain1[[x]][[8]][1:3,'Freqt'],datain1=datain1,domainn=3,linen=3,perct=c(0.5,0,0),position=8)
  ss2<-ptsinpop2(x=x,coutcome=datain1[[x]][[8]][4:6,'Freqt'],datain1=datain1,domainn=3,linen=4,perct=c(0.5,0,0),position=8)
  mm1<-ptsinpop2(x=x,coutcome=datain1[[x]][[9]][1:2,'Freqt'],datain1=datain1,domainn=2,linen=5,perct=c(0.5,0),position=11)
  mm2<-ptsinpop2(x=x,coutcome=datain1[[x]][[9]][3:4,'Freqt'],datain1=datain1,domainn=2,linen=6,perct=c(0.5,0),position=11)


  nn1<-max(length(aa1[[1]]),length(aa2[[1]]),length(ss1[[1]]),length(ss2[[1]]),
          length(mm1[[1]]),length(mm2[[1]]))
  length(aa1[[1]])<-nn1
  length(aa2[[1]])<-nn1
  length(ss1[[1]])<-nn1
  length(ss2[[1]])<-nn1
  length(mm1[[1]])<-nn1
  length(mm2[[1]])<-nn1
  patientintrialt<-rbind(aa1[[1]],aa2[[1]],ss1[[1]],ss2[[1]],mm1[[1]],mm2[[1]])

  nn2<-max(length(aa1[[2]]),length(aa2[[2]]),length(ss1[[2]]),length(ss2[[2]]),
           length(mm1[[2]]),length(mm2[[2]]))
  length(aa1[[2]])<-nn2
  length(aa2[[2]])<-nn2
  length(ss1[[2]])<-nn2
  length(ss2[[2]])<-nn2
  length(mm1[[2]])<-nn2
  length(mm2[[2]])<-nn2
  patientinpopulationt<-rbind(aa1[[2]],aa2[[2]],ss1[[2]],ss2[[2]],mm1[[2]],mm2[[2]])

  return(list(patientinpopulationt,patientintrialt,triallen,decision))
}



#############################
enrollratef<-0.9

repn<-1000
sim11<-vector("list",length(enrollratef))
jj=1
  seeds<-1:repn
  sim11[[jj]]<-lapply(1:repn,function(x) {
     set.seed(seeds[x])
    pop(vPatsPerMonth=119,nMaxQtyPats=50000,enrollrate1=0.9)})


#######################
set.seed(12345)
sim11a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(0,0,0,0,0),
           sss=c(0,0,0),mmm=c(0,log(1.2)),interaction1=c(0,0))},mc.cores = 60)

set.seed(12345)
sim11b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim11a,datapop=sim11,jj=1,x=x)},mc.cores = 60)


########################
set.seed(12345)
sim22a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(-0.2,-0.2,-0.2,-0.2,0),
           sss=c(0,-0.17,-0.17),mmm=c(0,-0.05),interaction1=c(-0.1,-0.1))},mc.cores = 60)


sim22b<-parallel::mclapply(1:1000,function(x){
  set.seed(12345)
  simbrarp(datain1=sim22a,datapop=sim11,jj=1,x=x)},mc.cores = 60)

########################
set.seed(12345)
sim33a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(-0.8,-0.8,-0.8,-0.8,0),
           sss=c(0,-0.8,-0.8),mmm=c(0,-0.8),interaction1=c(-0.8,-0.8))},mc.cores = 60)


set.seed(12345)
sim33b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim33a,datapop=sim11,jj=1,x=x)},mc.cores = 60)


########################
set.seed(12345)
sim44a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(-0.5,-0.5,-0.5,-0.5,0),
           sss=c(0,-0.5,-0.5),mmm=c(0,-0.5),interaction1=c(-0.5,-0.5))},mc.cores = 60)


set.seed(12345)
sim44b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim44a,datapop=sim11,jj=1,x=x)},mc.cores = 60)


########################
set.seed(12345)
sim55a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(-0.2,-0.4,-0.5,-0.7,0),
           sss=c(0,-0.2,-0.4),mmm=c(0,-0.3),interaction1=c(-0.4,-0.4))},mc.cores = 60)

set.seed(12345)
sim55b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim55a,datapop=sim11,jj=1,x=x)},mc.cores = 60)


########################
set.seed(12345)
sim66a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(-0.2,-0.3,-0.25,-0.3,0),
           sss=c(0,-0.2,-0.24),mmm=c(0,-0.3),interaction1=c(-0.25,-0.25))},mc.cores = 60)


set.seed(12345)
sim66b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim66a,datapop=sim11,jj=1,x=x)},mc.cores = 60)

########################
set.seed(12345)
sim77a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(-0.6,-0.5,-0.55,-0.62,0),
           sss=c(0,-0.4,-0.4),mmm=c(0,-0.2),interaction1=c(-0.3,-0.3))},mc.cores = 60)

set.seed(12345)
sim77b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim77a,datapop=sim11,jj=1,x=x)},mc.cores = 60)

########################
set.seed(12345)
sim88a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(0.2,0.2,0.2,0.2,0),
           sss=c(0,0.2,0.2),mmm=c(0,0.2),interaction1=c(0.2,0.2))},mc.cores = 60)

set.seed(12345)
sim88b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim88a,datapop=sim11,jj=1,x=x)},mc.cores = 60)

########################
set.seed(12345)
sim99a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(0.5,0.5,0.5,0.5,0),
           sss=c(0,0.5,0.5),mmm=c(0,0.5),interaction1=c(0.5,0.5))},mc.cores = 60)

set.seed(12345)
sim99b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim99a,datapop=sim11,jj=1,x=x)},mc.cores = 60)

#######################
set.seed(12345)
sim100a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(0.8,0.8,0.8,0.8,0),
           sss=c(0,0.8,0.8),mmm=c(0,0.8),interaction1=c(0.8,0.8))},mc.cores = 60)

set.seed(12345)
sim100b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim100a,datapop=sim11,jj=1,x=x)},mc.cores = 60)

#######################
set.seed(12345)
sim110a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(0.2,0.4,0.5,0.7,0),
           sss=c(0,0.2,0.4),mmm=c(0,0.3),interaction1=c(0.4,0.4))},mc.cores = 60)


set.seed(12345)
sim110b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim110a,datapop=sim11,jj=1,x=x)},mc.cores = 60)

########################
set.seed(12345)
sim120a<-parallel::mclapply(1:1000,function(x){
 snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(0.2,0.3,0.25,0.3,0),
           sss=c(0,0.2,0.24),mmm=c(0,0.3),interaction1=c(0.25,0.25))},mc.cores = 60)

set.seed(12345)
sim120b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim120a,datapop=sim11,jj=1,x=x)},mc.cores = 60)

########################
set.seed(12345)
sim130a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(0.6,0.5,0.55,0.62,0),
           sss=c(0,0.4,0.4),mmm=c(0,0.2),interaction1=c(0.3,0.3))},mc.cores = 60)

set.seed(12345)
sim130b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim130a,datapop=sim11,jj=1,x=x)},mc.cores = 60)

#######################
set.seed(12345)
sim140a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(0.3,-0.2,0.15,-0.1,0),
           sss=c(0,0.2,-0.1),mmm=c(0,0.4),interaction1=c(-0.35,0.35))},mc.cores = 60)

set.seed(12345)
sim140b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim140a,datapop=sim11,jj=1,x=x)},mc.cores = 60)


#######################
set.seed(12345)
sim150a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(0.6,-0.4,0.25,-0.8,0),
           sss=c(0,0.5,-0.25),mmm=c(0,0.3),interaction1=c(0.2,0.35))},mc.cores = 60)

set.seed(12345)
sim150b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim150a,datapop=sim11,jj=1,x=x)},mc.cores = 60)


########################
set.seed(12345)
sim160a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(0.33,-0.4,0.2,0.5,0),
           sss=c(0,0.24,0.25),mmm=c(0,-0.5),interaction1=c(0.2,-0.2))},mc.cores = 60)

set.seed(12345)
sim160b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim160a,datapop=sim11,jj=1,x=x)},mc.cores = 60)


########################
set.seed(12345)
sim17a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(-0.18,-0.16,-0.14,-0.12,0),
           sss=c(0,-0.1,-0.1),mmm=c(0,0.1),interaction1=c(-0.1,-0.1))},mc.cores = 60)

set.seed(12345)
sim17b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim17a,datapop=sim11,jj=1,x=x)},mc.cores = 60)


########################
set.seed(12345)
sim18a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(-0.08,-0.06,-0.04,-0.02,0),
           sss=c(0,-0.05,-0.05),mmm=c(0,0.05),interaction1=c(-0.05,-0.05))},mc.cores = 60)

set.seed(12345)
sim18b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim18a,datapop=sim11,jj=1,x=x)},mc.cores = 60)


########################
set.seed(12345)
sim19a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(-0.18,-0.16,-0.14,-0.12,0),
           sss=c(0,-0.05,-0.1),mmm=c(0,0.1),interaction1=c(-0.06,-0.08))},mc.cores = 60)

set.seed(12345)
sim19b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim19a,datapop=sim11,jj=1,x=x)},mc.cores = 60)


########################
set.seed(12345)
sim20a<-parallel::mclapply(1:1000,function(x){
  snap_fun(x=x,intercept1=c(log(0.2/(1-0.2)),log(0.35/(1-0.35))),aaa=c(-0.08,-0.06,-0.04,-0.02,0),
           sss=c(0,-0.08,-0.06),mmm=c(0,0.05),interaction1=c(0.04,-0.04))},mc.cores = 60)

set.seed(12345)
sim20b<-parallel::mclapply(1:1000,function(x){
  simbrarp(datain1=sim20a,datapop=sim11,jj=1,x=x)},mc.cores = 60)

