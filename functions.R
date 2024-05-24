library("stringr")  
superfun<-function(dat,category1,decisionf1,decisionf2,type,stopnum,stopp1,stopp2,jjj,
                   stopA1,stopA2,stopS1,stopS2,stopM1,stopM2){
  

   stop1<-data.frame(t(apply(dat, 1,exp)))#exp(dat)
   colnames(stop1)<-c("aS[1]", "aS[2]" ,"bA[1]","bA[2]" ,"bA[3]","bA[4]", "bS[2]", "bS[3]" , "bM[2]"  ,  
                  "aSBS[2,2]" ,"aSBS[2,3]")
  stop1$shock_s2<-stop1$'bS[2]'*stop1$'aSBS[2,2]'
  stop1$shock_s3<-stop1$'bS[3]'*stop1$'aSBS[2,3]'
  stop1<-stop1[,!(names(stop1) %in% c('aS[1]','aS[2]', 'aSBS[2,2]','aSBS[2,3]'))]
  
  stop11<-do.call(rbind,lapply(1:ncol(stop1), function(x) {sum(stop1[,x]<1)/20000}))
  rownames(stop11)<-c('bA[1]','bA[2]','bA[3]','bA[4]','bS[2]','bS[3]',
                      'bM[2]', 'shock_S2','shock_S3')

  
  stop12<-do.call(rbind,lapply(1:ncol(stop1), function(x) {sum(stop1[,x]<(1/1.2))/20000}))
  rownames(stop12)<-c('bA[1]','bA[2]','bA[3]','bA[4]','bS[2]','bS[3]',
                      'bM[2]', 'shock_S2','shock_S3')

  
  stop13<-do.call(rbind,lapply(1:ncol(stop1), function(x) {sum(stop1[,x]<1.2)/20000}))
  rownames(stop13)<-c('bA[1]','bA[2]','bA[3]','bA[4]','bS[2]','bS[3]',
                      'bM[2]', 'shock_S2','shock_S3')

  if (jjj>1){
    stopA11<- gsub(" ", "", paste('bA[',stopA1,']'))
    stopA22<- gsub(" ", "", paste('bA[',stopA2,']'))
    stopM11<-gsub(" ", "", paste('bM[',stopM2,']'))
    stopM22<-gsub(" ", "", paste('bM[',stopM1,']'))
    stopS22<-gsub(" ", "", paste("shock_S",stopS2))
    stopS11<-gsub(" ", "", paste('bS[',stopS1,']'))
    
    stop13<-stop13[rownames(stop13) %in% c(stopM11,stopM22),,drop=F]
    stop11<-stop11[rownames(stop11) %in% c(stopA11,stopS11,stopA22,stopS22),,drop=F]
    stop12<-stop12[rownames(stop12) %in% c(stopA11,stopS11,stopA22,stopS22),,drop=F]

  }

 if (type=='superiority'){
   if (category1!='S'){
     stop11<-stop11[str_detect(rownames(stop11), category1),,drop=F]
     superiority1<-stop11[ stop11[,1]>0.99 , ,drop=F]#stop11[ stop11[str_detect(rownames(stop11), category1),1]>0.99 , ,drop=F]#shock
     superiority2<-stop11[ stop11[,1]>0.99 , ,drop=F]#stop11[ stop11[str_detect(rownames(stop11), category1),1]>0.99 , ,drop=F ]#non-shock
     
     stop12<-stop12[str_detect(rownames(stop12), category1),,drop=F]
     inferiority1<-stop12[ !stop12[,1]<0.01 , ,drop=F]#stop12[( !stop12[str_detect(rownames(stop12), category1),1]<0.01), ,drop=F]
     inferiority2<-stop12[ !stop12[,1]<0.01 , ,drop=F]#stop12[( !stop12[str_detect(rownames(stop12), category1),1]<0.01), ,drop=F]
     inferioritydrop1<-stop12[ stop12[,1]<0.01 , ,drop=F]#stop12[stop12[str_detect(rownames(stop12), category1),1]<0.01, ,drop=F]
     inferioritydrop2<-stop12[ stop12[,1]<0.01 , ,drop=F]#stop12[stop12[str_detect(rownames(stop12), category1),1]<0.01, ,drop=F]
     
   }else if (category1=='S'){
 
     
     stop11<-stop11[str_detect(rownames(stop11), category1),,drop=F]
     stop112<-stop11[rownames(stop11) %in% c('shock_S2','shock_S3'),1 , drop=FALSE]
     stop111<-stop11[rownames(stop11) %in% c('bS[2]','bS[3]'),1 , drop=FALSE]
     superiority1<-stop111[ stop111[,1]>0.99 , ,drop=F]#stop11[ stop11[str_detect(rownames(stop11), category1),1]>0.99 , ,drop=F]#shock
     superiority2<-stop112[ stop112[,1]>0.99 , ,drop=F]#stop11[ stop11[str_detect(rownames(stop11), category1),1]>0.99 , ,drop=F ]#non-shock
     
     stop12<-stop12[str_detect(rownames(stop12), category1),,drop=F]
     stop122<-stop12[rownames(stop12) %in% c('shock_S2','shock_S3'),1 , drop=FALSE]
     stop121<-stop12[rownames(stop12) %in% c('bS[2]','bS[3]'),1 , drop=FALSE]
     inferiority1<-stop121[!stop121[,1]<0.01, ,drop=F]
     inferiority2<-stop122[!stop122[,1]<0.01, ,drop=F]
     inferioritydrop1<-stop121[stop121[,1]<0.01, ,drop=F]
     inferioritydrop2<-stop122[stop122[,1]<0.01, ,drop=F]

   }
   
   
 }else if (type=='non-inferiority'){
   if (category1!='S'){
     stop13<-stop13[str_detect(rownames(stop13), category1),,drop=F]
     superiority1<-stop13[ stop13[,1]>0.99 , ,drop=F]#stop13[ stop13[str_detect(rownames(stop13), category1),1]>0.99 , ,drop=F]#shock
     superiority2<-stop13[ stop13[,1]>0.99 , ,drop=F]#stop13[ stop13[str_detect(rownames(stop13), category1),1]>0.99 , ,drop=F ]#non-shock
     
     stop13<-stop13[str_detect(rownames(stop13), category1),,drop=F]
     inferiority1<-stop13[ !stop13[,1]<0.01 , ,drop=F]#stop13[( !stop13[str_detect(rownames(stop13), category1),1]<0.01), ,drop=F]
     inferiority2<-stop13[ !stop13[,1]<0.01 , ,drop=F]#stop13[( !stop13[str_detect(rownames(stop13), category1),1]<0.01), ,drop=F]
     inferioritydrop1<-stop13[ stop13[,1]<0.01 , ,drop=F]#stop13[stop13[str_detect(rownames(stop13), category1),1]<0.01, ,drop=F]
     inferioritydrop2<-stop13[ stop13[,1]<0.01 , ,drop=F]#stop13[stop13[str_detect(rownames(stop13), category1),1]<0.01, ,drop=F]
     
   }else if (category1=='S'){
    # stop13<-stop13[str_detect(rownames(stop13), category1),,drop=F]
       
     
     stop13<-stop13[str_detect(rownames(stop13), category1),,drop=F]
     stop132<-stop13[rownames(stop13) %in% c('shock_S2','shock_S3'),1 , drop=FALSE]#stop13[1:2,1 , drop=FALSE]
     stop131<-stop13[rownames(stop13) %in% c('bS[2]','bS[3]'),1 , drop=FALSE]#stop13[3:4,1 , drop=FALSE]
     superiority1<-stop131[ stop131[,1]>0.99 , ,drop=F]#stop13[ stop13[str_detect(rownames(stop13), category1),1]>0.99 , ,drop=F]#shock
     superiority2<-stop132[ stop132[,1]>0.99 , ,drop=F]#stop13[ stop13[str_detect(rownames(stop13), category1),1]>0.99 , ,drop=F ]#non-shock
     
     inferiority1<-stop131[!stop131[,1]<0.01, ,drop=F]
     inferiority2<-stop132[!stop132[,1]<0.01, ,drop=F]
     inferioritydrop1<-stop131[stop131[,1]<0.01, ,drop=F]
     inferioritydrop2<-stop132[stop132[,1]<0.01, ,drop=F]

   }
 }


  if (nrow(superiority1)>=1){
    superiority11<- as.numeric(str_extract_all(rownames(superiority1), '[0-9]+') )
    decisionf1[superiority11]<-'S'#success
    if (nrow(inferioritydrop1)>0){
      decisionf1[as.numeric(str_extract_all(rownames(inferioritydrop1), '[0-9]+') ) ]<-'I'
      stopp1[as.numeric(str_extract_all(rownames(inferioritydrop1), '[0-9]+') )]<-stopnum
    }
      stopp1[superiority11]<-stopnum
      stopp1[which(decisionf1=='')]<-stopnum
      decisionf1[which(decisionf1=='')]<-'ND' #no decision

  }else{
    superiority11<-as.numeric(str_extract_all(rownames(inferiority1), '[0-9]+') ) 
    if (nrow(inferioritydrop1)>0){
      decisionf1[as.numeric(str_extract_all(rownames(inferioritydrop1), '[0-9]+') ) ]<-'I' #inforiority
      stopp1[as.numeric(str_extract_all(rownames(inferioritydrop1), '[0-9]+') ) ]<-stopnum
    }   
  }
  
  if (nrow(superiority2)>=1){
    superiority22<- as.numeric(str_extract_all(rownames(superiority2), '[0-9]+') )
    decisionf2[superiority22]<-'S'#success
    if (nrow(inferioritydrop2)>0){
      decisionf2[as.numeric(str_extract_all(rownames(inferioritydrop2), '[0-9]+') ) ]<-'I'
      stopp2[as.numeric(str_extract_all(rownames(inferioritydrop2), '[0-9]+') )]<-stopnum
    }
      stopp2[superiority22]<-stopnum
      stopp2[which(decisionf2=='')]<-stopnum
      decisionf2[which(decisionf2=='')]<-'ND' #no decision
  }else{
    superiority22<-as.numeric(str_extract_all(rownames(inferiority2), '[0-9]+') ) 
    if (nrow(inferioritydrop2)>0){
      decisionf2[as.numeric(str_extract_all(rownames(inferioritydrop2), '[0-9]+') ) ]<-'I' #inforiority
      stopp2[as.numeric(str_extract_all(rownames(inferioritydrop2), '[0-9]+') ) ]<-stopnum
    }   
  }

  return(list(superiority11,superiority22,decisionf1,decisionf2,stopp1,stopp2))
}  




