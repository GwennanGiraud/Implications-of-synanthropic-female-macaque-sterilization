#######################Script R Data analysis Impact sterilization on females' social dynamics#######################################

#Book reference : Imbens & Rubin, 2015 Causal inference for statistics, social, and biomedical sciences: An introduction

rm(list = ls(all.names = TRUE)) 

durcomp <- read.csv("Dataset.csv", header=TRUE, sep=";") 
durtot<- read.csv("Sampling_effort.csv", header=TRUE, sep=";") 
attribut <- read.csv("Attributs.csv", header=TRUE, sep=";") 

durcomp2<-durcomp
#####################################
c<-1                  # variables (behaviors) to test: 1 to 14

# Test T1 = mean differences after sterilization
# Test T2= mean differences of after-before sterilization differences
# 1=AGORF,2=AgoPARTRF, 3=AGOEF, 4=AGOPARTEF, 5=SexProPart, 6=SexPro, 7=SexRecPart, 
# 8=SexRec, 9=SexAttPart, 10=SexAtt, 11=GroomRM, 12=GroomPartRM, 13=GroomEM, 
# 14=GroomPartEM)

stenum<-vector() 
durcomp3<-as.data.frame(matrix(nrow=35)) 
durcomp4<-as.data.frame(matrix(nrow=35))

for (i in 1:56){          # 56 females
  k<-1+c+(i-1)*14     
  print (k)
  if (grepl("Part",colnames(durcomp2)[k])){  
    division<-TRUE                           
    pc<-durcomp2[,k]/(durtot[,i+1]/3600)}  
  else{                                   
    division<-FALSE                          
    pc<-durcomp2[,k]}   
  
  t<-durtot[,i+1]                     
  
  durcomp3<-cbind(durcomp3,pc)   
  durcomp4<-cbind(durcomp4,t)     
  
  sta<-attribut[,i+1]                     
  if(sum(sta>0)){stenum<-append(stenum,i)} 
}

durcomp3<-durcomp3[,2:57]     
durcomp4<-durcomp4[,2:57]    

colnames(durcomp3)<-names(durtot)[2:57] 
colnames(durcomp4)<-names(durtot)[2:57] 

v=1
w=1
for (v in 1:35){       
  for (w in 1:56){     
    if (!is.na (durcomp4[v,w])&durcomp4[v,w]==0){durcomp3[v,w]<-NA} 
    if (is.na (durcomp3[v,w])){durcomp4[v,w]<-NA} 
  }
}

stenum
stetime<-attribut[,1+(stenum)] 


# Statistical test starts here
mf<-1000

T<-as.data.frame(matrix(rep(NA,8*mf),ncol=8)) 
colnames(T)<-c("T1","pT1","mT1","sdT1","T2","pT2","mT2","sdT2") 

av0<-matrix(rep(NA,39*mf),nrow=39) 
ap0<-matrix(rep(NA,39*mf),nrow=39) 

for (m in 1:mf){                   
  attribut2<-attribut
  steril<-rep(0,56)                                  
  
  for (i in 1:56){   
    if(sum(attribut[,i+1])>0){steril[i]<-1}         
  } 
  li<-round(runif(39,1,17),0)       
  n<-1                                       
  for (i in 1:56){                                      
    if(steril[i]==0){attribut2[,1+i]<-stetime[,li[n]]
    n<-n+1}                                            
  }                                          
  avant<-rep(NA,56) 
  apres<-rep(NA,56)
  totavant<-rep(NA,56)
  totapres<-rep(NA,56)
  
  i=1
  for (i in 1:56){                
    obs<-durcomp3[,i]
    totobs<-durcomp4[,i]
    if (division==FALSE){                              
      avant[i]<-sum(obs[attribut2[,i+1]==0],na.rm=TRUE) 
      apres[i]<-sum(obs[attribut2[,i+1]==1],na.rm=TRUE) 
      totavant[i]<-sum(totobs[attribut2[,i+1]==0],na.rm=TRUE) 
      if (totavant[i]==0){totavant[i]<-NA}         
      totapres[i]<-sum(totobs[attribut2[,i+1]==1],na.rm=TRUE)
      
      avant[i]<-avant[i]/totavant[1]*1000               #Correction for freq/h (*1000) or for duration (*100)
      apres[i]<-apres[i]/totapres[1]*1000}              #Correction for freq/h or for duration (*100)
    else {                                               
      avant[i]<-mean(obs[attribut2[,i+1]==0],na.rm=TRUE) 
      apres[i]<-mean(obs[attribut2[,i+1]==1],na.rm=TRUE) 
      totavant[i]<-sum(totobs[attribut2[,i+1]==0],na.rm=TRUE) 
      totapres[i]<-sum(totobs[attribut2[,i+1]==1],na.rm=TRUE)
    }}
  
  dataset1<-as.data.frame(cbind(steril,avant,totavant,apres,totapres)) 
  dataset1$totavant[dataset1$totavant==0]<-NA 
  dataset1$totapres[dataset1$totapres==0]<-NA
  write.table(dataset1,"dataset1.txt",sep="\t",row.names=F,col.names=TRUE,quote=F)
  
  T$T1[m]<-abs(mean(apres[steril==0],na.rm=TRUE)-mean(apres[steril==1],na.rm=TRUE))   
  
  T$T2[m]<-abs(mean(apres[steril==0]-avant[steril==0],na.rm=TRUE)-mean(apres[steril==1]-avant[steril==1],na.rm=TRUE))
  
  if (m==1){
    av1<-avant[steril==1] 
    ap1<-apres[steril==1] 
  }
  
  Trd<-matrix(rep(NA,2000),ncol=2)
  n=1
  for (n in 1:1000){  #Test Imbens & Rubin
    rd<-sample(steril, replace=FALSE) 
    # l'attribution du traitement
    Trd[n,1]<-abs(mean(apres[rd==0],na.rm=TRUE)-mean(apres[rd==1],na.rm=TRUE)) 
    Trd[n,2]<-abs(mean(apres[rd==0]-avant[rd==0],na.rm=TRUE)-mean(apres[rd==1]-avant[rd==1],na.rm=TRUE)) 
  }
  av0[,m]<-avant[steril==0] 
  ap0[,m]<-apres[steril==0] 
  
  #T1
  T$pT1[m]<-table(Trd[,1]<T$T1[m])[1]/1000 
  T$mT1[m]<-mean(Trd[,1], na.rm=TRUE) 
  T$sdT1[m]<-sd(Trd[,1],na.rm=TRUE)  
  
  #T2
  T$pT2[m]<-table(Trd[,2]<T$T2[m])[1]/1000
  T$mT2[m]<-mean(Trd[,2],na.rm=TRUE)
  T$sdT2[m]<-sd(Trd[,2],na.rm=TRUE)
}

write.table(T,paste("T-",c,".txt",sep=""),sep="\t",row.names=F,col.names=TRUE,quote=F) 
write.table(ap0,paste("ap0-",c,".txt",sep=""),sep="\t",row.names=F,col.names=TRUE,quote=F) 
write.table(av0,paste("av0-",c,".txt",sep=""),sep="\t",row.names=F,col.names=TRUE,quote=F)
write.table(ap1,paste("ap1-",c,".txt",sep=""),sep="\t",row.names=F,col.names=TRUE,quote=F)
write.table(av1,paste("av1-",c,".txt",sep=""),sep="\t",row.names=F,col.names=TRUE,quote=F)



###Results###

###T1: After sterilisation
Meanav1<-mean(av1,na.rm=TRUE)
Meanav1
Meanav0<-mean(av0,na.rm=TRUE) 
Meanav0
Meanap1<-mean(ap1,na.rm=TRUE) 
Meanap1
Meanap0<-mean(ap0,na.rm=TRUE) 
Meanap0

table(T$pT1<0.05) 
R.T1<-mean(T$T1, na.rm=TRUE) 
R.T1
R.mT1<-mean(T$mT1, na.rm=TRUE)  
R.mT1
R.pT1<-mean(T$pT1, na.rm=TRUE) 
R.pT1
R.sdT1<-mean(T$sdT1, na.rm=TRUE) 
R.sdT1


###T2: After-before sterilisation
Meanap1av1<-mean(ap1-av1,na.rm=TRUE) 
Meanap0av0<-mean(ap0-av0,na.rm=TRUE)
Meanap0av0

table(T$pT2<0.05)
R.T2<-mean(T$T2, na.rm=TRUE) 
R.T2
R.mT2<-mean(T$mT2, na.rm=TRUE)
R.mT2
R.pT2<-mean(T$pT2, na.rm=TRUE)
R.pT2
R.sdT2<-mean(T$sdT2, na.rm=TRUE)
R.sdT2
