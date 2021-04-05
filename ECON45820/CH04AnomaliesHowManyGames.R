setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyr)
library(ggplot2)
library(stringr)
library(dplyr)
D<- (read.csv("Data/HowManyGames.csv")
     %>% filter(Period==20)     
     %>% select(uid,contains("chi"),contains("invest"))
)
IDList<-unique(D$uid)
chi1<-c(2.25,2.75,2.35,5.50)
chi2<-c(2.2,2.5,3.0,4.0)

Part1<-(D 
        %>% dplyr::select(uid,contains("invest1"))
        %>% gather(x,invest,invest1_1:invest1_4)
        %>% mutate(instance = as.integer(str_sub(x,-1,-1)),
                   Part = "Part 1")
        %>% select(uid,invest,instance,Part)
        %>% mutate(chi = chi1[instance])
)


Part2<-(D 
        %>% select(uid,contains("invest2"))
        %>% mutate(
          invest2_1=round(0.5*(invest2_1_1+invest2_2_1)),
          invest2_2=round(0.5*(invest2_1_2+invest2_2_2)),
          invest2_3=round(0.5*(invest2_1_3+invest2_2_3)),
          invest2_4=round(0.5*(invest2_1_4+invest2_2_4)),
          
        )
        %>% gather(x,invest,invest2_1:invest2_4)
        %>% mutate(instance = as.integer(str_sub(x,-1,-1)),
                   Part = "Part 2"
        )
        %>% select(uid,invest,instance,Part)
        %>% mutate(chi = chi2[instance])
        
)

LotteryData<-rbind(Part1,Part2)

LotteryData %>% head()

(ggplot(LotteryData,aes(x=invest))
  +geom_histogram(bins=15)
  +facet_grid(Part ~ instance)
  +theme_bw()
  
)

Y<-seq(0,100,length=101)

r<-seq(0.1,1,length=91)

EstimateRN<-function(data) {
  SSR<-r*0
  for (rr in 1:length(r)) {
    for (ii in 1:dim(data)[1])
    U<-0.5*(100-Y)^r[rr]+0.5*(100-Y+data$chi[ii]*Y)^r[rr]
    Ystar<-Y[which.max(U)]
    SSR[rr]<-SSR[rr]+(Ystar-data$invest[ii])^2
  }
  r[which.min(SSR)]
}

WWNBD<-function(data) {
  # estimate r from Part 1
  r<-EstimateRN(data %>% filter(Part=="Part 1"))
  p2Data<-data %>% filter(Part == "Part 2")
  
  Ystar<-c()
  for (ii in 1:dim(p2Data)[1]) {
    U<-0.5*(100-Y)^r+0.5*(100-Y+p2Data$chi[ii]*Y)^r
    Ystar[ii]<-Y[which.max(U)]
  }
  Ystar
}
WWBBD<-function(data) {
  # estimate r from Part 1
  r<-EstimateRN(data %>% filter(Part=="Part 1"))
  p2Data<-data %>% filter(Part == "Part 2")
  
  Ystar<-c()
  for (ii in 1:dim(p2Data)[1]) {
    U<-0.25*(200-2*Y)^r+0.5*(200-2*Y+p2Data$chi[ii]*Y)^r+0.25*(200-2*Y+2*p2Data$chi[ii]*Y)^r
    #U<-0.5*(100-Y)^r+0.5*(100-Y+p2Data$chi[ii]*Y)^r
    Ystar[ii]<-Y[which.max(U)]
  }
  Ystar
}
  

EstimateRB<-function(data) {
  SSR<-r*0
  for (rr in 1:length(r)) {
    for (ii in 1:dim(data)[1])
      U<-0.25*(200-2*Y)^r[rr]+0.5*(200-2*Y+data$chi[ii]*Y)^r[rr]+0.25*(200-2*Y+2*data$chi[ii]*Y)^r[rr]
    Ystar<-Y[which.max(U)]
    SSR[rr]<-SSR[rr]+(Ystar-data$invest[ii])^2
  }
  r[which.min(SSR)]
}

RPart1<-c()
RPart2n<-c()
RPart2b<-c()

predictions<-data.frame()

for (ii in 1:length(IDList)) {
  RPart1[ii]<-EstimateRN(LotteryData %>% filter(Part == "Part 1",uid==IDList[ii]))
  RPart2n[ii]<-EstimateRN(LotteryData %>% filter(Part == "Part 2",uid==IDList[ii]))
  RPart2b[ii]<-EstimateRB(LotteryData %>% filter(Part == "Part 2",uid==IDList[ii]))
  
  NB<-WWNBD(LotteryData %>% filter(uid==IDList[ii]))
  BB<-WWBBD(LotteryData %>% filter(uid==IDList[ii]))
  
  Actual<-LotteryData %>% filter(uid==IDList[ii] & Part == "Part 2") %>% select(invest)
  
  
  tmp<-data.frame(NB,BB,Actual)
  tmp$uid<-IDList[ii]
  tmp$r<-RPart1[ii]
  predictions<-predictions %>% rbind(tmp)
  
}

(ggplot()
  +geom_histogram(bins=10,aes(x=RPart1,fill="Part 1"),alpha=0.5)
  +geom_histogram(bins=10,aes(x=RPart2n,fill="Part 2"),alpha=0.5)
  +theme_bw()
  )

(ggplot(data.frame(RPart1,RPart2n,RPart2b))
  +geom_point(aes(x=RPart1,y=RPart2n,color="Narrow"))
  +geom_point(aes(x=RPart1,y=RPart2b,color="Broad"))
  +theme_bw()
  +geom_abline(slope=1,intercept=0)
  +geom_smooth(method="lm",aes(x=RPart1,y=RPart2n,color="Narrow"))
  +geom_smooth(method="lm",aes(x=RPart1,y=RPart2b,color="Broad"))
)

estimates<-data.frame(RPart1,RPart2n,RPart2b)

lm1<-lm(data=estimates,RPart2b~RPart1)
lm2<-lm(data=estimates,RPart2n~RPart1)

stargazer::stargazer(lm1,lm2,type="text")

(
  ggplot(predictions)
  +geom_point(aes(x=invest,y=NB,color="narrow"))
  +geom_smooth(aes(x=invest,y=NB,color="narrow"))
  +geom_point(aes(x=invest,y=BB,color="broad"))
  +geom_smooth(aes(x=invest,y=BB,color="broad"))
  +theme_bw()
  +geom_abline(slope=1,intercept=0)
)

lm3<-lm(data=predictions,invest~NB)
lm4<-lm(data=predictions,invest~BB)

stargazer::stargazer(lm3,lm4,type="text")

(ggplot(
  data=(predictions 
        %>% group_by(uid) 
        %>% summarize(MSEnarrow = mean((invest-NB)^2),
                      MSEbroad = mean((invest-BB)^2 ),
                      r = mean(r))
))
  +geom_point(aes(x=MSEnarrow,y=MSEbroad,color=r))
  +theme_bw()
  +xlab("MSE - narrow predictions")
  +ylab("MSE - broad predictions")
  +geom_abline(slope=1,intercept=0)
  +scale_x_continuous(trans="log10")
  +scale_y_continuous(trans="log10")
  
  )
