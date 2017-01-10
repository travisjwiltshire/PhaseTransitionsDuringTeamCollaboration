# Loading packages, functions, and data -----------------------------------

require(tseriesChaos)
require(plyr)
require(entropy)
require(zoo)
require(ggplot2)
require(Hmisc)
require(plotrix)
require(dplyr)
require(tidyr)



#data should be in long format with a repeating ID variable for each team and a column for the communications code time series
data<-read.csv(file= "CommCodesData_Example.csv",header=TRUE, sep=",")
data$TeamIDvar<-as.numeric(as.factor(data$TeamID))

list<-unique(data$TeamIDvar)
for (j in 1:length(list)){
  case<-list[j]
  subdata <- subset(data, data$TeamIDvar==case)
  subdata$Order<-1:nrow(subdata)
  if(j==1){
    data_c<-subdata
  } else{
    data_c<-rbind(data_c,subdata)
  }
}




#create cbind function needed for windowed entropy loop
cbind.fill <- function(...) {                                                                                                                                                       
  transpoted <- lapply(list(...),t)                                                                                                                                                 
  transpoted_dataframe <- lapply(transpoted, as.data.frame)                                                                                                                         
  return (data.frame(t(rbind.fill(transpoted_dataframe))))                                                                                                                          
} 

#This function converts the input time series of categorical codes into the required observed counts vector, calculates the entropy, and returns the value 
CountAndEntropy <- function(x) {
  require(entropy)
  C1<-sum(x==1)
  C2<-sum(x==2)
  C3<-sum(x==3)
  C4<-sum(x==4)
  C5<-sum(x==5)
  C6<-sum(x==6)
  C7<-sum(x==7)
  C8<-sum(x==8)
  C9<-sum(x==9)
  C10<-sum(x==10)
  C11<-sum(x==11)
  C12<-sum(x==12)
  C13<-sum(x==13)
  C14<-sum(x==14)
  C15<-sum(x==15)
  temp<-c(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15)
  entout<-entropy(temp)
  return(entout)
}  




# Examining AMI for appropriate window size -------------------------------


#This loop calculates average mutual information for each time series
mutualdata<- c()
list<-unique(data$TeamIDvar)
for (i in 1:length(list)){
  rm(res)
  case<-list[i]
  subdata <- subset(data, data$TeamIDvar==case)
  res = mutual(subdata$Code, lag.max=100)  
  res = as(res,"numeric")
  if(i==1){
    mutualdata<-res
  } else{
    mutualdata<-cbind.fill(mutualdata,res)
  }
}

# This loop returns a vector of the lags for first local minima in each of the AMI series and then descriptives can be calculated 
amiMins<-c(1:NCOL(mutualdata))
for (j in 1:NCOL(mutualdata)){
  amiMins[j]<-which.min(mutualdata[,j])
  }
mean(amiMins)
range(amiMins)
sd(amiMins)


# Calculating sliding window entropy --------------------------------------


#here is an example way to calculate and plot the results of sliding window entropy on a single team
#width determines the window size and by determines the iteration interval (i.e., how far to slide the window)
Team1Comm<-subset(data,TeamIDvar==1)

windowentropy25<-rollapply(Team1Comm$Code, width=25,by=1, FUN=CountAndEntropy,align="left") #this one seems to look most appropriate
plot(windowentropy25, type='l', ylab="entropy", xlab="Window Number/Time")


#calculate sliding windpow entropy for each team
#return results in long data format where all values are in one column
#also adds in grouping variable and order number (Time)
winentropydatalong <- NULL
list<-unique(data$TeamIDvar)
for (j in 1:length(list)){
  rm(res)
  case<-list[j]
  subdata <- subset(data, data$TeamIDvar==case)
  res = as.data.frame(rollapply(subdata$Code, width=25,by=1, FUN=CountAndEntropy,align="left"))    
  res$Order<-1:nrow(res)
  temp<-cbind(res,rep(subdata$TeamIDvar[1],times=nrow(res)))
  if(j==1){
    winentropydatalong<-temp
  } else{
    winentropydatalong<-rbind(winentropydatalong,temp)
  }
}
names(winentropydatalong)<-c("entropy","Order","TeamID")

#write out the file
write.csv(winentropydatalong,file="WinEntropyDataLong.csv",row.names=FALSE,na="-999")





# Peak identification -----------------------------------------------------
#This uses the winentropydatalong data frame/file created above

#Generate a lead and a lag and remove rows with missing data
list<-unique(winentropydatalong$TeamID)
for (j in 1:length(list)){
  case<-list[j]
  subdata <- subset(winentropydatalong, winentropydatalong$TeamID==case)
  subdata$entropy_ld1<-Lag(subdata$entropy,-1)
  subdata$entropy_lg1<-Lag(subdata$entropy,1)
  if(j==1){
    data1<-subdata
  } else{
    data1<-rbind(data1,subdata)
  }
}
data1<-na.omit(data1)

#Example loop that creates a peak point variable based on current entropy value being higher than preceding and following entropy value (lead 1 and lag 1)
for (i in 1:length(data1$entropy)){
  if(data1$entropy[i]>data1$entropy_ld1[i] && data1$entropy[i]>data1$entropy_lg1[i]){
    data1$peakpoint[i]<-1
  } else {
    data1$peakpoint[i]<-0
  }
}

#example plot for a single team showing identified peak points
#The xlim needs to be manually adjusted here to the lenght of the time series
Team1<-subset(data1,TeamID==1)
ss<-subset(Team1,Team1$peakpoint==1)
p<-ggplot(data=Team1, aes(x=Order,y=entropy))+geom_line(linetype=1)+geom_point() + geom_vline(xintercept=ss$Order)+xlim(1,347)+labs(list(title="Entropy Time Series", x= "Time Unit", y = "Entropy Value"))


# Smoothing procedure----------------------------------------------------


#create smoothing varables using a moving average from 2-15 time points
#the values in parentheses can be modified if certain time points should be weighted differently
ma2<-c(1,1)/2
ma3<-c(1,1,1)/3
ma4<-c(1,1,1,1)/4
ma5<-c(1,1,1,1,1)/5
ma6<-c(1,1,1,1,1,1)/6
ma7<-c(1,1,1,1,1,1,1)/7
ma8<-c(1,1,1,1,1,1,1,1)/8
ma9<-c(1,1,1,1,1,1,1,1,1)/9
ma10<-c(1,1,1,1,1,1,1,1,1,1)/10
ma11<-c(1,1,1,1,1,1,1,1,1,1,1)/11
ma12<-c(1,1,1,1,1,1,1,1,1,1,1,1)/12
ma13<-c(1,1,1,1,1,1,1,1,1,1,1,1,1)/13
ma14<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1)/14
ma15<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/15

#Example smoothing for a single team using the moving average of 5
smEnt<-as.data.frame(na.omit(stats::filter(Team1$entropy, ma5)))

#apply this loop to orignal windowed entropy data frame
#be sure to change ma# in filter function if a different moving average window is desired
smoothdata <- NULL
list<-unique(winentropydatalong$TeamID)
for (j in 1:length(list)){
  rm(res)
  case<-list[j]
  subdata <- subset(winentropydatalong, winentropydatalong$TeamID==case)
  res = as.data.frame(as.numeric(na.omit(stats::filter(subdata$entropy, ma5))))  #modify this line for different moving averages
  res$Order<-1:nrow(res)
  temp<-as.data.frame(cbind(res,rep(subdata$TeamID[j],times=nrow(res))))
  if(j==1){
    smoothdata<-temp
  } else{
    smoothdata<-as.data.frame(rbind(smoothdata,temp))
  }
}
names(smoothdata)<-c("ent_smooth","Order","TeamID")

#lead and lag the smoothed data
list<-unique(smoothdata$TeamID)
for (j in 1:length(list)){
  case<-list[j]
  subdata <- subset(smoothdata, smoothdata$TeamID==case)
  subdata$ent_smooth_ld1<-Lag(subdata$ent_smooth,-1)
  subdata$ent_smooth_lg1<-Lag(subdata$ent_smooth,1)
  if(j==1){
    smoothdata1<-subdata
  } else{
    smoothdata1<-rbind(smoothdata1,subdata)
  }
}
smoothdata1<-na.omit(smoothdata1)



# Peak identification for smoothed data----------------------------------------------------

#peak picking on smoothed data
for (i in 1:length(smoothdata1$ent_smooth)){
  if(smoothdata1$ent_smooth[i]>smoothdata1$ent_smooth_ld1[i] && smoothdata1$ent_smooth[i]>smoothdata1$ent_smooth_lg1[i]){
    smoothdata1$peakpoint[i]<-1
  } else {
    smoothdata1$peakpoint[i]<-0
  }
}

ss<-subset(smoothdata1,smoothdata1$peakpoint==1)

#example plot for a single team showing identified peak points for smoothed data
#The xlim needs to be manually adjusted here to the length of the time series, but it should be kept at the same value of plot p if a multiplot is desired.
Team1<-subset(smoothdata1,TeamID==1)
ss<-subset(Team1,Team1$peakpoint==1)
s<-ggplot(data=Team1, aes(x=Order,y=ent_smooth))+geom_line(linetype=1)+geom_point() + geom_vline(xintercept=ss$Order)+xlim(1,347)+labs(list(title="Smoothed Entropy Time Series", x= "Time Unit", y = "Entropy Value"))


# Peak total and proportion measures -----------------------------------------------------------

#Create new variables that take the total number of peak points and proportion of peaks to total length of time series
##For original entropy data frame with peak point variable added
list<-unique(data1$TeamID)
for (j in 1:length(list)){
  case<-list[j]
  subdata <- subset(data1, data1$TeamID==case)
  TotalPeaks<-sum(subdata$peakpoint)
  PeakProp<-sum(subdata$peakpoint)/length(subdata$entropy)
  TeamID<-j
  TsLength<-length(subdata$entropy)
  temp<-as.data.frame(cbind(TeamID,TsLength, TotalPeaks,PeakProp))
  if(j==1){
    PeakRes<-temp
  } else{
    PeakRes<-rbind(PeakRes,temp)
  }
}
write.csv(PeakRes, file="PeakResOrig.csv",row.names=FALSE)

#Create new variables that takes the total number of peak points and proportion of peaks to total length of time series
##For Smoothed data
list<-unique(smoothdata1$TeamID)
for (j in 1:length(list)){
  case<-list[j]
  subdata <- subset(smoothdata1, smoothdata1$TeamID==case)
  TotalPeaks<-sum(subdata$peakpoint)
  PeakProp<-sum(subdata$peakpoint)/length(subdata$ent_smooth)
  TeamID<-j
  TsLength<-length(subdata$ent_smooth)
  temp<-as.data.frame(cbind(TeamID,TsLength, TotalPeaks,PeakProp))
  if(j==1){
    PeakResSm<-temp
  } else{
    PeakResSm<-rbind(PeakResSm,temp)
  }
}

write.csv(PeakResSm, file="PeakResSmMA5.csv",row.names=FALSE)



# Smoothing Measures ------------------------------------------------------
#This code uses a combined data (perfdata) file that utilizes the peak proportion data from the moving average window sizes 2-15
#calculate confidence intervals on means of each moving average window for peak proportions
#must not have any missing values
perfdata<--read.csv(file= "PeakResAll.csv",header=TRUE, sep=",")
perfdata<-abs(perfdata)

se_orig<-std.error(perfdata$OrigPeakProp,na.rm=TRUE) 
mn_orig<-mean(perfdata$OrigPeakProp,na.rm=TRUE)
orig_ci<-se_orig*1.96
orig_ci_up<-mn_orig+orig_ci
orig_ci_low<-mn_orig-orig_ci

se_MA2<-std.error(perfdata$MA2PeakProp,na.rm=TRUE) 
mn_MA2<-mean(perfdata$MA2PeakProp,na.rm=TRUE)
MA2_ci<-se_MA2*1.96
MA2_ci_up<-mn_MA2+MA2_ci
MA2_ci_low<-mn_MA2-MA2_ci

se_MA3<-std.error(perfdata$MA3PeakProp,na.rm=TRUE) 
mn_MA3<-mean(perfdata$MA3PeakProp,na.rm=TRUE)
MA3_ci<-se_MA3*1.96
MA3_ci_up<-mn_MA3+MA3_ci
MA3_ci_low<-mn_MA3-MA3_ci

se_MA3<-std.error(perfdata$MA3PeakProp,na.rm=TRUE) 
mn_MA3<-mean(perfdata$MA3PeakProp,na.rm=TRUE)
MA3_ci<-se_MA3*1.96
MA3_ci_up<-mn_MA3+MA3_ci
MA3_ci_low<-mn_MA3-MA3_ci

se_MA4<-std.error(perfdata$MA4PeakProp,na.rm=TRUE)
mn_MA4<-mean(perfdata$MA4PeakProp,na.rm=TRUE)
MA4_ci<-se_MA4*1.96
MA4_ci_up<-mn_MA4+MA4_ci
MA4_ci_low<-mn_MA4-MA4_ci

se_MA5<-std.error(perfdata$MA5PeakProp,na.rm=TRUE) 
mn_MA5<-mean(perfdata$MA5PeakProp,na.rm=TRUE)
MA5_ci<-se_MA5*1.96
MA5_ci_up<-mn_MA5+MA5_ci
MA5_ci_low<-mn_MA5-MA5_ci

se_MA6<-std.error(perfdata$MA6PeakProp,na.rm=TRUE)
mn_MA6<-mean(perfdata$MA6PeakProp,na.rm=TRUE)
MA6_ci<-se_MA6*1.96
MA6_ci_up<-mn_MA6+MA6_ci
MA6_ci_low<-mn_MA6-MA6_ci

se_MA7<-std.error(perfdata$MA7PeakProp,na.rm=TRUE)
mn_MA7<-mean(perfdata$MA7PeakProp,na.rm=TRUE)
MA7_ci<-se_MA7*1.96
MA7_ci_up<-mn_MA7+MA7_ci
MA7_ci_low<-mn_MA7-MA7_ci

se_MA8<-std.error(perfdata$MA7PeakProp,na.rm=TRUE)
mn_MA8<-mean(perfdata$MA8PeakProp,na.rm=TRUE)
MA8_ci<-se_MA8*1.96
MA8_ci_up<-mn_MA8+MA8_ci
MA8_ci_low<-mn_MA8-MA8_ci

se_MA9<-std.error(perfdata$MA9PeakProp,na.rm=TRUE)
mn_MA9<-mean(perfdata$MA9PeakProp,na.rm=TRUE)
MA9_ci<-se_MA9*1.96
MA9_ci_up<-mn_MA9+MA9_ci
MA9_ci_low<-mn_MA9-MA9_ci

se_MA10<-std.error(perfdata$MA10PeakProp,na.rm=TRUE)
mn_MA10<-mean(perfdata$MA10PeakProp,na.rm=TRUE)
MA10_ci<-se_MA10*1.96
MA10_ci_up<-mn_MA10+MA10_ci
MA10_ci_low<-mn_MA10-MA10_ci

se_MA11<-std.error(perfdata$MA11PeakProp,na.rm=TRUE) 
mn_MA11<-mean(perfdata$MA11PeakProp,na.rm=TRUE)
MA11_ci<-se_MA11*1.96
MA11_ci_up<-mn_MA11+MA11_ci
MA11_ci_low<-mn_MA11-MA11_ci

se_MA12<-std.error(perfdata$MA12PeakProp,na.rm=TRUE)
mn_MA12<-mean(perfdata$MA12PeakProp,na.rm=TRUE)
MA12_ci<-se_MA12*1.96
MA12_ci_up<-mn_MA12+MA12_ci
MA12_ci_low<-mn_MA12-MA12_ci

se_MA13<-std.error(perfdata$MA13PeakProp,na.rm=TRUE) 
mn_MA13<-mean(perfdata$MA13PeakProp,na.rm=TRUE)
MA13_ci<-se_MA13*1.96
MA13_ci_up<-mn_MA13+MA13_ci
MA13_ci_low<-mn_MA13-MA13_ci

se_MA14<-std.error(perfdata$MA14PeakProp,na.rm=TRUE)
mn_MA14<-mean(perfdata$MA14PeakProp,na.rm=TRUE)
MA14_ci<-se_MA14*1.96
MA14_ci_up<-mn_MA14+MA14_ci
MA14_ci_low<-mn_MA14-MA14_ci

se_MA15<-std.error(perfdata$MA15PeakProp,na.rm=TRUE) 
mn_MA15<-mean(perfdata$MA15PeakProp,na.rm=TRUE)
MA15_ci<-se_MA15*1.96
MA15_ci_up<-mn_MA15+MA15_ci
MA15_ci_low<-mn_MA15-MA15_ci

#create lines for each variable
mean_peak_prop<-c(mn_orig,mn_MA2,mn_MA3,mn_MA4,mn_MA5,mn_MA6,mn_MA7,mn_MA8,mn_MA9, mn_MA10, mn_MA11,mn_MA12,mn_MA13,mn_MA14,mn_MA15)
up_ci_peak<-c(orig_ci_up,MA2_ci_up,MA3_ci_up,MA4_ci_up,MA5_ci_up,MA6_ci_up,MA7_ci_up,MA8_ci_up,MA9_ci_up,MA10_ci_up,MA11_ci_up,MA12_ci_up,MA13_ci_up,MA14_ci_up,MA15_ci_up)
low_ci_peak<-c(orig_ci_low,MA2_ci_low,MA3_ci_low,MA4_ci_low,MA5_ci_low,MA6_ci_low,MA7_ci_low,MA8_ci_low,MA9_ci_low,MA10_ci_low,MA11_ci_low,MA12_ci_low,MA13_ci_low,MA14_ci_low,MA15_ci_low)
ci<-c(orig_ci,MA2_ci,MA3_ci,MA4_ci,MA5_ci,MA6_ci,MA7_ci,MA8_ci,MA9_ci,MA10_ci,MA11_ci,MA12_ci,MA13_ci,MA14_ci,MA15_ci)

peak_df<-as.data.frame(cbind(mean_peak_prop,ci))

#Use the results of these plots to determine the appropriate moving average to use
#The window where the confidence intervals overlap is a justifiable choice, but this example has two few data points for a reliable CI
plot(mean_peak_prop,type="b",main="Change in Peak Proportions", xlab="Window Size",ylab="Peak Proportion")
lines(up_ci_peak,lty=2)
lines(low_ci_peak,lty=2)

pd <- position_dodge(0.1)
ggplot(peak_df, aes(x=c(1:15), y=mean_peak_prop)) + 
  geom_errorbar(aes(ymin=mean_peak_prop-ci, ymax=mean_peak_prop+ci), width=.1) + geom_line() + geom_point()+
  labs(x="Window Size", y="Peak Proportion", title="Mean Peak Proportion Across Moving Average Window Sizes with 95% CIs" )



# Creating epoch data for MLM ---------------------------------------------


#Take peak data and then create epoch ID using smooth data
list<-unique(smoothdata1$TeamID)
for (j in 1:length(list)){
  case<-list[j]
  subdata <- subset(smoothdata1, smoothdata1$TeamID==case)
  epoc<-subdata$peakpoint
  subdata$epoch_id<-cumsum(epoc)+1
  if(j==1){
    data2<-subdata
  } else{
    data2<-rbind(data2,subdata)
  }
}

#take original code time series, subset by Team ID and Order 2:length(Order-29)
#add in epoch variable from previous data frame
#also adds in separate variables reprsenting the count of each individual code
#Names need to be modified dependoning on the number of codes
list<-unique(data_c$TeamID)
for (i in 1:length(list)){
  case<-list[i]
  subdata <- subset(data_c, data_c$TeamID==case)
  len<-length(subdata$Order)-29 #This must be modified based on window size and selected moving average
  subdata1<-subdata[2:len,]
  if(i==1){
    data3<-subdata1
  } else{
    data3<-rbind(data3,subdata1)
  }
}
data3$EpochID<-data2$epoch_id
ind<-model.matrix(~factor(data3$Code)-1)

data4<-cbind(data3,ind)
names(data4)<-c("TeamID","Code","TeamIDvar","Order","EpochID","C1","C2","C3","C4","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15")
#No C5 in the current examples so they are ommited in the above names 

#This adds the weight variable to the data frame based on the length of the epoch
#The output of this part is what can be used in MLM to test for differences in code distributions between epochs
list1<-unique(data4$TeamIDvar)
for (j in 1:length(list1)){
  case1<-list1[j]
  subdata1<-subset(data4,data4$TeamIDvar==case1)
  list<-unique(subdata1$EpochID)
  for (i in 1:length(list)){
    case<-list[i]
    subdata<-subset(subdata1,subdata1$EpochID==case)
    temp<-nrow(subdata)
    subdata$weight<-rep(1/temp, time=length(subdata$EpochID))
    if(i==1){
      data5<-subdata
    } else{
      data7<-rbind(data5,subdata)
    }
  }
  if(j==1){
    data6<-data5
  } else{
    data6<-rbind(data6,data5)
  }
}

write.csv(data6, file="CodeTs_w_EpochId.csv",row.names=FALSE)




# Output Entropy Descriptives from Peak Points Only -----------------------


#code to take the average entropy and variability in entropy from peaks points only
#This uses the data frame that includes the smoothed entropy, peakpoint, and epoc_id
teamvar<-read.csv(file= "teamvariables.csv",header=TRUE, sep=",")
#This extra data file is added in so that different team-level variables could be added into the data file for further analyses

i=1
list<-unique(data2$TeamID)
for (i in 1:length(list)){
  case<-list[i]
  subdata <- subset(data2, data2$TeamID==case)
  peaksubdata<-subset(subdata, subdata$peakpoint==1)
  avePent<-mean(peaksubdata$ent_smooth) #These are some examples of entropy properties at the peak points that could be extracted
  varPent<-var(peaksubdata$ent_smooth)
  minPent<-min(peaksubdata$ent_smooth)
  maxPent<-max(peaksubdata$ent_smooth)
  perf<-teamvar$Performance_Rescale_NoZeros[i]#This variable and the next two would need to be modified to use different variables
  gendercomp<-teamvar$GenderComp[i]
  minKnowAss<-teamvar$minKnowAss[i]
  temp<-as.data.frame(cbind(i,avePent, varPent, minPent,maxPent,perf, gendercomp,minKnowAss))
  if(i==1){
    entropyres<-temp
  } else{
    entropyres<-rbind(entropyres,temp)
  }
}
write.csv(entropyres,file="entropyPeakres2.csv",row.names=FALSE,na="-999")


# Random Entropy Generation -----------------------------------------------
#All the analyses can be rerun on the randomly generated data
#The sampling values should be changed based on the number of codes (1:15 in this case)

for (i in 1:2){
randent<-as.data.frame(sample(1:15,size=perfdata$OrigTsLength[i],replace=TRUE))
randent$TeamIDvar<-rep(i, perfdata$OrigTsLength[i])
randent$TeamID<-rep(i, perfdata$OrigTsLength[i])
names(randent)<-c("Code", "TeamIDvar", "TeamID")
if (i==1){
  randEnt<-randent
} else {
  randEnt<-rbind(randEnt,randent)
}
}

data<-randEnt

# Plotting Code Distributions By Phases ----------------------------------------------------------------

#This part of the code uses the epoch segmented data frame
#For the MLM the code below subsets by each team and for each epoch within a team
#It then provide frequency counts for each of the communication codes
#If a different number of codes are used then those lines summing the codes would need to be modified


epoc_freq<-NULL
list<-unique(data3$TeamIDvar)
for (j in 1:length(list)){
  case<-list[j]
  subdata <- subset(data3, data3$TeamIDvar==case)
  list2<-unique(subdata$EpochID)
  for (i in 1:length(list2)) {
    case2<-list2[i]
    subdata2<-subset(subdata,subdata$EpochID==case2)
    TeamIDvar<-subdata$TeamIDvar[j]
    EpochID<-i
    C1<-sum(subdata2$Code==1)
    C2<-sum(subdata2$Code==2)
    C3<-sum(subdata2$Code==3)
    C4<-sum(subdata2$Code==4)
    C5<-sum(subdata2$Code==5)
    C6<-sum(subdata2$Code==6)
    C7<-sum(subdata2$Code==7)
    C8<-sum(subdata2$Code==8)
    C9<-sum(subdata2$Code==9)
    C10<-sum(subdata2$Code==10)
    C11<-sum(subdata2$Code==11)
    C12<-sum(subdata2$Code==12)
    C13<-sum(subdata2$Code==13)
    C14<-sum(subdata2$Code==14)
    C15<-sum(subdata2$Code==15)
    temp<-as.data.frame(cbind(TeamIDvar,EpochID, C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15))
    if (i==1){
      epoc.1<-temp
    } else {
      epoc.1<-rbind(epoc.1,temp)
    }
  }  
  if(j==1){
    epoc_freq<-epoc.1
  } else{
    epoc_freq<-rbind(epoc_freq,epoc.1)
  }
}

write.csv(epoc_freq, file="Code_Epoch_freq.csv",row.names=FALSE)

#Generate Phase Barplots for an Example Team
Team2<-subset(epoc_freq,epoc_freq$TeamIDvar==2)
Team2<-Team2[3:17]
par(mfrow=c(5, 3))


barplot(unlist(Team2[1,]),main="Phase 1",ylim = c(0,7),ylab="Code Freq")
barplot(unlist(Team2[2,]),main="Phase 2",ylim = c(0,7),ylab="Code Freq")
barplot(unlist(Team2[3,]),main="Phase 3",ylim = c(0,7),ylab="Code Freq")
barplot(unlist(Team2[4,]),main="Phase 4",ylim = c(0,7),ylab="Code Freq")
barplot(unlist(Team2[5,]),main="Phase 5",ylim = c(0,7),ylab="Code Freq")
barplot(unlist(Team2[6,]),main="Phase 6",ylim = c(0,7),ylab="Code Freq")
barplot(unlist(Team2[7,]),main="Phase 7",ylim = c(0,7),ylab="Code Freq")
barplot(unlist(Team2[8,]),main="Phase 8",ylim = c(0,7),ylab="Code Freq")
barplot(unlist(Team2[9,]),main="Phase 9",ylim = c(0,7),ylab="Code Freq")
barplot(unlist(Team2[10,]),main="Phase 10",ylim = c(0,7),ylab="Code Freq")
barplot(unlist(Team2[11,]),main="Phase 11",ylim = c(0,7),ylab="Code Freq")
barplot(unlist(Team2[12,]),main="Phase 12",ylim = c(0,7),ylab="Code Freq")
barplot(unlist(Team2[13,]),main="Phase 13",ylim = c(0,7),ylab="Code Freq")