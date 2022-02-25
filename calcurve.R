library(ggplot2) #Needed for pretty plots
library(dplyr) #Needed for %>%
library(outliers) #Needed for grubbs.test
library(sjPlot) #For pretty tables
source("L:/Google drive/R Scripts/integratePeaks.R")

setwd("L:/OneDrive/OneDrive - Washington State University (email.wsu.edu)/CE - SO2/Cider_SO2/Raw data/03042021 cal curve 3")

##Calibration curve level 1 - 3ppm SO2

three_1<-read.csv("3PPM-4.CSV",header=F,sep=",",fileEncoding="UTF-16")
three_2<-read.csv("3PPM-5.CSV",header=F,sep=",",fileEncoding="UTF-16")
three_3<-read.csv("3PPM-6.CSV",header=F,sep=",",fileEncoding="UTF-16")
three_4<-read.csv("3PPM-1.CSV",header=F,sep=",",fileEncoding="UTF-16")
three_5<-read.csv("3PPM-2.CSV",header=F,sep=",",fileEncoding="UTF-16")
three_6<-read.csv("3PPM-3.CSV",header=F,sep=",",fileEncoding="UTF-16")
three_7<-read.csv("3PPM-7.CSV",header=F,sep=",",fileEncoding="UTF-16")
three_8<-read.csv("3PPM-8.CSV",header=F,sep=",",fileEncoding="UTF-16")

lvl1<-data.frame(matrix(ncol=2))
lvl1[1:3,1]<-3

lvl1[1,2]<-integratePeaks(three_1,1,2,"gauss")
lvl1[2,2]<-integratePeaks(three_2,1,2,"gauss")
lvl1[3,2]<-integratePeaks(three_3,1,2,"gauss")

LOD1<-lvl1[1,2]
LOD2<-lvl1[2,2]
LOD3<-lvl1[3,2]
LOD4<-integratePeaks(three_4,2,3,"gauss")
LOD5<-integratePeaks(three_5,1,2,"gauss")
LOD6<-integratePeaks(three_6,1,2,"gauss")  
LOD7<-integratePeaks(three_7,1,2,"gauss")
LOD8<-integratePeaks(three_8,1,2,"gauss")

##Calibration curve level 2 - 5ppm SO2

five_1<-read.csv("5PPM-1.CSV",header=F,sep=",",fileEncoding="UTF-16")
five_2<-read.csv("5PPM-2.CSV",header=F,sep=",",fileEncoding="UTF-16")
five_3<-read.csv("5PPM-3.CSV",header=F,sep=",",fileEncoding="UTF-16")

lvl2<-data.frame(matrix(ncol=2))
lvl2[1:3,1]<-5

lvl2[1,2]<-integratePeaks(five_1,1,2,"gauss")
lvl2[2,2]<-integratePeaks(five_2,1,3,"gauss")
lvl2[3,2]<-integratePeaks(five_3,1,5,"gauss")

##Calibration curve level 3 - 10ppm SO2

ten_1<-read.csv("10PPM-1.CSV",header=F,sep="",fileEncoding="UTF-16")
ten_2<-read.csv("10PPM-2.CSV",header=F,sep=",",fileEncoding="UTF-16")
ten_3<-read.csv("10PPM-3.CSV",header=F,sep=",",fileEncoding="UTF-16")

lvl3<-data.frame(matrix(ncol=2))
lvl3[1:3,1]<-10

lvl3[1,2]<-integratePeaks(ten_1,1,2,"gauss")
lvl3[2,2]<-integratePeaks(ten_2,1,2,"gauss")
lvl3[3,2]<-integratePeaks(ten_3,1,2,"gauss")

##Calibration curve level 4 - 50ppm SO2

fif_1<-read.csv("50PPM-1.CSV",header=F,sep=",",fileEncoding="UTF-16")
fif_2<-read.csv("50PPM-2.CSV",header=F,sep=",",fileEncoding="UTF-16")
fif_3<-read.csv("50PPM-3.CSV",header=F,sep=",",fileEncoding="UTF-16")

lvl4<-data.frame(matrix(ncol=2))
lvl4[1:3,1]<-50

lvl4[1,2]<-integratePeaks(fif_1,1,2,"gauss")
lvl4[2,2]<-integratePeaks(fif_2,1,2,"gauss")
lvl4[3,2]<-integratePeaks(fif_3,1,2,"gauss")

##Calibration curve level 5 - 100ppm SO2

hun_1<-read.csv("100PPM-1.CSV",header=F,sep=",",fileEncoding="UTF-16")
hun_2<-read.csv("100PPM-2.CSV",header=F,sep=",",fileEncoding="UTF-16")
hun_3<-read.csv("100PPM-3.CSV",header=F,sep=",",fileEncoding="UTF-16")

lvl5<-data.frame(matrix(ncol=2))
lvl5[1:3,1]<-100

lvl5[1,2]<-integratePeaks(hun_1,1,2,"gauss")
lvl5[2,2]<-integratePeaks(hun_2,1,2,"gauss")
lvl5[3,2]<-integratePeaks(hun_3,1,2,"gauss")

#Create table of all integration methods

cal<-rbind(lvl1,lvl2,lvl3,lvl4,lvl5)

colnames(cal)<-c("Std conc (ppm)","RF")

#Create linear model
model<-lm(cal[,2]~cal[,1],data=cal)

#Outlier test on absolute value of residuals to determine if points can be removed
grubbs.test(abs(summary(model)$residual))

#Grubbs' test determined point 13 was an outlier (G=3.25216, P=0.0004, n=18)
cal<-cal[-10,]

#Final calibration curve
model<-lm(cal[,2]~cal[,1],data=cal)

#Plot fit
ggplot(cal,aes(x=cal[,1],y=cal[,2])) + geom_point() + geom_abline(intercept=summary(model)$coef[[1]],slope=summary(model)$coef[[2]],colour="blue",size=1) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="white")) + xlab("Free SO2 Concentration (ppm)") + ylab("Response Factor") + labs(title=paste("R2 = ",signif(summary(model)$r.squared,5),"Intercept = ",signif(summary(model)$coef[[1]],5),"Slope = ",signif(summary(model)$coef[[2]],5)))

write.csv(cal,"curve.csv",row.names=F)

#Equation for calculating concentration from response factor
b<-summary(model)$coeff[[1]]
m<-summary(model)$coeff[[2]]
curve<-function(x) (x-b)/m
summary(model)
#For calculating error
sy<-sigma(model)
n<-as.numeric(nrow(cal))
devsq<-sum((cal[,1]-mean(cal[,1]))^2)
ybar<-mean(cal[,2])
sx<-function(z,k) (sy/abs(m))*sqrt(((1/k)+(1/n)+((z-ybar)^2)/((m^2)*devsq)))
LOD<-sd(c(LOD1,LOD2,LOD3,LOD4,LOD5,LOD6,LOD7,LOD8))*3.3
LOQ<-sd(c(LOD1,LOD2,LOD3,LOD4,LOD5,LOD6,LOD7,LOD8))*10

#Percent recovery of 75ppm standard
LFB_1<-read.csv("LFB_75PPM-1.CSV",header=F,sep=",",fileEncoding="UTF-16")
LFB_2<-read.csv("LFB_75PPM-2.CSV",header=F,sep=",",fileEncoding="UTF-16")
LFB_3<-read.csv("LFB_75PPM-3.CSV",header=F,sep=",",fileEncoding="UTF-16")

lfb<-data.frame(matrix(ncol=2))
lfb[1:3,1]<-75

#Use same peak approximation as curve
lfb[1,2]<-integratePeaks(LFB_1,1,2,"gauss")
lfb[2,2]<-integratePeaks(LFB_2,1,2,"gauss")
lfb[3,2]<-integratePeaks(LFB_3,1,2,"gauss")

rec1<-curve(lfb[1,2])
rec2<-curve(lfb[2,2])
rec3<-curve(lfb[3,2])
rec<-mean(rec1,rec2,rec3)
(rec/75)*100

##-----------------------------------------------Sample Data--------------------------------------

####Day two samples: 03/05/2021####

setwd("L:/OneDrive/OneDrive - Washington State University (email.wsu.edu)/CE - SO2/Cider_SO2/Raw data/03052021 second samples")

nmns1_1<-read.csv("NMNS1-1.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns1_2<-read.csv("NMNS1-2.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns1_3<-read.csv("NMNS1-3.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns2_1<-read.csv("NMNS2-1.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns2_2<-read.csv("NMNS2-2.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns2_3<-read.csv("NMNS2-3.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns3_1<-read.csv("NMNS3-1.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns3_2<-read.csv("NMNS3-2.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns3_3<-read.csv("NMNS3-3.CSV",header=F,sep="",fileEncoding="UTF-16")

nms1_1<-read.csv("NMS1-1.CSV",header=F,sep="",fileEncoding="UTF-16")
nms1_2<-read.csv("NMS1-2.CSV",header=F,sep="",fileEncoding="UTF-16")
nms1_3<-read.csv("NMS1-3.CSV",header=F,sep="",fileEncoding="UTF-16")
nms2_1<-read.csv("NMS2-1.CSV",header=F,sep="",fileEncoding="UTF-16")
nms2_2<-read.csv("NMS2-2.CSV",header=F,sep="",fileEncoding="UTF-16")
nms2_3<-read.csv("NMS2-3.CSV",header=F,sep="",fileEncoding="UTF-16")
nms3_1<-read.csv("NMS3-1.CSV",header=F,sep="",fileEncoding="UTF-16")
nms3_2<-read.csv("NMS3-2.CSV",header=F,sep="",fileEncoding="UTF-16")
nms3_3<-read.csv("NMS3-3.CSV",header=F,sep="",fileEncoding="UTF-16")

data<-data.frame(matrix(ncol=3))

data[1,1]<-"NMNS1"
data[1,2]<-as.numeric(mean(c(curve(integratePeaks(nmns1_1,2,3,"gauss")),curve(integratePeaks(nmns1_2,2,3,"gauss")),curve(integratePeaks(nmns1_3,2,3,"gauss")))))
data[1,3]<-sx(data[1,2],3)

data[2,1]<-"NMNS2"
data[2,2]<-0
data[2,3]<-0
#data[2,2]<-as.numeric(mean(curve(integratePeaks(nmns2_1,1,2,"gauss")),curve(integratePeaks(nmns2_2,2,3,"gauss")),curve(integratePeaks(nmns2_3,2,3,"gauss"))))
#data[2,3]<-sx(data[2,2],3)

data[3,1]<-"NMNS3"
data[3,2]<-as.numeric(mean(c(curve(integratePeaks(nmns3_1,2,3,"gauss")),curve(integratePeaks(nmns3_2,2,3,"gauss")),curve(integratePeaks(nmns3_3,2,3,"gauss")))))
data[3,3]<-sx(data[3,2],3)

data[4,1]<-"NMS1"
data[4,2]<-mean(c(curve(integratePeaks(nms1_1,2,3,"gauss")),curve(integratePeaks(nms1_2,2,3,"gauss")),curve(integratePeaks(nms1_3,2,3,"gauss"))))
data[4,3]<-sx(data[4,2],3)

data[5,1]<-"NMS2"
data[5,2]<-mean(c(curve(integratePeaks(nms2_1,2,3,"gauss")),curve(integratePeaks(nms2_2,2,3,"gauss")),curve(integratePeaks(nms2_3,2,3,"gauss"))))
data[5,3]<-sx(data[5,2],3)

data[6,1]<-"NMS3"
data[6,2]<-mean(c(curve(integratePeaks(nms3_1,2,3,"gauss")),curve(integratePeaks(nms3_2,2,3,"gauss")),curve(integratePeaks(nms3_3,2,3,"gauss"))))
data[6,3]<-sx(data[6,2],3)

colnames(data)<-c("Sample","Free SO2 (ppm)","Standard Error (ppm)")

write.csv(data,"data.csv",row.names=F)

####Day three samples: 03/12/2021####

setwd("L:/OneDrive/OneDrive - Washington State University (email.wsu.edu)/CE - SO2/Cider_SO2/Raw data/03122021 third samples")

mns1_1<-read.csv("MNS1-1.CSV",header=F,sep="",fileEncoding="UTF-16")
mns1_2<-read.csv("MNS1-2.CSV",header=F,sep="",fileEncoding="UTF-16")
mns1_3<-read.csv("MNS1-3.CSV",header=F,sep="",fileEncoding="UTF-16")
mns2_1<-read.csv("MNS2-1.CSV",header=F,sep="",fileEncoding="UTF-16")
mns2_2<-read.csv("MNS2-2.CSV",header=F,sep="",fileEncoding="UTF-16")
mns2_3<-read.csv("MNS2-3.CSV",header=F,sep="",fileEncoding="UTF-16")
mns3_1<-read.csv("MNS3-1.CSV",header=F,sep="",fileEncoding="UTF-16")
mns3_2<-read.csv("MNS3-2.CSV",header=F,sep="",fileEncoding="UTF-16")
mns3_3<-read.csv("MNS3-3.CSV",header=F,sep="",fileEncoding="UTF-16")

ms1_1<-read.csv("MS1-1.CSV",header=F,sep="",fileEncoding="UTF-16")
ms1_2<-read.csv("MS1-2.CSV",header=F,sep="",fileEncoding="UTF-16")
ms1_3<-read.csv("MS1-3.CSV",header=F,sep="",fileEncoding="UTF-16")
ms2_1<-read.csv("MS2-1.CSV",header=F,sep="",fileEncoding="UTF-16")
ms2_2<-read.csv("MS2-2.CSV",header=F,sep="",fileEncoding="UTF-16")
ms2_3<-read.csv("MS2-3.CSV",header=F,sep="",fileEncoding="UTF-16")
ms3_1<-read.csv("MS3-1.CSV",header=F,sep="",fileEncoding="UTF-16")
ms3_2<-read.csv("MS3-2.CSV",header=F,sep="",fileEncoding="UTF-16")
ms3_3<-read.csv("MS3-3.CSV",header=F,sep="",fileEncoding="UTF-16")

data<-data.frame(matrix(ncol=3))

data[1,1]<-"MNS1"
data[1,2]<-as.numeric(mean(c(curve(integratePeaks(mns1_1,4,5,"gauss")),curve(integratePeaks(mns1_2,4,7,"gauss")),curve(integratePeaks(mns1_3,3,4,"gauss")))))
data[1,3]<-sx(data[1,2],3)

data[2,1]<-"MNS2"
data[2,2]<-as.numeric(mean(c(curve(integratePeaks(mns2_1,3,4,"gauss")),curve(integratePeaks(mns2_2,2,3,"gauss")),curve(integratePeaks(mns2_3,2,3,"gauss")))))
data[2,3]<-sx(data[2,2],3)

data[3,1]<-"MNS3"
data[3,2]<-as.numeric(mean(c(curve(integratePeaks(mns3_1,2,3,"gauss")),curve(integratePeaks(mns3_2,2,3,"gauss")),curve(integratePeaks(mns3_3,2,3,"gauss")))))
data[3,3]<-sx(data[3,2],3)

data[4,1]<-"MS1"
data[4,2]<-as.numeric(mean(c(curve(integratePeaks(ms1_1,2,3,"gauss")),curve(integratePeaks(ms1_2,2,3,"gauss")),curve(integratePeaks(ms1_3,2,3,"gauss")))))
data[4,3]<-sx(data[4,2],3)

data[5,1]<-"MS2"
data[5,2]<-as.numeric(mean(c(curve(integratePeaks(ms2_1,2,3,"gauss")),curve(integratePeaks(ms2_2,1,2,"gauss")),curve(integratePeaks(ms2_3,2,3,"gauss")))))
data[5,3]<-sx(data[5,2],3)

data[6,1]<-"MS3"
data[6,2]<-as.numeric(mean(c(curve(integratePeaks(ms3_1,2,3,"gauss")),curve(integratePeaks(ms3_2,1,2,"gauss")),curve(integratePeaks(ms3_3,2,3,"gauss")))))
data[6,3]<-sx(data[6,2],3)

colnames(data)<-c("Sample","Free SO2 (ppm)","Standard Error (ppm)")

write.csv(data,"data.csv",row.names=F)

####Day four samples: 05/24/2021####

setwd("L:/OneDrive/OneDrive - Washington State University (email.wsu.edu)/CE - SO2/Cider_SO2/Raw data/05242021 fourth samples")

nmns1_1<-read.csv("NMNS1-1.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns1_2<-read.csv("NMNS1-2.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns1_3<-read.csv("NMNS1-3.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns2_1<-read.csv("NMNS2-1.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns2_2<-read.csv("NMNS2-2.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns2_3<-read.csv("NMNS2-3.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns3_1<-read.csv("NMNS3-1.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns3_2<-read.csv("NMNS3-2.CSV",header=F,sep="",fileEncoding="UTF-16")
nmns3_3<-read.csv("NMNS3-3.CSV",header=F,sep="",fileEncoding="UTF-16")

nms1_1<-read.csv("NMS1-1.CSV",header=F,sep="",fileEncoding="UTF-16")
nms1_2<-read.csv("NMS1-2.CSV",header=F,sep="",fileEncoding="UTF-16")
nms1_3<-read.csv("NMS1-3.CSV",header=F,sep="",fileEncoding="UTF-16")
nms2_1<-read.csv("NMS2-1.CSV",header=F,sep="",fileEncoding="UTF-16")
nms2_2<-read.csv("NMS2-2.CSV",header=F,sep="",fileEncoding="UTF-16")
nms2_3<-read.csv("NMS2-3.CSV",header=F,sep="",fileEncoding="UTF-16")
nms3_1<-read.csv("NMS3-1.CSV",header=F,sep="",fileEncoding="UTF-16")
nms3_2<-read.csv("NMS3-2.CSV",header=F,sep="",fileEncoding="UTF-16")
nms3_3<-read.csv("NMS3-3.CSV",header=F,sep="",fileEncoding="UTF-16")

mns1_1<-read.csv("MNS1-1.CSV",header=F,sep="",fileEncoding="UTF-16")
mns1_2<-read.csv("MNS1-2.CSV",header=F,sep="",fileEncoding="UTF-16")
mns1_3<-read.csv("MNS1-3.CSV",header=F,sep="",fileEncoding="UTF-16")
mns2_1<-read.csv("MNS2-1.CSV",header=F,sep="",fileEncoding="UTF-16")
mns2_2<-read.csv("MNS2-2.CSV",header=F,sep="",fileEncoding="UTF-16")
mns2_3<-read.csv("MNS2-3.CSV",header=F,sep="",fileEncoding="UTF-16")
mns3_1<-read.csv("MNS3-1.CSV",header=F,sep="",fileEncoding="UTF-16")
mns3_2<-read.csv("MNS3-2.CSV",header=F,sep="",fileEncoding="UTF-16")
mns3_3<-read.csv("MNS3-3.CSV",header=F,sep="",fileEncoding="UTF-16")

ms1_1<-read.csv("MS1-1.CSV",header=F,sep="",fileEncoding="UTF-16")
ms1_2<-read.csv("MS1-2.CSV",header=F,sep="",fileEncoding="UTF-16")
ms1_3<-read.csv("MS1-3.CSV",header=F,sep="",fileEncoding="UTF-16")
ms2_1<-read.csv("MS2-1.CSV",header=F,sep="",fileEncoding="UTF-16")
ms2_2<-read.csv("MS2-2.CSV",header=F,sep="",fileEncoding="UTF-16")
ms2_3<-read.csv("MS2-3.CSV",header=F,sep="",fileEncoding="UTF-16")
ms3_1<-read.csv("MS3-1.CSV",header=F,sep="",fileEncoding="UTF-16")
ms3_2<-read.csv("MS3-2.CSV",header=F,sep="",fileEncoding="UTF-16")
ms3_3<-read.csv("MS3-3.CSV",header=F,sep="",fileEncoding="UTF-16")

LFB1<-read.csv("LFM50-1.CSV",header=F,sep="",fileEncoding="UTF-16")
LFB2<-read.csv("LFM50-2.CSV",header=F,sep="",fileEncoding="UTF-16")
LFB3<-read.csv("LFM50-3.CSV",header=F,sep="",fileEncoding="UTF-16")

data<-data.frame(matrix(ncol=3))
data[1,1]<-"NMNS1"
data[1,2]<-mean(c(curve(integratePeaks(nmns1_1,2,3,"gauss")),curve(integratePeaks(nmns1_2,2,3,"gauss")),curve(integratePeaks(nmns1_3,2,3,"gauss"))))
data[1,3]<-sx(data[1,2],3)

data[2,1]<-"NMNS2"
data[2,2]<-mean(c(curve(integratePeaks(nmns2_1,2,3,"gauss")),curve(integratePeaks(nmns2_2,2,3,"gauss")),curve(integratePeaks(nmns2_3,2,3,"gauss"))))
data[2,3]<-sx(data[2,2],3)

data[3,1]<-"NMNS3"
data[3,2]<-mean(c(curve(integratePeaks(nmns3_1,2,3,"gauss")),curve(integratePeaks(nmns3_2,2,3,"gauss")),curve(integratePeaks(nmns3_3,2,3,"gauss"))))
data[3,3]<-sx(data[3,2],3)

data[4,1]<-"NMS1"
data[4,2]<-mean(c(curve(integratePeaks(nms1_1,2,3,"gauss")),curve(integratePeaks(nms1_2,2,3,"gauss")),curve(integratePeaks(nms1_3,2,3,"gauss"))))
data[4,3]<-sx(data[4,2],3)

data[5,1]<-"NMS2"
data[5,2]<-mean(c(curve(integratePeaks(nms2_1,2,3,"gauss")),curve(integratePeaks(nms2_2,2,3,"gauss")),curve(integratePeaks(nms2_3,2,3,"gauss"))))
data[5,3]<-sx(data[5,2],3)

data[6,1]<-"NMS3"
data[6,2]<-mean(c(curve(integratePeaks(nms3_1,2,3,"gauss")),curve(integratePeaks(nms3_2,2,3,"gauss")),curve(integratePeaks(nms3_3,2,3,"gauss"))))
data[6,3]<-sx(data[6,2],3)

data[7,1]<-"MNS1"
data[7,2]<-mean(c(curve(integratePeaks(mns1_1,2,3,"gauss")),curve(integratePeaks(mns1_2,2,3,"gauss")),curve(integratePeaks(mns1_3,2,3,"gauss"))))
data[7,3]<-sx(data[7,2],3)

data[8,1]<-"MNS2"
data[8,2]<-mean(c(curve(integratePeaks(mns2_1,2,3,"gauss")),curve(integratePeaks(mns2_2,2,3,"gauss")),curve(integratePeaks(mns2_3,2,3,"gauss"))))
data[8,3]<-sx(data[8,2],3)

data[9,1]<-"MNS3"
data[9,2]<-mean(c(curve(integratePeaks(mns3_1,2,3,"gauss")),curve(integratePeaks(mns3_2,2,3,"gauss")),curve(integratePeaks(mns3_3,2,3,"gauss"))))
data[9,3]<-sx(data[9,2],3)

data[10,1]<-"MS1"
data[10,2]<-mean(c(curve(integratePeaks(ms1_1,2,3,"gauss")),curve(integratePeaks(ms1_2,2,3,"gauss")),curve(integratePeaks(ms1_3,2,3,"gauss"))))
data[10,3]<-sx(data[10,2],3)

data[11,1]<-"MS2"
data[11,2]<-mean(c(curve(integratePeaks(ms2_1,2,3,"gauss")),curve(integratePeaks(ms2_2,2,3,"gauss")),curve(integratePeaks(ms2_3,2,3,"gauss"))))
data[11,3]<-sx(data[11,2],3)

data[12,1]<-"MS3"
data[12,2]<-mean(c(curve(integratePeaks(ms3_1,2,3,"gauss")),curve(integratePeaks(ms3_2,2,3,"gauss")),curve(integratePeaks(ms3_3,2,3,"gauss"))))
data[12,3]<-sx(data[12,2],3)

data[13,1]<-"LFB"
data[13,2]<-mean(c(curve(integratePeaks(LFB1,1,2,"gauss")),curve(integratePeaks(LFB2,1,2,"gauss")),curve(integratePeaks(LFB3,1,2,"gauss"))))
data[13,3]<-sx(data[13,2],3)

colnames(data)<-c("Sample","Free SO2 (ppm)", "Standard Error (ppm)")

write.csv(data,"data.csv",row.names=F)

nmns1<-read.csv("nmns1.csv",header=T,sep=",")
nmns2<-read.csv("nmns2.csv",header=T,sep=",")
nmns3<-read.csv("nmns3.csv",header=T,sep=",")
nms1<-read.csv("nms1.csv",header=T,sep=",")
nms2<-read.csv("nms2.csv",header=T,sep=",")
nms3<-read.csv("nms3.csv",header=T,sep=",")
mns1<-read.csv("mns1.csv",header=T,sep=",")
mns2<-read.csv("mns2.csv",header=T,sep=",")
mns3<-read.csv("mns3.csv",header=T,sep=",")
ms1<-read.csv("ms1.csv",header=T,sep=",")
ms2<-read.csv("ms2.csv",header=T,sep=",")
ms3<-read.csv("ms3.csv",header=T,sep=",")

big<-read.csv("table.csv",header=T,sep=",")
avg<-read.csv("averaged.csv",header=T,sep=",")

summary(aov(SO2~Method,data=big))
TukeyHSD(aov(SO2~Method,data=big))

summary(aov(SO2~Method+MLB,data=big))
TukeyHSD(aov(SO2~Method+MLB,data=big))

summary(aov(SO2~Method,data=big[1:3,]))

##Red Cider samples

setwd("L:/OneDrive/OneDrive - Washington State University (email.wsu.edu)/CE - SO2/Cider_SO2/Raw data/11052021 red samples")

blkb1<-read.csv("BLACKBERRY-1.CSV",header=F,sep="",fileEncoding="UTF-16")
blkb2<-read.csv("BLACKBERRY-2.CSV",header=F,sep="",fileEncoding="UTF-16")
blkb3<-read.csv("BLACKBERRY-3.CSV",header=F,sep="",fileEncoding="UTF-16")

cherry1<-read.csv("CHERRY-1.CSV",header=F,sep="",fileEncoding="UTF-16")
cherry2<-read.csv("CHERRY-2.CSV",header=F,sep="",fileEncoding="UTF-16")
cherry3<-read.csv("CHERRY-3.CSV",header=F,sep="",fileEncoding="UTF-16")

cran1<-read.csv("CRAN-1.CSV",header=F,sep="",fileEncoding="UTF-16")
cran2<-read.csv("CRAN-2.CSV",header=F,sep="",fileEncoding="UTF-16")
cran3<-read.csv("CRAN-3.CSV",header=F,sep="",fileEncoding="UTF-16")

model1<-read.csv("MODEL-1.CSV",header=F,sep="",fileEncoding="UTF-16")
model2<-read.csv("MODEL-2.CSV",header=F,sep="",fileEncoding="UTF-16")
model3<-read.csv("MODEL-3.CSV",header=F,sep="",fileEncoding="UTF-16")

rasp1<-read.csv("RASP-1.CSV",header=F,sep="",fileEncoding="UTF-16")
rasp2<-read.csv("RASP-2.CSV",header=F,sep="",fileEncoding="UTF-16")
rasp3<-read.csv("RASP-3.CSV",header=F,sep="",fileEncoding="UTF-16")

red1<-read.csv("RED-1.CSV",header=F,sep="",fileEncoding="UTF-16")
red2<-read.csv("RED-2.CSV",header=F,sep="",fileEncoding="UTF-16")
red3<-read.csv("RED-3.CSV",header=F,sep="",fileEncoding="UTF-16")

lfb1<-read.csv("LFM50-1.csv",header=F,sep="",fileEncoding="UTF-16")
lfb2<-read.csv("LFM50-2.csv",header=F,sep="",fileEncoding="UTF-16")
lfb3<-read.csv("LFM50-3.csv",header=F,sep="",fileEncoding="UTF-16")

lfb_1<-curve(integratePeaks(lfb1,1,2,"gauss"))
lfb_2<-curve(integratePeaks(lfb2,1,2,"gauss"))
lfb_3<-curve(integratePeaks(lfb3,1,2,"gauss"))
lfb/50
data<-data.frame(matrix(ncol=3))

data[1,1]<-"Blackberry"
data[1,2]<-mean(c(curve(integratePeaks(blkb1,2,5,"gauss")),curve(integratePeaks(blkb2,2,3,"gauss")),curve(integratePeaks(blkb3,3,6,"gauss"))))
data[1,3]<-sx(data[1,2],3)

data[2,1]<-"Cherry"
data[2,2]<-0
data[2,3]<-0

data[3,1]<-"Cranberry"
data[3,2]<-mean(c(curve(integratePeaks(cran1,2,4,"gauss")),curve(integratePeaks(cran2,2,4,"gauss")),curve(integratePeaks(cran3,8,11,"gauss"))))
data[3,3]<-sx(data[3,2],3)

data[4,1]<-"Model Cider"
data[4,2]<-mean(c(curve(integratePeaks(model1,1,2,"gauss")),curve(integratePeaks(model2,1,2,"gauss")),curve(integratePeaks(model3,1,2,"gauss"))))
data[4,3]<-sx(data[4,2],3)

data[5,1]<-"Raspberry"
data[5,2]<-mean(c(curve(integratePeaks(rasp1,2,3,"gauss")),curve(integratePeaks(rasp2,2,3,"gauss")),curve(integratePeaks(rasp3,3,4,"gauss"))))
data[5,3]<-sx(data[5,2],3)

data[6,1]<-"Red Flesh"
data[6,2]<-0
data[6,3]<-0

colnames(data)<-c("Sample","Free SO2 (ppm)","Standard Error (ppm)")

write.csv(data,"data.csv",row.names=F)

red<-read.csv("red.csv",header=T,sep=",")
summary(aov(SO2~Method,data=red))
TukeyHSD(aov(SO2~Method,data=red))

summary(aov(SO2~Method*Cider,data=red))
TukeyHSD(aov(SO2~Method*Cider,data=red))

#-----------------------------------------TOP/BOTTOM PLOTS-----------------------------------------------
install.packages("ggpubr")

library(ggplot2)
library(ggpubr)

setwd("L:/Google drive/School/graduate research/Dissertation/Free and total so2 chapter")
fif_2<-fif_2<-read.csv("50PPM-2.CSV",header=F,sep=",",fileEncoding="UTF-16")
nmns3_3<-read.csv("NMNS3-3.CSV",header=F,sep="",fileEncoding="UTF-16")

ggplot(x=V1,y=V2,data=fif_2) + geom_line(aes(x=V1,y=V2)) + xlab("Migration Time (min)") + ylab("Absorbance @ 192nm (uAU)") + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_rect(fill="white"),panel.border=element_rect(color="black",fill=NA)) + geom_label(x=1,y=7,label="IS (1.1min)") + geom_label(x=1.6,y=16.5,label="SO2 (1.7min)") + ylim(-2,18)
ggplot(x=V1,y=V2,data=nmns3_3) + geom_line(aes(x=V1,y=V2)) + xlab("Migration Time (min)") + ylab("Absorbance @ 192nm (uAU)") + theme(axis.text=element_text(size=10),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_rect(fill="white"),panel.border=element_rect(color="black",fill=NA)) + geom_label(x=1.15,y=6,label="IS") + geom_label(x=1.9,y=10,label="SO2") + ylim(-2,15)


model<-ggplot(x=V1,y=V2,data=fif_2) + geom_line(aes(x=V1,y=V2)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_rect(fill="white"),panel.border=element_rect(color="black",fill=NA)) + ylim(-2,16.5) + xlab(element_blank()) + ylab(element_blank())
cider<-ggplot(x=V1,y=V2,data=nmns3_3) + geom_line(aes(x=V1,y=V2)) + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_rect(fill="white"),panel.border=element_rect(color="black",fill=NA)) + ylim(-2,16.5) + xlab(element_blank()) + ylab(element_blank())

