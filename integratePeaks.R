library(OrgMassSpecR) #Needed for DrawChromatogram() function
library(RpeakChrom) #Needed for findpeaks() function
library(DescTools) #Needed for Area Under the Curve (AUC) function
library(ggplot2) #Needed for pretty plots
library(dplyr) #Needed for %>%

FWHM<-function(df,widthOnly)
{
  #Takes dataframe (df), interpolates points using linear approximation to determine full width at half maximum height
  #Returns either desired width or a dataframe containing the top half of the peak depending on boolean widthOnly input
  
  colnames(df)<-c("x","y")
  rownames(df)<-1:nrow(df)
  
  peakmax<-max(df[,2]) #Find the maximum point of the peak
  peak50<-df[df$y>peakmax/2,] #Keep only points from the peak where y is greater than half peak max
  
  #Need to interpolate points that don't exist on curve to create midway points
  
  i2<-as.numeric(rownames(peak50[1,]))
  i1<-i2-1
  
  j1<-as.numeric(rownames(peak50[nrow(peak50),]))
  j2<-j1+1
  
  xi1<-df[i1,1]
  xi2<-df[i2,1]
  yi1<-df[i1,2]
  yi2<-df[i2,2]
  slope1<-((yi2-yi1)/(xi2-xi1))
  b1<-yi2-slope1*xi2
  
  midy1<-peakmax/2
  midx1<-(midy1-b1)/slope1
  
  xj1<-df[j1,1]
  xj2<-df[j2,1]
  yj1<-df[j1,2]
  yj2<-df[j2,2]
  slope2<-((yj2-yj1)/(xj2-xj1))
  b2<-yj2-slope2*xj2
  
  midy2<-peakmax/2
  midx2<-(midy2-b2)/slope2
  
  mu<-df$x[df$y==peakmax]
  
  width50 <- midx2-midx1
  
  
  #plot(df) #Can be included for visual of where points get added
  #points(midx1,midy1,col="red") #Provides visual of where point is added
  left<-data.frame()
  left<-c(midx1,midy1)
  
  
  #points(midx2,midy2,col="red") #Visual of second point
  right<-data.frame()
  right<-c(midx2,midy2)
  
  results<-data.frame()
  results<-rbind(left,right)
  
  if (widthOnly==T)
    return(width50)
  else
    return(results)
}

integratePeaks<-function(df,peak1=1,peak2=2,method="gauss")
{
  
  #Integrates selected peaks at the width at 1/2 peak height
  #df = dataframe of electropherogram
  #peak1 = integer, number of internal standard peak
  #peak2 = integer, number of analyte peak
  #method = "half" for integral at half max peak height
  #         "gauss" for gaussian peak approximation
  #         "triangle" for triangle peak approximation
  #         "full" for full peak integration
  
  require(OrgMassSpecR)
  require(RpeakChrom)
  require(DescTools)
  
  colnames(df)<-c("MT","mAU") #Create dataframe header
  
  df$mAU[df$mAU<0] <- 0 #Set baseline = 0
  
  peaks<-findpeaks(df$mAU,minpeakheight=0) # Create table of peaks of minimum height of 0, includes RTmax as well as RT ranges
  
  ISTDstart<-as.numeric(peaks[peak1,3]) #Datapoint at beginning of ISTD peak
  ISTDend<-as.numeric(peaks[peak1,4]) #Datapoint at end of ISTD peak
  ISTDmax<-as.numeric(peaks[peak1,1]) #ISTD peak max intensity value
  ISTD<-df[ISTDstart:ISTDend,] #Isolate internal standard peak
  
  ANstart<-as.numeric(peaks[peak2,3]) #Datapoint at beginning of analyte peak
  ANend<-as.numeric(peaks[peak2,4]) #Datapoint at end of analyte peak
  ANmax<-as.numeric(peaks[peak2,1]) #Analyte peak max intesntiy value
  AN<-df[ANstart:ANend,] #Isolate analyte peak
  
  #Visual of electropherogram with ISTD peak in RED and Analyte peak in BLUE to ensure proper peak selection
  DrawChromatogram(df$MT,df$mAU,
                   range=(list(start=c(ISTD[1,1],AN[1,1]),stop=c(ISTD[nrow(ISTD),1],AN[nrow(AN),1]))),
                   color=c("red","blue"),
                   xlab="Migration Time (min)",ylab="Absorbance (mAU)") 
  if(method=="full")
  {
    ISTDarea<-AUC(ISTD[,1],ISTD[,2],method=c("trapezoid")) #Find area under entire ISTD peak
    
    ANarea<-AUC(AN[,1],AN[,2],method=c("trapezoid")) #Find area under entire analyte peak
    
    RF<-ANarea/ISTDarea #Calculate response factor
    
    return(RF)
  }
  
  if(method=="half")
  {
    ISTDwidth<-FWHM(ISTD,widthOnly=F)
    ANwidth<-FWHM(AN,widthOnly=F)
    
    ISTD50<-ISTD[ISTD$mAU>ISTDmax/2,] #ISTD peak at 50% max peak height
    ISTD50<-rbind(ISTDwidth[1,],ISTD50,ISTDwidth[2,]) #Add left and right points from FWMH
    
    AN50<-AN[AN$mAU>ANmax/2,] #Analyte peak at 50% max peak height
    AN50<-rbind(ANwidth[1,],AN50,ANwidth[2,])
    
    ISTDarea<-AUC(ISTD50[,1],ISTD50[,2],method=c("trapezoid")) #Integrate ISTD50
    ANarea<-AUC(AN50[,1],AN50[,2],method=c("trapezoid")) #Integrate Analyte50
    RF<-ANarea/ISTDarea #Final response factor
    
    return(RF)
  }
  
  if(method=="gauss")
  {
    ISTDmu<-ISTD$MT[ISTD$mAU==ISTDmax] #Find mu for gaussian approximation (RT at which ISTDmax occurs)
    ISTDwidth<-FWHM(ISTD,widthOnly=T) #Use FWHM with widthOnly=True to find full width at half max height
    ISTDsig<-ISTDwidth/(2*sqrt(2*log(2))) #Find stdev for gaussian approximation using FWHM/2*sqrt(2*ln(2))
    
    ANmu<-AN$MT[AN$mAU==ANmax] #Analyte mu
    ANwidth<-FWHM(AN,widthOnly=T) #Analyte FWHM
    ANsig<-ANwidth/(2*sqrt(2*log(2))) #Analyte sigma
    
    ISTDgaus<-function(x) ISTDmax*exp(-0.5*((x-ISTDmu)/ISTDsig)**2) #Gaussian approximation of ISTD peak
    ISTDarea<-AUC(ISTD[,1],ISTDgaus(ISTD[,1]),method=c("trapezoid")) #Find area under curve
    
    ANgaus<-function(x) ANmax*exp(-0.5*((x-ANmu)/ANsig)**2) #Gaussian approximation of analyte peak
    ANarea<-AUC(AN[,1],ANgaus(AN[,1]),method=c("trapezoid")) #Find area under the curve
    
    RF<-ANarea/ISTDarea #Final response factor
    
    return(RF)
  }
  
  if(method=="triangle")
  {
    ISTDwidth<-FWHM(ISTD,widthOnly=F)
    ANwidth<-FWHM(AN,widthOnly=F)
    
    ISTD50<-ISTD[ISTD$mAU>ISTDmax/2,] #ISTD peak at 50% max peak height
    ISTD50<-rbind(ISTDwidth[1,],ISTD50,ISTDwidth[2,]) #Add left and right points from FWMH
    rownames(ISTD50)<-1:nrow(ISTD50)
    
    AN50<-AN[AN$mAU>ANmax/2,] #Analyte peak at 50% max peak height
    AN50<-rbind(ANwidth[1,],AN50,ANwidth[2,])
    rownames(AN50)<-1:nrow(AN50)
    
    #Initial point for ISTD peak - point at max value
    ISTDmaxx<-ISTD50$MT[ISTD50$mAU==ISTDmax]
    ISTDmaxy<-ISTDmax
    
    #Second point - leftmost point of peak at 1/2 max height
    #Creates line between max point and this point
    ISTDx1<-ISTD50[1,1]
    ISTDy1<-ISTD50[1,2]
    ISTDm1<-(ISTDmaxy-ISTDy1)/(ISTDmaxx-ISTDx1)
    ISTDb1<-ISTDy1-ISTDm1*ISTDx1
    ISTDint1<-(-1*ISTDb1)/ISTDm1
    
    #Third point - Rightmost point of peak at 1/2 max height
    #creates second line from max point to this point
    ISTDx2<-ISTD50[nrow(ISTD50),1]
    ISTDy2<-ISTD50[nrow(ISTD50),2]
    ISTDm2<-(ISTDy2-ISTDmaxy)/(ISTDx2-ISTDmaxx)
    ISTDb2<-ISTDy2-ISTDm2*ISTDx2
    ISTDint2<-(-1*ISTDb2)/ISTDm2
    
    #Finds area of resulting triangle = (1/2)*base*height
    ISTDarea<-0.5*(ISTDint2-ISTDint1)*ISTDmaxy
    
    #Repeat process for analyte peak
    AN50<-AN[AN$mAU>ANmax/2,]
    AN50<-rbind(ANwidth[1,],AN,ANwidth[2,])
    
    ANmaxx<-AN50$MT[AN50$mAU==ANmax]
    ANmaxy<-ANmax
    
    ANx1<-AN50[1,1]
    ANy1<-AN50[1,2]
    ANm1<-(ANmaxy-ANy1)/(ANmaxx-ANx1)
    ANb1<-ANy1-ANm1*ANx1
    ANint1<-(-1*ANb1)/ANm1
    
    ANx2<-AN50[nrow(AN50),1]
    ANy2<-AN50[nrow(AN50),2]
    ANm2<-(ANy2-ANmaxy)/(ANx2-ANmaxx)
    ANb2<-ANy2-ANm2*ANx2
    ANint2<-(-1*ANb2)/ANm2
    
    ANarea<-0.5*(ANint2-ANint1)*ANmaxy
    
    RF<-ANarea/ISTDarea #Final response factor
    
    return(RF)
  }
  
  else
    stop("Error")
}