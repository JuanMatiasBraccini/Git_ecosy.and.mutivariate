
#SCRIPT FOR ANALYSING ECOSYSTEM PATTERNS IN THE SHARK GILLNET FISHERIES OF AUSTRALIA

#INDEX
# ----- DATA SECTION
#         1. Bring in WA shark observer data
#         2. Bring in MAFFRI shark gillnet data
#         3. Bring in WA Species names
#         4. Bring in Length coefficients
#         5. Bring in Commercial data

# ----- PROCEDURE SECTION
#         1. Manipulated WA shark observer data
#         2. Ecosystems indicators analysis
#           2.1 WA Fisheries observer data
#           2.2 WA Fisheries commercial data
#         3. Multivariate analysis

#missing: Commercial TDGDLF. Keep Monthly and Daily separate!!!
#         Do Sections 3 to 6 and 7.1
#         For PRIMER notes and recommendations, look in C:\Matias\Workshops\2019_Primer workshop\2019_PRIMER_workshop.docx

rm(list=ls(all=TRUE))

library(maps)
library(mapdata)
library(RColorBrewer)
library(vegan)
#library(gamlss)
library(reshape2)
#library(iNEXT) 
library(ggplot2)
library(caret)      #for all things data mining
library(plotrix)
library(ReporteRs)
library(mvtnorm)      #for multivariate normal pdf
library(lme4) #mixed models
library(MuMIn)  #model selection and pseudoR2 mixed effect models
library(car)    #to get ANOVA from lmer model
library(data.table)
library(tidyverse)

#Define working directory
handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')

setwd(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch"))
User="Matias"

#Source functions
source(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Git_ecosy.and.mutivariate/Ecosystem_functions.R"))


#choose if doing .jpeg or .tiff figures
Do.jpeg="YES"
Do.tiff="NO"
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R"))


# 1 Data section-----------------------------------------------------------------------

# 1. Bring in WA shark observer data
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R"))
rm(DATA)
DATA=DATA.ecosystems


# 2. Bring in WA Species names + PCS + FATE
setwd(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch"))
SPECIES_PCS_FATE=read.csv("SPECIES+PCS+FATE.csv",stringsAsFactors=F)
SPECIES_PCS_FATE_Com=read.csv("SPECIES+PCS+FATE_Com.csv",stringsAsFactors=F)

# 3. Bring in Length coefficients
Len.cof=read.csv("Raw coefficents table.csv")


#4. Bring in Commercial data
source(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Git_ecosy.and.mutivariate/Commercial_data_for_Ecosystem_Analysis.R"))



#CONTROL
  #Selection of recrods
Min.shts=10 #USe records with at least 10 shots per year-block
Min.recs=5 #minimum number of records
Min.individuals=5   #minimum number of individuals per shot to use

  #vessel used as mixed effect
MixedEff="BOAT"   

  #choose if doing data exploration
do.exploratory="NO"  

niter=100    #number of iterations for MC procedure for confidence intervals

#Percent.show.bycatch=0.8   #proportion of bycatch explained
#Prop.TC.bycatch=0.05       #proportion of catch explained

use.soak.time="YES"  #define if using soak time in the calculation of effort


ResVar="INDIVIDUALS"
MultiVar="SPECIES"
IDVAR=c("SHEET_NO","YEAR","MONTH","BOTDEPTH","BLOCK","ZONE","BOAT","SKIPPER","MESH_SIZE","BIOREGION","SEASON","EFFORT")   
Predictors=c("YEAR","BLOCK","BOAT","MONTH","BOTDEPTH")
Expl.varS=c("YEAR","BOAT","MONTH","BLOCK","BOTDEPTH")
FactoRS=Expl.varS[-match(c("BOTDEPTH"),Expl.varS)]
OFFSETT=NA

# 2 Functions-----------------------------------------------------------------------


#Deviance explained
Dsquared <- function(model, adjust = FALSE)
{
  # version 1.1 (13 Aug 2013)
  # calculates the explained deviance of a GLM
  # model: a model object of class "glm"
  # adjust: logical, whether or not to use the adjusted deviance taking into acount 
  #        the number of observations and parameters (Weisberg 1980; Guisan & Zimmermann 2000)
  d2 <- (model$null.deviance - model$deviance) / model$null.deviance
  if (adjust) {
    n <- length(model$fitted.values)
    p <- length(model$coefficients)
    d2 <- 1 - ((n - 1) / (n - p)) * (1 - d2)
  }
  return(d2)
}  

#function for exporting table
Export.tbl=function(WD,Tbl,Doc.nm,caption,paragph,HdR.col,HdR.bg,Hdr.fnt.sze,Hdr.bld,
                    body.fnt.sze,Zebra,Zebra.col,Grid.col,Fnt.hdr,Fnt.body,
                    HDR.names,HDR.span,HDR.2nd)
{
  mydoc = docx(Doc.nm)  #create r object
  mydoc = addSection( mydoc, landscape = T )   #landscape table
  # add title
  if(!is.na(caption))mydoc = addParagraph(mydoc, caption, stylename = "TitleDoc" )
  
  # add a paragraph
  if(!is.na(paragph))mydoc = addParagraph(mydoc , paragph, stylename="Citationintense")
  
  #add table
  MyFTable=FlexTable(Tbl,header.column=F,add.rownames =F,
                     header.cell.props = cellProperties(background.color=HdR.bg), 
                     header.text.props = textProperties(color=HdR.col,font.size=Hdr.fnt.sze,
                                                        font.weight="bold",font.family =Fnt.hdr), 
                     body.text.props = textProperties(font.size=body.fnt.sze,font.family =Fnt.body))
  
  #Add header
  MyFTable = addHeaderRow(MyFTable,text.properties=textBold(),value=HDR.names,colspan=HDR.span)
  
  #Add second header
  MyFTable = addHeaderRow(MyFTable, text.properties = textBold(),value =HDR.2nd)
  
  
  # zebra stripes - alternate colored backgrounds on table rows
  if(Zebra=="YES") MyFTable = setZebraStyle(MyFTable, odd = Zebra.col, even = "white" )
  
  # table borders
  MyFTable = setFlexTableBorders(MyFTable,
                                 inner.vertical = borderNone(),inner.horizontal = borderNone(),
                                 outer.vertical = borderNone(),
                                 outer.horizontal = borderProperties(color=Grid.col, style="solid", width=4))
  
  # set columns widths (in inches)
  #MyFTable = setFlexTableWidths( MyFTable, widths = Col.width)
  
  mydoc = addFlexTable( mydoc, MyFTable)   
  mydoc = addSection( mydoc, landscape = F ) 
  
  # write the doc 
  writeDoc( mydoc, file = paste(Doc.nm,".docx",sep=''))
}


#handy function for exporting figures
fn.fig=function(NAME,Width,Height)
{
  if(Do.tiff=="YES") tiff(file=paste(NAME,".tiff",sep=""),width=Width,height=Height,units="px",res=300,compression="lzw")
  if(Do.jpeg=="YES") jpeg(file=paste(NAME,".jpeg",sep=""),width=Width,height=Height,units="px",res=300)
}

fn.bx.plt=function(d,Var)  boxplot(d[,match(Var,names(d))]~d$YEAR,main=Var)
density.response=function(d,Var,Log) 
{
  d=d[,match(Var,names(d))]
  if(Log=="YES")d=log(d[!is.na(d)])
  a=round(range(d),2)
  plot(density(d,adjust=2),main=paste("log=",Log,"  ",Var," (range=",a[1],"-",a[2],")",sep=""),
       ylab="",xlab="",cex.main=.8)
}


#Look at predictors distribution
jplot=function(var)
{
  hs=hist(eval(parse(text=var)),main="",xlab=var,ylab="Frequency",
          col="grey90",ylim=c(0,90),border=T)
  dens=density(eval(parse(text=var)),na.rm=T)
  rs=max(hs$counts/max(dens$y))
  lines(dens$x,dens$y*rs,type='l',col=2)
  rug(eval(parse(text=var)))  #add a rug to plot
  title(main=paste("Distribution of",var))
}

#Model functions
Mst.cmn=function(D)
{
  D=table(D)
  return(names(sort(D)[length(D)]))
}

#glm approach
Mod.fn.glm=function(d,ResVar,Expl.vars,Predictrs,FactoRs,OFFSET,log.var,add.inter,MixedEff)
{
  d=d[,match(c(ResVar,Expl.vars),names(d))]
  d=d[!is.na(d[,match(ResVar,names(d))]),]
  d=d[d[,match(ResVar,names(d))]>0,]
  if(!is.na(OFFSET))d=subset(d,!is.na(EFFORT))
  if(log.var=="YES" & length(d$BOTDEPTH)>0)
  {
    d=subset(d,!is.na(BOTDEPTH))
    d=subset(d,BOTDEPTH>=3)
    d$log.BOTDEPTH=log(d$BOTDEPTH)
  }
  Y=d[,match(ResVar,names(d))]
  if(!is.na(OFFSET))d$log.EFFORT=log(d$EFFORT)
  
  #set factors
  for(x in 1:length(FactoRs)) d[,match(FactoRs[x],names(d))]=as.factor(d[,match(FactoRs[x],names(d))])
  
  #build formula
  if(add.inter=="YES")
  {
    IDS=match(c("YEAR","ZONE"),Predictrs)
    Pre.yr.zn=Predictrs[IDS]
    Pred.form=paste(Pre.yr.zn,collapse="*")
    if(!is.na(OFFSET))Pred.form=paste(c(Pred.form,Predictrs[-IDS],OFFSET),collapse="+")
    if(is.na(OFFSET))Pred.form=paste(c(Pred.form,Predictrs[-IDS]),collapse="+")
  }else
  {
    if(!is.na(OFFSET))Pred.form=paste(c(Predictrs,OFFSET),collapse="+")
    if(is.na(OFFSET))Pred.form=paste(c(Predictrs),collapse="+")
  }
  
  if(!is.na(MixedEff)) Pred.form=paste(paste(c(Predictrs[-match(MixedEff,Predictrs)]),collapse="+"),
                                       "+(1|",MixedEff,")",sep="")
  if(log.var=="YES")Formula=as.formula(paste("log(",ResVar,")", "~", Pred.form,collapse=NULL))
  if(log.var=="NO")Formula=as.formula(paste(ResVar, "~", Pred.form,collapse=NULL))
  
  
  #Fit model 
  if(is.na(MixedEff)) model <- glm(Formula, data=d)
  null.model=NULL
  if(!is.na(MixedEff))
  {
    model=lmer(Formula, data=d,REML = FALSE)
    null.model=lmer(as.formula(paste(ResVar, "~", paste("1+(1|",MixedEff,")",sep=""),collapse=NULL)), data=d,REML = FALSE)
  }
  
  
  return(list(model=model, null.model=null.model,data=d))
}

#data mining approach
Mod.fn.mining=function(d,ResVar,Predictrs,Y.type,Prop.train,ALL.models,nboot)
{
  d=d[,match(c(ResVar,Expl.vars),names(d))]
  d=d[!is.na(d[,match(ResVar,names(d))]),]
  d=subset(d,!is.na(EFFORT))
  if(length(d$BOTDEPTH)>0)
  {
    d=subset(d,!is.na(BOTDEPTH))
    d=subset(d,BOTDEPTH>=3)
    d$log.BOTDEPTH=log(d$BOTDEPTH)
  }
  Y=d[,match(ResVar,names(d))]
  d$log.EFFORT=log(d$EFFORT)
  
  #set factors
  for(x in 1:length(FactoRs)) d[,match(FactoRs[x],names(d))]=as.factor(d[,match(FactoRs[x],names(d))])
  Pred.form=paste(Predictrs,collapse="+")
  Formula=as.formula(paste(ResVar, "~", Pred.form,collapse=NULL)) 
  
  #A. Select best model
  
  #A.1. Select random training and testing data
  training_dp <- createDataPartition(Y, p = Prop.train, list = FALSE)
  training_partition <- d[training_dp, ]
  testing_partition <- d[-training_dp, ]
  
  #Reset factor levels
  for(x in 1:length(FactoRs))
  {
    Fac.lev=levels(d[,match(FactoRs[x],names(d))])
    training_partition[,match(FactoRs[x],names(d))]=factor(as.character(training_partition[,match(FactoRs[x],names(d))]),levels=Fac.lev)
    testing_partition[,match(FactoRs[x],names(d))]=factor(as.character(testing_partition[,match(FactoRs[x],names(d))]),levels=Fac.lev)
  }
  
  # set response var as factor if appropriate
  yy=match(ResVar,colnames(training_partition))
  if(Y.type=="Factor")
  {
    training_partition[,yy] <- as.factor(training_partition[,yy])
    testing_partition[,yy] <- as.factor(testing_partition[,yy])
  }
  
  # create an empty numeric vector to calculate out of sample error
  outOfSampleError <- numeric()
  
  
  # add some parameters for train control
  TC <- trainControl(method = "cv", number = 12, returnData=FALSE, 
                     returnResamp="none", savePredictions=FALSE, 
                     verboseIter=FALSE , preProcOptions="pca", allowParallel=TRUE)
  
  
  #2.1  (Neural Net)
  nnet <- train(Formula, method="nnet", data=training_partition, trControl= TC)
  nnetPrediction <- predict(nnet, testing_partition)
  if(Y.type=="Factor")
  {
    nnetAccuracy <- sum(nnetPrediction == testing_partition[,yy]) / length(nnetPrediction)
    nnetOutOfSampleError <- c(outOfSampleError, 1-nnetAccuracy)
  }else
  {
    nnetAccuracy=sqrt(mean((testing_partition[,yy]-nnetPrediction)^2))
    nnetOutOfSampleError <- c(outOfSampleError, nnetAccuracy)
  }
  
  #2.2 rf (Random Forest)
  rf <- train(Formula, method="rf", data=training_partition, trControl= TC)
  rfPrediction <- predict(rf, testing_partition)
  if(Y.type=="Factor")
  {
    rfAccuracy <- sum(rfPrediction == testing_partition[,yy]) / length(rfPrediction)
    rfOutOfSampleError <- c(outOfSampleError, 1-rfAccuracy)
  }else
  {
    rfAccuracy=sqrt(mean((testing_partition[,yy]-rfPrediction)^2))
    rfOutOfSampleError <- c(outOfSampleError, rfAccuracy)
  }
  
  #2.3  gbm (Generalized Boosted Regression)
  gbm <- train(Formula, method="gbm", data=training_partition, trControl= TC)
  gbmPrediction <- predict(gbm, testing_partition)
  if(Y.type=="Factor")
  {
    gbmAccuracy <- sum(gbmPrediction == testing_partition[,yy]) / length(gbmPrediction)
    gbmOutOfSampleError <- c(outOfSampleError, 1-gbmAccuracy)
  }else
  {
    gbmAccuracy=sqrt(mean((testing_partition[,yy]-gbmPrediction)^2))
    gbmOutOfSampleError <- c(outOfSampleError, gbmAccuracy)
  }
  
  #2.4 svmLinear (Support Vector Machines Linear)
  svml <- train(Formula, method="svmLinear", data=training_partition, trControl= TC)
  svmlPrediction <- predict(svml, testing_partition)
  if(Y.type=="Factor")
  {
    svmlAccuracy <- sum(svmlPrediction == testing_partition[,yy]) / length(svmlPrediction)
    svmlOutOfSampleError <- c(outOfSampleError, 1-svmlAccuracy)    
  }else
  {
    svmlAccuracy=sqrt(mean((testing_partition[,yy]-svmlPrediction)^2))
    svmlOutOfSampleError <- c(outOfSampleError, svmlAccuracy)
  }
  
  if(ALL.models=="YES")
  {
    
    #2.6 knn (K Nearest Neighbor)
    knn <- train(Formula, method="knn", data=training_partition, trControl= TC)
    knnPrediction <- predict(knn, testing_partition)
    if(Y.type=="Factor")
    {
      knnAccuracy <- sum(knnPrediction == testing_partition[,yy]) / length(knnPrediction)
      knnOutOfSampleError <- c(outOfSampleError, 1-knnAccuracy)
    }else
    {
      knnAccuracy=sqrt(mean((testing_partition[,yy]-knnPrediction)^2))
      knnOutOfSampleError <- c(outOfSampleError, knnAccuracy)
    }
    
    
    #2.8 rpart (Recursive Partitioning and Regression Trees)
    rpart <- train(Formula, method="rpart", data=training_partition, trControl= TC)
    rpartPrediction <- predict(rpart, testing_partition)
    if(Y.type=="Factor")
    {
      rpartAccuracy <- sum(rpartPrediction == testing_partition[,yy]) / length(rpartPrediction)
      rpartOutOfSampleError <- c(outOfSampleError, 1-rpartAccuracy)
    }else
    {
      rpartAccuracy=sqrt(mean((testing_partition[,yy]-rpartPrediction)^2))
      rpartOutOfSampleError <- c(outOfSampleError, rpartAccuracy)
    }
    
    #2.9 svmRadial (Support Vector Machines Radial)  
    svmr <- train(Formula, method="svmRadial", data=training_partition, trControl= TC)
    svmrPrediction <- predict(svmr, testing_partition)
    if(Y.type=="Factor")
    {
      svmrAccuracy <- sum(svmrPrediction == testing_partition[,yy]) / length(svmrPrediction)
      svmrOutOfSampleError <- c(outOfSampleError, 1-svmrAccuracy)
    }else
    {
      svmrAccuracy=sqrt(mean((testing_partition[,yy]-svmrPrediction)^2))
      svmrOutOfSampleError <- c(outOfSampleError, svmrAccuracy)
    }
    
    #2.10 treebag (Bagged Classification and Regression Trees)
    treebag <- train(Formula, method="treebag", data=training_partition, trControl= TC)
    treebagPrediction <- predict(treebag, testing_partition)
    if(Y.type=="Factor")
    {
      treebagAccuracy <- sum(treebagPrediction == testing_partition[,yy]) / length(treebagPrediction)
      treebagOutOfSampleError <- c(outOfSampleError, 1-treebagAccuracy)
    }else
    {
      treebagAccuracy=sqrt(mean((testing_partition[,yy]-treebagPrediction)^2))
      treebagOutOfSampleError <- c(outOfSampleError, treebagAccuracy)
    }
    
  }
  
  
  #A.3. Display values in table ranked by accuracy
  trainMethods <- c("Neural Net","Random Forest","Generalized Boosted Regression",
                    "Support Vector Machines Linear")
  accuracy <- c(nnetAccuracy,rfAccuracy,gbmAccuracy,svmlAccuracy)
  outOfSampleError <- c(nnetOutOfSampleError,rfOutOfSampleError,
                        gbmOutOfSampleError,svmlOutOfSampleError)
  
  if(ALL.models=="YES")
  {
    trainMethods <- c(trainMethods, 
                      "K Nearest Neighbor",
                      "Recursive Partitioning and Regression Trees",
                      "Support Vector Machines Radial",
                      "Bagged Classification and Regression Trees")
    accuracy <- c(accuracy,knnAccuracy,rpartAccuracy,
                  svmrAccuracy, treebagAccuracy)
    outOfSampleError=c(outOfSampleError,  
                       knnOutOfSampleError, 
                       rpartOutOfSampleError,svmrOutOfSampleError,
                       treebagOutOfSampleError)
  }
  
  results <- data.frame(trainMethods, accuracy, outOfSampleError)
  results=results[order(results$accuracy),]
  if(!Y.type=="Factor") colnames(results)[match("accuracy",colnames(results))]="RMSE"
  All.models.table=results
  
  
  #A.4 Choose best model
  if(Y.type=="Factor") Best.mod=as.character(results$trainMethods[which(results$accuracy==max(results$accuracy))])
  if(!Y.type=="Factor")Best.mod=as.character(results$trainMethods[which(results$RMSE==min(results$RMSE))])
  
  This.method=ifelse(Best.mod=="Bagged Classification and Regression Trees","treebag", 
                     ifelse(Best.mod=="Generalized Boosted Regression","gbm",             
                            ifelse(Best.mod=="K Nearest Neighbor","knn",                         
                                   ifelse(Best.mod=="Neural Net","nnet",                                 
                                          ifelse(Best.mod=="Random Forest","rf",                              
                                                 ifelse(Best.mod=="Recursive Partitioning and Regression Trees","rpart",
                                                        ifelse(Best.mod=="Support Vector Machines Linear","svmLinear",             
                                                               ifelse(Best.mod=="Support Vector Machines Radial","svmRadial",NA))))))))
  
  
  #Predict year effect using best model in bootstrap look to get CI
  
  #create new data
  NEWDATA=data.frame(YEAR=d[1:length(levels(d$YEAR)),match("YEAR",names(d))])                   
  NEWDATA$YEAR=factor(levels(d$YEAR),levels=levels(d$YEAR))
  Other.preds=Predictrs[-match("YEAR",Predictrs)]
  d.frm=t(matrix(rep(NA,length(Other.preds))))
  colnames(d.frm)=Other.preds
  d.frm=as.data.frame(d.frm)
  for(q in 1:length(Other.preds))  
  {
    Vec=d[,match(Other.preds[q],names(d))]
    CLs=class(Vec)
    if(CLs=='factor') Value=factor(Mst.cmn(Vec),levels=levels(Vec))
    if(CLs=='numeric') Value=mean(Vec)
    d.frm[,q]=Value
  }
  NEWDATA=cbind(NEWDATA,d.frm)
  
  Preds=matrix(NA,nrow=nboot,ncol=nrow(NEWDATA))
  for(nn in 1:nboot)
  {
    index=sample(1:nrow(d),nrow(d),replace=T)
    dboot=d[index,]
    
    #use best model to fit data
    model <- train(Formula, method=This.method, data=dboot, trControl= TC)
    Preds[nn,]=predict(model,newdata=NEWDATA)
  }
  return(list(All.models.table=All.models.table,Best.method=This.method, Preds=Preds))
}

#Get mean and CI function
mean.quan=function(r,CI=c(0.975,0.025))
{
  Mean=Quan1=Quan2=vector("numeric",ncol(r)) 
  for(i in 1:ncol(r)) 
  {
    X=r[,i]
    Mean[i]=mean(X)
    Quan1[i]=quantile(X,CI[1])
    Quan2[i]=quantile(X,CI[2])
  }
  return(list(MEAN=Mean,UQ=Quan1,LQ=Quan2))
}

#see blockx
fn.see=function(d)
{
  par(mfcol=c(4,3),mai=c(.15,.2,.15,.05),mgp=c(1,.5,0))
  Mn=1:12
  Ylim=c(-36,-19)
  Xlim=c(113,130)
  for(n in 1:length(Mn))
  {
    x=subset(d,MONTH==Mn[n],select=c(LATITUDE, LONGITUDE, N))
    if(nrow(x)==0)plot(1,xaxt='n',yaxt='n',col="transparent",main=Mn[n],ylab='',xlab='')else
    {
      x=x[order(x$LATITUDE),]
      Reshp=reshape(x,v.names = "N", idvar = "LONGITUDE",timevar = "LATITUDE", direction = "wide")
      Reshp=Reshp[order(Reshp$LONGITUDE),]
      z=as.matrix(Reshp[,-match("LONGITUDE",colnames(Reshp))])
      Breaks=range(unlist(z),na.rm=TRUE)
      Breaks=seq(Breaks[1],Breaks[2],1)
      if(length(Breaks)>1)
      {
        numberLab=2
        numInt=length(Breaks)-1
        couleurs=rev(heat.colors(numInt))
        
        image(Reshp$LONGITUDE,sort(unique(x$LATITUDE)),z,main=Mn[n],ylab="",xlab="",col =couleurs,breaks=Breaks,
              ylim=Ylim,xlim=Xlim)
        color.legend(129,-20,130,-34,Breaks,rect.col=couleurs,gradient="y",col=1,cex=0.95)
      }else
        plot(1,xaxt='n',yaxt='n',col="transparent",main=Mn[n],ylab='',xlab='') 
      
    }
  }
  mtext(yrs[i],3,-3,cex=1.5,outer=T)
}

#Function for expressing index in relative terms
fn.relative=function(D)
{
  Mn=mean(D$MEAN)
  D$LOW=D$LOW/Mn
  D$UP=D$UP/Mn
  D$MEAN=D$MEAN/Mn
  return(D)
}

#Plot predictions
fun.plt.pred=function(d,Show.pred,normalised,PredictorS,MAIN,log.var,Cx,YLIM,Cx.axs)
{
  #create new data
  Dat=d$data
  Other.preds=c(PredictorS[-match(Show.pred,PredictorS)])
  if(!is.na(OFFSETT))Other.preds=c(PredictorS[-match("year",PredictorS)],"log.EFFORT")
  
  idPred=match(Show.pred,names(Dat))
  if(is.factor(Dat[,idPred]))
  {
    NEWDATA=data.frame(x=sort(unique(as.character(Dat[,idPred]))))
    NEWDATA[,1]=factor(NEWDATA[,1],levels=levels(Dat[,idPred]))
    names(NEWDATA)=Show.pred
  }else
  {
    NEWDATA=data.frame(x=sort(seq(min(Dat[,idPred]),max(Dat[,idPred]),by=BY)))
    names(NEWDATA)=Show.pred
  }
  
  d.frm=t(matrix(rep(NA,length(Other.preds))))
  colnames(d.frm)=Other.preds
  d.frm=as.data.frame(d.frm)
  for(q in 1:length(Other.preds))  
  {
    Vec=Dat[,match(Other.preds[q],names(Dat))]
    CLs=class(Vec)
    if(CLs=='factor') Value=factor(Mst.cmn(Vec),levels=levels(Vec))
    if(CLs=='numeric') Value=mean(Vec)
    d.frm[,q]=Value
  }                  
  NEWDATA=cbind(NEWDATA,d.frm)
  Covar.pos=as.matrix(vcov(d$model))
  set.seed(999);Pos.pars.rand=rmvnorm(niter,mean=coef(d$model),sigma=Covar.pos)    
  MC.preds=matrix(nrow=niter,ncol=length(NEWDATA[,match(Show.pred,names(NEWDATA))]))
  
  for(n in 1:niter)
  {
    model=d$model
    model$coefficients=Pos.pars.rand[n,]
    a=predict(model,newdata=NEWDATA,type='response',se.fit=T)
    if(log.var=="YES") Pred=exp(a$fit+(a$se.fit^2)/2)  #apply bias correction for log transf
    if(log.var=="NO") Pred=a$fit
    MC.preds[n,]=Pred
  }
  
  PRD=data.frame(X=NEWDATA[,match(Show.pred,names(NEWDATA))],
                 MEAN=colMeans(MC.preds,na.rm=T),
                 LOW=apply(MC.preds, 2, function(x) quantile(x, 0.025,na.rm=T)),
                 UP=apply(MC.preds, 2, function(x) quantile(x, 0.975,na.rm=T)))
  names(PRD)[1]=Show.pred
  
  #standardise to a mean score of 1
  if(normalised=="YES") PRD=fn.relative(PRD)
  XX=as.numeric(as.character(PRD[,match(Show.pred,names(PRD))]))
  MeAn=PRD$MEAN
  UppCI=PRD$UP
  LowCI=PRD$LOW
  dat.plt=data.frame(yr=XX,MeAn=MeAn,UppCI=UppCI,LowCI=LowCI)
  
  #add missing years
  if(Show.pred=="year")
  {
    mis.yr=seq(ALL.yrs[1],ALL.yrs[length(ALL.yrs)])
    mis.yr=mis.yr[which(!mis.yr%in%dat.plt$yr)]
    if(length(mis.yr)>0)
    {
      dummy=dat.plt[1:length(mis.yr),]
      dummy[,]=NA
      dummy$yr=mis.yr
      dat.plt=rbind(dat.plt,dummy)
    }
    dat.plt=dat.plt[order(dat.plt$yr),]
  }
  
  if(is.null(YLIM)) YLIM=c(min(dat.plt$LowCI,na.rm=T),max(dat.plt$UppCI,na.rm=T))
  
  with(dat.plt,plot(yr,MeAn,pch=19,main="",xlab="",ylab="",
                    cex=1.25,cex.axis=1.25,ylim=YLIM))
  #with(dat.plt,plot(yr,MeAn,pch=19,main=MAIN,xlab="",ylab="",
  #                  cex=1.25,xaxt="n",cex.axis=1.25,ylim=YLIM))
  with(dat.plt,arrows(x0=yr, y0=LowCI, x1=yr, y1=UppCI,code = 3,angle=90,length=.025))
  #axis(1,dat.plt$yr,F,tck=-0.02)
  #with(dat.plt,axis(1,seq(yr[1],yr[length(yr)],5),F,tck=-0.05))
  #if(i %in%c(3,6)) with(dat.plt,axis(1,seq(yr[1],
  #     yr[length(yr)],5),seq(yr[1],yr[length(yr)],5),tck=-0.05,cex.axis=1.25))
  mtext(MAIN,3,cex=1.25)
}

#Plot year predictions
fun.plt.yr.pred= function(d,normalised,PredictorS,MAIN,log.var,Cx,YLIM,Cx.axs)
{
  #create new data
  Dat=d$data
  NEWDATA=data.frame(YEAR=sort(unique(as.character(Dat$YEAR))))                   
  NEWDATA$YEAR=factor(NEWDATA$YEAR,levels=levels(d$data$YEAR))
  Other.preds=c(PredictorS[-match(c("YEAR"),PredictorS)])
  #Other.preds=c(PredictorS[-match(c("YEAR"),PredictorS)],"log.EFFORT")
  d.frm=t(matrix(rep(NA,length(Other.preds))))
  colnames(d.frm)=Other.preds
  d.frm=as.data.frame(d.frm)
  for(q in 1:length(Other.preds))  
  {
    Vec=Dat[,match(Other.preds[q],names(Dat))]
    CLs=class(Vec)
    if(CLs=='factor') Value=factor(Mst.cmn(Vec),levels=levels(Vec))
    if(CLs=='numeric') Value=mean(Vec)
    d.frm[,q]=Value
  }
  NEWDATA=cbind(NEWDATA,d.frm)
  
  # Preds=predict(d$model,newdata=NEWDATA,type='response',se.fit=T)
  # 
  # if(log.var=="NO")PRD=data.frame(Year=NEWDATA$YEAR,Mean=Preds$fit,SE=Preds$se.fit)
  # if(log.var=="YES")PRD=data.frame(Year=NEWDATA$YEAR,Mean=exp(Preds$fit),SE=exp(Preds$se.fit))
  
  Covar.pos=as.matrix(vcov(d$model))
  set.seed(999);Pos.pars.rand=rmvnorm(niter,mean=coef(d$model),sigma=Covar.pos)    
  MC.preds=matrix(nrow=niter,ncol=length(NEWDATA$YEAR))
  
  for(n in 1:niter)
  {
    model=d$model
    model$coefficients=Pos.pars.rand[n,]
    a=predict(model,newdata=NEWDATA,type='response',se.fit=T)
    
    if(log.var=="YES") Pred=exp(a$fit+(a$se.fit^2)/2)  #apply bias correction for log transf
    if(log.var=="NO") Pred=a$fit
    
    MC.preds[n,]=Pred
  }
  
  PRD=data.frame(Year=NEWDATA$YEAR,
                 MEAN=colMeans(MC.preds,na.rm=T),
                 LOW=apply(MC.preds, 2, function(x) quantile(x, 0.025,na.rm=T)),
                 UP=apply(MC.preds, 2, function(x) quantile(x, 0.975,na.rm=T)))
  
  #standardise to a mean score of 1
  if(normalised=="YES") PRD=fn.relative(PRD)
  
  yr=as.numeric(as.character(PRD$Year))
  MeAn=PRD$MEAN
  UppCI=PRD$UP
  LowCI=PRD$LOW
  
  
  dat.plt=data.frame(yr=yr,MeAn=MeAn,UppCI=UppCI,LowCI=LowCI)
  
  
  #add missing years
  mis.yr=seq(ALL.yrs[1],ALL.yrs[length(ALL.yrs)])
  mis.yr=mis.yr[which(!mis.yr%in%yr)]
  if(length(mis.yr)>0)
  {
    dummy=dat.plt[1:length(mis.yr),]
    dummy[,]=NA
    dummy$yr=mis.yr
    dat.plt=rbind(dat.plt,dummy)
  }
  dat.plt=dat.plt[order(dat.plt$yr),]
  
  if(is.null(YLIM)) YLIM=c(min(dat.plt$LowCI,na.rm=T),max(dat.plt$UppCI,na.rm=T))
  with(dat.plt,plot(yr,MeAn,pch=19,main=MAIN,xlab="",ylab="",
                    cex=Cx,xaxt="n",cex.axis=1.25,cex.lab=1.5,cex.main=1.75,ylim=YLIM))
  with(dat.plt,arrows(x0=yr, y0=LowCI, x1=yr, y1=UppCI,code = 3,angle=90,length=.025))
  with(dat.plt,axis(1,yr,F,tck=-0.025))
  with(dat.plt,axis(1,seq(yr[1],yr[length(yr)],5),F,tck=-0.05))
  with(dat.plt,axis(1,seq(yr[1],yr[length(yr)],5),seq(yr[1],yr[length(yr)],5),
                    tck=-0.05,cex.axis=Cx.axs))
}

fn.reshp=function(d,Y,TimeVar,IdVAR)
{
  Form=formula(paste(Y,paste(c(TimeVar,IdVAR),collapse="+"),sep="~"))
  d=aggregate(Form,d,sum)
  DATA.wide=reshape(d[,match(c(ResVar,IdVAR,TimeVar),names(d))],v.names=ResVar,
                    idvar=IdVAR,timevar=TimeVar,direction="wide")
  DATA.wide[is.na(DATA.wide)]=0
  colnames(DATA.wide)=gsub(paste(ResVar,".",sep=""), "", names(DATA.wide))
  return(DATA.wide)
}

#Function for looping over different observer data sets
fn.loop.over.obsrvr.data=function(DaTA,dat.nm,normalised,Drop.yrs)
{
  #Define years to use
  if(Drop.yrs=="YES")DaTA=subset(DaTA,!YEAR%in%as.character(1993:1999))
  
  #select shots with a minimum number of individuals
  N.ind.shot=aggregate(INDIVIDUALS~SHEET_NO,DaTA,sum)
  N.ind.shot=subset(N.ind.shot,INDIVIDUALS>=Min.individuals)
  DaTA=subset(DaTA,SHEET_NO%in%unique(N.ind.shot$SHEET_NO))
  
  # 3.1.2 Diversity and ecosystem indicators
  Dat=fn.reshp(d=DaTA,Y=ResVar,TimeVar=MultiVar,IdVAR=IDVAR) 
  DATA.shots.diversity=Dat[,match(IDVAR,names(Dat))]
  Dat.y=Dat[,-match(IDVAR,names(Dat))]
  
  DATA.shots.diversity$Shannon = diversity(Dat.y, index = "shannon")
  DATA.shots.diversity$Simpson = diversity(Dat.y, index = "simpson")
  sp_names=colnames(Dat.y)
  Evenness = function(H, S) {J = H / log(S)}  #where H = Shannon index for each community; S = Number of species
  DATA.shots.diversity$Pielou = Evenness(H = DATA.shots.diversity$Shannon, S = length(sp_names))
  DATA.shots.diversity$Pielou[DATA.shots.diversity$Pielou == 0] = NA
  
  #ecosystems indicators   
  Numbers=aggregate(INDIVIDUALS~SHEET_NO+SPECIES,DaTA,sum,na.rm=T)
  MTL=aggregate(TROPHIC_LEVEL~SHEET_NO+SPECIES,DaTA,mean,na.rm=T) 
  MTL=merge(MTL,Numbers,by=c("SHEET_NO","SPECIES"),all=T)
  MTL <- data.table(MTL)
  MTL=as.data.frame(MTL[,list(MTL = weighted.mean(TROPHIC_LEVEL,INDIVIDUALS,na.rm=T)),by=SHEET_NO])
  
  
  DaTA$FL=with(DaTA,ifelse(is.na(FL),TL,FL))
  MML=aggregate(FL~SHEET_NO+SPECIES,DaTA,mean,na.rm=T)
  MML=merge(MML,Numbers,by=c("SHEET_NO","SPECIES"),all.y=T)
  MML <- data.table(MML)
  MML=as.data.frame(MML[,list(MML = weighted.mean(FL,INDIVIDUALS,na.rm=T)),by=SHEET_NO])
  
  
  DATA.shots.diversity=merge(DATA.shots.diversity,MTL,by=c("SHEET_NO"),all.x=T)
  DATA.shots.diversity=merge(DATA.shots.diversity,MML,by=c("SHEET_NO"),all.x=T)
  
  
  DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(Shannon))
  #DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(Margalef))
  #DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(Pielou))
  DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(Simpson))
  DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(MTL))
  DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(MML))
  DATA.shots.diversity=subset(DATA.shots.diversity,Simpson>0)
  DATA.shots.diversity=subset(DATA.shots.diversity,MML>0)
  
  
  
  #Exploratory analyses
  if (do.exploratory=="YES")
  {
    #Look at annual averages
    fn.fig("Preliminary_bxplt_year",2000,2000) 
    smart.par(n.plots=length(resp.vars),MAR=c(2,2,1,1),OMA=c(1,1.5,.1,.1),MGP=c(2.5,.7,0))
    for(ix in 1:length(resp.vars))
    {
      boxplot(DATA.shots.diversity[,match(resp.vars[ix],names(DATA.shots.diversity))]~DATA.shots.diversity$YEAR)
      mtext(resp.vars[ix],3)
    }
    dev.off()
    
    #look at number of skippers
    Skip.tab=sort(table(DATA.shots.diversity$SKIPPER))
    #Drop.skip=c(names(Skip.tab)[Skip.tab<10])
    #note: SKIPPER removed from Predictors because data don't allow estimation of coefficients
    
    
    Table.yrs.blk=with(DATA.shots.diversity,table(YEAR,BLOCK,useNA='ifany'))
    
    # par(mfcol=c(2,2),las=1,mar=c(3,1.5,1,1.2),oma=c(1,1,.1,1))
    # for(i in 1:n.rv)
    # {
    #   for(a in 1:length(AREAS))
    #   {
    #     fn.bx.plt(d=subset(DATA.shots.diversity,AREA==AREAS[a]),Var=resp.vars[i])
    #     legend("top",AREAS[a],bty='n')
    #   }
    # }
    
    fn.fig("Preliminary_error_structure_Obs",2000,2000) 
    par(mfrow=c(n.rv,2),las=1,mar=c(1,1.5,1,1.2),oma=c(1,1,.1,1),mgp=c(1,.5,0))
    for(i in 1:n.rv)
    {
      density.response(d=DATA.shots.diversity,Var=resp.vars[i],Log="NO")
      density.response(d=DATA.shots.diversity,Var=resp.vars[i],Log="YES")
    }
    dev.off()
    
    #Summary
    summary(DATA.shots.diversity)
    
    #apply function to all predictors
    attach(DATA.shots.diversity)
    par(mfcol=c(2,1),mai=c(1,1,.1,.1))
    v=c("BOTDEPTH")
    lapply(v,jplot)
    detach(DATA.shots.diversity)
  }
  
  #Remove boats with less than Min.recs
  AA=sort(table(DATA.shots.diversity$BOAT))
  AA=AA[AA>Min.recs]
  DATA.shots.diversity=subset(DATA.shots.diversity,BOAT%in%names(AA))
  DATA.shots.diversity=subset(DATA.shots.diversity,!BOAT%in%c("E67","E34","E10"))  #cannot estimate coef
  #Remove years with less than Min.recs
  AA=sort(table(DATA.shots.diversity$YEAR))
  AA=AA[AA>Min.recs]
  DATA.shots.diversity=subset(DATA.shots.diversity,YEAR%in%names(AA))
  
  #Remove blocks with less than Min.recs
  AA=sort(table(DATA.shots.diversity$BLOCK))
  AA=AA[AA>Min.recs]
  DATA.shots.diversity=subset(DATA.shots.diversity,BLOCK%in%names(AA))
  
  #Check if correlation between effort and indices
  fn.fig("Correlation_effort_indices_Observer",2000,2000) 
  smart.par(n.plots=length(resp.vars),MAR=c(1.75,3.5,1.5,.1),OMA=c(3,1,.3,2),MGP=c(1,.8,0))
  for(i in 1:n.rv)
  {
    dd=DATA.shots.diversity[,match(c("EFFORT",resp.vars[i]),names(DATA.shots.diversity))]
    plot(dd$EFFORT,dd[,2])
    mtext(paste(resp.vars[i]," cor=",round(cor(dd[,1],dd[,2]),2),sep=""),3,col=2)
    rm(dd)
  }
  dev.off()
  
  #run model
  Store.mod.out.observer=vector('list',length(resp.vars))
  names(Store.mod.out.observer)=resp.vars
  
  #glm approach
  for(i in 1:n.rv)
  {
    Store.mod.out.observer[[i]]=Mod.fn.glm(d=DATA.shots.diversity,
                                           ResVar=resp.vars[i],Expl.vars=Expl.varS,
                                           Predictrs=Predictors,FactoRs=FactoRS,
                                           OFFSET=OFFSETT,
                                           log.var=Res.var.in.log[i],add.inter="NO",
                                           MixedEff=MixedEff)
  }
  
  #data mining
  #system.time(for(i in 1:n.rv)Store.mod.out.observer[[i]]=Mod.fn.mining(d=DATA.shots.diversity,ResVar=resp.vars[i],
  #                   Predictrs=Predictors,Y.type="Continuous",Prop.train=.7,ALL.models="NO",nboot=1))
  
  #Anova tables
  TABL=vector('list',length(resp.vars))
  names(TABL)=resp.vars
  
  for(i in 1:n.rv)
  {
    Modl=Store.mod.out.observer[[i]]$model
    if(is.na(MixedEff))
    {
      #each term
      Anova.tab=anova(Modl, test = "Chisq")
      n=2:length(Anova.tab$Deviance)
      Term.dev.exp=100*(Anova.tab$Deviance[n]/Modl$null.deviance)
      names(Term.dev.exp)=rownames(Anova.tab)[n]
      
      #nice table
      Anov.tab=as.data.frame.matrix(Anova.tab)
      Term.tab=data.frame(Percent.dev.exp=Term.dev.exp)
      Anova.tab=Anova.tab[-1,match(c("Deviance","Pr(>Chi)"),names(Anova.tab))]
      Anova.tab=cbind(Anova.tab,Term.tab)
      Anova.tab=Anova.tab[,-match("Deviance",names(Anova.tab))]
      Anova.tab$"Pr(>Chi)"=ifelse(Anova.tab$"Pr(>Chi)"<0.001,"<0.001",round(Anova.tab$"Pr(>Chi)",3))
      Total=Anova.tab[1,]
      Total$"Pr(>Chi)"=""
      Total$Percent.dev.exp=sum(Anova.tab$Percent.dev.exp)
      rownames(Total)="Total"
      Anova.tab=rbind(Anova.tab,Total)
      Anova.tab$Percent.dev.exp=round(Anova.tab$Percent.dev.exp,2)
    }
    if(!is.na(MixedEff))
    {
      Dev.ex.fixed.and.fixed_random=r.squaredGLMM(Modl)
      Anova.tab=as.data.frame.matrix(Anova(Modl))
      dd=Anova.tab[1:3,]
      dd[,]=""
      row.names(dd)=c("","r2_fixed_effects","r2_fixed&random_effects")
      dd[2:3,1]=Dev.ex.fixed.and.fixed_random
      Anova.tab=rbind(Anova.tab,dd)
    }
    TABL[[i]]=Anova.tab
  }
  
  TABL=do.call(cbind,TABL)
  TABL=cbind(TERM=row.names(TABL),TABL)
  
  if(is.na(MixedEff))
  {
    LBL_1st=c(1,rep(2,length(resp.vars)))
    LBL_2nd=c("",rep(c("P","% dev. expl."),length(resp.vars)))
  }
  
  if(!is.na(MixedEff))
  {
    LBL_1st=c(1,rep(3,length(resp.vars)))
    LBL_2nd=c("",rep(c("Chisq","Df","P"),length(resp.vars)))
  }
  #export anova tables as word doc
  Export.tbl(WD=getwd(),Tbl=TABL,Doc.nm=paste("Anova.table.obs",dat.nm,sep="_"),caption=NA,paragph=NA,
             HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
             Zebra='NO',Zebra.col='grey60',Grid.col='black',
             Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
             HDR.names=c('TERM', resp.vars),
             HDR.span=LBL_1st,
             HDR.2nd=LBL_2nd)
  
  #Plot model fit
  for(i in 1:n.rv)
  {
    fn.fig(paste("GLM_Fit_Observer",names(Store.mod.out.observer)[i],sep="_"),2000,2000) 
    par(mfcol=c(2,2))
    plot(Store.mod.out.observer[[i]]$model)
    dev.off()
  }
  
  
  #Plot diversity and ecosystem indices
  #note: present in relative terms as we are interested in the trend
  fn.fig(paste("Div & Eco indices.observer_normalised",dat.nm,sep="_"),1800,2000)  
  smart.par(n.plots=length(resp.vars),MAR=c(1.75,3.5,1.5,.1),OMA=c(3,1,.3,2),MGP=c(1,.8,0))
  #loop over each index
  for(i in 1:length(Store.mod.out.observer))
  {
    fun.plt.pred(d=Store.mod.out.observer[[i]],Show.pred="YEAR",normalised=normalised,
                 PredictorS=Predictors,MAIN=Main.title[i],
                 log.var=Res.var.in.log[i],Cx=1.125,YLIM=NULL,Cx.axs=1.35)
  }
  mtext("Year",1,outer=T,cex=1.5,line=1.75)
  mtext("Relative value",2,-.75,outer=T,cex=1.5,las=3)
  dev.off()
  
  
  #Observer effort used
  # if(dat.nm=="BC")
  # {
  #   Annual.eff=aggregate(EFFORT~YEAR,DATA.shots.diversity,sum,na.rm=T)
  #   Yrs.all=as.numeric(sort(Annual.eff$YEAR))
  #   Yrs.all=seq(Yrs.all[1],Yrs.all[length(Yrs.all)])
  #   misn.yr=Yrs.all[which(!Yrs.all%in%as.numeric(Annual.eff$YEAR))]
  #   ADD=matrix(NA,nrow=length(misn.yr),ncol=ncol(Annual.eff))
  #   colnames(ADD)=colnames(Annual.eff)
  #   ADD=as.data.frame(ADD)
  #   ADD$YEAR=misn.yr
  #   Annual.eff=rbind(Annual.eff,ADD)  
  #   Annual.eff=Annual.eff[order(Annual.eff$YEAR),]  
  #   
  #   fn.fig(paste("Observers Effort",dat.nm,sep="_"),1500,1500)  
  #   par(las=1,mar=c(3.5,3.5,1.5,.1),oma=c(1,1,.3,2),mgp=c(2.1,.7,0))
  #   plot(Annual.eff$YEAR,Annual.eff$EFFORT,type='l',
  #        ylab='Effort',xlab='Year',cex.lab=1.6)
  #   dev.off()
  # }
  
  
  #3.1.3 Bycatch patterns     
  #note: for consistency, use same records as for indices
  Bycatch=subset(DaTA,FATE%in%c("C","D") &SHEET_NO%in%unique(DATA.shots.diversity$SHEET_NO) ) 
  Bycatch$Discard=with(Bycatch,ifelse(FATE=="D",1,0))
  
  
  #Is there a change in what prortion is discarded by species?  
  Disc=aggregate(Discard~YEAR+SPECIES,Bycatch,sum)
  Tot= aggregate(INDIVIDUALS~YEAR+SPECIES,Bycatch,sum)
  Disc=merge(Disc,Tot,by=c("YEAR","SPECIES"),all=T)
  Disc$Prop=Disc$Discard/Disc$INDIVIDUALS
  dummy=as.numeric(sort(unique(Disc$YEAR)))
  all.yrs=dummy[1]:dummy[length(dummy)]
  msn.yrs=all.yrs[which(!all.yrs%in%dummy)]
  dummy=Disc[1:length(msn.yrs),]
  dummy$YEAR=msn.yrs
  dummy$Prop=NA
  Disc=rbind(Disc,dummy)
  Disc=Disc[order(Disc$YEAR),]
  
  wide <- reshape(Disc[,-match(c("Discard","INDIVIDUALS"),names(Disc))], v.names = "Prop", idvar = "SPECIES",
                  timevar = "YEAR", direction = "wide")
  colnames(wide)=gsub(paste("Prop",".",sep=""), "", names(wide))
  wide[is.na(wide)]="transparent"
  wide[wide=="1"]="black"
  wide[wide=="0"]="grey60"
  xx=as.numeric(colnames(wide)[-1])
  wide=wide[order(wide$SPECIES),]
  
  fn.fig("Discard_species_by_yr",1200,2400) 
  par(las=1,mar=c(1,2,0,0),oma=c(1.5,1,0,0),mgp=c(1,.6,0),cex.lab=1.25)
  plot(xx,rep(0,length(xx)),pch=19,col=unlist(wide[1,2:ncol(wide)]),cex=.85,ylab="",xlab="",yaxt='n',ylim=c(0,nrow(wide)+1))
  for(n in 2:nrow(wide)) points(xx,rep(n-1,length(xx)),pch=19,cex=.85,col=unlist(wide[n,2:ncol(wide)]))
  axis(2,0:(nrow(wide)-1),wide$SPECIES,cex.axis=.4)
  mtext("Species",2,1.8,cex=1.5,las=3)
  mtext("Year",1,1.5,cex=1.5)
  legend('top',c("Discarded","Retained"),bty='n',pch=19,col=c("black","grey60"),horiz=T,cex=1.25)
  dev.off()
  
  
  #set factors
  for(x in 1:length(FactoRS)) Bycatch[,match(FactoRS[x],names(Bycatch))]=as.factor(Bycatch[,match(FactoRS[x],names(Bycatch))])
  
  #Binominal GLM
  bycatch.model=glm(Discard ~ YEAR + BLOCK + BOAT + MONTH + BOTDEPTH,data=Bycatch,family=binomial)
  
  #Binominal GLM  
  Bycatch.glm=fn.reshp(d=Bycatch,Y=ResVar,TimeVar="FATE",IdVAR=IDVAR) 
  Bycatch.glm$Tot=Bycatch.glm$D+Bycatch.glm$C
  Bycatch.glm$Prop.disc=Bycatch.glm$D/(Bycatch.glm$Tot)
  Bycatch.glm[is.na(Bycatch.glm)]=0
  
  
  #Anova tables
  Anova.tab=anova(bycatch.model, test = "Chisq")
  n=2:length(Anova.tab$Deviance)
  Term.dev.exp=100*(Anova.tab$Deviance[n]/bycatch.model$null.deviance)
  names(Term.dev.exp)=rownames(Anova.tab)[n]
  Anov.tab=as.data.frame.matrix(Anova.tab)
  Term.tab=data.frame(Percent.dev.exp=Term.dev.exp)
  Anova.tab=Anova.tab[-1,match(c("Deviance","Pr(>Chi)"),names(Anova.tab))]
  Anova.tab=cbind(Anova.tab,Term.tab)
  Anova.tab=Anova.tab[,-match("Deviance",names(Anova.tab))]
  Anova.tab$"Pr(>Chi)"=ifelse(Anova.tab$"Pr(>Chi)"<0.001,"<0.001",round(Anova.tab$"Pr(>Chi)",3))
  Total=Anova.tab[1,]
  Total$"Pr(>Chi)"=""
  Total$Percent.dev.exp=sum(Anova.tab$Percent.dev.exp)
  rownames(Total)="Total"
  Anova.tab=rbind(Anova.tab,Total)
  Anova.tab$Percent.dev.exp=round(Anova.tab$Percent.dev.exp,2)
  Anova.tab=cbind(rownames(Anova.tab),Anova.tab)
  
  #export anova tables as word doc
  Export.tbl(WD=getwd(),Tbl=Anova.tab,Doc.nm=paste("Anova.table_Discard_observer.table",dat.nm,sep="_"),caption=NA,paragph=NA,
             HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
             Zebra='NO',Zebra.col='grey60',Grid.col='black',
             Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
             HDR.names=c('TERM', "P","% dev. expl."),
             HDR.span=c(1,1,1),
             HDR.2nd=c("","",""))
  
  #show annual predictions of proportion of discards
  ALL.yrs.bycatch=as.numeric(levels(Bycatch.glm$YEAR))
  
  fun.plt.yr.pred.bycatch=function(d,normalised,PredictorS,MAIN,log.var,Cx,YLIM,Cx.axs)
  {
    #create new data
    Dat=d$data
    NEWDATA=data.frame(YEAR=sort(unique(as.character(Dat$YEAR))))                   
    NEWDATA$YEAR=factor(NEWDATA$YEAR,levels=levels(Dat$YEAR))
    Other.preds=c(PredictorS[-match("YEAR",PredictorS)])
    d.frm=t(matrix(rep(NA,length(Other.preds))))
    colnames(d.frm)=Other.preds
    d.frm=as.data.frame(d.frm)
    for(q in 1:length(Other.preds))  
    {
      Vec=Dat[,match(Other.preds[q],names(Dat))]
      CLs=class(Vec)
      if(CLs=='factor') Value=factor(Mst.cmn(Vec),levels=levels(Vec))
      if(CLs=='numeric') Value=mean(Vec)
      d.frm[,q]=Value
    }
    NEWDATA=cbind(NEWDATA,d.frm)
    
    # Preds=predict(d$model,newdata=NEWDATA,type='response',se.fit=T)
    # 
    # if(log.var=="NO")PRD=data.frame(Year=NEWDATA$YEAR,Mean=Preds$fit,SE=Preds$se.fit)
    # if(log.var=="YES")PRD=data.frame(Year=NEWDATA$YEAR,Mean=exp(Preds$fit),SE=exp(Preds$se.fit))
    
    Covar.pos=as.matrix(vcov(d$model))
    set.seed(999);Pos.pars.rand=rmvnorm(niter,mean=coef(d$model),sigma=Covar.pos)    
    MC.preds=matrix(nrow=niter,ncol=length(NEWDATA$YEAR))
    
    for(n in 1:niter)
    {
      model=d$model
      model$coefficients=Pos.pars.rand[n,]
      a=predict(model,newdata=NEWDATA,type='response',se.fit=T)
      
      if(log.var=="YES") Pred=exp(a$fit+(a$se.fit^2)/2)  #apply bias correction for log transf
      if(log.var=="NO") Pred=a$fit
      
      MC.preds[n,]=Pred
    }
    
    PRD=data.frame(Year=NEWDATA$YEAR,
                   MEAN=colMeans(MC.preds,na.rm=T),
                   LOW=apply(MC.preds, 2, function(x) quantile(x, 0.025,na.rm=T)),
                   UP=apply(MC.preds, 2, function(x) quantile(x, 0.975,na.rm=T)))
    
    #standardise to a mean score of 1
    if(normalised=="YES") PRD=fn.relative(PRD)
    
    yr=as.numeric(as.character(PRD$Year))
    MeAn=PRD$MEAN
    UppCI=PRD$UP
    LowCI=PRD$LOW
    
    
    dat.plt=data.frame(yr=yr,MeAn=MeAn,UppCI=UppCI,LowCI=LowCI)
    
    
    #add missing years
    mis.yr=seq(ALL.yrs.bycatch[1],ALL.yrs.bycatch[length(ALL.yrs.bycatch)])
    mis.yr=mis.yr[which(!mis.yr%in%yr)]
    if(length(mis.yr)>0)
    {
      dummy=dat.plt[1:length(mis.yr),]
      dummy[,]=NA
      dummy$yr=mis.yr
      dat.plt=rbind(dat.plt,dummy)
    }
    dat.plt=dat.plt[order(dat.plt$yr),]
    
    if(is.null(YLIM)) YLIM=c(min(dat.plt$LowCI,na.rm=T),max(dat.plt$UppCI,na.rm=T))
    with(dat.plt,plot(yr,MeAn,pch=19,main=MAIN,xlab="",ylab="",
                      cex=Cx,xaxt="n",cex.axis=1.25,cex.lab=1.5,cex.main=1.75,ylim=YLIM))
    with(dat.plt,arrows(x0=yr, y0=LowCI, x1=yr, y1=UppCI,code = 3,angle=90,length=.025))
    with(dat.plt,axis(1,yr,F,tck=-0.05))
    # with(dat.plt,axis(1,seq(yr[1],yr[length(yr)],5),F,tck=-0.05))
    # with(dat.plt,axis(1,seq(yr[1],yr[length(yr)],5),seq(yr[1],yr[length(yr)],5),
    #                                    tck=-0.05,cex.axis=Cx.axs))
  }
  
  fn.fig(paste("Discard_observer_predicted_prop",dat.nm,sep="_"),2000,2400)  
  par(mfrow=c(2,1),mai=c(.2,.1,.1,1.1),oma=c(2,4,.1,.1),xpd=TRUE,mgp=c(2.5,.75,0),las=1)
  
    #Species composition
  Bycatch=subset(Bycatch,SHEET_NO%in%unique(Bycatch.glm$SHEET_NO))
  YrS=sort(unique(as.numeric(as.character(Bycatch$YEAR))))
  YrS=seq(as.numeric(YrS[1]),as.numeric(YrS[length(YrS)]))
  d=Bycatch
  d$SPECIES=as.character(d$SPECIES)
  d$SP.fate=with(d,ifelse(SPECIES%in%
                            c("BB.T","DM.T","ER","PJ"),SPECIES,
                          ifelse(FATE=="C","Comm.","Other")))
  
  Tot=aggregate(INDIVIDUALS~SP.fate+YEAR,d,sum)
  Tot.reshaped=reshape(Tot,v.names = "INDIVIDUALS", idvar =c("YEAR"),
                       timevar = "SP.fate", direction = "wide")
  ID=2:ncol(Tot.reshaped)
  colnames(Tot.reshaped)=gsub(paste("INDIVIDUALS",".",sep=""), "", names(Tot.reshaped))
  Tot.reshaped$YEAR=as.numeric(as.character(Tot.reshaped$YEAR))
  Tot.reshaped[,ID]=Tot.reshaped[,ID]/rowSums(Tot.reshaped[,ID],na.rm=T)
  missn=YrS[which(!YrS%in%Tot.reshaped$YEAR)]
  ADD=matrix(NA,nrow=length(missn),ncol=ncol(Tot.reshaped))
  colnames(ADD)=colnames(Tot.reshaped)
  ADD=as.data.frame(ADD)
  ADD$YEAR=missn
  Tot.reshaped=rbind(Tot.reshaped,ADD)
  Tot.reshaped=Tot.reshaped[order(Tot.reshaped$YEAR),]
  DD=as.matrix(Tot.reshaped[,2:ncol(Tot.reshaped)])
  like.this=c("Comm.","Other","PJ","ER","BB.T","DM.T")
  DD=DD[,match(like.this,colnames(DD))]
  DD=t(DD)
  Xmax=ncol(DD)+5
  a=barplot(DD,plot=F)
  CLSs=c("black","grey90","grey30","grey80","grey50","grey70")
  #CLSs=grey.colors(nrow(DD))
  a=barplot(DD,xaxt='n',cex.axis=1.5,col=CLSs,xlim=c(0,a[length(a)]+1))
  LEN.nms=row.names(DD)
  LEN.nms=ifelse(LEN.nms=="BB.T","BB",ifelse(LEN.nms=="DM.T","DM",LEN.nms))
  legend(a[length(a)]+1,1,rev(LEN.nms), cex=1.15,
         fill=rev(CLSs),bty='n',title="A)")
  axis(1,a,F)
  box()
  
  #Predicted discard
  List=list(model=bycatch.model,data=Bycatch.glm)
  fun.plt.yr.pred.bycatch(d=List,normalised="NO",PredictorS=Predictors,MAIN="",
                          log.var="NO",Cx=1.125,YLIM=c(0,1),Cx.axs=1.35)
  legend(YrS[length(YrS)]+1,1,"", bty='n',title="B)")
  axis(1,a,YrS,cex.axis=1.5)
  
  mtext("Proportion",2,outer=T,cex=1.75,las=3,line=2.2)
  mtext("Year",1,outer=T,line=1.1,cex=1.75)
  dev.off()
  
  #Species table
  d$SPEC.tabl=with(d,ifelse(FATE=="C","Comm.",SPECIES))
  d$SPEC.tabl=with(d,ifelse(SPEC.tabl=="OS","BS",SPEC.tabl))
  TLB.sp=table(d$SPEC.tabl)
  TLB.sp=TLB.sp/sum(TLB.sp)
  TLB.sp=data.frame(SPECIES=names(TLB.sp),PROP=c(TLB.sp))
  #a=subset(SPECIES.names,Species%in%unique(TLB.sp$SPECIES))
  a=subset(SPECIES_PCS_FATE,SPECIES%in%unique(TLB.sp$SPECIES))
  TLB.sp=merge(TLB.sp,a,by="SPECIES",all.x=T)
  TLB.sp=TLB.sp[order(-TLB.sp$PROP),]
  TLB.sp$PROP=round(TLB.sp$PROP,3)
  TLB.sp$PROP=with(TLB.sp,ifelse(PROP<0.001,"<0.001",PROP))
  TLB.sp=TLB.sp[,match(c("SPECIES","PROP","COMMON_NAME","SCIENTIFIC_NAME"),names(TLB.sp))]
  write.csv(TLB.sp,paste("Table_species_disc_com",dat.nm,"csv",sep="."),row.names=F)
}

#Function for anova commercial
fn.anova.com=function(MODEL)
{
  TABL=vector('list',length(resp.vars))
  names(TABL)=resp.vars
  
  for(i in 1:n.rv)       #takes 60 seconds
  {
    Modl=MODEL[[i]]$model
    
    #each term
    Anova.tab=anova(Modl, test = "Chisq")
    n=2:length(Anova.tab$Deviance)
    Term.dev.exp=100*(Anova.tab$Deviance[n]/Modl$null.deviance)
    names(Term.dev.exp)=rownames(Anova.tab)[n]
    
    #nice table
    Anov.tab=as.data.frame.matrix(Anova.tab)
    Term.tab=data.frame(Percent.dev.exp=Term.dev.exp)
    Anova.tab=Anova.tab[-1,match(c("Deviance","Pr(>Chi)"),names(Anova.tab))]
    Anova.tab=cbind(Anova.tab,Term.tab)
    Anova.tab=Anova.tab[,-match("Deviance",names(Anova.tab))]
    Anova.tab$"Pr(>Chi)"=ifelse(Anova.tab$"Pr(>Chi)"<0.001,"<0.001",round(Anova.tab$"Pr(>Chi)",3))
    Total=Anova.tab[1,]
    Total$"Pr(>Chi)"=""
    Total$Percent.dev.exp=sum(Anova.tab$Percent.dev.exp)
    rownames(Total)="Total"
    Anova.tab=rbind(Anova.tab,Total)
    Anova.tab$Percent.dev.exp=round(Anova.tab$Percent.dev.exp,2)
    TABL[[i]]=Anova.tab
  }
  TABL=do.call(cbind,TABL)
  TABL=cbind(TERM=row.names(TABL),TABL)
  
  return(TABL)
}

#function for predicting commercial indices
fun.plt.yr.pred.com=function(d,normalised,PredictorS,log.var,ALL.YRS)
{
  #create new data
  Dat=d$data
  NEWDATA=data.frame(YEAR=sort(unique(as.character(Dat$YEAR))))                   
  NEWDATA$YEAR=factor(NEWDATA$YEAR,levels=levels(d$data$YEAR))
  Other.preds=PredictorS[-match(c("YEAR"),PredictorS)]
  d.frm=t(matrix(rep(NA,length(Other.preds))))
  colnames(d.frm)=Other.preds
  d.frm=as.data.frame(d.frm)
  for(q in 1:length(Other.preds))  
  {
    Vec=Dat[,match(Other.preds[q],names(Dat))]
    CLs=class(Vec)
    if(CLs=='factor') Value=factor(Mst.cmn(Vec),levels=levels(Vec))
    if(CLs=='numeric') Value=mean(Vec)
    d.frm[,q]=Value
  }
  NEWDATA=cbind(NEWDATA,d.frm)
  
  # Preds=predict(d$model,newdata=NEWDATA,type='response',se.fit=T)
  # 
  # if(log.var=="NO")PRD=data.frame(Year=NEWDATA$YEAR,Mean=Preds$fit,SE=Preds$se.fit)
  # if(log.var=="YES")PRD=data.frame(Year=NEWDATA$YEAR,Mean=exp(Preds$fit),SE=exp(Preds$se.fit))
  
  Covar.pos=as.matrix(vcov(d$model))
  set.seed(999);Pos.pars.rand=rmvnorm(niter,mean=coef(d$model),sigma=Covar.pos)    
  MC.preds=matrix(nrow=niter,ncol=length(NEWDATA$YEAR))
  
  for(n in 1:niter)
  {
    model=d$model
    model$coefficients=Pos.pars.rand[n,]
    a=predict(model,newdata=NEWDATA,type='response',se.fit=T)
    
    if(log.var=="YES") Pred=exp(a$fit+(a$se.fit^2)/2)  #apply bias correction for log transf
    if(log.var=="NO") Pred=a$fit
    
    MC.preds[n,]=Pred
  }
  
  PRD=data.frame(Year=NEWDATA$YEAR,
                 MEAN=colMeans(MC.preds,na.rm=T),
                 LOW=apply(MC.preds, 2, function(x) quantile(x, 0.025,na.rm=T)),
                 UP=apply(MC.preds, 2, function(x) quantile(x, 0.975,na.rm=T)))
  
  #standardise to a mean score of 1
  if(normalised=="YES") PRD=fn.relative(PRD)
  
  yr=as.numeric(as.character(PRD$Year))
  MeAn=PRD$MEAN
  UppCI=PRD$UP
  LowCI=PRD$LOW
  
  dat.plt=data.frame(yr=yr,MeAn=MeAn,UppCI=UppCI,LowCI=LowCI)
  
  #add missing years
  mis.yr=seq(ALL.YRS[1],ALL.YRS[length(ALL.YRS)])
  mis.yr=mis.yr[which(!mis.yr%in%yr)]
  if(length(mis.yr)>0)
  {
    dummy=dat.plt[1:length(mis.yr),]
    dummy[,]=NA
    dummy$yr=mis.yr
    dat.plt=rbind(dat.plt,dummy)
  }
  dat.plt=dat.plt[order(dat.plt$yr),]
  return(dat.plt)
}

#function for plotting commercial indices
plot.comm=function(dat.plt,MAIN,Cx,YLIM,Cx.axs)
{
  if(is.null(YLIM)) YLIM=c(min(dat.plt$LowCI,na.rm=T),max(dat.plt$UppCI,na.rm=T))
  with(dat.plt,plot(yr,MeAn,pch=19,main=MAIN,xlab="",ylab="",
                    cex=Cx,xaxt="n",cex.axis=1.25,cex.lab=1.5,cex.main=1.75,ylim=YLIM))
  with(dat.plt,arrows(x0=yr, y0=LowCI, x1=yr, y1=UppCI,code = 3,angle=90,length=.025))
  with(dat.plt,axis(1,yr,F,tck=-0.025))
  with(dat.plt,axis(1,seq(yr[1],yr[length(yr)],5),F,tck=-0.05))
}



# 3 Precedure section-----------------------------------------------------------------------

#3.1. Manipulate trophic levels
  #add SD from FishBase to Cortes'; use FishBase for Teleosts
SPECIES_PCS_FATE$TL_SD=with(SPECIES_PCS_FATE,ifelse(is.na(TL_SD)& NATURE%in%c("T"),TL_SD2,
            ifelse(is.na(TL_SD)& NATURE%in%c("S","R"),TROPHIC_LEVEL*TL_SD2/TROPHIC_LEVEL2,
            TL_SD)))
SPECIES_PCS_FATE$TROPHIC_LEVEL=with(SPECIES_PCS_FATE,
            ifelse(is.na(TROPHIC_LEVEL),TROPHIC_LEVEL2,TROPHIC_LEVEL))
SPECIES_PCS_FATE_Com$TL_SD=with(SPECIES_PCS_FATE_Com,ifelse(is.na(TL_SD)& NATURE%in%c("T"),TL_SD2,
            ifelse(is.na(TL_SD)& NATURE%in%c("S","R"),TROPHIC_LEVEL*TL_SD2/TROPHIC_LEVEL2,
            TL_SD)))            
SPECIES_PCS_FATE_Com$TROPHIC_LEVEL=with(SPECIES_PCS_FATE_Com,
            ifelse(is.na(TROPHIC_LEVEL),TROPHIC_LEVEL2,TROPHIC_LEVEL))

  #Manually fixing inconsistencies in trophic level
SPECIES_PCS_FATE$TROPHIC_LEVEL=with(SPECIES_PCS_FATE,ifelse(SPECIES=="TG",TROPHIC_LEVEL2,
            ifelse(SPECIES=="PD",4.1,TROPHIC_LEVEL)))
SPECIES_PCS_FATE_Com$TROPHIC_LEVEL=with(SPECIES_PCS_FATE_Com,
            ifelse(SPECIES==18022,4.5,TROPHIC_LEVEL))


#3.2. Manipulated WA shark observer data

  #Extract year and month
DATA$DATE=as.Date(DATA$date,format="%d/%m/%Y")
DATA$YEAR=as.numeric(strftime(DATA$DATE, format="%Y"))
DATA$MONTH=as.numeric(strftime(DATA$DATE, format="%m"))

  #Fix method and mesh size issues
DATA$Method=with(DATA,ifelse(is.na(Method) & BOAT%in%c("B67","F244","F517","E35"),"GN",Method))
DATA$MESH_SIZE=with(DATA,ifelse(is.na(MESH_SIZE) & BOAT=="B67" & YEAR>1997,"7",
                                ifelse(is.na(MESH_SIZE) & BOAT=="F517" & YEAR>2002,"7",MESH_SIZE)))

  #Select commercial gillnet and mesh size
DATA=subset(DATA,Method=="GN" & MESH_SIZE%in%c("6","6.5","7"))     
Research.vess=c("HAM","HOU","RV GANNET","RV BREAKSEA","NATT","NAT","FLIN","RV SNIPE 2")
DATA=subset(DATA,!BOAT%in%Research.vess)
ALL.yrs=sort(as.numeric(unique(DATA$YEAR)))

  #Get latitude and longitude from de first end of the net for location
DATA$LATITUDE=as.numeric(paste(DATA$END1LATD,".",ifelse(trunc(DATA$END1LATM*100/60)>=10,
                                                        paste(substr(DATA$END1LATM*100/60,1,2),substr(DATA$END1LATM*100/60,4,5),sep=""),
                                                        paste(0,substr(DATA$END1LATM*100/60,1,1),substr(DATA$END1LATM*100/60,3,4),sep="")),sep=""))
DATA$LONGITUDE=as.numeric(paste(DATA$END1LNGD,".",ifelse(trunc(DATA$END1LNGM*100/60)>=10,
                                                         paste(substr(DATA$END1LNGM*100/60,1,2),substr(DATA$END1LNGM*100/60,4,5),sep=""),
                                                         paste(0,substr(DATA$END1LNGM*100/60,1,1),substr(DATA$END1LNGM*100/60,3,4),sep="")),sep=""))
  #Data range
DATA=subset(DATA, LATITUDE<(-26) | LATITUDE==0) #zero to include the dodgy sheet numbers F00001, F00002 and F00003 with zero position

  #Fix sex
DATA$SEX=as.character(DATA$SEX)
DATA$SEX=with(DATA,ifelse(SEX%in%c("f","F"),"F",ifelse(SEX%in%c("m","M"),"M","U")))

  #Get TEPS interactions from Comments in boat.hdr  
Comments=subset(DATA,select=c(SHEET_NO,COMMENTS.hdr))
Comments=Comments[!is.na(Comments$COMMENTS.hdr),]
Comments=Comments[!duplicated(Comments$SHEET_NO),]

  #Add common and scientific names, nature, fate, PCS and trophic level
DATA=merge(DATA,subset(SPECIES_PCS_FATE,SPECIES%in%unique(DATA$SPECIES)),by="SPECIES",all.x=T)

  #Special treatment for White shark (WP) and Grey nurse shark (GN) as they became protected throughout the period
DATA$FATE=ifelse(DATA$SPECIES=="WP",ifelse(DATA$YEAR<1997,"C","D"), 
            ifelse(DATA$SPECIES=="GN",ifelse(DATA$YEAR<2001,"C","D"),as.character(DATA$FATE)))
              
DATA$NATURE=ifelse(DATA$SPECIES=="WP" | DATA$SPECIES=="GN",ifelse(DATA$FATE=="C","S","TEPS"),as.character(DATA$NATURE))           
                                               
  #remove unkwnown species codes
a=subset(DATA,is.na(FATE),selec=c(SPECIES,COMMON_NAME));a=a[!duplicated(a$SPECIES),]
NN.sp=a$SPECIES
DATA=subset(DATA,!SPECIES%in%NN.sp)

  #Add fishing zones (from Department of Fisheries WA)
DATA$LATITUDE=with(DATA,ifelse(LATITUDE==0,-as.numeric(substr(BLOCK,1,2)),LATITUDE))
DATA$LONGITUDE=with(DATA,ifelse(LONGITUDE==0,100+as.numeric(substr(BLOCK,3,4)),LONGITUDE))

DATA=subset(DATA,!BLOCK==0)
DATA$ZONE=as.character(with(DATA,ifelse(LONGITUDE>=116.5 & LATITUDE<=(-26),"Zone2",
        ifelse(LONGITUDE<116.5 & LATITUDE<=(-33),"Zone1",
        ifelse(LATITUDE>(-33) & LATITUDE<=(-26) & LONGITUDE<116.5,"West",
         ifelse(LATITUDE>(-26) & LONGITUDE<114,"Closed",
         ifelse(LATITUDE>(-26) & LONGITUDE>=114 & LONGITUDE<123.75,"North",
        ifelse(LATITUDE>(-26) & LONGITUDE>=123.75,"Joint",NA))))))))

  #Add regions (from Hall & Wise 2011)
DATA$REGION=as.character(with(DATA,ifelse(LONGITUDE>=124 & LONGITUDE<=129,"Region1",
       ifelse(LONGITUDE>=119 & LONGITUDE<124,"Region2",
       ifelse(LONGITUDE>=116 & LONGITUDE<119,"Region3",
       ifelse(LONGITUDE<116 & LATITUDE<=(-33),"Region4",
       ifelse(LATITUDE>(-33) & LATITUDE<=(-30),"Region5",
       ifelse(LATITUDE>(-30) & LATITUDE<=(-27),"Region6","Out.of.region"))))))))

  #Add area variable for temporal comparison
DATA$AREA=as.character(with(DATA,ifelse(BLOCK%in%c(2813,2814,2914),"Area1",
        ifelse(BLOCK%in%c(3114,3115,3214,3215),"Area2",
        ifelse(BLOCK%in%c(3314,3315,3414,3415),"Area3",
        ifelse(BLOCK%in%c(3322,3323,3324),"Area4",NA))))))

  #Add period variable for spatial comparison
DATA$SEASON=as.character(with(DATA,ifelse(MONTH==12 | MONTH<3,"Summer",
         ifelse(MONTH>=3 & MONTH<6,"Autumn",
         ifelse(MONTH>=6 & MONTH<9,"Winter",
         ifelse(MONTH>=9 & MONTH<12,"Spring",NA))))))
YR.selected=c(1995:1998,2002:2003)
DATA$Yr.dummy=with(DATA,ifelse(YEAR%in%YR.selected,YEAR,NA))
DATA$PER=as.character(with(DATA,paste(Yr.dummy,SEASON)))
DATA$PERIOD=with(DATA,ifelse(substr(PER,1,2)=="NA",NA,PER))
DATA=DATA[,-match(c("Yr.dummy","PER"),colnames(DATA))]

  #Add gillnet effort
if(use.soak.time=="YES") DATA$EFFORT=DATA$NET_LENGTH*DATA$SOAK_TIME
if(use.soak.time=="NO") DATA$EFFORT=DATA$NET_LENGTH

  #Remove tropical species that belong to the Northern Territory and the "other sharks and scale fish" group
DATA=subset(DATA, FATE!="A" & SPECIES!="XX")

  #Add number caught
DATA$INDIVIDUALS=1

  #Fix fork length by total length
DATA$FL=with(DATA,{ifelse(NATURE=="T",TL,FL)})#for teleosts copy TL to FL as is the only measure and we want the indeces with the TL 

  #Sharks and rays
Sp.len.dat=as.character(Len.cof$Species)
ID=with(DATA,which(COMMON_NAME%in%Sp.len.dat & is.na(FL) & !is.na(TL)))
DATA.fix=DATA[ID,]
DATA=DATA[-ID,]
DATA.fix=merge(DATA.fix,Len.cof[,-14],by.x="COMMON_NAME",by.y="Species",all.x=T)
DATA.fix$FL=with(DATA.fix,(TL-Intercept)/Slope)
DATA.fix=DATA.fix[,match(names(DATA),names(DATA.fix))]
DATA=rbind(DATA,DATA.fix)

  #Select relevant variables
DATA=subset(DATA,select=c(SHEET_NO,YEAR,MONTH,BLOCK,BOAT,SKIPPER,SPECIES,FATE,
                          COMMON_NAME,NATURE,TL,FL,SEX,TROPHIC_LEVEL,TL_SD,BOTDEPTH,
                          MESH_SIZE,NET_LENGTH,SOAK_TIME,EFFORT,LATITUDE,LONGITUDE,
                          ZONE,REGION,BIOREGION,AREA,SEASON,PERIOD,INDIVIDUALS))

#Define attributes to variables
DATA$YEAR=as.character(DATA$YEAR)
DATA$MONTH=as.character(DATA$MONTH)
DATA$BLOCK=as.character(DATA$BLOCK)


# 4 Ecosystems indicators analysis-----------------------------------------------------------------------

  #--- 4.1 WA Fisheries observer data---

# 4.1.1 Preliminary analyses
#No of sheet number per year
surveys=data.frame(YEAR=as.factor(""), SURVEYS=as.numeric(0))#dummy
for(i in unique(DATA$YEAR))
   {
      year=subset(DATA, YEAR==i)
      n=length(unique(year$SHEET_NO))
      surveys=rbind(surveys, data.frame(YEAR=as.factor(i), SURVEYS=n))
   }
surveys=surveys[-1,]#deleting dummy row
surveys=surveys[order(as.numeric(as.character(surveys$YEAR))),]#sorting by year

barplot(surveys$SURVEYS,ylim=c(0,500),main="No of sheet # per year",space=0.5)
text(seq(1,25.5,1.5),surveys$SURVEYS+20,surveys$YEAR,srt=45,cex=0.75)


#No of sheet # by Area-Year-Season
AREAS=unique(DATA$AREA)
AREAS=sort(AREAS)
AREAS=AREAS[!is.na(AREAS)]
Num.sheet.area.yr.sea=with(subset(DATA,!is.na(AREA)),table(AREA,YEAR,SEASON))

Yrs_by_Area=vector('list',length(AREAS))
names(Yrs_by_Area)=AREAS
for(a in 1:length(AREAS))
{
  dd=subset(DATA,AREA==AREAS[a])
  Yy=sort(unique(dd$YEAR))
  stor=rep(NA,length(Yy))
  names(stor)=Yy
  for(i in 1:length(Yy))
  {
    year=subset(dd, YEAR==Yy[i])
    stor[i]=length(unique(year$SHEET_NO))
  }
  Yrs_by_Area[[a]]=stor
}


#No of sheet # per block
surveys=data.frame(BLOCK=as.factor(""), SURVEYS=as.numeric(0))#dummy
for(i in unique(DATA$BLOCK))
   {
      block=subset(DATA, BLOCK==i)
      n=length(unique(block$SHEET_NO))
      surveys=rbind(surveys, data.frame(BLOCK=as.factor(i), SURVEYS=n))
   }
surveys=surveys[-1,]#deleting dummy row

surveys$LAT_BLOCK=-as.numeric(substr(surveys$BLOCK,1,2))
surveys$LONG_BLOCK=as.numeric(substr(surveys$BLOCK,3,4))+100

ocean.pal=colorRampPalette(brewer.pal(n=9,'YlOrRd'))(max(surveys$SURVEYS))#colors
map("worldHires",xlim=c(112.95, 129.5),ylim=c(-36, -26.5),col="grey80",fill=F)
for(i in seq(surveys$SURVEYS))
   {
      rect(surveys$LONG_BLOCK[i],surveys$LAT_BLOCK[i]-1,surveys$LONG_BLOCK[i]+1,surveys$LAT_BLOCK[i],border="blue",col=ocean.pal[surveys$SURVEYS[i]])
   }
map("worldHires",xlim=c(112.95, 129.5),ylim=c(-36, -26.5),col="grey80",fill=T,border="grey80",add=T)
rect(surveys$LONG_BLOCK,surveys$LAT_BLOCK-1,surveys$LONG_BLOCK+1,surveys$LAT_BLOCK,border="blue")
mtext("No of sheet # per block",line=1,cex=1.5)
text(surveys$LONG_BLOCK+0.5,surveys$LAT_BLOCK-0.5,surveys$SURVEYS,cex=0.8)
text(127.5,-27,paste("n =",sum(surveys[2])))

#boats and skippers
A=table(as.character(DATA$BOAT),as.character(DATA$SKIPPER),useNA = "ifany")
A[A>0]=1
# a=subset(DATA,is.na(SKIPPER))
# table(a$BOAT,a$YEAR)
# b=subset(DATA,(BOAT=="F505"))
# table(b$YEAR,b$SKIPPER)

DATA$SKIPPER=with(DATA,ifelse(is.na(SKIPPER) & YEAR%in%c(2001:2002),"C. GULLOTTI",SKIPPER))


#not using marine mammals
DATA=subset(DATA,!COMMON_NAME=="Whale") 


#Define Species as factor
DATA$SPECIES=as.factor(DATA$SPECIES)


aa=DATA
aa=aa[!duplicated(aa$SHEET_NO),]
YR.ARE.sht=table(aa$YEAR,aa$AREA)

YR.ARE.sht[YR.ARE.sht<Min.shts]=0
YR.ARE.sht[YR.ARE.sht>0]=1
YR.ARE.sht=as.numeric(rownames(YR.ARE.sht))*YR.ARE.sht
YR.ARE.sht=as.data.frame.matrix(YR.ARE.sht)

#keep records from shots according to criteria
d=DATA[!duplicated(DATA$SHEET_NO),]
d$N=1
d1=aggregate(N~BLOCK+YEAR,d,sum)
d1=subset(d1,!BLOCK=="0")
N.base.case=subset(d1,N>Min.shts)
N.sens.1=subset(d1,N>Min.shts.sens[1])
N.sens.2=subset(d1,N>Min.shts.sens[2])

N.base.case$Yr.blk=with(N.base.case,paste(YEAR,BLOCK,sep="_"))
N.sens.1$Yr.blk=with(N.sens.1,paste(YEAR,BLOCK,sep="_"))
N.sens.2$Yr.blk=with(N.sens.2,paste(YEAR,BLOCK,sep="_"))

DATA$Yr.blk=with(DATA,paste(YEAR,BLOCK,sep="_"))


DATA=subset(DATA,!BOAT=="E7")   #cannot estimate coef for this, few obs

DATA.base.case=subset(DATA,Yr.blk%in%unique(N.base.case$Yr.blk))
DATA.sens.1=subset(DATA,Yr.blk%in%unique(N.sens.1$Yr.blk))
DATA.sens.2=subset(DATA,Yr.blk%in%unique(N.sens.2$Yr.blk))


#Loop over each data seta
setwd(paste(getwd(),"Outputs",sep="/"))
WD=getwd()

resp.vars=c("Shannon","Simpson","MTL","MML")
#resp.vars=c("Shannon","Margalef","Pielou","Simpson","MTL","MML")
Res.var.in.log=rep("NO",length(resp.vars))    #fit response var in log space or not?
n.rv=length(resp.vars)
Main.title=resp.vars


#each data sets takes .3 seconds per iteration
system.time({OUT.base.case=fn.loop.over.obsrvr.data(DaTA=DATA.base.case,dat.nm="BC",normalised="YES",Drop.yrs="NO")})

  #sensitivity tests
setwd(paste(getwd(),"/Sensitivity",sep="/"))
OUT.sens.1=fn.loop.over.obsrvr.data(DaTA=DATA.sens.1,dat.nm="S1",normalised="YES",Drop.yrs="NO")
OUT.sens.2=fn.loop.over.obsrvr.data(DaTA=DATA.sens.2,dat.nm="S2",normalised="YES",Drop.yrs="NO")

setwd(WD)


  # 4.2 WA Fisheries commercial data

# 4.2.1 Preliminary analyses
if(do.exploratory=="YES")
{
  #No of shots per year-block-month
  a=Data.monthly[!duplicated(Data.monthly$Same.return),]
  a$N=1
  a=aggregate(N~YEAR+MONTH+BLOCK+LATITUDE+LONGITUDE,a,sum)
  
  yrs=sort(unique(a$YEAR))
  for(i in 1:length(yrs)) fn.see(d=subset(a,YEAR==yrs[i]))
  rm(a)
}

#remove records from shots with less than Min.shts per year-block
d=Data.monthly[!duplicated(Data.monthly$Same.return),]
d$N=1
d1=aggregate(N~BLOCKX+FINYEAR,d,sum)
N.base.case=subset(d1,N>Min.shts)

N.base.case$Yr.blk=with(N.base.case,paste(FINYEAR,BLOCKX,sep="_"))
Data.monthly$Yr.blk=with(Data.monthly,paste(FINYEAR,BLOCKX,sep="_"))

Data.monthly=subset(Data.monthly,Yr.blk%in%N.base.case$Yr.blk)


#Remove boats with less than 10 records
AA=sort(table(Data.monthly$BOAT))
AA=AA[AA>Min.recs]
Data.monthly=subset(Data.monthly,BOAT%in%names(AA))

#Remove years with less than 10 records
AA=sort(table(Data.monthly$FINYEAR))
AA=AA[AA>Min.recs]
Data.monthly=subset(Data.monthly,FINYEAR%in%names(AA))


#add species trophic level,etc
Data.Mon=merge(Data.monthly,SPECIES_PCS_FATE_Com,by="SPECIES",all.x=T,all.y=F)

Data.monthly=subset(Data.Mon,select=c(SHEET_NO,YEAR,MONTH,BLOCK,BOAT,TYPE.DATA,SPECIES,
                                      SCIENTIFIC_NAME,NATURE,LIVEWT.c,TROPHIC_LEVEL,TL_SD2,EFFORT,
                                      LATITUDE,LONGITUDE,ZONE,BIOREGION))
Data.monthly=subset(Data.monthly,!is.na(Data.monthly$TROPHIC_LEVEL))
Data.monthly$SPECIES=as.factor(Data.monthly$SPECIES)
Data.monthly$INDIVIDUALS=1
Data.monthly$FINYEAR=Data.monthly$YEAR
Data.monthly$YEAR=as.numeric(substr(Data.monthly$YEAR,1,4))

#Export species by year
TBLA=table(paste(Data.monthly$SCIENTIFIC_NAME,Data.monthly$TROPHIC_LEVEL),Data.monthly$YEAR)
write.csv(TBLA,"TBLA_species.TL_year.csv",row.names=T)

TBLA=table(Data.monthly$SCIENTIFIC_NAME,Data.monthly$YEAR)
write.csv(TBLA,"TBLA_species_year.csv",row.names=T)

TBLA=log(TBLA+0.0001)
CL=heat.colors(nrow(TBLA))
plot(as.numeric(colnames(TBLA)),TBLA[1,],ylab="Number of records (logged)",pch=19,col=CL[1],ylim=c(0,max(TBLA)))
for(n in 2:nrow(TBLA)) points(as.numeric(colnames(TBLA)),TBLA[n,],pch=19,col=CL[n])

TBLA=table(Data.monthly$TROPHIC_LEVEL,Data.monthly$YEAR)
#TBLA=log(TBLA)
colfunc <- colorRampPalette(c("red", "yellow"))
#CL <-colfunc(nrow(TBLA))
CL=rep("grey40",nrow(TBLA))
xx=as.numeric(colnames(TBLA))
mltplr=5

fn.fig("TrophicLevel_by_yr",2000,2000) 
par(las=1,mar=c(1.75,3.5,1.5,.1),oma=c(1.5,1,.1,.1),mgp=c(1,.8,0),cex.lab=1.25)
plot(xx,rep(1,length(xx)),cex=(TBLA[1,]/max(TBLA))*mltplr,ylab="",pch=19,col=CL[1],xlab="",yaxt='n',ylim=c(0,nrow(TBLA)))
for(n in 2:nrow(TBLA)) points(xx,rep(n,length(xx)),cex=(TBLA[n,]/max(TBLA))*mltplr,pch=19,col=CL[n])
axis(2,1:nrow(TBLA),row.names(TBLA))
mtext("Trophic level",2,2.5,cex=1.5,las=3)
mtext("Year",1,2,cex=1.5)
dev.off()

Start.yr=1989   #before 1989 most species not reported in commercial logbooks


#calculate indices for each shot
d.anlsys=subset(Data.monthly,YEAR>=Start.yr)
SHOTS=unique(d.anlsys$SHEET_NO)
Store.shots=vector('list',length(SHOTS))
names(Store.shots)=SHOTS
system.time(for(i in 1:length(SHOTS))       #takes 102 seconds
{
  dat=subset(d.anlsys,SHEET_NO==SHOTS[i])
  dat.indx=dat[!duplicated(dat$SHEET_NO),]
  dat.indx=subset(dat.indx,select=c(SHEET_NO,YEAR,MONTH,ZONE,BLOCK,BOAT,
                                    EFFORT,BIOREGION))
  #d=table(dat$SPECIES)                     #use number of species
  d=tapply(dat$LIVEWT.c, dat$SPECIES, sum)  #use catch weight     
  d[is.na(d)]=0
  
  #Diversity indicators
  Div.InDX=data.frame(
    Shannon = diversity(d, index = "shannon"),
    Simpson = diversity(d, index = "simpson")
  )
  dat.indx=cbind(dat.indx,Div.InDX)
  
  #Ecosystem indicators 
  dd=d
  dd=subset(dd,dd>0)
  t_level=aggregate(TROPHIC_LEVEL ~ SPECIES, FUN=mean, data=dat)
  fl_max=NULL
  Eco.InDX=Eco.ind.shot(data=dd,TROPHIC.LEVEL=t_level,MAX.BODY.LENGTH=NULL,TROPHIC.EFFICIENCY = 0.1)
  Eco.InDX=as.data.frame(do.call(cbind,Eco.InDX))
  dat.indx=cbind(dat.indx,Eco.InDX)
  Store.shots[[i]]=dat.indx
})

DATA.shots.diversity.commercial=do.call(rbind,Store.shots)


DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,!is.na(Shannon))
#DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,!is.na(Margalef))
#DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,!is.na(Pielou))
DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,!is.na(Simpson))
DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,!is.na(MTL))

DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,Shannon>0)
#DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,Margalef>0)
#DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,Pielou>0)
DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,Simpson>0)
DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,MTL>0)


# Q999=quantile(DATA.shots.diversity.commercial$Shannon,probs=0.9999)
# DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,Shannon<=Q999)
# 
# Q999=quantile(DATA.shots.diversity.commercial$Margalef,probs=0.9999)
# DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,Margalef<=Q999)
# 
# Q999=quantile(DATA.shots.diversity.commercial$Pielou,probs=0.9999)
# DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,Pielou<=Q999)
# 
# Q999=quantile(DATA.shots.diversity.commercial$Simpson,probs=0.9999)
# DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,Simpson<=Q999)
# 
# Q999=quantile(DATA.shots.diversity.commercial$MTL,probs=0.9999)
# DATA.shots.diversity.commercial=subset(DATA.shots.diversity.commercial,MTL<=Q999)



# 3.2.3 Modelling
#note: Margalef removed because not compatible with weight data
resp.vars=c("Shannon","Simpson","MTL")
n.rv=length(resp.vars)
Predictors=c("YEAR","BOAT","BLOCK","MONTH")
Expl.varS=Predictors
FactoRS=Expl.varS

Store.mod.out.commercial=vector('list',length(resp.vars))
names(Store.mod.out.commercial)=resp.vars

Res.var.in.log=rep("NO",length(resp.vars))
Yrs_commercial=sort(unique(DATA.shots.diversity.commercial$YEAR))
Separate.monthly_daily="NO"
if(Separate.monthly_daily=="YES")
{
  Yrs_monthly=Yrs_commercial[1]:2005
  Yrs_daily=2006:Yrs_commercial[length(Yrs_commercial)]
  Store.mod.out.commercial.daily=Store.mod.out.commercial
}



# 5 glm approach-------------------------------------------------------------------------
#takes 11 seconds
               
if(Separate.monthly_daily=="NO")
{
  system.time(for(i in 1:n.rv)
  {
    #split data into montlhy and daily       
    DDD=DATA.shots.diversity.commercial
    
    DDD$YEAR=as.factor(DDD$YEAR)
    DDD$BLOCK=as.factor(DDD$BLOCK)
    DDD$MONTH=as.factor(DDD$MONTH)
    DDD$BOAT=as.factor(DDD$BOAT)
    DDD$ZONE=as.factor(DDD$ZONE)

    Store.mod.out.commercial[[i]]=Mod.fn.glm(d=DDD,
                                             ResVar=resp.vars[i],Expl.vars=Expl.varS,
                                             Predictrs=Predictors,
                                             FactoRs=FactoRS,OFFSET=OFFSETT,
                                             log.var=Res.var.in.log[i],add.inter="NO",
                                             MixedEff=MixedEff)
    
  })
}
if(Separate.monthly_daily=="YES")
{
  system.time(for(i in 1:n.rv)
  {
    #split data into montlhy and daily       
    DDD_comm=subset(DATA.shots.diversity.commercial,YEAR%in%Yrs_monthly)
    
    DDD_comm$YEAR=as.factor(DDD_comm$YEAR)
    DDD_comm$MONTH=as.factor(DDD_comm$MONTH)
    DDD_comm$BOAT=as.factor(DDD_comm$BOAT)
    DDD_comm$ZONE=as.factor(DDD_comm$ZONE)
    
    Store.mod.out.commercial[[i]]=Mod.fn.glm(d=DDD_comm,
                                             ResVar=resp.vars[i],Expl.vars=Expl.varS,
                                             Predictrs=Predictors,
                                             FactoRs=FactoRS,OFFSET="offset(log.EFFORT)",
                                             log.var=Res.var.in.log[i],add.inter="NO")
    
    
    DDD_daily=subset(DATA.shots.diversity.commercial,YEAR%in%Yrs_daily)
    
    DDD_daily$YEAR=as.factor(DDD_daily$YEAR)
    DDD_daily$MONTH=as.factor(DDD_daily$MONTH)
    DDD_daily$BOAT=as.factor(DDD_daily$BOAT)
    DDD_daily$ZONE=as.factor(DDD_daily$ZONE)
    
    Store.mod.out.commercial.daily[[i]]=Mod.fn.glm(d=DDD_daily,
                                                   ResVar=resp.vars[i],Expl.vars=Expl.varS,
                                                   Predictrs=Predictors,
                                                   FactoRs=FactoRS,OFFSET="offset(log.EFFORT)",
                                                   log.var=Res.var.in.log[i],add.inter="NO")
    
    
  }) 
}

# 6 data mining-------------------------------------------------------------------------
# system.time(for(i in 1:n.rv)Store.mod.out.commercial[[i]]=Mod.fn.mining(d=DATA.shots.diversity.commercial,
#             ResVar=resp.vars[i],Predictrs=Predictors,Y.type="Continuous",Prop.train=.7,ALL.models="NO",nboot=1))

  #Get anova table
TABL.monthly=fn.anova.com(MODEL=Store.mod.out.commercial)
if(Separate.monthly_daily=="YES")TABL.daily=fn.anova.com(MODEL=Store.mod.out.commercial.daily)


#export anova tables as word doc
Export.tbl(WD=getwd(),Tbl=TABL.monthly,Doc.nm="Anova.table.commercial",caption=NA,paragph=NA,
           HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
           Zebra='NO',Zebra.col='grey60',Grid.col='black',
           Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
           HDR.names=c('TERM', resp.vars),
           HDR.span=c(1,rep(2,length(resp.vars))),
           HDR.2nd=c("",rep(c("P","% dev. expl."),length(resp.vars))))

if(Separate.monthly_daily=="YES")Export.tbl(WD=getwd(),Tbl=TABL.daily,Doc.nm="Anova.table.commercial_daily",caption=NA,paragph=NA,
           HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
           Zebra='NO',Zebra.col='grey60',Grid.col='black',
           Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
           HDR.names=c('TERM', resp.vars),
           HDR.span=c(1,rep(2,length(resp.vars))),
           HDR.2nd=c("",rep(c("P","% dev. expl."),length(resp.vars))))


#Plot diversity and ecosystem indices by zone
zoN=sort(unique(DATA.shots.diversity.commercial$ZONE))
YrSs=sort(unique(DATA.shots.diversity.commercial$YEAR))

# normalisation done on each data set separately   
fn.fig("Div & Eco indices.commercial_normalised",1200,2400)  #takes 5 sec per iteration
smart.par(n.plots=length(resp.vars),MAR=c(1.75,3.5,1.5,.1),OMA=c(3,1,.1,2),MGP=c(1,.8,0))
Main.title=resp.vars 

  #loop over each index   
system.time(for(i in 1:length(Store.mod.out.commercial))    #takes 1 sec per niter
{

  if(Separate.monthly_daily=="NO")
  {
    MN=fun.plt.yr.pred.com(d=Store.mod.out.commercial[[i]],
                           normalised="YES",PredictorS=Predictors,log.var=Res.var.in.log[i],
                           ALL.YRS=YrSs)
    
    plot.comm(dat.plt=MN,MAIN=Main.title[i],Cx=1.125,YLIM=NULL,Cx.axs=1.35)
    
     axis(1,seq(Yrs_commercial[1],Yrs_commercial[length(Yrs_commercial)],5)
                        ,seq(Yrs_commercial[1],Yrs_commercial[length(Yrs_commercial)],5),tck=-0.05,cex.axis=1.15)
  }
  if(Separate.monthly_daily=="YES")
  {
    MN=fun.plt.yr.pred.com(d=Store.mod.out.commercial[[i]],
                           normalised="YES",PredictorS=Predictors,log.var=Res.var.in.log[i],
                           ALL.YRS=Yrs_monthly)
    DAY=fun.plt.yr.pred.com(d=Store.mod.out.commercial.daily[[i]],
                            normalised="YES",PredictorS=Predictors,log.var=Res.var.in.log[i],
                            ALL.YRS=Yrs_monthly)
    plot.comm(dat.plt=rbind(MN,DAY),MAIN=Main.title[i],
              Cx=1.125,YLIM=NULL,Cx.axs=1.35)
    
    #highlight daily records
    polygon(x=c(2006,YrSs[length(YrSs)]+1,YrSs[length(YrSs)]+1,2006),
            y=c(-10,-10,10,10),col=rgb(.1,.1,.1,alpha=0.25),border=rgb(.1,.1,.1,alpha=0.25))
    
    if(i%in%c(3,4))axis(1,seq(Yrs_commercial[1],Yrs_commercial[length(Yrs_commercial)],5)
                        ,seq(Yrs_commercial[1],Yrs_commercial[length(Yrs_commercial)],5),tck=-0.05,cex.axis=1.15)
    
  }
  })
mtext("Year",1,outer=T,cex=1.5,line=1.5)
mtext("Relative value",2,-.6,outer=T,cex=1.5,las=3)
dev.off()

#ACA
# 7 Multivariate analysis-----------------------------------------------------------------------
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Multivariate_statistics.R"))

  #7.1  Observer data
DataSets=c("proportion")   #response variables
Predictors=c("YEAR","BLOCK","BOAT","MONTH","BOTDEPTH")
IDVAR=Predictors
Prop.sp=table(DATA.base.case$SPECIES)
Prop.sp=Prop.sp/sum(Prop.sp)
Grp.sp=names(100*Prop.sp[100*Prop.sp<0.1]) #group species occurring less than 0.1%, otherwise multivars fail
Numbers.block.year.mn.vesl=DATA.base.case
Numbers.block.year.mn.vesl$SPECIES=as.character(Numbers.block.year.mn.vesl$SPECIES)
Numbers.block.year.mn.vesl$SPECIES=with(Numbers.block.year.mn.vesl,
  ifelse(SPECIES%in%Grp.sp,"Other",SPECIES))   
Numbers.block.year.mn.vesl$SPECIES=factor(Numbers.block.year.mn.vesl$SPECIES)

Numbers.block.year.mn.vesl=aggregate(formula(paste(ResVar,paste(c(MultiVar,Predictors),collapse="+"),sep="~")),
                                     Numbers.block.year.mn.vesl,sum)
Numbers.block.year.mn.vesl$YEAR=factor(Numbers.block.year.mn.vesl$YEAR)
Numbers.block.year.mn.vesl$BLOCK=factor(Numbers.block.year.mn.vesl$BLOCK)
Numbers.block.year.mn.vesl$BOAT=factor(Numbers.block.year.mn.vesl$BOAT)
Numbers.block.year.mn.vesl$MONTH=factor(Numbers.block.year.mn.vesl$MONTH)

STore.multi.var.observer=Multivar.fn(DATA=Numbers.block.year.mn.vesl,ResVar=ResVar,MultiVar=MultiVar,
                                   Predictors=Predictors,IDVAR=IDVAR,
                                   Formula=formula("d.res.var~."),DataSets=DataSets)

Hndl=handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Multivariate")
for(s in 1:length(STore.multi.var.observer)) 
{
  with(STore.multi.var.observer[[s]],
       {
         fn.display.multivar(d=d,IDVAR=IDVAR,Predictors=Predictors,MDS=MDS,
                             Permanova.table=Permanova.table,
                             permanova.pairwise=permanova.pairwise,
                             Simper=Simper,NM=names(STore.multi.var.observer)[s],
                             hndl=paste(Hndl,"Observer/",sep="/"),cexMDS=.75)
       })
}

  #7.2 Commercial
fn.min.obs.ef=function(d,subset.sp)
{
  N.samples=length(unique(d$SHEET_NO))
  d=d%>%mutate(ktch=1,
               Eff.breaks=cut(EFFORT,breaks=50))
  Sp.occ=d%>%group_by(SPECIES)%>%
              summarise(Tot=sum(ktch))%>%
              mutate(Occ=100*Tot/N.samples)%>%
              filter(Occ>subset.sp)%>%pull(SPECIES)
  
  b1=d%>%filter(SPECIES%in%Sp.occ)%>%
    group_by(SHEET_NO,Eff.breaks)%>%
    summarise(N.species=sum(ktch))
  
  boxplot(N.species~Eff.breaks,b1)
  
  #ACA: Need to extract the min effort when it asymtotes.. can I automate this??
    # min.effort=
  
  return(min.effort)
  
}
Min.occ=0  #in 100%, if set to 0, then all species used

    #7.2.1 Monthly                                    AM I USING MONTHLY????????????

    #7.2.2 Daily
Min.obs.effort=fn.min.obs.ef(d=Data.daily,subset.sp=Min.occ)

#remove records with too low effort          #note: this should be move up front and    MISSING!!!!
                                             # all analysis done on this data subset
Data.daily=Data.daily%>%filter(EFFORT>=Min.obs.effort)


#ACA implement PRIMER stuff....

names(Data.daily)



Predictors=c("YEAR","BLOCK","BOAT","MONTH")
IDVAR=Predictors
Prop.sp=table(d.anlsys$SPECIES)
Prop.sp=Prop.sp/sum(Prop.sp)
Grp.sp=names(100*Prop.sp[100*Prop.sp<0.1])  #grouped species occurring less than 0.1%

Prop.sp=table(DATA.base.case$SPECIES)
Prop.sp=Prop.sp/sum(Prop.sp)
Grp.sp=names(100*Prop.sp[100*Prop.sp<0.1]) #group species occurring less than 0.1%, otherwise multivars fail
Kg.block.year.mn.vesl=d.anlsys
Kg.block.year.mn.vesl$SPECIES=as.character(Kg.block.year.mn.vesl$SPECIES)
Kg.block.year.mn.vesl$SPECIES=with(Kg.block.year.mn.vesl,
                                        ifelse(SPECIES%in%Grp.sp,"Other",SPECIES))   
Kg.block.year.mn.vesl$SPECIES=factor(Kg.block.year.mn.vesl$SPECIES)

Kg.block.year.mn.vesl=aggregate(formula(paste(ResVar,paste(c(MultiVar,Predictors),collapse="+"),sep="~")),
                                Kg.block.year.mn.vesl,sum)
Kg.block.year.mn.vesl$YEAR=factor(Kg.block.year.mn.vesl$YEAR)
Kg.block.year.mn.vesl$BLOCK=factor(Kg.block.year.mn.vesl$BLOCK)
Kg.block.year.mn.vesl$BOAT=factor(Kg.block.year.mn.vesl$BOAT)
Kg.block.year.mn.vesl$MONTH=factor(Kg.block.year.mn.vesl$MONTH)

STore.multi.var.commercial=Multivar.fn(DATA=Kg.block.year.mn.vesl,ResVar="LIVEWT.c",MultiVar=MultiVar,
                                     Predictors=Predictors,IDVAR=IDVAR,
                                     Formula=formula("d.res.var~."),DataSets=DataSets)

for(s in 1:length(STore.multi.var.observer)) 
{
  with(STore.multi.var.commercial[[s]],
       {
         fn.display.multivar(d=d,IDVAR=IDVAR,Predictors=Predictors,MDS=MDS,
                             Permanova.table=Permanova.table,
                             permanova.pairwise=permanova.pairwise,
                             Simper=Simper,NM=names(STore.multi.var.observer)[s],
                             hndl=paste(Hndl,"Commercial/",sep="/"),cexMDS=.75)
       })
}