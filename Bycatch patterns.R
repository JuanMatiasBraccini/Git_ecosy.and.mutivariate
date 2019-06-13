
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

#Define working directory

Folder="Matias"
#Folder="Agustin"

if(Folder=="Matias") setwd("C:/Matias/Analyses/Ecosystem indices/Shark-bycatch")
if(Folder=="Agustin") setwd("C:/Users/Agustín/Desktop/Shark-bycatch")

if(Folder=="Matias") User="Matias"

#Source functions
source("Git_bycatch_TDGDLF/Ecosystem_functions.R")


#choose if doing .jpeg or .tiff figures
Do.jpeg="YES"
Do.tiff="NO"

if(User=="Matias")source("C:/Matias/Analyses/SOURCE_SCRIPTS/Smart_par.R")

##################
###DATA SECTION###
##################

# 1. Bring in WA shark observer data
if(User=="Matias")
{
  source("C:/Matias/Analyses/SOURCE_SCRIPTS/Source_Shark_bio.R")
  rm(DATA)
  DATA=DATA.ecosystems
}

if(Folder=="Agustin") DATA=read.csv("WA_Agustin.csv",stringsAsFactors=F)

if(Folder=="Matias") setwd("C:/Matias/Analyses/Ecosystem indices/Shark-bycatch")

# 2. Bring in MAFFRI shark gillnet data
#DATA.MAFFRI=read.csv("MAFFRI_data.csv")
#DATA_TEPS.MAFFRI=read.csv("MAFFRI_TEPS_data.csv")


# 3. Bring in WA Species names + PCS + FATE
SPECIES_PCS_FATE=read.csv("SPECIES+PCS+FATE.csv",stringsAsFactors=F)
SPECIES_PCS_FATE_Com=read.csv("SPECIES+PCS+FATE_Com.csv",stringsAsFactors=F)

# 4. Bring in Length coefficients
Len.cof=read.csv("Raw coefficents table.csv")


#5. Bring in Commercial data
if(Folder=="Matias") source("C:/Matias/Analyses/Ecosystem indices/Shark-bycatch/Git_bycatch_TDGDLF/Commercial_data_for_Ecosystem_Analysis.R")
if(Folder=="Agustin") Data.monthly=read.csv("Data.monthly.csv",stringsAsFactors=F)


#CONTROL

Min.shts=10 #USe records with at least 10 shots per year-block
Min.shts.sens=c(5,20)

Min.recs=5 #minimum number of records

Min.individuals=5   #minimum number of individuals per shot to use

MixedEff=NA
#MixedEff="BOAT"

do.exploratory="NO"  #choose if doing data exploration

niter=100    #number of iterations for MC procedure for confidence intervals

#Percent.show.bycatch=0.8   #proportion of bycatch explained
#Prop.TC.bycatch=0.05       #proportion of catch explained

use.soak.time="NO"  #define if using soak time in the calculation of effort


ResVar="INDIVIDUALS"
MultiVar="SPECIES"
IDVAR=c("SHEET_NO","YEAR","MONTH","BOTDEPTH","BLOCK","ZONE","BOAT","SKIPPER","MESH_SIZE","BIOREGION","SEASON","EFFORT")   
Predictors=c("YEAR","BLOCK","BOAT","MONTH","BOTDEPTH")
Expl.varS=c("YEAR","BOAT","MONTH","BLOCK","BOTDEPTH")
FactoRS=Expl.varS[-match(c("BOTDEPTH"),Expl.varS)]
OFFSETT=NA
#OFFSETT="offset(log.EFFORT)"

##########
#Functions
##########

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



#######################
###PROCEDURE SECTION###
#######################

#1. Manipulated WA shark observer data

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

#Manipulate trophic levels
#note: add SD from FishBase to Cortes'; use FishBase for Teleosts
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

#Add bioregions (from Department of Fisheries WA)
DATA$BIOREGION=ifelse(DATA$REGION=="Region1"|DATA$REGION=="Region2"|DATA$REGION=="Region3","SouthCoast","WestCoast")

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


#Add gillnet absolute effort
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

#Output
#write.csv(DATA, "DATA.csv", row.names=F)


# 2. Ecosystems indicators analysis

  #--- 2.1 WA Fisheries observer data---


# 2.1.1 Preliminary analyses
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


  # 2.2 WA Fisheries commercial data

# 2.2.1 Preliminary analyses

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



#glm approach               #takes 11 seconds
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


#data mining
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



# 3. Multivariate analysis          #ACA
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Multivariate_statistics.R")
  #3.1  Logbook
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

Hndl="C:/Matias/Analyses/Ecosystem indices/Shark-bycatch/Outputs/Multivariate"
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

  #3.2 Commercial
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





#NOTE USED
# if(Folder=="Agustin")
# {
#   #Communities diversity throughout years
#   years=table(DATA1$YEAR,DATA1$SPECIES)
#   Div.InDX.yrs=Div.ind(data=years)
#   par(mfcol=c(2,2))
#   with(Div.InDX.yrs,{
#     fun.plt.Indx(Var_names,Shannon,"Shannon-Wiener diversity index","YEAR",Var_names,c(1, max(Shannon) + 0.5))
#     fun.plt.Indx(Var_names,Pielou,"Pielou's evenness","YEAR",Var_names,c(0.2, max(Pielou) + 0.2))
#     fun.plt.Indx(Var_names,Simpson,"Simpson's index","YEAR",Var_names,c(0.45, max(Simpson) + 0.2))
#     fun.plt.Indx(Var_names,Margalef,"Margalef species richness","YEAR",Var_names,c(9, max(Margalef) + 1))  
#   })
#   
#   #########################################
#   #Communities diversity throughout regions
#   regions=table(DATA1$REGION,DATA1$SPECIES)
#   Div.InDX.reg=Div.ind(data=regions)
#   par(mfcol=c(2,2))
#   with(Div.InDX.reg,{
#     fun.plt.Indx(1:length(Var_names),Shannon,"Shannon-Wiener diversity index","REGION",Var_names,c(1, max(Shannon) + 1))
#     fun.plt.Indx(1:length(Var_names),Pielou,"Pielou's evenness","REGION",Var_names,c(0, max(Pielou) + 0.3))
#     fun.plt.Indx(1:length(Var_names),Simpson,"Simpson's index","REGION",Var_names,c(0.3, max(Simpson) + 0.3))
#     fun.plt.Indx(1:length(Var_names),Margalef,"Margalef species richness","REGION",Var_names,c(8, max(Margalef) + 1.2))  
#   })
#   
#   #######################################
#   #Communities diversity throughout zones
#   zones=table(DATA1$ZONE,DATA1$SPECIES)
#   Div.InDX.zone=Div.ind(data=zones)
#   par(mfcol=c(2,2))
#   with(Div.InDX.zone,{
#     fun.plt.Indx(1:length(Var_names),Shannon,"Shannon-Wiener diversity index","ZONE",Var_names,c(1, max(Shannon) + 1))
#     fun.plt.Indx(1:length(Var_names),Pielou,"Pielou's evenness","ZONE",Var_names,c(0, max(Pielou) + 0.3))
#     fun.plt.Indx(1:length(Var_names),Simpson,"Simpson's index","ZONE",Var_names,c(0.4, max(Simpson) + 0.3))
#     fun.plt.Indx(1:length(Var_names),Margalef,"Margalef species richness","ZONE",Var_names,c(8, max(Margalef) + 1.5))  
#   })
#   
#   ########################################
#   #Communities diversity throughout blocks
#   blocks=table(DATA1$BLOCK,DATA1$SPECIES)
#   Div.InDX.blo=Div.ind(data=blocks)
#   par(mfcol=c(2,2))
#   with(Div.InDX.blo,{
#     fun.plt.Indx(1:length(Var_names),Shannon,"Shannon-Wiener diversity index","BLOCK",Var_names,c(0, max(Shannon, na.rm=T) + 1))
#     fun.plt.Indx(1:length(Var_names),Pielou,"Pielou's evenness","BLOCK",Var_names,c(0, max(Pielou, na.rm=T) + 0.2))
#     fun.plt.Indx(1:length(Var_names),Simpson,"Simpson's index","BLOCK",Var_names,c(0, max(Simpson) + 0.2))
#     fun.plt.Indx(1:length(Var_names),Margalef,"Margalef species richness","BLOCK",Var_names,c(0, max(Margalef, na.rm=T) + 4))  
#   })
#   
#   ############################################################################
#   #Communities diversity throughout years and regions (space-time interaction)
#   interac_list=vector(mode="list")                         #empty list of length zero
#   for(i in sort(unique(DATA1$REGION)))
#   {
#     subset_region=subset(DATA1,REGION==i)
#     years=table(subset_region$YEAR,subset_region$SPECIES)
#     Div.InDX.interac=Div.ind(data=years)
#     interac_list[[i]]=Div.InDX.interac                    #list of lists with the indeces of the communities per year for each region
#   }
#   
#   par(mfcol=c(2,2), xpd=T)
#   years_range=range(as.numeric(as.character(DATA1$YEAR)))  #for XLIM
#   
#   #Shannon
#   index_range=max(c(interac_list[[1]][[1]],interac_list[[2]][[1]],interac_list[[3]][[1]],interac_list[[4]][[1]],interac_list[[5]][[1]],interac_list[[6]][[1]]),na.rm=T)                        #for YLIM   
#   with(interac_list[[1]],{                                 #plot the first Region
#     fun.plt.Indx(Var_names,Shannon,"Shannon-Wiener diversity index","YEAR",Var_names,c(0.5,index_range + 0.5),years_range)
#   })
#   for(i in 2:length(interac_list))                         #Adding the other regions points to the existing plot
#   {
#     with(interac_list[[i]],{
#       points(Var_names, Shannon, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = i)
#     })
#   }
#   
#   #Pielou
#   index_range=max(c(interac_list[[1]][[3]],interac_list[[2]][[3]],interac_list[[3]][[3]],interac_list[[4]][[3]],interac_list[[5]][[3]],interac_list[[6]][[3]]),na.rm=T)                
#   with(interac_list[[1]],{                                 
#     fun.plt.Indx(Var_names,Pielou,"Pielou's evenness","YEAR",Var_names,c(0,index_range + 0.3),years_range)
#   })
#   for(i in 2:length(interac_list))                         
#   {
#     with(interac_list[[i]],{
#       points(Var_names, Pielou, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = i)
#     })
#   }
#   
#   #Simpson
#   index_range=max(c(interac_list[[1]][[4]],interac_list[[2]][[4]],interac_list[[3]][[4]],interac_list[[4]][[4]],interac_list[[5]][[4]],interac_list[[6]][[4]]),na.rm=T)               
#   with(interac_list[[1]],{                                 
#     fun.plt.Indx(Var_names,Simpson,"Simpson's index","YEAR",Var_names,c(0.2,index_range + 0.2),years_range)
#   })
#   for(i in 2:length(interac_list))                        
#   {
#     with(interac_list[[i]],{
#       points(Var_names, Simpson, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = i)
#     })
#   }
#   
#   #Margalef
#   index_range=max(c(interac_list[[1]][[2]],interac_list[[2]][[2]],interac_list[[3]][[2]],interac_list[[4]][[2]],interac_list[[5]][[2]],interac_list[[6]][[2]]),na.rm=T)            
#   with(interac_list[[1]],{                                
#     fun.plt.Indx(Var_names,Margalef,"Margalef species richness","YEAR",Var_names,c(8,index_range + 4),years_range)
#   })
#   for(i in 2:length(interac_list))                       
#   {
#     with(interac_list[[i]],{
#       points(Var_names, Margalef, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = i)
#     })
#   }
#   
#   mtext(names(interac_list), side=3, line=2.5, at=c(1992,1996,2000,2004,2008,2012), cex=1)          #points references
#   points(c(1990.5,1994.5,1998.5,2002.5,2006.5,2010.5), rep(39,6), pch=rep(20, 6), col=1:6, cex=3)
#   
#   ##############################################################################
#   #Communities diversity throughout years and bioregions (space-time interaction)
#   interac_list=vector(mode="list")                         #empty list of length zero
#   for(i in sort(unique(DATA1$BIOREGION)))
#   {
#     subset_bioregion=subset(DATA1,BIOREGION==i)
#     years=table(subset_bioregion$YEAR,subset_bioregion$SPECIES)
#     Div.InDX.interac2=Div.ind(data=years)
#     interac_list[[i]]=Div.InDX.interac2                   #list of lists with the indeces of the communities per year for each bioregion
#   }
#   
#   par(mfcol=c(2,2), xpd=T)
#   years_range=range(as.numeric(as.character(DATA1$YEAR)))  #for XLIM
#   
#   #Shannon
#   index_range=max(c(interac_list[[1]][[1]],interac_list[[2]][[1]]),na.rm=T)                        #for YLIM   
#   with(interac_list[[1]],{                                 #plot the first bioregion
#     fun.plt.Indx(Var_names,Shannon,"Shannon-Wiener diversity index","YEAR",Var_names,c(1,index_range + 0.2),years_range)
#   })
#   with(interac_list[[2]],{                                 #plot the second bioregion
#     points(Var_names, Shannon, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = 2)
#   })
#   
#   #Pielou
#   index_range=max(c(interac_list[[1]][[3]],interac_list[[2]][[3]]),na.rm=T)                
#   with(interac_list[[1]],{                                 
#     fun.plt.Indx(Var_names,Pielou,"Pielou's evenness","YEAR",Var_names,c(0,index_range + 0.3),years_range)
#   })
#   with(interac_list[[2]],{
#     points(Var_names, Pielou, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = 2)
#   })
#   
#   #Simpson
#   index_range=max(c(interac_list[[1]][[4]],interac_list[[2]][[4]]),na.rm=T)               
#   with(interac_list[[1]],{                                 
#     fun.plt.Indx(Var_names,Simpson,"Simpson's index","YEAR",Var_names,c(0.3,index_range + 0.2),years_range)
#   })
#   with(interac_list[[2]],{
#     points(Var_names, Simpson, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = 2)
#   })
#   
#   #Margalef
#   index_range=max(c(interac_list[[1]][[2]],interac_list[[2]][[2]]),na.rm=T)            
#   with(interac_list[[1]],{                                
#     fun.plt.Indx(Var_names,Margalef,"Margalef species richness","YEAR",Var_names,c(8,index_range + 4),years_range)
#   })
#   with(interac_list[[2]],{
#     points(Var_names, Margalef, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = 2)
#   })
#   
#   mtext(names(interac_list), side=3, line=2.5, at=c(1997,2007), cex=1)          #points references
#   points(c(1995,2005), rep(28.3,2), pch=rep(20, 6), col=1:2, cex=3)
#   
#   
#   ############################
#   #Ecosystem-based indicators#
#   ############################
#   
#   #Check shots with few records
#   A=sort(table(DATA$SHEET_NO))
#   a=subset(A,A<2)
#   StoReE=vector('list',length(a))
#   names(StoReE)=a
#   for(p in 1:length(a)) StoReE[[p]]=with(subset(DATA,SHEET_NO==names(a)[p]),table(as.character(COMMON_NAME)))
#   
#   
#   DATA2=subset(DATA,NATURE!="TEPS" & REGION!="Out.of.region",select=c("SHEET_NO","YEAR","BLOCK","SPECIES","TL","FL","TROPHIC_LEVEL","ZONE","REGION","BIOREGION"))  #not using TEPS and out of region data
#   DATA2$YEAR=as.factor(DATA2$YEAR)
#   DATA2$SPECIES=as.factor(DATA2$SPECIES)
#   
#   t_level=aggregate(TROPHIC_LEVEL ~ SPECIES, FUN=mean, data=DATA2) #trophic level by species table
#   dummy0=NULL
#   species=for(i in levels(DATA2$SPECIES)){dummy0=c(dummy0,rep(i,17))} #17 times each species
#   
#   ################################################################################
#   #Ecosystem-based indicators throughout years and regions (space-time interaction)
#   
#   N.boots=5
#   Reg.temp=1:4
#   ECO.I.R=vector('list',length(Reg.temp))
#   names(ECO.I.R)=Reg.temp
#   Yr.list=ECO.I.R
#   for(r in 1:length(Reg.temp))
#   {
#     D.boot=Boot.eco.indx(Boot.what="MIN",D=subset(DATA2,Reg.Temp==Reg.temp[r]))
#     Eco.InDX=vector('list',N.boots)                            
#     for(n in 1:N.boots)
#     {
#       d=lapply(X = D.boot, FUN = `[[`, as.character(n))
#       subset_region=do.call(rbind,d)
#       subset_region$SPECIES=as.character(subset_region$SPECIES)
#       years=table(subset_region$YEAR,subset_region$SPECIES)         #table of the number of individuals for all species and years
#       fl_max=aggregate(FL~YEAR+SPECIES,FUN=max,data=subset_region) #maximum length of all species and years   
#       fl_max_table=acast(fl_max,YEAR~SPECIES,value.var="FL")              
#       
#       Eco.InDX[[n]]=Eco.ind(data=years,t_level,fl_max_table)
#     }
#     Yr.list[[r]]=sort(as.numeric(names(D.boot)))  
#     MTL=do.call(rbind,lapply(X = Eco.InDX, FUN = `[[`, "MTL"))
#     MML=do.call(rbind,lapply(X = Eco.InDX, FUN = `[[`, "MML"))
#     FIB=do.call(rbind,lapply(X = Eco.InDX, FUN = `[[`, "FIB"))  
#     ECO.I.R[[r]]=list(MTL=MTL,MML=MML,FIB=FIB)                        
#   }
#   
#   
#   #Plot indices
#   STATS=ECO.I.R
#   for(r in 1:length(Reg.temp))
#   {
#     STATS[[r]]$MTL=Stats(Yr=Yr.list[[r]],DD=ECO.I.R[[r]]$MTL,lw=0.025,up=0.975)
#     STATS[[r]]$MML=Stats(Yr=Yr.list[[r]],DD=ECO.I.R[[r]]$MML,lw=0.025,up=0.975)
#     STATS[[r]]$FIB=Stats(Yr=Yr.list[[r]],DD=ECO.I.R[[r]]$FIB,lw=0.025,up=0.975)
#   }
#   
#   YEARS=as.numeric(as.character(sort(unique(DATA2$YEAR))))
#   
#   
#   fn.plt(d=STATS$MTL) 
#   
#   # dummy=data.frame(YEAR=rep(levels(DATA2$YEAR),99),BLOCK=NA,SPECIES=dummy0,TL=NA,FL=NA,TROPHIC_LEVEL=NA,ZONE=NA,REGION=NA,BIOREGION=NA)#data.frame whit the combination of all species and years                          
#   # dummy[is.na(dummy)]=0
#   # interac_list2=vector(mode="list")                                #empty list of length zero
#   # for(i in sort(unique(DATA2$REGION)))
#   #    {
#   #    subset_region=subset(DATA2,REGION==i)                   
#   #    years=table(subset_region$YEAR,subset_region$SPECIES)         #table of the number of individuals for all species and years
#   #    subset_region2=rbind(subset_region,dummy)
#   #    fl_max=aggregate(FL~YEAR+SPECIES,FUN=max,data=subset_region2) #maximum length of all species and years   
#   #    fl_max_table=acast(fl_max,YEAR~SPECIES,value.var="FL")              
#   #    Eco.InDX.interac=Eco.ind(data=years,t_level,fl_max_table)
#   #    interac_list2[[i]]=Eco.InDX.interac                           #list of lists with the indeces of the communities per year for each region
#   #    }
#   
#   par(mfcol=c(2,2), xpd=T)
#   years_range=range(as.numeric(as.character(DATA2$YEAR)))  #for XLIM
#   
#   #Mean Trophic Level
#   index_range=max(c(interac_list2[[1]][[1]],interac_list2[[2]][[1]],interac_list2[[3]][[1]],interac_list2[[4]][[1]],interac_list2[[5]][[1]],interac_list2[[6]][[1]]),na.rm=T)                       #for YLIM   
#   with(interac_list2[[1]],{                                #plot the first Region
#     fun.plt.Indx(Var_names,Mean_trophic_level,"Mean Trophic Level","YEAR",Var_names,c(3.5,index_range + 0.2),years_range)
#   })
#   for(i in 2:length(interac_list2))                        #Adding the other regions points to the existing plot
#   {
#     with(interac_list2[[i]],{                             
#       points(Var_names, Mean_trophic_level, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = i)
#     })
#   }
#   
#   #Mean Maximum Length
#   index_range=max(c(interac_list2[[1]][[2]],interac_list2[[2]][[2]],interac_list2[[3]][[2]],interac_list2[[4]][[2]],interac_list2[[5]][[2]],interac_list2[[6]][[2]]),na.rm=T)                       #for YLIM   
#   with(interac_list2[[1]],{                                #plot the first Region
#     fun.plt.Indx(Var_names,Max_length,"Mean Maximum Length","YEAR",Var_names,c(20,index_range + 20),years_range)
#   })
#   for(i in 2:length(interac_list2))                        #Adding the other regions points to the existing plot
#   {
#     with(interac_list2[[i]],{
#       points(Var_names, Max_length, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = i)
#     })
#   }
#   
#   #Fishery in Balance
#   index_range=range(c(interac_list2[[1]][[3]],interac_list2[[2]][[3]],interac_list2[[3]][[3]],interac_list2[[4]][[3]],interac_list2[[5]][[3]],interac_list2[[6]][[3]]),na.rm=T)                       #for YLIM   
#   with(interac_list2[[1]],{                                #plot the first Region
#     fun.plt.Indx(Var_names,FIB,"Fishery in Balance","YEAR",Var_names,c(index_range[1] - 1,index_range[2] + 1),years_range)
#   })
#   for(i in 2:length(interac_list2))                        #Adding the other regions points to the existing plot
#   {
#     with(interac_list2[[i]],{
#       points(Var_names, FIB, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = i)
#     })
#   }
#   
#   mtext(names(interac_list2), side=1, line=4.25, at=c(1992,1996,2000,2004,2008,2012), cex=1)          #points references
#   points(c(1990.5,1994.5,1998.5,2002.5,2006.5,2010.5), rep(-9.75,6), pch=rep(20, 6), col=1:6, cex=3)
#   
#   ###################################################################################
#   #Ecosystem-based indicators throughout years and bioregions (space-time interaction)
#   dummy=data.frame(YEAR=rep(levels(DATA2$YEAR),99),BLOCK=NA,SPECIES=dummy0,TL=NA,FL=NA,TROPHIC_LEVEL=NA,ZONE=NA,REGION=NA,BIOREGION=NA)#data.frame whit the combination of all species and years
#   dummy[is.na(dummy)]=0
#   interac_list2=vector(mode="list")                                #empty list of length zero
#   for(i in sort(unique(DATA2$BIOREGION)))
#   {
#     subset_bioregion=subset(DATA2,BIOREGION==i)
#     years=table(subset_bioregion$YEAR,subset_bioregion$SPECIES)   #table of the number of individuals for all species and years
#     subset_bioregion2=rbind(subset_bioregion,dummy)
#     fl_max=aggregate(FL~YEAR+SPECIES,FUN=max,data=subset_bioregion2)#maximum length of all species and years
#     fl_max_table=acast(fl_max,YEAR~SPECIES,value.var="FL")
#     Eco.InDX.interac2=Eco.ind(data=years,t_level,fl_max_table)
#     interac_list2[[i]]=Eco.InDX.interac2                          #list of lists with the indeces of the communities per year for each bioregion
#   }
#   
#   par(mfcol=c(2,2), xpd=T)
#   years_range=range(as.numeric(as.character(DATA2$YEAR)))  #for XLIM
#   
#   #Mean Trophic Level
#   index_range=max(c(interac_list2[[1]][[1]],interac_list2[[2]][[1]]),na.rm=T)                       #for YLIM
#   with(interac_list2[[1]],{                                #plot the first bioregion
#     fun.plt.Indx(Var_names,Mean_trophic_level,"Mean Trophic Level","YEAR",Var_names,c(3.5,index_range + 0.2),years_range)
#   })
#   with(interac_list2[[2]],{                             #plot the second bioregion
#     points(Var_names, Mean_trophic_level, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = 2)
#   })
#   
#   #Mean Maximum Length
#   index_range=max(c(interac_list2[[1]][[2]],interac_list2[[2]][[2]]),na.rm=T)                       #for YLIM
#   with(interac_list2[[1]],{                                #plot the first Region
#     fun.plt.Indx(Var_names,Max_length,"Mean Maximum Length","YEAR",Var_names,c(20,index_range + 20),years_range)
#   })
#   with(interac_list2[[2]],{
#     points(Var_names, Max_length, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = 2)
#   })
#   
#   #Fishery in Balance
#   index_range=range(c(interac_list2[[1]][[3]],interac_list2[[2]][[3]]),na.rm=T)                       #for YLIM
#   with(interac_list2[[1]],{                                #plot the first Region
#     fun.plt.Indx(Var_names,FIB,"Fishery in Balance","YEAR",Var_names,c(index_range[1] - 1,index_range[2] + 1),years_range)
#   })
#   with(interac_list2[[2]],{
#     points(Var_names, FIB, type = "o", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = 2)
#   })
#   
#   mtext(names(interac_list2), side=1, line=4.5, at=c(1997,2007), cex=1)          #points references
#   points(c(1995,2005), rep(-4,2), pch=rep(20, 2), col=1:2, cex=3)
#   
#   
#   
#   #######
#   #INEXT#
#   #######
#   
#   ####################
#   #by-Year assamblages
#   DATA3=subset(DATA,NATURE!="TEPS" & REGION!="Out.of.region")
#   sample_vector=NULL
#   for(i in sort(unique(DATA3$YEAR)))                           #calculates number of surveys per region
#   {
#     sub_year=subset(DATA3,YEAR==i)
#     sample_number=length(unique(sub_year$SHEET_NO))
#     sample_vector=c(sample_vector,sample_number)  
#   }                                                            
#   freq_list=list("1993"=NULL,"1994"=NULL,"1995"=NULL,"1996"=NULL,"1997"=NULL,"1998"=NULL,"1999"=NULL,"2000"=NULL,"2001"=NULL,"2002"=NULL,"2003"=NULL,"2004"=NULL,"2005"=NULL,"2006"=NULL,"2007"=NULL,"2012"=NULL,"2013"=NULL)
#   for(i in sort(unique(DATA3$YEAR)))                          #calculates the list of species frequencies per region
#   {
#     sub_year=subset(DATA3,YEAR==i)  
#     sp_matrix=table(sub_year$SHEET_NO,sub_year$SPECIES)
#     sp_matrix=ifelse(sp_matrix>0,1,0)
#     freq_vector=NULL
#     for(j in colnames(sp_matrix))
#     {
#       freq_value=sum(sp_matrix[,j])
#       freq_vector=c(freq_vector,freq_value)
#     }
#     freq_list[[i]]=sort(freq_vector,decreasing=T)
#   }
#   for(i in 1:length(freq_list))                                 #introduce the number of surveys before the frequencies in each vector
#   {
#     freq_list[[i]]=c(sample_vector[i],freq_list[[i]]) 
#   }
#   
#   year_iNEXT=iNEXT(freq_list, q=c(0,1,2), datatype="incidence_freq", size=seq(1,500,20))
#   ggiNEXT(year_iNEXT, type=1, facet.var="site")
#   
#   ###################### 
#   #by-Region assamblages
#   sample_vector=NULL
#   for(i in sort(unique(DATA3$REGION)))                           #calculates number of surveys per region
#   {
#     sub_region=subset(DATA3,REGION==i)
#     sample_number=length(unique(sub_region$SHEET_NO))
#     sample_vector=c(sample_vector,sample_number)  
#   }                                                            
#   freq_list=list(Region1=NULL,Region2=NULL,Region3=NULL,Region4=NULL,Region5=NULL,Region6=NULL)
#   for(i in sort(unique(DATA3$REGION)))                          #calculates the list of species frequencies per region
#   {
#     sub_region=subset(DATA3,REGION==i)  
#     sp_matrix=table(sub_region$SHEET_NO,sub_region$SPECIES)
#     sp_matrix=ifelse(sp_matrix>0,1,0)
#     freq_vector=NULL
#     for(j in colnames(sp_matrix))
#     {
#       freq_value=sum(sp_matrix[,j])
#       freq_vector=c(freq_vector,freq_value)
#     }
#     freq_list[[i]]=sort(freq_vector,decreasing=T)
#   }
#   for(i in 1:length(freq_list))                                 #introduce the number of surveys before the frequencies in each vector
#   {
#     freq_list[[i]]=c(sample_vector[i],freq_list[[i]]) 
#   }
#   
#   region_iNEXT=iNEXT(freq_list, q=c(0,1,2), datatype="incidence_freq", size=seq(1,1200,20))
#   ggiNEXT(region_iNEXT, type=1, facet.var="site")
#   
#   ########################
#   #by-Bioegion assamblages
#   sample_vector=NULL
#   for(i in sort(unique(DATA3$BIOREGION)))                           #calculates number of surveys per region
#   {
#     sub_bioregion=subset(DATA3,BIOREGION==i)
#     sample_number=length(unique(sub_bioregion$SHEET_NO))
#     sample_vector=c(sample_vector,sample_number)  
#   }                                                            
#   freq_list=list(SouthCoast=NULL,WestCoast=NULL)
#   for(i in sort(unique(DATA3$BIOREGION)))                          #calculates the list of species frequencies per region
#   {
#     sub_bioregion=subset(DATA3,BIOREGION==i)  
#     sp_matrix=table(sub_bioregion$SHEET_NO,sub_bioregion$SPECIES)
#     sp_matrix=ifelse(sp_matrix>0,1,0)
#     freq_vector=NULL
#     for(j in colnames(sp_matrix))
#     {
#       freq_value=sum(sp_matrix[,j])
#       freq_vector=c(freq_vector,freq_value)
#     }
#     freq_list[[i]]=sort(freq_vector,decreasing=T)
#   }
#   for(i in 1:length(freq_list))                                 #introduce the number of surveys before the frequencies in each vector
#   {
#     freq_list[[i]]=c(sample_vector[i],freq_list[[i]]) 
#   }
#   
#   bioregion_iNEXT=iNEXT(freq_list, q=c(0,1,2), datatype="incidence_freq", size=seq(1,2100,100))
#   ggiNEXT(bioregion_iNEXT, type=1, facet.var="site")
#   
#   ##############################
#   #by-Year-and-Bioregion assamblages
#   #West Coast Bioregion
#   DATA4=subset(DATA3,BIOREGION=="WestCoast")
#   sample_vector=NULL
#   for(i in sort(unique(DATA4$YEAR)))                           #calculates number of surveys per region
#   {
#     sub_year=subset(DATA4,YEAR==i)
#     sample_number=length(unique(sub_year$SHEET_NO))
#     sample_vector=c(sample_vector,sample_number)  
#   }                                                            
#   freq_list=list("1993"=NULL,"1994"=NULL,"1995"=NULL,"1996"=NULL,"1997"=NULL,"1998"=NULL,"1999"=NULL,"2000"=NULL,"2001"=NULL,"2002"=NULL,"2003"=NULL,"2004"=NULL,"2005"=NULL,"2006"=NULL,"2007"=NULL,"2012"=NULL,"2013"=NULL)
#   for(i in sort(unique(DATA4$YEAR)))                          #calculates the list of species frequencies per region
#   {
#     sub_year=subset(DATA4,YEAR==i)  
#     sp_matrix=table(sub_year$SHEET_NO,sub_year$SPECIES)
#     sp_matrix=ifelse(sp_matrix>0,1,0)
#     freq_vector=NULL
#     for(j in colnames(sp_matrix))
#     {
#       freq_value=sum(sp_matrix[,j])
#       freq_vector=c(freq_vector,freq_value)
#     }
#     freq_list[[i]]=sort(freq_vector,decreasing=T)
#   }
#   for(i in 1:length(freq_list))                                 #introduce the number of surveys before the frequencies in each vector
#   {
#     freq_list[[i]]=c(sample_vector[i],freq_list[[i]]) 
#   }
#   
#   w.coast_iNEXT=iNEXT(freq_list, q=c(0,1,2), datatype="incidence_freq", size=seq(1,280,12))
#   ggiNEXT(w.coast_iNEXT, type=1, facet.var="site")
#   
#   ###################### 
#   #South Coast Bioregion
#   DATA5=subset(DATA3,BIOREGION=="SouthCoast")
#   sample_vector=NULL
#   for(i in sort(unique(DATA5$YEAR)))                           #calculates number of surveys per region
#   {
#     sub_year=subset(DATA5,YEAR==i)
#     sample_number=length(unique(sub_year$SHEET_NO))
#     sample_vector=c(sample_vector,sample_number)  
#   } 
#   freq_list=list("1993"=NULL,"1994"=NULL,"1995"=NULL,"1996"=NULL,"1997"=NULL,"1998"=NULL,"1999"=NULL,"2012"=NULL,"2013"=NULL)
#   for(i in sort(unique(DATA5$YEAR)))                          #calculates the list of species frequencies per region
#   {
#     sub_year=subset(DATA5,YEAR==i)  
#     sp_matrix=table(sub_year$SHEET_NO,sub_year$SPECIES)
#     sp_matrix=ifelse(sp_matrix>0,1,0)
#     freq_vector=NULL
#     for(j in colnames(sp_matrix))
#     {
#       freq_value=sum(sp_matrix[,j])
#       freq_vector=c(freq_vector,freq_value)
#     }
#     freq_list[[i]]=sort(freq_vector,decreasing=T)
#   }
#   for(i in 1:length(freq_list))                                 #introduce the number of surveys before the frequencies in each vector
#   {
#     freq_list[[i]]=c(sample_vector[i],freq_list[[i]]) 
#   }
#   
#   s.coast_iNEXT=iNEXT(freq_list, q=c(0,1,2), datatype="incidence_freq", size=seq(1,250,10))
#   ggiNEXT(s.coast_iNEXT, type=1, facet.var="site")
#   
#   
#   ######################################
#   #nMDS analysis of species composition#
#   ######################################
#   
#   DATA6=subset(DATA,NATURE!="TEPS" & REGION!="Out.of.region" & !is.na(EFFORT))#not using TEPS, out of region and no effort data
#   
#   ############
#   #All species
#   all_sp=aggregate(INDIVIDUALS~SHEET_NO+SPECIES+YEAR+REGION+BIOREGION+ZONE+EFFORT,FUN=sum,data=DATA6)#Individuals per sheet number and species
#   all_sp$CPUE=all_sp$INDIVIDUALS/all_sp$EFFORT         #CPUE calculation
#   all_sp$SPECIES=as.factor(all_sp$SPECIES)
#   
#   #YEARS
#   years=aggregate(CPUE~YEAR+SPECIES,all_sp,mean)                #aggregate the CPUE by YEAR and SPECIES
#   years1=acast(years,YEAR~SPECIES,value.var="CPUE")             #reshape de data frame into matrix-like for the metaMDS function
#   years1[is.na(years1)]=0                                       #replace NA with 0 values as it means no capture 
#   NMDS1=metaMDS(years1,wascores=T)                              #nMDS
#   plot(NMDS1,type="t",main=paste("nMDS-Bray, stress=",round(NMDS1$stress,3)))#nMDS biplot of a Bray–Curtis dissimilarity matrix of the species CPUE
#   
#   par(mfrow=c(1,2))
#   stressplot(NMDS1,main="Sheppard plot")                        #Sheppard plot shows the appropriateness of the nMDS + linear and non-linear fit
#   g.o.f=goodness(NMDS1)                                         #goodness of fit
#   plot(NMDS1,type="t",main="Goodness of fit")                   #nMDS biplot + goodness of fit
#   points(NMDS1,cex=g.o.f*200)                                   #poorly fitted years have larger bubbles
#   
#   #REGIONS
#   regions=aggregate(CPUE~REGION+SPECIES,all_sp,mean)            #aggregate the CPUE by REGION and SPECIES
#   regions1=acast(regions,REGION~SPECIES,value.var="CPUE")       #reshape de data frame into matrix-like for the metaMDS function
#   regions1[is.na(regions1)]=0                                   #replace NA with 0 values as it means no capture 
#   NMDS2=metaMDS(regions1,wascores=T)                            #nMDS
#   plot(NMDS2,type="t",main=paste("nMDS-Bray, stress=",round(NMDS2$stress,3)))#nMDS biplot of a Bray–Curtis dissimilarity matrix of the species CPUE
#   
#   par(mfrow=c(1,2))
#   stressplot(NMDS2,main="Sheppard plot")                        #Sheppard plot shows the appropriateness of the nMDS + linear and non-linear fit
#   g.o.f=goodness(NMDS2)                                         #goodness of fit
#   plot(NMDS2,type="t",main="Goodness of fit")                   #nMDS biplot + goodness of fit
#   points(NMDS2,cex=g.o.f*200000)                                #poorly fitted years have larger bubbles
#   
#   #REGIONS AND YEARS
#   years_regions=NULL
#   plot_col=NULL
#   for(i in sort(unique(all_sp$REGION)))
#   {
#     subset_region=subset(all_sp,REGION==i)                   
#     years=aggregate(CPUE~YEAR+SPECIES,subset_region,mean)      #aggregate the CPUE by YEAR and SPECIES
#     dummy.df=data.frame(YEAR=unique(years$YEAR)[1],SPECIES=unique(all_sp$SPECIES),CPUE=NA)#so all the species appear in the disimilarity matrix
#     years1=rbind(years,dummy.df)  
#     years2=acast(years1,YEAR~SPECIES,value.var="CPUE",fun.aggregate=mean,na.rm=T)#reshape de data frame into matrix-like for the metaMDS function
#     years_regions=rbind(years_regions,years2)                  #growing table of years and species
#     plot_col=c(plot_col,length(unique(years$YEAR)))            #for the plotting
#   }
#   years_regions[is.na(years_regions)]=0
#   NMDS3=metaMDS(years_regions^(1/4),wascores=F,autotransform=F) #nMDS + fourth-root transformation to reduce variance
#   plot(NMDS3,type="t",main=paste("nMDS-Bray, stress=",round(NMDS3$stress,3)))#nMDS biplot of a Bray–Curtis dissimilarity matrix of the species CPUE
#   points(NMDS3,pch=20,col=c(rep(1,plot_col[1]),rep(2,plot_col[2]),rep(3,plot_col[3]),rep(4,plot_col[4]),rep(5,plot_col[5]),rep(6,plot_col[6])))
#   legend("topright",sort(unique(all_sp$REGION)),pch=20,pt.cex=2,col=1:6)
#   
#   par(mfrow=c(1,2))
#   stressplot(NMDS3,main="Sheppard plot")                        #Sheppard plot shows the appropriateness of the nMDS + linear and non-linear fit
#   g.o.f=goodness(NMDS3)                                         #goodness of fit
#   plot(NMDS3,type="t",main="Goodness of fit")                   #nMDS biplot + goodness of fit
#   points(NMDS3,cex=g.o.f*200)                                   #poorly fitted years have larger bubbles
#   
#   #BIOREGIONS AND YEARS
#   years_bioregions=NULL
#   plot_col=NULL
#   for(i in sort(unique(all_sp$BIOREGION)))
#   {
#     subset_bioregion=subset(all_sp,BIOREGION==i)                   
#     years=aggregate(CPUE~YEAR+SPECIES,subset_bioregion,mean)   #aggregate the CPUE by YEAR and SPECIES
#     dummy.df=data.frame(YEAR=unique(years$YEAR)[1],SPECIES=unique(all_sp$SPECIES),CPUE=NA)#so all the species appear in the disimilarity matrix
#     years1=rbind(years,dummy.df)  
#     years2=acast(years1,YEAR~SPECIES,value.var="CPUE",fun.aggregate=mean,na.rm=T)#reshape de data frame into matrix-like for the metaMDS function
#     years_bioregions=rbind(years_bioregions,years2)            #growing table of years and species
#     plot_col=c(plot_col,length(unique(years$YEAR)))            #for the plotting
#   }
#   years_bioregions[is.na(years_bioregions)]=0
#   NMDS4=metaMDS(years_bioregions^(1/4),wascores=F,autotransform=F)#nMDS + fourth-root transformation to reduce variance
#   plot(NMDS4,type="t",main=paste("nMDS-Bray, stress=",round(NMDS4$stress,3)))#nMDS biplot of a Bray–Curtis dissimilarity matrix of the species CPUE
#   points(NMDS4,pch=20,col=c(rep(1,plot_col[1]),rep(2,plot_col[2])))
#   legend("topright",sort(unique(all_sp$BIOREGION)),pch=20,pt.cex=2,col=1:2)
#   
#   par(mfrow=c(1,2))
#   stressplot(NMDS4,main="Sheppard plot")                        #Sheppard plot shows the appropriateness of the nMDS + linear and non-linear fit
#   g.o.f=goodness(NMDS4)                                         #goodness of fit
#   plot(NMDS4,type="t",main="Goodness of fit")                   #nMDS biplot + goodness of fit
#   points(NMDS4,cex=g.o.f*200)                                   #poorly fitted years have larger bubbles
#   
#   
#   #######################
#   #Only discarded species
#   
#   discarded=subset(DATA6,FATE=="D",select=c("SHEET_NO","YEAR","REGION","BIOREGION","ZONE","EFFORT","SPECIES","INDIVIDUALS"))#subset of dicarded records    
#   discarded2=aggregate(INDIVIDUALS~SHEET_NO+SPECIES+YEAR+REGION+BIOREGION+ZONE+EFFORT,FUN=sum,data=discarded)#Individuals per sheet number and species
#   discarded2$CPUE=discarded2$INDIVIDUALS/discarded2$EFFORT#CPUE calculation
#   discarded2$SPECIES=as.factor(discarded2$SPECIES)
#   
#   #YEARS
#   years=aggregate(CPUE~YEAR+SPECIES,discarded2,mean)           #aggregate the CPUE by YEAR and SPECIES
#   years1=acast(years,YEAR~SPECIES,value.var="CPUE")            #reshape de data frame into matrix-like for the metaMDS function
#   years1[is.na(years1)]=0                                      #replace NA with 0 values as it means no capture 
#   NMDS5=metaMDS(years1^(1/4),wascores=T,autotransform=F)       #nMDS + fourth-root transformation to reduce variance
#   plot(NMDS5,type="t",main=paste("nMDS-Bray, stress=",round(NMDS5$stress,3)))#nMDS biplot of a Bray–Curtis dissimilarity matrix of the species CPUE
#   
#   par(mfrow=c(1,2))
#   stressplot(NMDS5,main="Sheppard plot")                       #Sheppard plot shows the appropriateness of the nMDS + linear and non-linear fit
#   g.o.f=goodness(NMDS5)                                        #goodness of fit
#   plot(NMDS5,type="t",main="Goodness of fit")                  #nMDS biplot + goodness of fit
#   points(NMDS5,cex=g.o.f*200)                                  #poorly fitted years have larger bubbles
#   
#   #REGIONS
#   regions=aggregate(CPUE~REGION+SPECIES,discarded2,mean)       #aggregate the CPUE by REGION and SPECIES
#   regions1=acast(regions,REGION~SPECIES,value.var="CPUE")      #reshape de data frame into matrix-like for the metaMDS function
#   regions1[is.na(regions1)]=0                                  #replace NA with 0 values as it means no capture 
#   NMDS6=metaMDS(regions1^(1/4),wascores=T,autotransform=F)     #nMDS + fourth-root transformation to reduce variance
#   plot(NMDS6,type="t",main=paste("nMDS-Bray, stress=",round(NMDS6$stress,3)))#nMDS biplot of a Bray–Curtis dissimilarity matrix of the species CPUE
#   
#   par(mfrow=c(1,2))
#   stressplot(NMDS6,main="Sheppard plot")                       #Sheppard plot shows the appropriateness of the nMDS + linear and non-linear fit
#   g.o.f=goodness(NMDS6)                                        #goodness of fit
#   plot(NMDS6,type="t",main="Goodness of fit")                  #nMDS biplot + goodness of fit
#   points(NMDS6,cex=g.o.f*200000)                               #poorly fitted years have larger bubbles
#   
#   #REGIONS AND YEARS
#   years_regions=NULL
#   plot_col=NULL
#   for(i in sort(unique(discarded2$REGION)))
#   {
#     subset_region=subset(discarded2,REGION==i)                   
#     years=aggregate(CPUE~YEAR+SPECIES,subset_region,mean)     #aggregate the CPUE by YEAR and SPECIES
#     dummy.df=data.frame(YEAR=unique(years$YEAR)[1],SPECIES=unique(discarded2$SPECIES),CPUE=NA)#so all the species appear in the matrix
#     years1=rbind(years,dummy.df)  
#     years2=acast(years1,YEAR~SPECIES,value.var="CPUE",fun.aggregate=mean,na.rm=T)#reshape de data frame into matrix-like for the metaMDS function
#     years_regions=rbind(years_regions,years2)                 #growing table of years and species
#     plot_col=c(plot_col,length(unique(years$YEAR)))           #for the plotting
#   }
#   years_regions[is.na(years_regions)]=0
#   NMDS7=metaMDS(years_regions^(1/4),wascores=F,autotransform=F)#nMDS + fourth-root transformation to reduce variance
#   plot(NMDS7,type="t",main=paste("nMDS-Bray, stress=",round(NMDS7$stress,3)))#nMDS biplot of a Bray–Curtis dissimilarity matrix of the species CPUE
#   points(NMDS7,pch=20,col=c(rep(1,plot_col[1]),rep(2,plot_col[2]),rep(3,plot_col[3]),rep(4,plot_col[4]),rep(5,plot_col[5]),rep(6,plot_col[6])))
#   legend("topright",sort(unique(discarded2$REGION)),pch=20,pt.cex=2,col=1:6)
#   
#   par(mfrow=c(1,2))
#   stressplot(NMDS7,main="Sheppard plot")                       #Sheppard plot shows the appropriateness of the nMDS + linear and non-linear fit
#   g.o.f=goodness(NMDS7)                                        #goodness of fit
#   plot(NMDS7,type="t",main="Goodness of fit")                  #nMDS biplot + goodness of fit
#   points(NMDS7,cex=g.o.f*200)                                  #poorly fitted years have larger bubbles
#   
#   #BIOREGIONS AND YEARS
#   years_bioregions=NULL
#   plot_col=NULL
#   for(i in sort(unique(discarded2$BIOREGION)))
#   {
#     subset_bioregion=subset(discarded2,BIOREGION==i)                   
#     years=aggregate(CPUE~YEAR+SPECIES,subset_bioregion,mean)  #aggregate the CPUE by YEAR and SPECIES
#     dummy.df=data.frame(YEAR=unique(years$YEAR)[1],SPECIES=unique(discarded2$SPECIES),CPUE=NA)#so all the species appear in the matrix
#     years1=rbind(years,dummy.df)  
#     years2=acast(years1,YEAR~SPECIES,value.var="CPUE",fun.aggregate=mean,na.rm=T)#reshape de data frame into matrix-like for the metaMDS function
#     years_bioregions=rbind(years_bioregions,years2)           #growing table of years and species
#     plot_col=c(plot_col,length(unique(years$YEAR)))           #for the plotting
#   }
#   years_bioregions[is.na(years_bioregions)]=0
#   NMDS8=metaMDS(years_bioregions^(1/4),wascores=T,autotransform=F)#nMDS + fourth-root transformation to reduce variance
#   plot(NMDS8,type="t",main=paste("nMDS-Bray, stress=",round(NMDS8$stress,3)))#nMDS biplot of a Bray–Curtis dissimilarity matrix of the species CPUE
#   points(NMDS8,pch=20,col=c(rep(1,plot_col[1]),rep(2,plot_col[2])))
#   legend("topright",sort(unique(discarded2$BIOREGION)),pch=20,pt.cex=2,col=1:2)
#   
#   ########
#   #ANOSIM#
#   ########
#   
#   ANOSIM=unique(discarded2[,c("SHEET_NO","YEAR","ZONE","REGION","BIOREGION")])#data frame with the factors for each sample (i.e. sheet number)
#   ano.cpue=acast(discarded2[,c("SHEET_NO","SPECIES","CPUE")],SHEET_NO~SPECIES,value.var="CPUE")#reshape de data frame into matrix-like for the anosim function
#   ano.cpue[is.na(ano.cpue)]=0                                 #replace NA with 0 values as it means no capture 
#   ano.cpue2=ano.cpue^(1/4)                                    #fourth-root transformation to reduce variance
#   ano.cpue3=vegdist(ano.cpue2)
#   
#   #YEARS
#   ANOSIM.YEARS=anosim(dat=ano.cpue3,grouping=ANOSIM$YEAR,permutations=999)#takes ~4 minutes to compute
#   
#   #REGIONS
#   ANOSIM.REGIONS=anosim(dat=ano.cpue3,grouping=ANOSIM$REGION,permutations=999)#takes ~4 minutes to compute
#   
#   #BIOREGIONS
#   ANOSIM.BIOREGIONS=anosim(dat=ano.cpue3,grouping=ANOSIM$BIOREGION,permutations=999)#takes ~4 minutes to compute
#   
#   ########
#   #SIMPER#
#   ########
#   
#   #YEARS
#   SIMPER.YEARS=simper(comm=ano.cpue2,group=ANOSIM$YEAR)
#   summary(SIMPER.YEARS)    
#   
#   #REGIONS
#   SIMPER.REGIONS=simper(comm=ano.cpue2,group=ANOSIM$REGION)
#   summary(SIMPER.REGIONS)
#   
#   #REGIONS
#   SIMPER.BIOREGIONS=simper(comm=ano.cpue2,group=ANOSIM$BIOREGION)
#   summary(SIMPER.BIOREGIONS)
#   
#   
#   #############################
#   #Bycatch proportion analyses#
#   #############################
#   
#   DATA7=subset(DATA,NATURE!="TEPS" & REGION!="Out.of.region" & !is.na(EFFORT))#not using TEPS, out of region and no effort data
#   
#   #Number of discarded sharks, rays and TEPS in total surveys
#   no.surveys=length(unique(DATA7$SHEET_NO))                    #number of total surveys
#   discarded=subset(DATA7, FATE=="D")                           #subset of dicarded records
#   no.surveys2=length(unique(discarded$SHEET_NO))               #number of surveys with discards
#   no.surveys2*100/no.surveys                                   #percentage of surveys with discards from total surveys
#   dim(discarded)[1]                                            #total discarded individuals
#   dim(discarded)[1]*100/dim(DATA7)[1]                          #percentage of discarded individuals from total individuals
#   range(table(discarded$SHEET_NO))                             #range of discarded individuals per survey
#   hist(table(discarded$SHEET_NO),breaks=200)                   #histogram of discarded individuals per survey
#   
#   #Discard proportion calculation
#   discards=subset(DATA7, FATE=="D" | FATE=="C")                #only discarded or commercial species
#   vector=NULL                                                  #dummy
#   for(i in unique(discards$SHEET_NO)){
#     sheet_no=subset(discards, SHEET_NO==i)
#     proportion=sum(sheet_no$FATE=="D")/dim(sheet_no)[1]
#     vector=c(vector,proportion)
#   }                                                         #this loop takes ~3 minutes to compute 
#   glm_dataframe=data.frame(SHEET_NO=unique(discards$SHEET_NO),PROPORTION=vector)#new data frame with discarded individual proportions
#   
#   #Proportion distribution
#   hist(glm_dataframe$PROPORTION)
#   boxplot(glm_dataframe$PROPORTION)                               
#   
#   #GLM modelling
#   glm_data_frame=merge(glm_dataframe,discards,by="SHEET_NO")   #data frame for the glm
#   merging=aggregate(PROPORTION~SHEET_NO+YEAR+MONTH+REGION+BLOCK+BOTDEPTH+NET_LENGTH+SOAK_TIME,glm_data_frame,mean)
#   glm1=glm(PROPORTION~YEAR+BLOCK+BOTDEPTH+NET_LENGTH+SOAK_TIME+REGION,family=quasibinomial,data=merging) 
#   
#   #GAMLSS modelling fitted with a zero-inflated beta distribution (BEZI)
#   model=gamlss(PROPORTION~YEAR+BOTDEPTH+NET_LENGTH+SOAK_TIME+REGION,family=BEINF,data=na.omit(merging))#use beta inflated distribution
#   summary(model)
#   plot(model)
#   
#   #BRT/Random forest
#   
#   #??????
#   
#   
#   ################################################################################################################################################
#   ################################################################################################################################################
#   ################################################################################################################################################
#   
#   #General plotting
#   
#   #Data plotting (data+blocks+towns+bioregions)
#   DATA$LAT_BLOCK=-as.numeric(substr(DATA$BLOCK,1,2))
#   DATA$LONG_BLOCK=as.numeric(substr(DATA$BLOCK,3,4))+100
#   
#   map("worldHires",xlim=c(112.95,129.5),ylim=c(-36,-26.5),type="n")
#   polygon(c(105,117,117,115.5,115.5,105),c(-27,-27,-34.15,-34.15,-40,-40),col="lightblue",border=NA)#Source:Department of Fisheries of WA (website)
#   polygon(c(115.5,117,117,130,130,115.5),c(-34.15,-34.15,-29,-29,-40,-40),col="lightgreen",border=NA)
#   map("worldHires",xlim=c(112.95,129.5),ylim=c(-36,-26.5),col="grey80",interior=T,fill=T,border="black",add=T)
#   points(DATA$LONGITUDE,DATA$LATITUDE,col=2,pch=20,cex=0.1)
#   rect(DATA$LONG_BLOCK,DATA$LAT_BLOCK-1,DATA$LONG_BLOCK+1,DATA$LAT_BLOCK,border="blue")
#   text(c(116.35,115.41,115.86,115.8,117.884,121.9,128.848),c(-31.95,-28.773,-30.5,-34.314,-34.8,-33.6,-31.5),labels=c("Perth","Geraldton","Cervantes","Augusta","Albany","Esperance","Eucla"),cex=c(rep(0.5,7),0.7))
#   legend(c(124,124),c(-27,-27.5),c("West Coast","South Coast"),fill=c("lightblue","lightgreen"),bty="n")
#   box()
#   
#   #Zones plotting
#   map("worldHires",xlim=c(112,129),ylim=c(-36,-12.5),type="n")
#   polygon(c(116.5,116.5,130,130),c(-26,-40,-40,-26),col="blue",border=NA)#Zone1
#   polygon(c(116.5,116.5,110,110),c(-33,-40,-40,-33),col="red",border=NA)#Zone2
#   polygon(c(116.5,116.5,110,110),c(-26,-33,-33,-26),col="green",border=NA)#West
#   polygon(c(114,114,110,110),c(-26,-10,-10,-26),col="black",border=NA)#Closed
#   polygon(c(114,114,123.75,123.75),c(-26,-10,-10,-26),col="orange",border=NA)#North
#   polygon(c(123.75,123.75,130,130),c(-26,-10,-10,-26),col="yellow",border=NA)#Joint
#   map("worldHires",xlim=c(112,129),ylim=c(-36,-12.5),col="grey80",interior=T,fill=T,border=NA,add=T)
#   legend(x=122,y=-21,legend=c("Zone 1","Zone 2","West","Closed","North","Joint"),fill=c("blue","red","green","black","orange","yellow"),bty="n")
#   box()
#   
#   #Regions plotting
#   map("worldHires",xlim=c(108,129),ylim=c(-39,-27),type="n")
#   polygon(c(129,124,124,129,129),c(-30,-30,-36,-35,-30),col="lightblue")#Region1
#   polygon(c(124,119,119,124,124),c(-33,-33,-38,-37,-33),col="lightblue")#Region2
#   polygon(c(119,116,116,119,119),c(-34,-34,-38,-38,-34),col="lightblue")#Region3
#   polygon(c(116,112,112,116,116),c(-33,-33,-37,-37,-33),col="lightblue")#Region4
#   polygon(c(116,110,112,116,116),c(-30,-30,-33,-33,-30),col="lightblue")#Region5
#   polygon(c(115,109,110,115,115),c(-27,-27,-30,-30,-27),col="lightblue")#Region6
#   map("worldHires",xlim=c(108,129),ylim=c(-39,-27),col="grey80",interior=T,fill=T,add=T)
#   text(c(126.5,121.5,117.5,114,113,112),c(-34,-35.5,-36,-35.5,-31.5,-28.5),c("Region 1","Region 2","Region 3","Region 4","Region 5","Region 6"),cex=0.8)
#   box()
#   
#   #Plot for spatio-temporal decisions
#   DATA$LAT_BLOCK=-as.numeric(substr(DATA$BLOCK,1,2))
#   DATA$LONG_BLOCK=as.numeric(substr(DATA$BLOCK,3,4))+100
#   DATA$MONTH=as.numeric(DATA$MONTH)
#   DATA99=subset(DATA,REGION!="Out.of.region")
#   DATA99$N=1
#   DATA99=DATA99[!duplicated(DATA99$SHEET_NO),]
#   aggr=aggregate(N~YEAR+SEASON+BLOCK+LAT_BLOCK+LONG_BLOCK,FUN=sum,data=DATA99)
#   for(i in sort(unique(aggr$YEAR)))
#   {
#     sub.set=subset(aggr,YEAR==i)
#     #pdf(file=paste(i,".pdf",sep=""))
#     tiff(file=paste(getwd(),"/MAPS/",i,".tiff",sep=""),width=2400,height=2000,units="px",res=300,compression="lzw")
#     par(mfrow=c(2,2),oma=c(.1,.1,.1,.1))
#     for(j in sort(unique(aggr$SEASON)))
#     {
#       sub.set2=subset(sub.set,SEASON==j)
#       if(dim(sub.set2)[1]==0) 
#       {plot(1,type="n",main=paste(i,j,sep=" "),axes=F,xlab=NA,ylab=NA)
#         text(1,1,"NO SURVEYS")
#       }else
#       {
#         map("worldHires",xlim=c(112.95,129.5),ylim=c(-36,-26.5),main=j,col="grey80",interior=T,fill=T,border=F)
#         rect(DATA$LONG_BLOCK,DATA$LAT_BLOCK-1,DATA$LONG_BLOCK+1,DATA$LAT_BLOCK,border="blue")
#         numbers=as.character(sub.set2$N)
#         text(x=sub.set2$LONG_BLOCK+0.5,y=sub.set2$LAT_BLOCK-0.5,labels=numbers,cex=0.75)
#         title(paste(i,j,sep=" "))
#       }
#     }
#     dev.off()
#   }
#   
#   
#   #E.g. LOOp
#   
#   #Define spatial groups
#   DATA$Group=
#     
#     Grupos=unique(DATA$Group)
#   
#   Lista=vector('list',length(Grupos))
#   names(Lista)=Grupos
#   
#   for (i in 1:length(Grupos))
#   {
#     Ll=list(Diversity=NULL,Ecosystem=NULL)
#     
#     D=subset(DATA,Group==Grupos[i])
#     
#     #Diversity indices
#     Ll$Diversity=  #hills function outputs
#       
#       #Ecosystem indices
#       Ll$Ecosystem= # ecosystem function output
#       
#       
#       
#       Lista[[i]]=Ll
#   }
#   
#   
#   ################################################################################
#   ################################################################################
#   ################################################################################
#   
#   ####################################################################################################
#   #2. Manipulated MAFFRI shark gillnet data
#   DATA.MAFFRI$SHEET_NO=with(DATA.MAFFRI,paste(Cruise,Station))
#   
#   DATA.MAFFRI$Mid.Lat=-DATA.MAFFRI$StartLat/100
#   DATA.MAFFRI$Mid.Long=DATA.MAFFRI$StartLong/100
#   
#   DATA_TEPS.MAFFRI$Mid.Lat=-DATA_TEPS.MAFFRI$StartLat/100
#   DATA_TEPS.MAFFRI$Mid.Long=DATA_TEPS.MAFFRI$StartLong/100
#   
#   DATA.MAFFRI=DATA.MAFFRI[,-match(c("StartLat","StartLong","Dummy"),names(DATA.MAFFRI))]
#   DATA_TEPS.MAFFRI=DATA_TEPS.MAFFRI[,-match(c("StartLat","StartLong","Dummy"),names(DATA_TEPS.MAFFRI))]
#   ####################################################################################################
#   
#   #3.2 MAFFRI
#   DATA.MAFFRI
#   
#   
#   #TEPS analysis
#   DATA_TEPS.MAFFRI
#   
# }
