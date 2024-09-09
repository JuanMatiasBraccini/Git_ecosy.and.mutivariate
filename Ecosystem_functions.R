
########### SCRIPT FOR COMPUTING ECOSYSTEM ANALYSIS OF FISHING IMPACTS ###################


#Reference: Hall & Wise (2011).

# "fun.plt.Indx"  --> for plotting indices
# "fn.plt"        --> completing missing years for plots
# "Stats"         --> stats summary
# "Boot.eco.indx" --> to generate uncertainty in eco indices through boostrapping
# "Eco.ind.shot"  --> calculates Mean trophic level (MTL) and Mean maximum length (MML)
# "Div.ind.shot"  --> calculates Shannon's and Simpson's diversity, Pielou's evenness and Margalef's richness



#require(BEQI2)
require(asbio)
require(marindicators)   #see fishingInBalance() margalef() meanTLLandings() etc etc, lots of functions
library(FD)   #functional diversity
library(flextable)
require(ecodist)
require(mgcv)
require(mvabund)
library(statmod)
require(tweedie)
require(pairwiseAdonis)
require(vegan)
library(boral)  #HMSC package is another good alternative
library(corrplot)

#--FUNCTIONS---

#Function for plotting
fun.plt.Indx = function(var.names, Indx, YLAB, XLAB, LABELS, YLIM, XLIM = NULL)       
{
  plot(var.names, Indx, ylim = YLIM, xlim = XLIM, ylab = YLAB, xlab = XLAB, type = "o", xaxt = "n", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = 1) #temporal trend
  axis(side = 1, at = var.names, labels = LABELS, las = 3)    
}

#Function for calculating Diversity indices
Div.ind.shot = function(data)
{
  dat.pos=subset(data,data>0)
  sp_names = names(dat.pos) 
  
  #Shannon-Wiener diversity index
  Shannon = unname(diversity(data, index = "shannon"))
  Shannon[Shannon == 0] = NA
      
  #Simpson's index
  Simpson = unname(diversity(data, index = "simpson"))
  Simpson[Simpson == 1] = NA
  
  #Margalef species richness
  Margalef= margalef(sp_names, dat.pos)
  Margalef[Margalef == 0] = NA
  
  #Pielou's evenness                        
  Evenness = function(H, S) {J = H / log(S)}  #where H = Shannon index for each community; S = Number of species
  Pielou = Evenness(H = Shannon, S = length(sp_names))
  Pielou[Pielou == 0] = NA
    
  #Summary
  return(list(Shannon = round(Shannon, 2), Margalef = round(Margalef, 2), Pielou = round(Pielou, 2), Simpson = round(Simpson, 2)))
}
Evenness = function(H, S) {J = H / log(S)}  #where H = Shannon index for each community; S = Number of species


#Function for calculating Ecosystem indicators
Eco.ind.shot = function(data, TROPHIC.LEVEL, MAX.BODY.LENGTH, TROPHIC.EFFICIENCY = 0.1)
{
  Trophic_level=Max_length=NA
  
  if(length(data)>1)
  {
    dat.pos=subset(data,data>0)
    sp_names = names(dat.pos) 
    
    ddd=data.frame(SPECIES=sp_names,Number=as.vector(data))
    TROPHIC.LEVEL=merge(TROPHIC.LEVEL,ddd,by="SPECIES")
   
    #Mean Trophic Level
    Trophic_level=with(TROPHIC.LEVEL,(sum(Number*TROPHIC_LEVEL)/sum(Number)))
    Trophic_level=round(Trophic_level, 3)
    
    #Mean Maximum Length
    if(!is.null(nrow(MAX.BODY.LENGTH)))
    {
      MAX.BODY.LENGTH=merge(MAX.BODY.LENGTH,ddd,by="SPECIES")
      Max_length=with(MAX.BODY.LENGTH,(sum(Number*TL)/sum(Number)))
      Max_length=round(Max_length, 3)
    }
  }
 
  return(list(MTL =Trophic_level , MML = Max_length))
}

#Bootstrapping for uncertainty in ecosystem indices
Boot.eco.indx=function(Boot.what,D)
{
  #loop over years 
  Yrs=sort(as.numeric(as.character(unique(D$YEAR))))
  
  dummy=D[!duplicated(D$SHEET_NO),]
  ShtN=sort(table(dummy$YEAR))
  if(Boot.what=="MIN")N.boot=ShtN[1]
  if(Boot.what=="MAX")N.boot=ShtN[length(ShtN)]
  
  boot.yr=vector('list',length(Yrs))
  names(boot.yr)=Yrs
  for(y in 1:length(Yrs))
  {
    d=subset(D,YEAR==Yrs[y])
    ShNy=unique(d$SHEET_NO)
    
    #boostrap
    bot.dat=vector('list',N.boots)
    names(bot.dat)=1:N.boots
    for(b in 1:N.boots)
    {
      ID=sample(ShNy, N.boot, replace = T)
      ID.list=vector('list',length(ID))
      for(x in 1:N.boot) ID.list[[x]]=subset(d,SHEET_NO==ID[x])
      bot.dat[[b]]=do.call(rbind,ID.list)
    }
    boot.yr[[y]]=bot.dat
  }
  return(boot.yr)
}

#Summary stats
Stats=function(Yr,DD,lw,up)
{
  MEAN=colMeans(DD,na.rm=T)
  LOW=apply(DD, 2, function(x) quantile(x, lw,na.rm=T))
  UP=apply(DD, 2, function(x) quantile(x, up,na.rm=T))
  names(MEAN)=names(LOW)=names(UP)=Yr
  return(list(MEAN=MEAN,LOW=LOW,UP=UP))
}

#Function for filling in missing year
fn.plt=function(d)
{
  #add missing year
  for(y in 1:length(d)) 
  {
    if(length(names(d[[y]]$MEAN))<length(YEARS))
    {
      id=which(!YEARS%in%names(d[[y]]$MEAN))
      addY=YEARS[id]
      dummy=rep(NA,length(addY))
      names(dummy)=addY
      mn=c(d[[y]]$MEAN,dummy)
      lw1=c(d[[y]]$LOW,dummy)
      up1=c(d[[y]]$UP,dummy)
      id=sort(names(mn))
      d[[y]]$MEAN=mn[id]
      d[[y]]$LOW=lw1[id]
      d[[y]]$UP=up1[id]
    }
  }
}

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
Export.tbl=function(WD,Tbl,Doc.nm)
{
  aid=match(unique(Tbl$Indicator), Tbl$Indicator)
  t.first <- Tbl$Indicator[aid]
  Tbl$Indicator=''
  Tbl$Indicator[aid]=t.first
  Tbl%>%
    mutate(p.value=formatC(as.numeric(p.value), format = "e", digits = 2))%>%
    flextable()%>%
    fontsize(size=8, part='body')%>%
    fontsize(size=9, part='header')%>%
    width(width=2,unit='cm')%>%
    width(j='Term',width=4,unit='cm')%>%
    align(align="left", part = "all")%>%
    valign(valign = "top", part = "all")%>%
    bg(bg = 'grey85', part = "header")%>%
    bold(part = "header")%>%
    save_as_docx( path = paste(WD,paste0(Doc.nm,".docx"),sep='/'))
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
  plot(density(d,adjust=1,na.rm=T),main=paste("log=",Log,"  ",Var," (range=",a[1],"-",a[2],")",sep=""),
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

#gamm approach
dev.expl=function(MOD) 100*(MOD$null.deviance-MOD$deviance)/MOD$null.deviance
Mod.fn.gamm=function(d,ResVar,Expl.vars,Predictrs,FactoRs,OFFSET,log.var,add.inter,MixedEff,temporal.autocorr)
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
  IDS=match(c("LATITUDE","LONGITUDE"),Predictrs)
  other.preds=Predictrs[-IDS]
  Pre.inter=Predictrs[IDS]
  if("MONTH" %in%other.preds)
  {
    d$MONTH=as.numeric(d$MONTH)
    other.preds[match('MONTH',other.preds)]="s(MONTH,k=12,bs='cc')"
  }
  if("BOTDEPTH"%in%other.preds) other.preds[match('BOTDEPTH',other.preds)]="s(BOTDEPTH)"
  if(!is.na(MixedEff))
  {
    Pred.1=paste(paste(c(other.preds[-match(MixedEff,other.preds)]),collapse="+"),
                 paste0("+s(",MixedEff,",bs='re')"),sep="")
    
  }
  if(is.na(MixedEff))
  {
    Pred.1=paste(other.preds,collapse="+")
  }
  if(add.inter=="YES")  Pred.form=paste0('s(',paste(Pre.inter,collapse=","),')') 
  if(!add.inter=="YES") Pred.form=paste(Pre.inter,collapse="+")
  if(!is.na(OFFSET)) Pred.form=paste(c(Pred.1,Pred.form,OFFSET),collapse="+")
  if(is.na(OFFSET))  Pred.form=paste(c(Pred.1,Pred.form),collapse="+")
  
  if(log.var=="YES")Formula=as.formula(paste("log(",ResVar,")", "~", Pred.form,collapse=NULL))
  if(log.var=="NO")Formula=as.formula(paste(ResVar, "~", Pred.form,collapse=NULL))
  
  
  #Fit model
    #full
  if(temporal.autocorr) 
  {
    model=gamm(Formula, correlation = corARMA(form = ~ 1|YEAR, p = 1), data=d)
  }else
  {
    #model=gamm(Formula, data=d)
    model=gam(Formula,data=d, method = 'REML')
  }
    #terms deviance
  Full.model.dev.expl=dev.expl(model)
  each.term=strsplit(Pred.form, "\\+")[[1]]
  mods.terms.dev=vector('list',length(each.term))
  for(m in 1:length(mods.terms.dev))
  {
    if(m==1)dummiies=gam(as.formula(paste(ResVar, "~", paste(each.term[m],collapse='+'),collapse=NULL)),data=d, method = 'REML')
    if(m>1)dummiies=gam(as.formula(paste(ResVar, "~", paste(each.term[1:m],collapse='+'),collapse=NULL)),data=d, method = 'REML')
    
    mods.terms.dev[[m]]=data.frame(Term=gsub("[()]","",str_remove(gsub("\\,.*", "", each.term[m]),'s')),
                                   Dev.explained=dev.expl(dummiies))
  } 
  mods.terms.dev=do.call(rbind,mods.terms.dev)%>%
                    mutate(lagged=lag(Dev.explained),
                           lagged=ifelse(is.na(lagged),0,lagged),
                           Dev.explained=Dev.explained-lagged)%>%
                    dplyr::select(-lagged)
  
  
  return(list(model=model, data=d,mods.terms.dev=mods.terms.dev,Full.model.dev.expl=Full.model.dev.expl))
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
fun.pred=function(d,Show.pred,normalised,PredictorS,log.var,MDL,BY=1,
                  add.dis.long.to.pred=FALSE,add.dis.lat.to.pred=FALSE)
{
  #create new data
  Dat=d$data
  Other.preds=c(PredictorS[-match(Show.pred,PredictorS)])
  if(!is.na(OFFSETT))Other.preds=c(PredictorS[-match("year",PredictorS)],"log.EFFORT")
  
  idPred=match(Show.pred,names(Dat))
  if(any(is.factor(Dat[,idPred])))
  {
    NEWDATA=data.frame(x=sort(unique(as.character(Dat[,idPred]))))
    NEWDATA[,1]=factor(NEWDATA[,1],levels=levels(Dat[,idPred]))
    names(NEWDATA)=Show.pred
  }else
  {
    if(length(idPred)==1)NEWDATA=data.frame(x=sort(seq(min(Dat[,idPred]),max(Dat[,idPred]),by=BY)))
    if(length(idPred)>1)
    {
      NEWDATA=Dat[,idPred]%>%mutate(dummy=paste(round(LATITUDE,2),round(LONGITUDE,2)))%>%
                                      distinct(dummy,.keep_all = T)%>%
                                      dplyr::select(-dummy)
    }
      
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
  if(MDL=='gamm')
  {
    mdl=d$model
    if("gamm"%in%class(mdl)) mdl=d$model$gam
      
    Covar.pos=as.matrix(vcov(mdl))
    Koef=coef(mdl)
    model=mdl
  }
    
  if(MDL=='glm')
  {
    Covar.pos=as.matrix(vcov(d$model))
    Koef=coef(d$model)
    model=d$model
  }
    
  set.seed(999);Pos.pars.rand=rmvnorm(niter,mean=Koef,sigma=Covar.pos)    
  
  if(any(!add.dis.long.to.pred==FALSE))
  {
    ddummi=vector('list',length(add.dis.long.to.pred))
    for(kk in 1:length(ddummi))
    {
      ddummi[[kk]]=NEWDATA%>%mutate(LATITUDE=add.dis.lat.to.pred[kk],LONGITUDE=add.dis.long.to.pred[kk])
    }
    NEWDATA=do.call(rbind,ddummi)
  }
  MC.preds=matrix(nrow=niter,ncol=length(NEWDATA[,match(Show.pred[1],names(NEWDATA))]))
  
  for(n in 1:niter)
  {
    
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
  names(PRD)[1:length(Show.pred)]=Show.pred
  
  #standardise to a mean score of 1
  if(normalised=="YES") PRD=fn.relative(PRD)
  if(length(Show.pred)==1) XX=as.numeric(as.character(PRD[,match(Show.pred,names(PRD))]))
  if(length(Show.pred)==2)
  {
    XX=PRD[,match(Show.pred,names(PRD))]
    XX=paste(XX[,1],XX[,2])
    Show.pred=paste(Show.pred,collapse=' ')
  }
    
  MeAn=PRD$MEAN
  UppCI=PRD$UP
  LowCI=PRD$LOW
  dat.plt=data.frame(Variable=Show.pred,Value=XX,MeAn=MeAn,UppCI=UppCI,LowCI=LowCI)
  
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
  
  #add location
  dat.plt=dat.plt%>%mutate(LATITUDE1=NEWDATA$LATITUDE, LONGITUDE1=NEWDATA$LONGITUDE)

  return(dat.plt)
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

fn.reshp=function(d,Y,TimeVar,IdVAR,do.reshape=FALSE)
{
  Form=formula(paste(Y,paste(c(TimeVar,IdVAR),collapse="+"),sep="~"))
  d=aggregate(Form,d,sum)
  if(do.reshape)
  {
    DATA.wide=reshape(d[,match(c(ResVar,IdVAR,TimeVar),names(d))],v.names=ResVar,
                      idvar=IdVAR,timevar=TimeVar,direction="wide")
    colnames(DATA.wide)=gsub(paste(ResVar,".",sep=""), "", names(DATA.wide))
  }else
  {
    DATA.wide=d%>%spread(TimeVar,ResVar)
  }
  DATA.wide[is.na(DATA.wide)]=0
  return(DATA.wide)
}
fn.viol.box=function(d,byzone,filcol=NULL,resp.vars)
{
  a=d%>%
    dplyr::select(c(resp.vars,YEAR,ZONE))%>%
    mutate(YEAR=as.character(YEAR))%>%
    gather(Indicator,Value,-c(YEAR,ZONE))%>%
    mutate(Indicator=case_when(Indicator=="FnRich_morph"~"Functional richness (morph.)",
                               Indicator=="FnRich_ecol"~"Functional richness (ecol.)",
                               Indicator=="MTL"~"Mean trophic level",
                               Indicator=="MML"~"Mean maximum length",
                               Indicator=="MeanML"~"Mean length",
                               Indicator=="Prop.Disc"~"Proportion of discards",
                               Indicator=="MaxAge"~"Maximum age",
                               Indicator=="Age.mat.prop"~"Proportion mature",
                               Indicator=="K"~"Growth coefficient",
                               Indicator=="MaxLen"~"Maximum length",
                               TRUE~Indicator))
  
  Yr.lev=seq(min(d$YEAR),max(d$YEAR))
  misn=Yr.lev[which(!Yr.lev%in%unique(d$YEAR))]
  a=rbind(a,
          expand.grid(YEAR=misn,
                      ZONE=unique(d$ZONE),
                      Indicator=unique(a$Indicator))%>%
            mutate(Value=NA))
  if(!byzone)
  {
    p=a%>%ggplot(aes(YEAR,Value))+
      geom_violin(width=1.5,fill='white')+
      geom_boxplot(width=0.5,fill=filcol)
  }
  if(byzone)
  {
    p=a%>%ggplot(aes(YEAR,Value,fill=ZONE))+
      geom_violin(width=1.5)+
      geom_boxplot(width=0.25)
  }
  p=p+
    facet_wrap(~Indicator,scales = 'free',ncol=2)+
    theme_PA()+
    theme(legend.position = 'top',
          axis.text.x = element_text(angle = 90, size=8,vjust = 0.5, hjust=1))+
    geom_smooth(method = "loess", se=TRUE,aes(group=1))+
    ylab('Ecological indicator value')+xlab('Year')
  
  return(p)
}

#Function for calculating ecological indicators
fn.calc.ecol.ind=function(DaTA,normalised,Drop.yrs,idvarS,resp.vars,TE=0.1,check.each.fun.rich.var=TRUE)
{
    #Define years to use
    if(Drop.yrs=="YES")DaTA=subset(DaTA,!YEAR%in%as.character(1993:1999))
    
    #select shots with a minimum number of individuals
    DaTA$ind.dummy=1
    N.ind.shot=aggregate(ind.dummy~SHEET_NO,DaTA,sum)
    N.ind.shot=subset(N.ind.shot,ind.dummy>=Min.individuals)
    DaTA=subset(DaTA,SHEET_NO%in%unique(N.ind.shot$SHEET_NO))
    
    # 3.1.2 Diversity and ecosystem indicators
    Dat=fn.reshp(d=DaTA,Y=ResVar,TimeVar=MultiVar,IdVAR=idvarS) 
    DATA.shots.diversity=Dat[,match(idvarS,names(Dat))]
    Dat.y=Dat[,-match(idvarS,names(Dat))]
    
    if("Shannon"%in%resp.vars)
    {
      DATA.shots.diversity$Shannon = diversity(Dat.y, index = "shannon")
      DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(Shannon))
    }
      
    if("Simpson"%in%resp.vars)
    {
      DATA.shots.diversity$Simpson = diversity(Dat.y, index = "simpson")
      DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(Simpson))
      DATA.shots.diversity=subset(DATA.shots.diversity,Simpson>0)
    }
    
    sp_names=colnames(Dat.y)
    
    if("Pielou"%in%resp.vars)
    {
      DATA.shots.diversity$Pielou = Evenness(H = DATA.shots.diversity$Shannon, S = length(sp_names))
      DATA.shots.diversity$Pielou[DATA.shots.diversity$Pielou == 0] = NA
      DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(Pielou))
    }

    Numbers=aggregate(INDIVIDUALS~SHEET_NO+SPECIES,DaTA,sum,na.rm=T)
    
    if("MTL"%in%resp.vars)
    {
      MTL=aggregate(TROPHIC_LEVEL~SHEET_NO+SPECIES,DaTA,mean,na.rm=T) 
      MTL=merge(MTL,Numbers,by=c("SHEET_NO","SPECIES"),all=T)
      MTL <- data.table(MTL)
      MTL=as.data.frame(MTL[,list(MTL = weighted.mean(TROPHIC_LEVEL,INDIVIDUALS,na.rm=T)),by=SHEET_NO])
      DATA.shots.diversity=merge(DATA.shots.diversity,MTL,by=c("SHEET_NO"),all.x=T)
      DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(MTL))
    }
    
    if("MML"%in%resp.vars)
    {
      MML=aggregate(TL~SHEET_NO+SPECIES,DaTA,max,na.rm=T)
      MML=merge(MML,Numbers,by=c("SHEET_NO","SPECIES"),all.y=T)
      MML <- data.table(MML)
      MML=as.data.frame(MML[,list(MML = weighted.mean(TL,INDIVIDUALS,na.rm=T)),by=SHEET_NO])
      DATA.shots.diversity=merge(DATA.shots.diversity,MML,by=c("SHEET_NO"),all.x=T)
      DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(MML))
      DATA.shots.diversity=subset(DATA.shots.diversity,MML>0)
    }
    
    if("MeanML"%in%resp.vars)
    {
      MeanML=aggregate(TL~SHEET_NO+SPECIES,DaTA,mean,na.rm=T)
      MeanML=merge(MeanML,Numbers,by=c("SHEET_NO","SPECIES"),all.y=T)
      MeanML <- data.table(MeanML)
      MeanML=as.data.frame(MeanML[,list(MeanML = weighted.mean(TL,INDIVIDUALS,na.rm=T)),by=SHEET_NO])
      DATA.shots.diversity=merge(DATA.shots.diversity,MeanML,by=c("SHEET_NO"),all.x=T)
      DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(MeanML))
      DATA.shots.diversity=subset(DATA.shots.diversity,MeanML>0)
    }
    
    if("FIB"%in%resp.vars)
    {
      Reference.year=1976  
      Ktch=aggregate(INDIVIDUALS~SHEET_NO+YEAR,DaTA,sum,na.rm=T)
      Ktch.0=median(with(subset(Ktch,YEAR==Reference.year),INDIVIDUALS))
      MTL.0=median(MTL%>%filter(SHEET_NO%in%unique(DaTA%>%filter(YEAR==Reference.year)%>%pull(SHEET_NO)))%>%pull(MTL))
      FIB=MTL%>%
        left_join(Ktch,by='SHEET_NO')%>%
        mutate(FIB=log(INDIVIDUALS * (1/TE)^MTL) - log(Ktch.0 * (1/TE)^MTL.0))%>%
        dplyr::select(SHEET_NO,FIB)
      DATA.shots.diversity=merge(DATA.shots.diversity,FIB,by=c("SHEET_NO"),all.x=T)
      DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(FIB))

    }
      
    if("Prop.Disc"%in%resp.vars)
    {
      Prop.Dis=DaTA%>%
              filter(!is.na(FATE))%>%
              group_by(SHEET_NO,FATE)%>%
              tally()%>%
              spread(FATE,n)%>%
              mutate(N=D+C,
                     Prop.Disc=D/N)%>%
              dplyr::select(SHEET_NO,Prop.Disc)
      DATA.shots.diversity=merge(DATA.shots.diversity,Prop.Dis,by=c("SHEET_NO"),all.x=T)
    }
    
    if("MaxAge"%in%resp.vars)
    {
      MaxAge=aggregate(tmax~SHEET_NO+SPECIES,DaTA,max,na.rm=T)
      MaxAge=merge(MaxAge,Numbers,by=c("SHEET_NO","SPECIES"),all.y=T)
      MaxAge <- data.table(MaxAge)
      MaxAge=as.data.frame(MaxAge[,list(MaxAge = weighted.mean(tmax,INDIVIDUALS,na.rm=T)),by=SHEET_NO])
      DATA.shots.diversity=merge(DATA.shots.diversity,MaxAge,by=c("SHEET_NO"),all.x=T)
      DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(MaxAge))
      DATA.shots.diversity=subset(DATA.shots.diversity,MaxAge>0)
    }
    
    if("Age.mat.prop"%in%resp.vars)
    {
      DaTA$Age.mat.prop=DaTA$tm/DaTA$tmax
      Age.mat.prop=aggregate(Age.mat.prop~SHEET_NO+SPECIES,DaTA,mean,na.rm=T)
      Age.mat.prop=merge(Age.mat.prop,Numbers,by=c("SHEET_NO","SPECIES"),all.y=T)
      Age.mat.prop <- data.table(Age.mat.prop)
      Age.mat.prop=as.data.frame(Age.mat.prop[,list(Age.mat.prop = weighted.mean(Age.mat.prop,INDIVIDUALS,na.rm=T)),by=SHEET_NO])
      DATA.shots.diversity=merge(DATA.shots.diversity,Age.mat.prop,by=c("SHEET_NO"),all.x=T)
      DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(Age.mat.prop))
      DATA.shots.diversity=subset(DATA.shots.diversity,Age.mat.prop>0)
    }
    
    if("K"%in%resp.vars)
    {
      K=aggregate(K~SHEET_NO+SPECIES,DaTA,mean,na.rm=T)
      K=merge(K,Numbers,by=c("SHEET_NO","SPECIES"),all.y=T)
      K <- data.table(K)
      K=as.data.frame(K[,list(K = weighted.mean(K,INDIVIDUALS,na.rm=T)),by=SHEET_NO])
      DATA.shots.diversity=merge(DATA.shots.diversity,K,by=c("SHEET_NO"),all.x=T)
      DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(K))
      DATA.shots.diversity=subset(DATA.shots.diversity,K>0)
    }
    
    if("MaxLen"%in%resp.vars)
    {
      MaxLen=aggregate(Lm~SHEET_NO+SPECIES,DaTA,max,na.rm=T)
      MaxLen=merge(MaxLen,Numbers,by=c("SHEET_NO","SPECIES"),all.y=T)
      MaxLen <- data.table(MaxLen)
      MaxLen=as.data.frame(MaxLen[,list(MaxLen = weighted.mean(Lm,INDIVIDUALS,na.rm=T)),by=SHEET_NO])
      DATA.shots.diversity=merge(DATA.shots.diversity,MaxLen,by=c("SHEET_NO"),all.x=T)
      DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(MaxLen))
      DATA.shots.diversity=subset(DATA.shots.diversity,MaxLen>0)
    }
    
    if("FnRich_ecol"%in%resp.vars)  #takes 6 mins for Logbook
    {
      #traits data set
      dd=DaTA[,match(c('SPECIES',traits_ecol),names(DaTA))]%>%distinct()%>%filter(SPECIES%in%colnames(Dat.y))
      DrOp=which(is.na(dd), arr.ind=TRUE)
      DrOp=as.character(dd[DrOp[,1],'SPECIES'])
      traits_sp=dd%>%filter(!SPECIES%in%DrOp)
      row.names(traits_sp)=traits_sp$SPECIES
      traits_sp=traits_sp%>%    
        dplyr::select(-SPECIES)
      
      needed=FALSE #convert to all numeric (needed for fd_fric from fundiversity package)
      if(needed)
      {
        gg=colnames(traits_sp)
        for(g in 1:length(gg))
        {
          if(!is.numeric(traits_sp[,g]))
          {
            UU=traits_sp[,g]
            U=unique(UU)  
            Riplase=data.frame(U,1:length(U))
            Nm=colnames(traits_sp)[g]
            names(Riplase)=c(Nm,'Rep')
            UU=data.frame(UU)
            names(UU)=Nm
            UU=UU%>%left_join(Riplase,by=Nm)
            traits_sp[,g]=UU$Rep
          }
        }
        traits_sp=traits_sp%>%  
          as.matrix()
      }

      #abundance data set
      site_sp=Dat.y[,-match(DrOp,names(Dat.y))]
      site_sp=as.matrix(site_sp)
      rownames(site_sp)=paste('Site',1:nrow(site_sp))
      
      #use only species accounting for 99% of catch as very rare species stuff up calculations
      CumS=rev(sort(colSums(site_sp)))
      CumS1=cumsum(CumS/sum(CumS))
      dis.sp=names(which(CumS1<0.99))
      dis.sp=subset(dis.sp,!dis.sp%in%rownames(traits_sp)[(apply(traits_sp, 1, function(r) any(r %in% '')))])
      dis.sp=subset(dis.sp,!dis.sp%in%rownames(traits_sp)[(apply(traits_sp, 1, function(r) any(is.na(r))))])
      traits_sp1=traits_sp[match(dis.sp,rownames(traits_sp)),]
      site_sp1=site_sp[,match(dis.sp,colnames(site_sp))]
      
      #calculate functional diversity
      for(g in 1:ncol(traits_sp1)) if(!is.numeric(traits_sp1[,g])) traits_sp1[,g]=as.factor(traits_sp1[,g])
      ex1 <- dbFD(traits_sp1, site_sp1) 
      Dat$FnRich_ecol=ex1$FRic
      
      #Add to data set  
      DATA.shots.diversity=DATA.shots.diversity%>%left_join(Dat%>%dplyr::select(SHEET_NO,FnRich_ecol),by=c("SHEET_NO"))
      #DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(FnRich_ecol))
      
    }
    
    if("FnRich_morph"%in%resp.vars)
    {
      #traits data set
      dd=DaTA[,match(c('SPECIES',traits_morph),names(DaTA))]%>%distinct()%>%filter(SPECIES%in%colnames(Dat.y))
      DrOp=which(is.na(dd), arr.ind=TRUE)
      DrOp=as.character(dd[DrOp[,1],'SPECIES'])
      traits_sp=dd%>%filter(!SPECIES%in%DrOp)
      row.names(traits_sp)=traits_sp$SPECIES
      traits_sp=traits_sp%>%    
        dplyr::select(-SPECIES)
      
      needed=FALSE  #convert to all numeric (needed for fd_fric from fundiversity package)
      if(needed)
      {      gg=colnames(traits_sp)
      for(g in 1:length(gg))
      {
        if(!is.numeric(traits_sp[,g]))
        {
          UU=traits_sp[,g]
          U=unique(UU)  
          Riplase=data.frame(U,1:length(U))
          Nm=colnames(traits_sp)[g]
          names(Riplase)=c(Nm,'Rep')
          UU=data.frame(UU)
          names(UU)=Nm
          UU=UU%>%left_join(Riplase,by=Nm)
          traits_sp[,g]=UU$Rep
        }
      }
      traits_sp=traits_sp%>%  
        as.matrix()}
      
      #abundance data set
      site_sp=Dat.y[,-match(DrOp,names(Dat.y))]
      site_sp=as.matrix(site_sp)
      rownames(site_sp)=paste('Site',1:nrow(site_sp))
      
      #use only species accounting for 95% of catch as very rare species stuff up calculations
      CumS=rev(sort(colSums(site_sp)))
      CumS1=cumsum(CumS/sum(CumS))
      dis.sp=names(which(CumS1<0.99))
      dis.sp=subset(dis.sp,!dis.sp%in%rownames(traits_sp)[(apply(traits_sp, 1, function(r) any(r %in% '')))])
      dis.sp=subset(dis.sp,!dis.sp%in%rownames(traits_sp)[(apply(traits_sp, 1, function(r) any(is.na(r))))])
      traits_sp1=traits_sp[match(dis.sp,rownames(traits_sp)),]
      site_sp1=site_sp[,match(dis.sp,colnames(site_sp))]
      
      #remove 0 species sites
      s=site_sp1
      s[s>0]=1
      b=rowSums(s)
      id=which(b==0)
      if(length(id)>0)site_sp1=site_sp1[-id,]
      
      #calculate functional diversity
      for(g in 1:ncol(traits_sp1)) if(!is.numeric(traits_sp1[,g])) traits_sp1[,g]=as.factor(traits_sp1[,g])
      ex1 <- dbFD(traits_sp1, site_sp1) 
      dummy=Dat
      if(length(id)>0)dummy=dummy[-id,]
      dummy$FnRich_morph=ex1$FRic
      
      #Add to data set  
      DATA.shots.diversity=DATA.shots.diversity%>%left_join(dummy%>%dplyr::select(SHEET_NO,FnRich_morph),by=c("SHEET_NO"))
      #DATA.shots.diversity=subset(DATA.shots.diversity,!is.na(FnRich_morph))
    }
    
    if(check.each.fun.rich.var)
    {
      fun.rich.vars=DaTA%>%
        dplyr::select( all_of(c('SHEET_NO','SPECIES',traits_ecol,traits_morph)))%>%
        filter(SPECIES%in%colnames(site_sp1))
      
      #categorical to numeric
      gg=colnames(fun.rich.vars)
      for(g in 3:length(gg))
      {
        if(!is.numeric(fun.rich.vars[,g]))
        {
          UU=fun.rich.vars[,g]
          U=unique(UU)  
          Riplase=data.frame(U,1:length(U))
          Nm=colnames(fun.rich.vars)[g]
          names(Riplase)=c(Nm,'Rep')
          UU=data.frame(UU)
          names(UU)=Nm
          UU=UU%>%left_join(Riplase,by=Nm)
          fun.rich.vars[,g]=UU$Rep
        }
      }
      
      fun.rich.vars=fun.rich.vars%>%
        distinct(SHEET_NO,SPECIES,.keep_all = T)%>%
        left_join(Numbers,by=c("SHEET_NO","SPECIES"))
      
      fun.rich.vars=fun.rich.vars%>%
        dplyr::select(-SPECIES)%>%
        group_by(SHEET_NO)%>%
        summarise(across(everything(), ~weighted.mean(., w=INDIVIDUALS, na.rm=TRUE)))%>%
        dplyr::select(-INDIVIDUALS)
      
      DATA.shots.diversity=merge(DATA.shots.diversity,fun.rich.vars,by=c("SHEET_NO"),all.x=T)
      
    }
    
    #Remove boats with less than Min.recs
    AA=sort(table(DATA.shots.diversity$BOAT))
    AA=AA[AA>Min.recs]
    DATA.shots.diversity=subset(DATA.shots.diversity,BOAT%in%names(AA))

    #Remove years with less than Min.recs
    AA=sort(table(DATA.shots.diversity$YEAR))
    AA=AA[AA>Min.recs]
    DATA.shots.diversity=subset(DATA.shots.diversity,YEAR%in%names(AA))
    
    #Remove blocks with less than Min.recs
    AA=sort(table(DATA.shots.diversity$BLOCK))
    AA=AA[AA>Min.recs]
    DATA.shots.diversity=subset(DATA.shots.diversity,BLOCK%in%names(AA))
    
    
    return(DATA.shots.diversity)
}

#Function for looping over different observer data sets
fn.apply.model=function(DaTA,dat.nm,normalised,Drop.yrs,idvarS,resp.vars,
                        do.data.mining=FALSE,do.glm=FALSE,do.gamm=TRUE,ADD.INTER)
{
  #1. run model
  Store.mod.out=vector('list',length(resp.vars))
  names(Store.mod.out)=resp.vars
  n.rv=length(resp.vars)
  Res.var.in.log=rep("NO",n.rv)
    #gamm
  if(do.gamm)
  {
    for(i in 1:n.rv)
    {
      Store.mod.out[[i]]=Mod.fn.gamm(d=DaTA,
                                     ResVar=resp.vars[i],
                                     Expl.vars=subset(Expl.varS,Expl.varS%in%idvarS),
                                     Predictrs=subset(Predictors,Predictors%in%idvarS),
                                     FactoRs=FactoRS,
                                     OFFSET=OFFSETT,
                                     log.var=Res.var.in.log[i],
                                     add.inter=ADD.INTER,
                                     MixedEff=MixedEff,
                                     temporal.autocorr=FALSE)  #not needed, no strong pattern, see '_autocorrelation' in fit folder
    }
  }
    #glm
  if(do.glm)
  {
    for(i in 1:n.rv)
    {
      Store.mod.out[[i]]=Mod.fn.glm(d=DaTA,
                                     ResVar=resp.vars[i],
                                     Expl.vars=Expl.varS,
                                     Predictrs=Predictors,
                                     FactoRs=FactoRS,
                                     OFFSET=OFFSETT,
                                     log.var=Res.var.in.log[i],
                                     add.inter="YES",
                                     MixedEff=MixedEff)  
    }
  }
    #data mining
  if(do.data.mining)
  {
    for(i in 1:n.rv)Store.mod.out[[i]]=Mod.fn.mining(d=DATA.shots.diversity,ResVar=resp.vars[i],
                      Predictrs=Predictors,Y.type="Continuous",Prop.train=.7,ALL.models="NO",nboot=1)
    
  }
  
  #2. Fit diagnostics
  if(do.gamm)
  {
    for(i in 1:n.rv)
    {
      Da=Store.mod.out[[i]]$data
      fn.fig(paste("Fit/GAMM_",dat.nm,names(Store.mod.out)[i],"distribution response var",sep="_"),2000,2000) 
      hist(Da[,match(resp.vars[i],names(Da))],main=resp.vars[i])
      dev.off()
      
      Modl=Store.mod.out[[i]]$model
      
      if("gamm"%in%class(Modl))
      {
        fn.fig(paste("Fit/GAMM_",dat.nm,names(Store.mod.out)[i],"gamm_lme",sep="_"),2000,2000)
        fv <- exp(fitted(Modl$lme)) ## predicted values (including re)
        rsd <- (Modl$gam$y - fv)/sqrt(fv) ## Pearson residuals (Poisson case)
        op <- par(mfrow=c(1,2))
        qqnorm(rsd)
        plot(fv^.5,rsd)
        dev.off()
        
        fn.fig(paste("Fit/GAMM_",dat.nm,names(Store.mod.out)[i],"autocorrelation",sep="_"),2000,2000)
        par(mfcol=c(2,1))
        acf(residuals(Modl$gam),main="raw residual ACF (gamm)") 
        acf(residuals(Modl$lme),main="raw residual ACF (lme)")
        dev.off()
        
        Modl=Modl$gam
      }
      
      fn.fig(paste("Fit/GAMM_",dat.nm,names(Store.mod.out)[i],"continuous",sep="_"),2000,2000) 
      par(mfcol=c(2,2))
      plot(Modl,pages=1,all.terms=TRUE)
      dev.off()
      
      fn.fig(paste("Fit/GAMM_",dat.nm,names(Store.mod.out)[i],"factor",sep="_"),2000,2000)
      vis.gam(Modl)
      dev.off()
      
      fn.fig(paste("Fit/GAMM_",dat.nm,names(Store.mod.out)[i],"fit",sep="_"),2000,2000) 
      par(mfcol=c(2,2))
      gam.check(Modl,pages=1)
      dev.off()
    }
  }
  if(do.glm)
  {
    for(i in 1:n.rv)
    {
      fn.fig(paste("GLM_Fit_Observer",names(Store.mod.out)[i],sep="_"),2000,2000) 
      par(mfcol=c(2,2))
      plot(Store.mod.out[[i]]$model)
      dev.off()
    }
  }

  #3. Anova tables
  TABL=vector('list',length(resp.vars))
  names(TABL)=resp.vars
  for(i in 1:n.rv)
  {
    Modl=Store.mod.out[[i]]$model
    if("gamm"%in%class(Modl)) Modl=Modl$gam
    if(do.gamm)
    {
      #each term
      Anova.tab=anova(Modl)
      
      Anova.smoothers=Anova.tab$s.table%>%
        data.frame()%>%
        tibble::rownames_to_column(var = "Term")%>%
        relocate(Term)%>%
        mutate(p.value=format(p.value, scientific = T, big.mark = ","),
               across(where(is.numeric), round, 3))
      Anova.factors=Anova.tab$pTerms.table%>%
        data.frame()%>%
        rename(edf=df)%>%
        mutate(Ref.df=NA)%>%
        tibble::rownames_to_column(var = "Term")%>%
        relocate(names(Anova.smoothers))%>%
        mutate(p.value=format(p.value, scientific = T, big.mark = ","),
               across(where(is.numeric), round, 3))
      
      Anova.tab=rbind(Anova.factors,Anova.smoothers)
    }
    if(do.glm)
    {
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
    }
    #deviance.expl= round(dev.expl(Modl),1) 
    qq=Anova.tab%>%
                mutate(Term=gsub("[()]","",str_remove(gsub("\\,.*", "", Term),'s')),
                       Indicator=names(TABL)[i])%>%
                left_join(Store.mod.out[[i]]$mods.terms.dev,by='Term')%>%
                relocate(Indicator,Term,Dev.explained)
    qq1=qq[1,]
    qq1[,-1]=''
    qq1=qq1%>%
      mutate(Term='Total',
             Dev.explained=Store.mod.out[[i]]$Full.model.dev.expl)
    TABL[[i]]=rbind(qq,qq1)%>%mutate(Term=ifelse(Term=='LATITUDE','LATITUDE x LONGITUDE',Term))
    
  }
  TABL=do.call(rbind,TABL)
  if(do.gamm)
  {
    #export anova tables as word doc   
    Export.tbl(WD=getwd(),Tbl=TABL,Doc.nm=paste("Anova.table",dat.nm,sep="_"))
    write.csv(TABL,paste0(paste(paste(getwd(),"Anova.table",sep='/'),dat.nm,sep="_"),'.csv'),row.names = F)
  }
  if(do.glm)
  {
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
    Export.tbl(WD=getwd(),Tbl=TABL,Doc.nm=paste("Anova.table",dat.nm,sep="_"),caption=NA,paragph=NA,
               HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
               Zebra='NO',Zebra.col='grey60',Grid.col='black',
               Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
               HDR.names=c('TERM', resp.vars),
               HDR.span=LBL_1st,
               HDR.2nd=LBL_2nd)
  }
    
  #4. Predictions
  if(do.glm) MDL='glm'
  if(do.gamm) MDL='gamm'
  Store.preds=vector('list',length(resp.vars))
  names(Store.preds)=resp.vars
  for(i in 1:n.rv)
  {
    Year=fun.pred(d=Store.mod.out[[i]],  
                  Show.pred="YEAR",
                  normalised="NO",
                  PredictorS=subset(Predictors,Predictors%in%idvarS),
                  log.var=Res.var.in.log[i],
                  MDL=MDL)        #add.dis.long.to.pred=Pred.long,
                                    #add.dis.lat.to.pred=Pred.lat  
    
    Year.R=fun.pred(d=Store.mod.out[[i]],  
                  Show.pred="YEAR",
                  normalised="YES",
                  PredictorS=subset(Predictors,Predictors%in%idvarS),
                  log.var=Res.var.in.log[i],
                  MDL=MDL)
      
    Month=fun.pred(d=Store.mod.out[[i]],
                   Show.pred="MONTH",
                   normalised=normalised,
                   PredictorS=subset(Predictors,Predictors%in%idvarS),
                   log.var=Res.var.in.log[i],
                   MDL=MDL,
                   BY=1)
    
    Lat.Long=fun.pred(d=Store.mod.out[[i]],
                          Show.pred=c("LATITUDE","LONGITUDE"),
                          normalised=normalised,
                          PredictorS=subset(Predictors,Predictors%in%idvarS),
                          log.var=Res.var.in.log[i],
                          MDL=MDL,
                          BY=0.5)    
        
    aa=rbind(Year%>%mutate(Relative='NO'),
             Year.R%>%mutate(Relative='YES'),
             Month%>%mutate(Relative=normalised),
             Lat.Long%>%mutate(Relative=normalised))
    
    Store.preds[[i]]=aa%>%mutate(Indicator=resp.vars[i],
                                 LowCI=ifelse(LowCI<0,0,LowCI),
                                 MeAn=ifelse(MeAn<0,0,MeAn))
    rm(aa)
  }
  
  return(do.call(rbind,Store.preds))
}

#Function for plotting predictions
fn.plot.preds=function(d,scale.y.low=NULL,rotate.axis=FALSE,add.smoother=FALSE,YLAB,add.zone=FALSE)
{
  if(add.zone)
  {
    p=d%>%
      ggplot(aes(Value,MeAn))+
      geom_point(aes(color=Zone))+
      geom_errorbar(aes(color=Zone,ymin = LowCI, ymax = UppCI))+
      facet_wrap(~Indicator,scales='free_y',ncol=2)+
      geom_line(aes(color=Zone),alpha=.6,linetype='dotted')+
      theme_PA(strx.siz=10)+ylab(YLAB)+xlab('Year')+
      theme(legend.position = 'top',
            legend.title = element_blank())
    
    if(add.smoother)
    {
      p=p+
        geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95,alpha=0.35)+
        geom_point(aes(color=Zone))
    }
      
  }else
  {
    p=d%>%
      ggplot(aes(Value,MeAn))+
      geom_point()+
      geom_errorbar(aes(ymin = LowCI, ymax = UppCI))+
      facet_wrap(~Indicator,scales='free_y',ncol=2)+
      geom_line(alpha=.6,linetype='dotted')+
      theme_PA(strx.siz=10)+ylab(YLAB)+xlab('Year')
    if(add.smoother)
    {
      p=p+
        geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95,alpha=0.35)+
        geom_point()
    }
  }
  if(!is.null(scale.y.low)) p=p+scale_y_continuous(limits = c(scale.y.low, NA))
  if(rotate.axis) p=p+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  return(p)
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

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#multivariate analysis 
manyany.fixed = function(formula, fn, family="negative.binomial", data=NULL, composition = FALSE, block = NULL, get.what="details", var.power=NA, na.action = "na.exclude", ...)
{
  #MANYANY - applies a function of your choice to each column of YMAT and computes logLik by taxon.
  # FORMULA is the formula to use in the call to FM.
  # FN is a character vector giving the name of the function to be applied to each taxon.  e.g. "glm"
  # a FAMILY argument needs to be specified. It can be a list with different families for different variables (with length matching the number of columns in YMAT)
  # COMPOSITION is a logical switch indicating whether or not to do a compositional analysis (i.e., include a
  #  row effect in the model to account for changes in total abundance across samples).
  # VAR.POWER is needed for tweedie distributions - the power parameter needs to be specified here as well as in the family argument.
  # ... any further arguments required by FN.
  #
  # Examples:
  # require(mvabund)
  # data(spider)
  # abund=spider$abund
  # X=data.frame(spider$x)
  # 
  # a manygam:
  # library(mgcv)
  # ft=manyany(abund~s(soil.dry),"gam",data=X,family="poisson")
  # 
  # a manyglmm:
  # library(lme4)
  # gr = rep(1:2,each=14)
  # ft=manyany(abund~soil.dry+1|gr,"lmer",data=X,family="poisson")
  #
  ## A manyglm:
  # ft=manyany(abund~soil.dry,"glm",data=X,family="poisson")
  ## note this gives the same answer as:  
  # ft=manyglm(mvabund(abund)~X$soil,family="poisson")
  
  #Now compatible with ordinal data, but note that linear predictor is for the lowest
  #category only, whereas fitted is for the observed category.
  
  
  tol=1.e-8 #for truncation of linear predictor
  
  #get the formula for extracting the data frame
  if(fn=="gam")
    mFormula = mgcv::interpret.gam(formula)$fake.formula
  else
    mFormula = formula
  
  if(missing(data)) # Only coerce to model frame if not specified by the user.
  {
    mf =model.frame(mFormula, parent.frame())
  }
  else
    mf = model.frame(mFormula,data=data)
  
  # get response and its dimensions
  yMat = model.response(mf)
  yMat = as.matrix(yMat)
  n.rows = dim(yMat)[1]
  n.vars = dim(yMat)[2]
  
  # for clm, need to change family argument and make each column of response a factor
  if(fn=="clm")
  {
    allargs <- match.call(expand.dots = FALSE)
    dots <- allargs$...
    if( "link" %in% names(dots) )
      link <- dots$link
    else
      link="logit"
    if(link=="loglog")
      fam.i = binomial("cloglog") #to avoid errors since "loglog" is not defined
    else
      fam.i = binomial(link) #although not binomial
    fam.i$family = "ordinal"
    fam = vector(mode="list",length=n.vars)
    family=fam
    # convert data to a dataframe of factors:
    yMat = data.frame(yMat, stringsAsFactors=TRUE) #converting to data frame so factor input is read as factors
    for(iVar in 1:n.vars)
      yMat[,iVar] = as.factor(yMat[,iVar])
  }
  
  # get names for response, or assign if empty
  yNames = dimnames(yMat)
  if(is.null(yNames[[1]]))
    yNames[[1]] = 1:n.rows
  if(length(yNames)==1)
    yNames[[2]] = paste("y",1:n.vars,sep="")
  #    yNames[[2]] = "y" #to avoid issues later.
  
  call=match.call()
  
  # fix this bit later  
  if(composition==TRUE)
  {
    yVec    = as.vector(yMat)
    rows     = factor(rep(1:n.rows,n.vars))
    cols     = factor(rep(1:n.vars,each=n.rows))
    mf    = data.frame(mf[rows,],rows,cols)
    formula = formula(paste(formula[2],"~rows+cols+cols:(",formula[3],")",sep=""))
    n.rows.orig = n.rows #save for later
    n.vars.orig = n.vars #save for later
    n.rows  = length(yVec)
    n.vars  = 1
    names(yVec) = paste( yNames[[2]][cols], ".", yNames[[1]][rows], sep="")
    yMat    = as.matrix(yVec)
    if(inherits(family,"family")==FALSE & length(family)>1)
      stop("when using composition=TRUE, family argument must have length one.")
    if(is.null(block))  #to make sure resampling is by row of original data, not of vectorised data.
      block = rows
    else
      block = block[rows]
  }
  
  #If family is specified once, turn it into a list of n.vars family objects
  if(inherits(family,"family") || length(family)==1)
  {
    fam = family #temporary store to slot into a big list
    family = vector("list",n.vars)
    for(i.var in 1:n.vars)
      family[[i.var]] = fam
  }
  if(length(family)!=n.vars)
    stop("family argument has length more than one but not equal to the number of columns in yMat (!?)")
  
  if(length(var.power)==1)
    var.power=rep(var.power,n.vars)
  
  fam = family
  
  # now ensure each family is a proper family function
  for(i.var in 1:n.vars)
  {
    if (is.character(family[[i.var]])) 
    {
      if (family[[i.var]] == "negbinomial" || family[[i.var]]=="negative.binomial")
      {
        fam[[i.var]] = negative.binomial(10^6)
        fam[[i.var]]$family = family[[i.var]]
      }
      else if (family[[i.var]] == "binomial(link=logit)")
      {
        fam[[i.var]] = binomial()
        fam[[i.var]]$family = family[[i.var]]
      }
      else if (family[[i.var]] == "binomial(link=cloglog)" || family[[i.var]] == "cloglog")
      {
        fam[[i.var]] = binomial("cloglog")
        fam[[i.var]]$family = family[[i.var]]
      }
      else
      {
        fam.fn = get(fam[[i.var]], mode = "function", envir = parent.frame())        
        fam[[i.var]] = fam.fn()
      }  
    }
    if(fn=="clm")
      fam[[i.var]] = family[[i.var]] = fam.i
    if(fam[[i.var]]$family=="binomial")
      warning("The binomial option of manyany currently assumes you have binary (presence/absence) response")
  }
  
  # find response variable in mf (should be first but better safe than sorry)
  nameOfResponse  = as.character(formula[[2]])
  whichIsResponse = which(names(mf)==nameOfResponse)
  
  # set up empty objects
  manyfit = vector(mode="list",length=n.vars)
  fits = matrix(NA,n.rows,n.vars)
  etas = matrix(NA,n.rows,n.vars)
  params = manyfit
  logL = rep(NA,n.vars)
  
  # fit model sequentially for each variable
  for(i.var in 1:n.vars)
  {
    # change response to just column iVar
    mf[[1]] = yMat[,i.var]
    
    # refit model via do.call
    manyfit[[i.var]] = do.call(fn, list(formula=formula, family=family[[i.var]], data=mf, na.action=na.action, ...)) #note use of family argument as originally specified
    
    # store logL, or get from dviance if undefined
    logL[i.var]  = logLik(manyfit[[i.var]])
    if(is.na(logL[i.var]))
      logL[i.var] = -0.5*deviance(manyfit[[i.var]]) #just in case logL function is undefined, e.g. tweedie 
    
    # only get extra stuff if get.what says to... this is skipped by anova.manyany
    if(get.what=="details"||get.what=="models")
    {
      fits[,i.var] = fitted(manyfit[[i.var]])
      
      etas[,i.var] = switch(fn,
                            "lmer"=manyfit[[i.var]]@eta,
                            "clm"=predict(manyfit[[i.var]],type="linear.predictor",newdata=mf)$eta1,
                            predict(manyfit[[i.var]])
      )
      #need to then truncate as if on logit scale...
      if(substr(fam[[i.var]]$family,1,3)=="bin" || fam[[i.var]]$family=="ordinal") #truncate linear predictor to more reasonable range
      {
        etas[,i.var] = pmax(etas[,i.var], fam[[i.var]]$linkfun(tol)/2)
        etas[,i.var] = pmin(etas[,i.var], fam[[i.var]]$linkfun(1-tol)/2)
      }
      if(fam[[i.var]]$link=="log"||fam[[i.var]]$link=="mu^0") #truncate linear predictor to more reasonable range
        etas[,i.var] = pmax(etas[,i.var], log(tol)/2)
      if(i.var==1)
      {
        cf = try(coef(manyfit[[i.var]]),silent=TRUE) #don't know if this function is defined
        if(inherits(cf, "try-error"))
        {
          do.coef   = FALSE
          coefs     = NULL
        } 
        else
        {
          coefs     = vector(mode="list",n.vars)
          coefs[[1]] = cf
          names(coefs[[1]])=dimnames(cf)[[1]]
          if(composition==FALSE & n.vars>1) #only name y variable if not compositional model
            names(coefs)=yNames[[2]]
          do.coef   = TRUE
        }
      }
      else
      {
        if(do.coef==TRUE)
          coefs[[i.var]] = coef(manyfit[[i.var]])
      }
      if(fam[[i.var]]$family=="poisson")
        params[[i.var]] = list(q=yMat[,i.var],lambda=fits[,i.var])
      if(substr(fam[[i.var]]$family,1,3)=="bin")
        params[[i.var]] = list(q=yMat[,i.var],prob=fits[,i.var],size=1)
      if(fam[[i.var]]$family=="Tweedie")
        params[[i.var]] = list(q=yMat[,i.var], power=var.power[i.var], mu=fits[,i.var], phi=summary(manyfit[[i.var]])$disp)
      if(fam[[i.var]]$family=="ordinal")
        params[[i.var]] = list(q=yMat[,i.var], mu=predict(manyfit[[i.var]], type="cum.prob"), muAll=predict(manyfit[[i.var]],type="cum.prob",newdata=mf[-whichIsResponse])$cprob2)
      if(grepl("egative",fam[[i.var]]$family) || fam[[i.var]]$family == "negbinomial")
      {
        if(any(names(manyfit[[i.var]])=="theta"))
          theta=manyfit[[i.var]]$theta
        else
        {
          if(any(names(manyfit[[i.var]])=="phi"))
            theta = 1/manyfit[[i.var]]$phi
          else # otherwise it must be fixed and tied up in the family argument 
            theta = 1/(fam[[i.var]]$var(1)-1)
        }
        params[[i.var]] = list(q=yMat[,i.var],mu=fits[,i.var],size=theta)
      }
      if(fam[[i.var]]$family=="gaussian")
      {
        s.ft=summary(manyfit[[i.var]])
        if(any(names(s.ft)=="sigma"))
          sd=s.ft$sigma
        else
          sd=s.ft$scale
        params[[i.var]] = list(q=yMat[,i.var],mean=fits[,i.var],sd=sd)
      }
    } #end get.what if statement
  } #end i.var loop
  object=list(logL=logL, get.what=get.what)
  
  #now format predictions and get residuals, if required
  if(get.what=="details"||get.what=="models")
  {
    if(composition==TRUE) #reshape to original data size if required
    {
      fits   = matrix(fits,n.rows.orig,n.vars.orig)
      etas   = matrix(etas,n.rows.orig,n.vars.orig)
    }    
    else #only name y variables for logL and params if not compositional model
    {
      if(n.vars>1)
      {
        names(logL) = yNames[[2]]
        names(params) = yNames[[2]]
      }
    }
    attributes(logL)$df = attributes(logLik(manyfit[[i.var]]))$df
    attributes(logL)$nobs = n.rows
    class(logL) = "logLik"
    resids = residuals.manyany(list(params=params, family=fam, composition=composition, fitted.values=fits, get.what=get.what))
    dimnames(resids) = yNames
    dimnames(fits)   = yNames
    dimnames(etas)   = yNames
    mf[[1]] = yMat #DW, 3/2/22 change: return full response
    object=list(logL=logL,fitted.values=fits,residuals=resids,linear.predictor=etas,family=fam, coefficients = coefs, call=call,params=params,model=mf, terms = terms(manyfit[[i.var]]), formula=formula, block=block, composition=composition, get.what=get.what)
    #    object=list(logL=logL,fitted.values=fits,residuals=resids,linear.predictor=etas,family=fam, coefficients = coefs, call=call,params=params,model=model.frame(manyfit[[i.var]]), terms = terms(manyfit[[i.var]]), formula=formula, block=block, composition=composition, get.what=get.what)
  }
  if(get.what=="models") #also return the model fits, if requested
  {
    object$fits = manyfit
    names(object$fits) = yNames[[2]]
  }
  class(object)=c("manyany", class(manyfit[[i.var]]) )
  return(object)
} 
residuals.manyany<- function(object, ...)
{
  if(object$get.what!="details" & object$get.what!="models")
    stop("To compute residuals, set get.what='details' in your manyany call")
  tol=1.e-6
  params = object$params
  n.rows = length(params[[1]]$q)
  n.vars = length(params)
  if(length(object$family)==1)
    family = rep(object$family,n.vars)
  else
    family=object$family
  resids=matrix(NA,n.rows,n.vars)
  dimnames(resids)[[1]] = names(params[[1]]$yMat)
  dimnames(resids)[[2]] = names(params)
  for(i.var in 1:n.vars)
  {
    if(family[[i.var]]$family=="ordinal")
    {
      u = runif(n.rows)
      resids[,i.var] = u*params[[i.var]]$mu$cprob1 + (1-u)*params[[i.var]]$mu$cprob2      
    }
    else
    {
      param.minus = params[[i.var]]
      param.minus$q = params[[i.var]]$q - tol
      if(grepl("egative",family[[i.var]]$family) || family[[i.var]]$family == "negbinomial")
        pfn = "pnbinom"
      if(family[[i.var]]$family=="poisson")
        pfn = "ppois"
      if(substr(family[[i.var]]$family,1,3)=="bin")
        pfn = "pbinom"
      if(family[[i.var]]$family=="gaussian")
      {
        pfn = "pnorm"
        param.minus$q = params[[i.var]]$q
      } 
      if(family[[i.var]]$family=="Tweedie")
        pfn = "ptweedie"
      u = runif(n.rows)
      #to avoid any values identically 1:
      pMinus = pmin(do.call(pfn, param.minus), 1-tol)
      resids[,i.var] = u*do.call(pfn, params[[i.var]]) + (1-u)*pMinus
    }
  }
  if(object$composition==TRUE) #reshape to original data size if required
    resids = matrix(resids, dim(object$fitted)[1], dim(object$fitted)[2])
  resids=qnorm(resids)
  return(resids)
}
multivariate.fn=function(d,Terms,Def.sp.term,Transf,Show.term,Group,hndl,simper.cumsum=0.85,LGSIZE=12,
                         All.species.names,dat.name,do.MDS=TRUE,do.pcoa=FALSE,
                         group.ordination=TRUE,Group.term.ordination=c('year','latitude','longitude'),
                         do.boral=FALSE,do.mvabund=TRUE,do.permanova=FALSE,do.simper=FALSE,
                         do.gam=FALSE,do.tweedie=TRUE,do.lm=FALSE,aggregate.monthly=TRUE,
                         term.form.model,term.form.permanova,use.Other=TRUE,dis.lat,dis.long,do.anova=FALSE)
{
  names(d)=tolower(names(d))
  
  if(aggregate.monthly)
  {
    Terms=Group.term.ordination
    d=d%>%
      mutate(latitude=round(latitude),
             longitude=round(longitude))%>%
      group_by_at(c('species',Terms))%>%
      summarise(individuals=mean(individuals,na.rm=T),
                ktch=sum(individuals*effort,na.rm=T))
  }
  
  #Define what species are grouped
  what.species=d%>%
    group_by_at(Def.sp.term)%>%
    summarise(n=sum(ktch))%>%
    spread(year,n,fill=0)%>%
    ungroup()%>%data.frame
  what.species$n=rowSums(what.species[,-1])

  what.species=what.species%>%    
    arrange(-n)%>%
    mutate(Row=1:n(),
           Cumktch=cumsum(n),
           Cumktch=Cumktch/sum(n))
  
  if(Group=='90') what.species=what.species%>%mutate(species2=ifelse(Cumktch<=0.9,species,'Other'))
  if(Group=='95') what.species=what.species%>%mutate(species2=ifelse(Cumktch<=0.95,species,'Other'))
  if(Group=='Top20') what.species=what.species%>%mutate(species2=ifelse(Row<=20,species,'Other'))
  what.species=what.species%>%
    distinct(species,species2)
  if(!use.Other) what.species=what.species%>%filter(!species2=='Other')
  
  
  d=d%>%
    left_join(what.species,by=c('species'))%>%
    filter(!is.na(species2))%>%
    group_by_at(c(Terms,'species2'))%>%
    summarise(individuals=sum(individuals,na.rm=T))%>%
    spread(species2,individuals,fill=0)%>%
    ungroup()
  
  d$ColSum=rowSums(d[-match(Terms,names(d))])
  d=d%>%filter(ColSum<100)
  d=d%>%filter(!is.na(ColSum))
  d=d%>%filter(ColSum>0)

  if(Transf=='proportion')
  {
    d[-match(c(Terms,'ColSum'),names(d))]=d[-match(c(Terms,'ColSum'),names(d))]/d$ColSum
  }
    
  d=d%>%
    dplyr::select(-ColSum)
  
  dd = d[,-match(Terms,names(d))]
  drop.this=which(rowSums(dd)==0)
  if(length(drop.this)>0) d=d[-drop.this,]
  rm(dd)
  rownames(d)=paste0('id.',1:nrow(d))
  Community = d[,-match(Terms,names(d))]
  if(Transf=='sqrt')Community = sqrt(Community)
  
  Terms.part=d[,match(Terms,names(d))]

  #1. Ordination
    #1.1 MDS 
  if(do.MDS)
  {
    Comm=Community
    
    if(group.ordination)
    {
      Comm <-cbind(Comm,d[,match(Group.term.ordination,names(d))])
     if('latitude'%in%Group.term.ordination) 
     {
       Comm=Comm%>%
         mutate(latitude=round(latitude),
                longitude=round(longitude))
     }
      Comm=Comm%>%
        group_by_at(Group.term.ordination)%>%
        summarise_all("mean")
      Terms.part.agg=Comm[,match(Group.term.ordination,names(Comm))]
      Comm = Comm[,-match(Group.term.ordination,names(Comm))]
    }
    
    MDS <- metaMDS(comm = Comm, distance = "bray",k=2,trymax=50, trace = T, autotransform = FALSE)
    mds.envfit <- envfit(MDS, Terms.part.agg, permutations = 999) # this fits environmental vectors
    mds.spp.fit <- envfit(MDS, Comm, permutations = 999) # this fits species vectors
    
    #create dataframes
    site.scrs <- as.data.frame(scores(MDS, display = "sites")) #save NMDS results into dataframe
    site.scrs <- cbind(site.scrs,Terms.part.agg) #add grouping variable "ecosystem" to dataframe
    
    #plot it up
    Kls=colfunc(length(levels(Terms.part.agg$year)))
    names(Kls)=levels(Terms.part.agg$year)
    nmds <- site.scrs%>%
            ggplot(aes(x=NMDS1, y=NMDS2,colour = year))+ #sets up the plot
            #stat_ellipse(aes(x=NMDS1,y=NMDS2,colour=year),level = 0.50,linewidth=1.2)+
            geom_point(aes(NMDS1, NMDS2, colour = year),alpha = 0.6, size = 4)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
            coord_fixed()+
            annotate(geom="text", x=0.85*max(site.scrs$NMDS1,na.rm=T), y=min(site.scrs$NMDS2,na.rm=T), 
                     label=paste("Stress=",round(MDS$stress,3)))+
            theme_PA(leg.siz=14,axs.t.siz=12,axs.T.siz=14)+
            theme(legend.title = element_blank(),
                  legend.position = 'top')+
            xlab('')+ylab('')+
            scale_color_manual(values=Kls)
    print(nmds)
    ggsave(paste(hndl,'/MDS_',dat.name,'.tiff',sep=""),width = 6,height = 6,compression = "lzw")
    
  }
  
    #1.2 PCOA
  if(do.pcoa)
  {
    Comm=Community
    Comm$dummy <- d[,match(Show.term,names(d))]
    if(group.ordination)
    {
      Comm=Comm%>%
        group_by(dummy) %>%
        summarise_all("mean")
    }
      
    
    Community.bray <- vegan::vegdist(Comm%>%dplyr::select(-dummy), method = "bray") 
    pcoaVS <- pco(Community.bray, negvals = "zero", dround = 0) 
    
    PCOA_xy <-  pcoaVS$vectors[,1:2]
    PCOA_xy <- as.data.frame(PCOA_xy)
    PCOA_xy$dummy <- Comm$dummy
    PCOA_xy$kolor=colorRampPalette(c("cyan4","cornflowerblue", "blue4"))(nrow(PCOA_xy))
    
    if(group.ordination)
    {
      p=PCOA_xy%>%                 
        ggplot(aes(x = V1, y = V2)) +
        geom_point(color = "transparent") +
        geom_segment(aes(xend = after_stat(lead(x)), yend = after_stat(lead(y)),colour = kolor),
                     arrow = arrow(length = unit(4, "mm")),linewidth=1.25) +
        geom_text(aes(label = dummy), size = 5, fontface = 2,alpha=0.6) +
        labs(x = "PCOA1", y = "PCOA2")+
        theme_PA(leg.siz=14,axs.t.siz=12,axs.T.siz=14)+
        theme(legend.position = "none")+scale_colour_identity() 
    }else
    {
      p=PCOA_xy%>%
        mutate(dummy=as.numeric(as.character(dummy)))%>%
        ggplot(aes(x = V1, y = V2)) +
        geom_point(aes(color = dummy),alpha = 0.4, size = 4) +
        labs(x = "PCOA1", y = "PCOA2")+
        theme_PA(leg.siz=14,axs.t.siz=12,axs.T.siz=14)
    }
    
    print(p)
    ggsave(paste(hndl,'/PCOA_',dat.name,'.tiff',sep=""),width = 6,height = 6,compression = "lzw")
    
  }
  
  #2. Multi stats 
  NEWDATA=data.frame(unique(Terms.part[,Show.term]))
  colnames(NEWDATA)=Show.term
  add.dis1=Terms.part[,-match(Show.term,names(Terms.part))]%>%data.frame
  add.dis=add.dis1[1,]
  for(o in 1:ncol(add.dis))
  {
    if(is.character(add.dis1[,o])|is.factor(add.dis1[,o])) sss=Mode(add.dis1[,o])
    if(is.numeric(add.dis1[,o])) sss=mean(add.dis1[,o],na.rm=T)
    add.dis[,o]=sss
  }
  if(length(dis.lat)>1)
  {
    add.dis=rbind(add.dis, add.dis[rep(1, length(dis.lat)-1), ])%>%
      mutate(latitude=dis.lat,
             longitude=dis.long)
  }else
  {
    add.dis$latitude=dis.lat 
    add.dis$longitude=dis.long
  }
  NEWDATA=tidyr::crossing(NEWDATA, add.dis)
  
  
  #2.1 boral   #takes 0.36 sec per row   #not working properly, cannot predict 'Newdata'
  if(do.boral)
  {
    example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, n.thin = 1)
    testpath <- file.path(tempdir(), "jagsboralmodel.txt")
    spp_mat <- as.matrix(Community)
    
    if(Transf=='proportion')
    {
      spp_mat[spp_mat == 0] <- 1e-4   #beta cannot take exact 0 or 1
      spp_mat[spp_mat == 1] <- 0.99999
      mod <- boral(spp_mat, X = Terms.part, formula.X = paste('~',term.form.model),
                        family = "beta", lv.control = list(num.lv = 2),save.model = T,
                        mcmc.control = example_mcmc_control,model.name = testpath)
    }else
    {
      mod <- boral(spp_mat, X = Terms.part, formula.X = paste('~',term.form.model),
                   family = "tweedie", lv.control = list(num.lv = 2),save.model = T, 
                   mcmc.control = example_mcmc_control, model.name = testpath)
      
    }

    #Predict year effect
    Preds=predict(object=mod, predict.type = "conditional",scale = "response")
    Preds=as.data.frame(Preds$linpred)
    names(Preds)=colnames(spp_mat)
    Preds=cbind(Terms.part,Preds)
 
    a=Preds%>%
      dplyr::select(colnames(spp_mat),year,latitude,longitude)%>%
      gather(species,Median,-year,-latitude, -longitude)%>%
      filter(latitude%in%dis.lat & longitude%in%dis.long)%>%
      mutate(species=str_remove(species, 'X'),
             latitude=round(latitude),
             longitude=round(longitude))%>%
      group_by(year,latitude,longitude,species)%>%
      summarise(Median=mean(Median))

    return(list(Preds=Preds))
  }
    #2.2 mvabund
  if(do.mvabund)
  {
    names(All.species.names)=tolower(names(All.species.names))
    
    spp_mat <- mvabund(as.matrix(Community))
    
    explore.dis=FALSE
    if(explore.dis)
    {
      par(mar = c(2, 10, 2, 2)) # adjusts the margins
      boxplot(Community, horizontal = TRUE, las = 2, main = "Abundance")
      meanvar.plot(spp_mat)  #check mean variance relationship
      plot(spp_mat ~ as.factor(d$year), cex.axis = 0.8, cex = 0.8)  
      
      left_join(d %>%
                  dplyr::select(-c(month,latitude,longitude))%>%
                  group_by(year) %>%
                  summarise_all(c(mean))%>%gather(species,mean,-year),
                d %>%
                  dplyr::select(-c(month,latitude,longitude))%>%
                  group_by(year) %>%
                  summarise_all(c(sd))%>%gather(species,sd,-year),
                by=c('year','species'))%>%
        mutate(lower95=mean -1.96*sd,
               upper95=mean +1.96*sd)%>%
        ggplot(aes(year,mean))+
        geom_bar(stat="identity")+
        facet_wrap(~species,scales='free_y')
    }
    
    null.mod=NULL
    if(do.gam) 
    {
       mod <- manyany.fixed(formula=formula(paste('spp_mat',term.form.model,sep='~')),fn='gam',    
                      family = gaussian(link = "identity"),data=Terms.part)
      #mod <- manyany.fixed(formula=formula(paste('spp_mat',term.form.model,sep='~')),fn='gam',    #doesn't work
      #                      family = betar(link="logit"),data=Terms.part)
       null.mod<-manyany.fixed(formula=spp_mat ~1,fn='gam',    
                              family = gaussian(link = "identity"),data=Terms.part)
    }
    if(do.tweedie) 
    {
       mod <- manyany(formula=formula(paste('spp_mat',term.form.model,sep='~')),fn="glm",data=Terms.part, 
                      family=tweedie(var.power=1.2, link.power=0), var.power=1.2) 
       
       if(do.anova)
       {
         null.mod<-manyany(formula=formula(paste('spp_mat',1,sep='~')),fn="glm",data=Terms.part,
                           family=tweedie(var.power=1.2, link.power=0), var.power=1.2)
         
         mod.alt<-manyany(formula=formula(spp_mat ~ year),fn="glm",data=Terms.part, 
                          family=tweedie(var.power=1.2, link.power=0), var.power=1.2)
         a=anova(null.mod, mod.alt, p.uni = "unadjusted", nBoot = 99)   
         anovt=data.frame(species=str_remove(names(a[[4]]), 'X'),LR=c(a[[3]]),p=a[[4]])%>%  #ACA
           left_join(All.species.names,by='species')%>%filter(!species=='Other')%>%
           relocate(scientific_name)%>%dplyr::select(-species)%>%
           mutate(p=ifelse(p==0.01,'<0.01',p))
         write.csv(anovt, paste0(hndl,"/Anova_year_",dat.name,".csv"),row.names=F)
         
         mod.alt<-manyany(formula=formula(spp_mat ~ latitude * longitude),fn="glm",data=Terms.part, 
                          family=tweedie(var.power=1.2, link.power=0), var.power=1.2)
         a=anova(null.mod, mod.alt, p.uni = "unadjusted", nBoot = 99)
         anovt=data.frame(species=str_remove(names(a[[4]]), 'X'),LR=c(a[[3]]),p=a[[4]])%>%
           left_join(All.species.names,by='species')%>%filter(!species=='Other')%>%
           relocate(scientific_name)%>%dplyr::select(-species)%>%
           mutate(p=ifelse(p==0.01,'<0.01',p))
         write.csv(anovt, paste0(hndl,"/Anova_lat.long_",dat.name,".csv"),row.names=F)
         
         rm(a)
       }

    }
    if(do.lm)
    {
      spp_mat.log=log(spp_mat+1e-6)
      mod <- manylm(formula=formula(paste('spp_mat.log',term.form.model,sep='~')),data=Terms.part)  
    }
    
    
    #Predict year effect 
    if(do.lm) Preds=predict(object=mod, newdata=NEWDATA, se.fit = TRUE,type = "response")
    if(do.tweedie) Preds=predict.manyglmTweedie(object=mod, newdata=NEWDATA, se.fit = TRUE,type = "response")
    if(do.gam) Preds=predict.manyany(object=mod,newdata=NEWDATA,type="response", se.fit = TRUE)
    
    Preds.median=Preds$fit%>%data.frame()
    Preds.se=Preds$se.fit%>%data.frame()
    names(Preds.se)=names(Preds.median)
    Preds.median=Preds.median%>%
      mutate(year=NEWDATA$year)
    if(length(dis.lat)>1)
    {
      Preds.median=Preds.median%>%
        mutate(latitude=NEWDATA$latitude,
               longitude=NEWDATA$longitude)%>%  
        gather(species,Median,-year,-latitude,-longitude)
    }else
    {
      Preds.median=Preds.median%>%
        gather(species,Median,-year)
    }
    
    Preds.se=Preds.se%>%
      mutate(year=NEWDATA$year)
    if(length(dis.lat)>1)
    {
      Preds.se=Preds.se%>%
        mutate(latitude=NEWDATA$latitude,
               longitude=NEWDATA$longitude)%>%  
        gather(species,SE,-year,-latitude,-longitude)
    }else
    {
      Preds.se=Preds.se%>%
        gather(species,SE,-year)
    }

    
    if(length(dis.lat)>1)
    {
      Preds.year=Preds.median%>%left_join(Preds.se,by=c('year','species','latitude','longitude'))
    }else
    {
      Preds.year=Preds.median%>%left_join(Preds.se,by=c('year','species'))
    }
    Preds.year=Preds.year%>%
      mutate(lower95=Median-1.96*SE,
             upper95=Median+1.96*SE,
             species=str_remove(species, 'X'))%>%
      left_join(All.species.names,by='species')
    
    #Predict latitude and longitude
    Preds.lat.lon=NULL
    do.dis=FALSE
    if(do.dis)
    {
      NEWDATA=Terms.part%>%
        mutate(latitude=round(latitude),
               longitude=round(longitude))%>%
        distinct(latitude,longitude)%>%
        mutate(latitude=latitude-0.5,
               longitude=longitude+0.5)
      
      add.dis1=Terms.part[,-match(c('latitude','longitude'),names(Terms.part))]
      add.dis=add.dis1[1,]
      for(o in 1:ncol(add.dis))
      {
        if(is.character(add.dis1[,o])|is.factor(add.dis1[,o])) sss=Mode(add.dis1[,o])
        if(is.numeric(add.dis1[,o])) sss=mean(add.dis1[,o],na.rm=T)
        add.dis[,o]=sss
      }
      suppressWarnings({NEWDATA=cbind(NEWDATA,add.dis)})
      if(do.lm) Preds=predict(mod, newdata=NEWDATA, se.fit = TRUE,type = "response")
      if(do.gam) Preds=predict.manyany(object=mod,newdata=NEWDATA,type="response", se.fit = TRUE)
      
      Preds.median=Preds$fit%>%data.frame()
      Preds.se=Preds$se.fit%>%data.frame()
      names(Preds.se)=names(Preds.median)
      Preds.median=Preds.median%>%
        mutate(longitude=NEWDATA$longitude,
               latitude=NEWDATA$latitude)%>%
        gather(species,Median,-longitude,-latitude)
      Preds.se=Preds.se%>%
        mutate(longitude=NEWDATA$longitude,
               latitude=NEWDATA$latitude)%>%
        gather(species,SE,-longitude,-latitude)
      Preds.lat.lon=Preds.median%>%left_join(Preds.se,by=c('longitude','latitude','species'))%>%
        mutate(lower95=Median-1.96*SE,
               upper95=Median+1.96*SE,
               species=str_remove(species, 'X'))%>%
        left_join(All.species.names,by='species')
      
    }
      
    
     return(list(null.mod=null.mod,mod=mod,Preds=Preds.year,Preds.lat.lon=Preds.lat.lon))
  }
  
    #2.3 Permanova  
  if(do.permanova)
  {
    # 2.1. overall significance test
    adon.results<-adonis2(formula(paste('Community',term.form.permanova,sep='~')),data=d, method="bray",perm=1e3,parallel=6)
    write.csv(as.data.frame(adon.results),paste(hndl,'/Permanova table_',dat.name,'.csv',sep=""))
    
    # 2.2. multilevel pairwise comparison with adjusted p-values
   do.this=TRUE
   if(do.this)
   {
     dummy=pairwise.adonis2(Community~year,data=d)
     adonis.pairwise=vector('list',(length(dummy)-1))
     for(qq in 2:length(dummy))
     {
       adonis.pairwise[[qq]]=data.frame(Pairs=names(dummy)[[qq]],
                                        P=dummy[[qq]]$`Pr(>F)`[1])
     }
     adonis.pairwise=do.call(rbind,adonis.pairwise)
     write.csv(adonis.pairwise,paste(hndl,'/Permanova table_pairwise_',dat.name,'.csv',sep=""),row.names = F)
     
   }
     
  }
  
    #2.4 Simper analysis to identify species that discriminate among groups
  if(do.simper)
  {
    SIMPER <- summary(simper(Community, d%>%pull(Show.term),parallel=7))
    
    #1. display species accounting for group differences
    Get=as.data.frame(str_split(names(SIMPER), "_", simplify = TRUE))%>%
      mutate(V1.method=sub("\\ .*", "", V1),
             V2.method=sub("\\ .*", "", V2),
             V1.zone=sub(".* ", "", V1),
             V2.zone=sub(".* ", "", V2),
             id=1:n())%>%
      filter(!V1.method==V2.method & V1.zone==V2.zone)
    SIMPER=SIMPER[Get$id]
    
    disp.simp=vector('list',length(SIMPER))
    for(n in 1:length(disp.simp))
    {
      disp.simp[[n]]=cbind(d[Show.term],
                           Community[row.names(SIMPER[[n]]%>%filter(cumsum<=simper.cumsum))])%>%
        filter(!!sym(Show.term)%in%str_split(names(SIMPER)[n], "_", simplify = TRUE))%>%
        gather(species,prop,-method.zone,)%>%
        mutate(groups=names(SIMPER)[n])
    }
    disp.simp=do.call(rbind,disp.simp)
    
    dis.cls=unique(disp.simp$species)
    dis.cls=All.species.names%>%
      filter(COMMON_NAME%in%unique(disp.simp$species))
    
    colfunc <- colorRampPalette(Shark.palette)
    n.col.elasmos=colfunc(length(dis.cls$CAES_Code[dis.cls$CAES_Code<50000]))
    names(n.col.elasmos)=dis.cls%>%filter(CAES_Code<50000)%>%pull(COMMON_NAME)
    
    colfunc <- colorRampPalette(Teleost.palette)
    n.col.teleos=colfunc(length(dis.cls$CAES_Code[dis.cls$CAES_Code>=50000]))
    names(n.col.teleos)=dis.cls%>%filter(CAES_Code>=50000)%>%pull(COMMON_NAME)
    
    p=disp.simp%>%
      group_by(method.zone,groups,species)%>%
      summarise(prop=mean(prop))%>%
      mutate(method=sub("\\ .*", "", method.zone),
             zone=sub(".* ", "", method.zone))%>%
      ggplot(aes(x=method,y=prop, fill=species))+
      geom_bar(stat="identity", width = 0.5)+
      facet_wrap(~zone,scales='free_y')+
      ylab("Average proportion")+xlab("Method")+
      scale_fill_manual(values=c(n.col.elasmos,n.col.teleos))+
      theme_PA(str.siz=14,leg.siz=LGSIZE,axs.t.siz=12,axs.T.siz=16)+
      theme(legend.position = "top",
            legend.title = element_blank(),
            plot.margin=unit(c(.1,.5,.1,.1),"cm"))+
      guides(fill = guide_legend(nrow = 4))
    print(p)
    ggsave(paste(hndl,'/Simper_',dat.name,'.tiff',sep=""),width = 6,height = 6,compression = "lzw")
    
    
    SIMPER.out=SIMPER
    for(s in 1:length(SIMPER.out))
    {
      x=SIMPER[[s]]%>%
        mutate(group=names(SIMPER)[s],
               species=row.names(SIMPER[[s]]))%>%
        relocate(group,species, .before=average)
      x=x%>%
        filter(cumsum<=simper.cumsum)
      SIMPER.out[[s]]=x 
    }
    write.csv(do.call(rbind,SIMPER.out),paste(hndl,'/Simper table_',dat.name,'.csv',sep=""),row.names = F) 
  }
  
}

predict.manyany=function (object, newdata = NULL, type = c("link", "response","terms"), se.fit = FALSE,
                          dispersion = NULL, terms = NULL, na.action = na.pass, ...) 
{
  object$family=object$family[[1]]
  if (any(object$family == "gaussian")) {
    if (type == "link") {
      stop("Possible type of predict.manylm is 'response' or 'term'.")
    }
    else {
      return(predict.manylm1(object, newdata = newdata, 
                             se.fit = se.fit, type = type, terms = terms, 
                             na.action = na.pass, ...))
    }
  }
  nVar <- NCOL(object$fitted.values)
  nObs <- NROW(object$fitted.values)
  if (is.null(newdata) == F) 
    nObs = NROW(newdata)
  ses <- fts <- matrix(NA, nObs, nVar)
  if (is.null(newdata)) 
    dimnames(fts)[[1]] = dimnames(object$fitted.values)[[1]]
  else dimnames(fts)[[1]] = rownames(newdata)
  dimnames(fts)[[2]] = dimnames(object$fitted.values)[[2]]
  type <- match.arg(type)
  na.act <- object$na.action
  object$na.action <- NULL
  fm <- formula(object)
  if ("K" %in% names(object$call)) {
    K = eval(object$call$K)
    if (K > 1 & object$family != "binomial") {
      warning("Argument K should only be specified when family is binomial, K reset to 1")
      K = 1
    }
  }
  else K = 1
  for (iVar in 1:nVar) {
    fam = switch(object$family, `binomial(link=logit)` = binomial(), 
                 `binomial(link=cloglog)` = binomial("cloglog"), 
                 poisson = poisson(), gaussian = gaussian(), gamma = Gamma(link = "log"), 
                 negative.binomial = negative.binomial(theta = object$theta[iVar]))
    if (is.null(object$data)) 
      dat.i = model.frame(object)
    else dat.i = data.frame(object$y[, iVar], object$data)
    if (K > 1) {
      dat.i$K = K
      form <- as.formula(paste("object$y[ ,", iVar, 
                               "]/K ~ ", fm[3]))
      object.i = glm(form, family = fam, data = dat.i, 
                     weights = K, start = as.vector(object$coef[, 
                                                                iVar]))
      object.i$coefficients = coef(object)[, iVar]
    }
    else {
      form <- as.formula(paste("object$y[ ,", iVar, 
                               "] ~ ", fm[3]))
      object.i = glm(form, family = fam, data = dat.i, 
                     start = as.vector(object$coef[, iVar]))
      object.i$coefficients = coef(object)[, iVar]
    }
    ft.i <- predict.glm(object.i, newdata = newdata, se.fit = se.fit, 
                        type = type, terms = terms, na.action = na.action)
    if (se.fit == T) {
      fts[, iVar] = ft.i$fit
      ses[, iVar] = ft.i$se
    }
    else fts[, iVar] = ft.i
  }
  if (se.fit) 
    out = list(fit = fts, se.fit = ses)
  else out = fts
  return(out)
}
predict.manyglmTweedie=function (object, newdata = NULL, type = c("link", "response", "terms"),
                                 se.fit = FALSE, dispersion = NULL, terms = NULL, na.action = na.pass, ...) 
{
  nVar <- NCOL(object$fitted.values)
  nObs <- NROW(object$fitted.values)
  if (is.null(newdata) == F) nObs = NROW(newdata)
  ses <- fts <- matrix(NA, nObs, nVar)
  if (is.null(newdata)) 
  {
    dimnames(fts)[[1]] = dimnames(object$fitted.values)[[1]]
  }else
  {
    dimnames(fts)[[1]] = rownames(newdata)
  }
  dimnames(fts)[[2]] = dimnames(object$fitted.values)[[2]]
  type <- match.arg(type)
  na.act <- object$na.action
  object$na.action <- NULL
  fm <- formula(object)
  if ("K" %in% names(object$call))
  {
    K = eval(object$call$K)
    if (K > 1 & object$family != "binomial") {
      warning("Argument K should only be specified when family is binomial, K reset to 1")
      K = 1
    }
  }else
  {
    K = 1
  }
   
  for (iVar in 1:nVar)
  {
    # fam = switch(object$family[[iVar]]$family,
    #               `binomial(link=logit)` = binomial(), 
    #               `binomial(link=cloglog)` = binomial("cloglog"), 
    #               poisson = poisson(),
    #               gaussian = gaussian(),
    #               gamma = Gamma(link = "log"),
    #               'Tweedie' = tweedie(var.power=1.2, link.power=0),
    #               negative.binomial = negative.binomial(theta = object$theta[iVar]))
    if (is.null(object$data)) 
    {
      dat.i = model.frame(object)%>%
        data.frame()%>%dplyr::select(all_of(Terms))
    }else
    {
      dat.i = data.frame(object$y[, iVar], object$data)
    }
    
    if(K > 1)
    {
      dat.i$K = K
      form <- as.formula(paste("object$fitted.values[ ,", iVar, 
                               "]/K ~ ", fm[3]))
      object.i = glm(form, family = fam, data = dat.i, 
                     weights = K, start = as.vector(object$coef[,iVar]))
      object.i$coefficients = coef(object)[iVar]
    }else
    {
      dd=dat.i%>%mutate(Response=object$fitted.values[, iVar])
      object.i = glm(paste("Response ~ ", fm[3]), family = tweedie(var.power=1.2, link.power=0), data = dd,
                     start = as.vector(object$coef[[iVar]]))
                     #start = as.vector(object$coef[,iVar]))
      object.i$coefficients = coef(object)[[iVar]]
    }
    
    ft.i <- predict.glm(object.i, newdata = newdata, se.fit = se.fit, type = type, terms = terms)
    if (se.fit == T)
    {
      fts[, iVar] = ft.i$fit
      ses[, iVar] = ft.i$se
    }
    else fts[, iVar] = ft.i
  }
  
  if (se.fit) 
    out = list(fit = fts, se.fit = ses)
  else out = fts
  return(out)
}
predict.manylm1=function (object, newdata = NULL, se.fit = FALSE, type = c("response","terms"),
                          terms = NULL, na.action = na.pass, ...) 
{
  if(!"y"%in%names(object))
  {
    object$y=object$model[,grepl(all.vars(formula(object))[1],names(object$model))]
  }
  
  
  nVar <- NCOL(object$fitted.values)
  nObs <- NROW(object$fitted.values)
  if (is.null(newdata) == F) 
    nObs = NROW(newdata)
  ses <- fts <- matrix(NA, nObs, nVar)
  dimnames(fts) = vector(2, mode = "list")
  dimnames(fts)[[2]] = dimnames(object$fitted.values)[[2]]
  if (is.null(newdata)) 
    dimnames(fts)[[1]] = dimnames(object$fitted.values)[[1]]
  else dimnames(fts)[[1]] = rownames(newdata)
  type <- match.arg(type)
  na.act <- object$na.action
  object$na.action <- NULL
  fm <- formula(object)
  for (iVar in 1:nVar)
  {
    form <- as.formula(paste("object$y[ ,", iVar, "] ~ ",fm[3]))
    if (is.null(object$data)) 
      dat.i = model.frame(object)
    else dat.i = data.frame(object$y[, iVar], object$data)
    object.i = gam(form, data = dat.i)
    ft.i <- predict(object.i, newdata = newdata, se.fit = se.fit, 
                    type = type, terms = terms, na.action = na.action)
    if (se.fit == T) {
      fts[, iVar] = ft.i$fit
      ses[, iVar] = ft.i$se
    }
    else fts[, iVar] = ft.i
  }
  if (se.fit) 
    out = list(fit = fts, se.fit = ses)
  else out = fts
  return(out)
}

#Catch composition thru time
Catch.comp=function(ddd,All.sp,Display)
{
  Main=unique(ddd$Dataset)
  dumi=ddd%>%
    group_by(SCIENTIFIC_NAME)%>%
    summarise(Prop=sum(Prop))%>%
    ungroup()%>%
    arrange(-Prop)%>%
    mutate(Cumktch=cumsum(Prop),
           Tot=sum(Prop),
           Cumktch2=Cumktch/Tot,
           SCIENTIFIC_NAME2=ifelse(Cumktch2<=0.95,SCIENTIFIC_NAME,'Other'))%>%
    left_join(All.sp,by='SCIENTIFIC_NAME')
  N.other=dumi%>%group_by(SCIENTIFIC_NAME2,Group)%>%tally()%>%filter(SCIENTIFIC_NAME2=='Other')
  dumi=dumi%>%
    mutate(SCIENTIFIC_NAME2=case_when(SCIENTIFIC_NAME2=='Other' & Group=='Elasmobranch'~paste0('Other elasmobranchs',' (',N.other%>%filter(Group=='Elasmobranch')%>%pull(n),' species)'),
                                      SCIENTIFIC_NAME2=='Other' & Group=='Teleost'~paste0('Other teleosts',' (',N.other%>%filter(Group=='Teleost')%>%pull(n),' species)'),
                                      TRUE~SCIENTIFIC_NAME2))%>%
    dplyr::select(SCIENTIFIC_NAME,SCIENTIFIC_NAME2,Group)
  
  LVLs=dumi%>%arrange(Group,SCIENTIFIC_NAME2)%>%distinct(Group,SCIENTIFIC_NAME2)%>%pull(SCIENTIFIC_NAME2)
  ddd=ddd%>%
    left_join(dumi,by='SCIENTIFIC_NAME')%>%
    mutate(color=ifelse(Group=='Teleost','dodgerblue4','firebrick3'))%>%
    group_by(SCIENTIFIC_NAME2,YEAR,color,Group)%>%
    summarise(Prop=sum(Prop,na.rm=T))%>%
    ungroup()%>%
    mutate(SCIENTIFIC_NAME2=factor(SCIENTIFIC_NAME2,levels=LVLs))
  if(Display=='tile')
  {
    a=ddd%>%distinct(SCIENTIFIC_NAME2,color)%>%arrange(factor(SCIENTIFIC_NAME2,levels=LVLs))
    p=ddd%>%
      rename(Proportion=Prop)%>%
      ggplot(aes(YEAR, SCIENTIFIC_NAME2 , fill= Proportion)) + 
      geom_tile()+
      scale_fill_gradient2(low="lightgoldenrodyellow",mid="darkgoldenrod2", high="brown4",midpoint = mean(range(ddd$Prop)))+
      ylab('')+xlab('Financial year')+
      theme_PA(axs.t.siz=10,Ttl.siz=13)+
      theme(legend.position = 'top',
            plot.title.position = "plot",
            axis.text.y = element_text(colour = a%>%pull(color)))
  }
  
  if(Display=='barplot')
  {
    p=ddd%>%
      ggplot(aes(YEAR,Prop,fill=SCIENTIFIC_NAME2))+
      geom_bar(stat='identity',position="stack")+
      theme_PA()+
      theme(legend.position = 'top',
            legend.title = element_blank())
  }
  
  p=p+ggtitle(Main)
  return(p)
}

# Cluster analysis for metiers
Cluster.fn=function(d,Terms,n.recent.years,percent.ktch.explained,var)
{
  #select years
  all.yrs=sort(unique(d$YEAR))
  an.yrs=all.yrs[(length(all.yrs)-4):length(all.yrs)]
  
  #set colnames to tolower
  d=d%>%
    filter(YEAR%in%an.yrs)%>%
    rename_with(tolower)
  
  #find top species
  top.sp=d%>%
    group_by(name)%>%
    summarise(Sum=sum(var))%>%
    arrange(-Sum)%>%
    ungroup()%>%
    mutate(CumSum=cumsum(Sum),
           Percent=CumSum/sum(Sum))%>%
    filter(Percent<=percent.ktch.explained)%>%
    pull(name)
  d1=d%>%
    data.frame%>%
    mutate(name=ifelse(!name%in%top.sp,'other',name))%>%
    group_by_at(Terms)%>%
    summarise(var=sum(var,na.rm=T))%>%
    ungroup()%>%
    pivot_wider(names_from = name, values_from = var,values_fill=0)
  
  #wide table 
  dt <- d1%>% dplyr::select(-subset(Terms,!Terms=='name'))
  if(var=='proportion') dt <- dt%>%mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) 
  dt=dt%>% as.matrix() 
  
  #find number of clusters
  b=fviz_nbclust(dt, clara, method = "silhouette",print.summary=T)
  #print(b+theme_classic())
  num.clus=as.numeric(as.character(b$data$clusters[match(max(b$data$y),b$data$y)]))
  
  
  #apply cluster analysis
  clara.res <- clara(dt, num.clus, samples = 50, pamLike = TRUE)
  dd <- cbind(d1, cluster = clara.res$cluster)

  #visualize clusters
  p=fviz_cluster(clara.res, 
               palette = rainbow(num.clus), # color palette
               ellipse.type = "t", # Concentration ellipse
               geom = "point", pointsize = 1,
               ggtheme = theme_classic())
  #print(p)
  
  return(list(num.clus=b+theme_classic(),cluster=p,d=dd,sp=c(top.sp,'other'),
              dt=cbind(d1%>% dplyr::select(subset(Terms,!Terms=='name')),dt,cluster = clara.res$cluster)))
}


####################################################################
#Old stuff

# {
#   #3.1.3 Bycatch patterns     
#   #note: for consistency, use same records as for indices
#   Bycatch=subset(DaTA,FATE%in%c("C","D") &SHEET_NO%in%unique(DATA.shots.diversity$SHEET_NO) ) 
#   Bycatch$Discard=with(Bycatch,ifelse(FATE=="D",1,0))
#   
#   
#   #Is there a change in what prortion is discarded by species?  
#   Disc=aggregate(Discard~YEAR+SPECIES,Bycatch,sum)
#   Tot= aggregate(INDIVIDUALS~YEAR+SPECIES,Bycatch,sum)
#   Disc=merge(Disc,Tot,by=c("YEAR","SPECIES"),all=T)
#   Disc$Prop=Disc$Discard/Disc$INDIVIDUALS
#   dummy=as.numeric(sort(unique(Disc$YEAR)))
#   all.yrs=dummy[1]:dummy[length(dummy)]
#   msn.yrs=all.yrs[which(!all.yrs%in%dummy)]
#   dummy=Disc[1:length(msn.yrs),]
#   dummy$YEAR=msn.yrs
#   dummy$Prop=NA
#   Disc=rbind(Disc,dummy)
#   Disc=Disc[order(Disc$YEAR),]
#   
#   wide <- reshape(Disc[,-match(c("Discard","INDIVIDUALS"),names(Disc))], v.names = "Prop", idvar = "SPECIES",
#                   timevar = "YEAR", direction = "wide")
#   colnames(wide)=gsub(paste("Prop",".",sep=""), "", names(wide))
#   wide[is.na(wide)]="transparent"
#   wide[wide=="1"]="black"
#   wide[wide=="0"]="grey60"
#   xx=as.numeric(colnames(wide)[-1])
#   wide=wide[order(wide$SPECIES),]
#   
#   fn.fig("Discard_species_by_yr",1200,2400) 
#   par(las=1,mar=c(1,2,0,0),oma=c(1.5,1,0,0),mgp=c(1,.6,0),cex.lab=1.25)
#   plot(xx,rep(0,length(xx)),pch=19,col=unlist(wide[1,2:ncol(wide)]),cex=.85,ylab="",xlab="",yaxt='n',ylim=c(0,nrow(wide)+1))
#   for(n in 2:nrow(wide)) points(xx,rep(n-1,length(xx)),pch=19,cex=.85,col=unlist(wide[n,2:ncol(wide)]))
#   axis(2,0:(nrow(wide)-1),wide$SPECIES,cex.axis=.4)
#   mtext("Species",2,1.8,cex=1.5,las=3)
#   mtext("Year",1,1.5,cex=1.5)
#   legend('top',c("Discarded","Retained"),bty='n',pch=19,col=c("black","grey60"),horiz=T,cex=1.25)
#   dev.off()
#   
#   
#   #set factors
#   for(x in 1:length(FactoRS)) Bycatch[,match(FactoRS[x],names(Bycatch))]=as.factor(Bycatch[,match(FactoRS[x],names(Bycatch))])
#   
#   #Binominal GLM
#   bycatch.model=glm(Discard ~ YEAR + BLOCK + BOAT + MONTH + BOTDEPTH,data=Bycatch,family=binomial)
#   
#   #Binominal GLM  
#   Bycatch.glm=fn.reshp(d=Bycatch,Y=ResVar,TimeVar="FATE",IdVAR=idvarS) 
#   Bycatch.glm$Tot=Bycatch.glm$D+Bycatch.glm$C
#   Bycatch.glm$Prop.disc=Bycatch.glm$D/(Bycatch.glm$Tot)
#   Bycatch.glm[is.na(Bycatch.glm)]=0
#   
#   
#   #Anova tables
#   Anova.tab=anova(bycatch.model, test = "Chisq")
#   n=2:length(Anova.tab$Deviance)
#   Term.dev.exp=100*(Anova.tab$Deviance[n]/bycatch.model$null.deviance)
#   names(Term.dev.exp)=rownames(Anova.tab)[n]
#   Anov.tab=as.data.frame.matrix(Anova.tab)
#   Term.tab=data.frame(Percent.dev.exp=Term.dev.exp)
#   Anova.tab=Anova.tab[-1,match(c("Deviance","Pr(>Chi)"),names(Anova.tab))]
#   Anova.tab=cbind(Anova.tab,Term.tab)
#   Anova.tab=Anova.tab[,-match("Deviance",names(Anova.tab))]
#   Anova.tab$"Pr(>Chi)"=ifelse(Anova.tab$"Pr(>Chi)"<0.001,"<0.001",round(Anova.tab$"Pr(>Chi)",3))
#   Total=Anova.tab[1,]
#   Total$"Pr(>Chi)"=""
#   Total$Percent.dev.exp=sum(Anova.tab$Percent.dev.exp)
#   rownames(Total)="Total"
#   Anova.tab=rbind(Anova.tab,Total)
#   Anova.tab$Percent.dev.exp=round(Anova.tab$Percent.dev.exp,2)
#   Anova.tab=cbind(rownames(Anova.tab),Anova.tab)
#   
#   #export anova tables as word doc
#   Export.tbl(WD=getwd(),Tbl=Anova.tab,Doc.nm=paste("Anova.table_Discard_observer.table",dat.nm,sep="_"),caption=NA,paragph=NA,
#              HdR.col='black',HdR.bg='white',Hdr.fnt.sze=10,Hdr.bld='normal',body.fnt.sze=10,
#              Zebra='NO',Zebra.col='grey60',Grid.col='black',
#              Fnt.hdr= "Times New Roman",Fnt.body= "Times New Roman",
#              HDR.names=c('TERM', "P","% dev. expl."),
#              HDR.span=c(1,1,1),
#              HDR.2nd=c("","",""))
#   
#   #show annual predictions of proportion of discards
#   ALL.yrs.bycatch=as.numeric(levels(Bycatch.glm$YEAR))
#   
#   fn.fig(paste("Discard_observer_predicted_prop",dat.nm,sep="_"),2000,2400)  
#   par(mfrow=c(2,1),mai=c(.2,.1,.1,1.1),oma=c(2,4,.1,.1),xpd=TRUE,mgp=c(2.5,.75,0),las=1)
#   
#   #Species composition
#   Bycatch=subset(Bycatch,SHEET_NO%in%unique(Bycatch.glm$SHEET_NO))
#   YrS=sort(unique(as.numeric(as.character(Bycatch$YEAR))))
#   YrS=seq(as.numeric(YrS[1]),as.numeric(YrS[length(YrS)]))
#   d=Bycatch
#   d$SPECIES=as.character(d$SPECIES)
#   d$SP.fate=with(d,ifelse(SPECIES%in%
#                             c("BB.T","DM.T","ER","PJ"),SPECIES,
#                           ifelse(FATE=="C","Comm.","Other")))
#   
#   Tot=aggregate(INDIVIDUALS~SP.fate+YEAR,d,sum)
#   Tot.reshaped=reshape(Tot,v.names = "INDIVIDUALS", idvar =c("YEAR"),
#                        timevar = "SP.fate", direction = "wide")
#   ID=2:ncol(Tot.reshaped)
#   colnames(Tot.reshaped)=gsub(paste("INDIVIDUALS",".",sep=""), "", names(Tot.reshaped))
#   Tot.reshaped$YEAR=as.numeric(as.character(Tot.reshaped$YEAR))
#   Tot.reshaped[,ID]=Tot.reshaped[,ID]/rowSums(Tot.reshaped[,ID],na.rm=T)
#   missn=YrS[which(!YrS%in%Tot.reshaped$YEAR)]
#   ADD=matrix(NA,nrow=length(missn),ncol=ncol(Tot.reshaped))
#   colnames(ADD)=colnames(Tot.reshaped)
#   ADD=as.data.frame(ADD)
#   ADD$YEAR=missn
#   Tot.reshaped=rbind(Tot.reshaped,ADD)
#   Tot.reshaped=Tot.reshaped[order(Tot.reshaped$YEAR),]
#   DD=as.matrix(Tot.reshaped[,2:ncol(Tot.reshaped)])
#   like.this=c("Comm.","Other","PJ","ER","BB.T","DM.T")
#   DD=DD[,match(like.this,colnames(DD))]
#   DD=t(DD)
#   Xmax=ncol(DD)+5
#   a=barplot(DD,plot=F)
#   CLSs=c("black","grey90","grey30","grey80","grey50","grey70")
#   #CLSs=grey.colors(nrow(DD))
#   a=barplot(DD,xaxt='n',cex.axis=1.5,col=CLSs,xlim=c(0,a[length(a)]+1))
#   LEN.nms=row.names(DD)
#   LEN.nms=ifelse(LEN.nms=="BB.T","BB",ifelse(LEN.nms=="DM.T","DM",LEN.nms))
#   legend(a[length(a)]+1,1,rev(LEN.nms), cex=1.15,
#          fill=rev(CLSs),bty='n',title="A)")
#   axis(1,a,F)
#   box()
#   
#   #Predicted discard
#   List=list(model=bycatch.model,data=Bycatch.glm)
#   fun.plt.yr.pred.bycatch(d=List,normalised="NO",PredictorS=Predictors,MAIN="",
#                           log.var="NO",Cx=1.125,YLIM=c(0,1),Cx.axs=1.35)
#   legend(YrS[length(YrS)]+1,1,"", bty='n',title="B)")
#   axis(1,a,YrS,cex.axis=1.5)
#   
#   mtext("Proportion",2,outer=T,cex=1.75,las=3,line=2.2)
#   mtext("Year",1,outer=T,line=1.1,cex=1.75)
#   dev.off()
#   
#   #Species table
#   d$SPEC.tabl=with(d,ifelse(FATE=="C","Comm.",SPECIES))
#   d$SPEC.tabl=with(d,ifelse(SPEC.tabl=="OS","BS",SPEC.tabl))
#   TLB.sp=table(d$SPEC.tabl)
#   TLB.sp=TLB.sp/sum(TLB.sp)
#   TLB.sp=data.frame(SPECIES=names(TLB.sp),PROP=c(TLB.sp))
#   #a=subset(SPECIES.names,Species%in%unique(TLB.sp$SPECIES))
#   a=subset(SPECIES_PCS_FATE,SPECIES%in%unique(TLB.sp$SPECIES))
#   TLB.sp=merge(TLB.sp,a,by="SPECIES",all.x=T)
#   TLB.sp=TLB.sp[order(-TLB.sp$PROP),]
#   TLB.sp$PROP=round(TLB.sp$PROP,3)
#   TLB.sp$PROP=with(TLB.sp,ifelse(PROP<0.001,"<0.001",PROP))
#   TLB.sp=TLB.sp[,match(c("SPECIES","PROP","COMMON_NAME","SCIENTIFIC_NAME"),names(TLB.sp))]
#   write.csv(TLB.sp,paste("Table_species_disc_com",dat.nm,"csv",sep="."),row.names=F)
#   
# }


# 
# #Function for calculating Diversity indices
# Div.ind = function(data)
# {
#   var_names = row.names(data)  
#   N = length(var_names)
#   sp_names = dimnames(data)[[2]]
#   
#   #Shannon-Wiener diversity index
#   Shannon = unname(diversity(data, index = "shannon"))
#   Shannon[Shannon == 0] = NA
#   
#   #Simpson's index
#   Simpson = unname(diversity(data, index = "simpson"))
#   Simpson[Simpson == 1] = NA
#   
#   #Margalef species richness
#   Margalef = rep(NA, N) 
#   for(i in 1:N) Margalef[i] = margalef(sp_names, data[i, ])
#   Margalef[Margalef == 0] = NA
#   
#   #Pielou's evenness                        
#   Evenness = function(H, S) {J = H / log(S)}  #where H = Shannon index for each community; S = Number of species
#   Pielou = rep(NA, N)
#   for(i in 1:N) Pielou[i] = Evenness(H = Shannon[[i]], S = length(sp_names))
#   Pielou[Pielou == 0] = NA
#   
#   #Summary
#   return(list(Shannon = round(Shannon, 2), Margalef = round(Margalef, 2), Pielou = round(Pielou, 2), Simpson = round(Simpson, 2), Var_names = var_names))
# }
# 
# 
# #Function for calculating Ecosystem indicators
# Eco.ind = function(data, TROPHIC.LEVEL, MAX.BODY.LENGTH, TROPHIC.EFFICIENCY = 0.1)
# {
#   var_names = row.names(data)  
#   N = length(var_names)
#   sp_names = dimnames(data)[[2]]
#   data[data == 0] = NA
#   
#   TROPHIC.LEVEL=subset(TROPHIC.LEVEL,SPECIES%in%colnames(data))
#   TROPHIC.LEVEL=TROPHIC.LEVEL[order(TROPHIC.LEVEL$SPECIES),]
#   
#   #Mean Trophic Level
#   Trophic_level = rep(NA,N)
#   for(i in 1:N)
#   {
#     a = data[i, ] * TROPHIC.LEVEL[, 2]
#     Trophic_level[i] =  sum(a, na.rm = T) / sum(data[i, ], na.rm = T)
#   }
#   
#   #Mean Maximum Length
#   Max_length = rep(NA,N)
#   a = data[,match(colnames(MAX.BODY.LENGTH),colnames(data))] * MAX.BODY.LENGTH
#   for(i in 1:N)
#   {
#     Max_length[i] = sum(a[i,], na.rm = T) / sum(data[i, ], na.rm = T)
#   }
#   
#   #Fishery in Balance
#   FIB = rep(NA,N)
#   for(i in 1:N) if(sum(data[i,], na.rm = T) > 0) {break}
#   
#   for(j in 1:N)
#   {
#     FIB[j] = log(sum(data[j,], na.rm = T) * (1 / TROPHIC.EFFICIENCY) ^ Trophic_level[j]) - log(sum(data[i,], na.rm = T) * (1 / TROPHIC.EFFICIENCY) ^ Trophic_level[i])
#   }
#   
#   #Summary
#   return(list(MTL = round(Trophic_level, 2), MML = round(Max_length, 2), FIB = round(FIB, 2)))
# }
# 
# 
