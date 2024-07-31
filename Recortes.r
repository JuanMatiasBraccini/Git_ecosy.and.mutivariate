#keep records from shots according to criteria  
# d=DATA[!duplicated(DATA$SHEET_NO),]
# d$N=1
# d1=aggregate(N~BLOCK+YEAR,d,sum)
# d1=subset(d1,!BLOCK=="0")
# N.base.case=subset(d1,N>Min.shts)
# N.sens.1=subset(d1,N>Min.shts.sens[1])


# N.base.case$Yr.blk=with(N.base.case,paste(YEAR,BLOCK,sep="_"))
# N.sens.1$Yr.blk=with(N.sens.1,paste(YEAR,BLOCK,sep="_"))


# DATA$Yr.blk=with(DATA,paste(YEAR,BLOCK,sep="_"))




# DATA.base.case=subset(DATA,Yr.blk%in%unique(N.base.case$Yr.blk))
# DATA.sens.1=subset(DATA,Yr.blk%in%unique(N.sens.1$Yr.blk))





# 4.2 WA Fisheries commercial data
if(do.commercial)
{
  
  
  
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
  
  
  
  # 4.2.3 Modelling
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
  
}


# 5 glm approach-------------------------------------------------------------------------
#takes 11 seconds
if(do.commercial)
{
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
}


# 6 data mining-------------------------------------------------------------------------
if(do.commercial)
{
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
  
}



# 7 Multivariate analysis-----------------------------------------------------------------------
#note: Consider the Multivariate stats used for Parks Australia 2019!!!

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
if(do.commercial)
{
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
}
