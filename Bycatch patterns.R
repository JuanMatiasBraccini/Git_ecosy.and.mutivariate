
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

#missing: Deal with imbalance data set   VIP!!
#         Models (ACA). 
#             Could combine indicators into one (see notes in Manuscript)
#             autocorrelation in time-space?
#             is Year the covariate or the fishery effort??
#             usesampling effort as offset()
#         For Multivariate, implement Parks_2019 and look in C:\Matias\Workshops\2019_Primer workshop\2019_PRIMER_workshop.docx

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
library(Hmisc)

#Define working directory
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

setwd(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch"))
User="Matias"


# 1 Data section-----------------------------------------------------------------------

# 1. Bring in WA shark observer data
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R"))
rm(DATA)
DATA=DATA.ecosystems%>%
      filter(!COMMON_NAME%in%c('Unidentified','','Other sharks'))


# 2. Bring in WA Species names-PCS-FATE
setwd(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch"))
SPECIES_PCS_FATE=read.csv("SPECIES+PCS+FATE.csv",stringsAsFactors=F)
SPECIES_PCS_FATE_Com=read.csv("SPECIES+PCS+FATE_Com.csv",stringsAsFactors=F)

# 3. Bring in Length coefficients
Len.cof=read.csv("Raw coefficents table.csv")


#4. Bring in Commercial landings data
source(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Git_ecosy.and.mutivariate/Commercial_data_for_Ecosystem_Analysis.R"))

#Source functions
source(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Git_ecosy.and.mutivariate/Ecosystem_functions.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_Population.dynamics/fn.fig.R"))
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Smart_par.R"))
source(handl_OneDrive('Analyses/SOURCE_SCRIPTS/Git_other/ggplot.themes.R'))


# 2 CONTROL section-----------------------------------------------------------------------

  #choose if doing .jpeg or .tiff figures
Do.jpeg="NO"
Do.tiff="YES"

  #Selection of records
Min.shts=5 #USe records with at least 5 shots per year-block
Min.shts.sens=c(Min.shts*2)
Min.recs=10        #minimum number of records per boat to be selected
Min.individuals=5   #minimum number of individuals per shot to use

  #vessel used as mixed effect
MixedEff="BOAT"   

  #choose if using commercial data
do.commercial=TRUE
combine.daily.monthly=TRUE
get.traits=FALSE  #get life history traits from FishLife
Start.yr=1975   #before 1989 some species not reported in commercial logbooks. Interpret results within this caveat

  #choose if doing data exploration
do.exploratory="NO"  

niter=100    #number of iterations for MC procedure for confidence intervals

#Percent.show.bycatch=0.8   #proportion of bycatch explained
#Prop.TC.bycatch=0.05       #proportion of catch explained


  #choose the variables used in models 
ResVar="INDIVIDUALS"
MultiVar="SPECIES"
IDVAR=c("SHEET_NO","YEAR","MONTH","BOTDEPTH","BLOCK","ZONE","BOAT","SKIPPER","MESH_SIZE",
        "BIOREGION","SEASON","EFFORT","LATITUDE","LONGITUDE")   
Predictors=c("YEAR","BLOCK","BOAT","MONTH","BOTDEPTH","LATITUDE","LONGITUDE")
Expl.varS=c("YEAR","BOAT","MONTH","BLOCK","BOTDEPTH","LATITUDE","LONGITUDE")
FactoRS=Expl.varS[-match(c("BOTDEPTH","LATITUDE","LONGITUDE"),Expl.varS)]
OFFSETT=NA

  #choose indicators 
Ecol.Indicators=c("Shannon","Pielou","Simpson","MTL","MML","MeanML","Prop.Disc")
# FIB not applicable to Effort managed fishery (if TL is maintained by catches reduced due to Management
#                                               then FIB is <0)
traits=c("MaxAge","Age.mat.prop","K","MaxLen")



# 3 Procedure section-----------------------------------------------------------------------

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
DATA=DATA%>%
  left_join(SPECIES_PCS_FATE%>%
            filter(SPECIES%in%unique(DATA$SPECIES))%>%
              dplyr::select(-c(COMMON_NAME,SCIENTIFIC_NAME)),
            by="SPECIES")

  #Special treatment for White shark (WP) and Grey nurse shark (GN) as they became protected throughout the period
DATA$FATE=ifelse(DATA$SPECIES=="WP",ifelse(DATA$YEAR<1997,"C","D"), 
            ifelse(DATA$SPECIES=="GN",ifelse(DATA$YEAR<2001,"C","D"),as.character(DATA$FATE)))
              
DATA$NATURE=ifelse(DATA$SPECIES=="WP" | DATA$SPECIES=="GN",ifelse(DATA$FATE=="C","S","TEPS"),as.character(DATA$NATURE))           
                                               
  #remove unknown species codes
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

  #Add bioregion
DATA=DATA%>%
  mutate(BIOREGION=case_when(LONGITUDE>=115.5 & LONGITUDE<=129 & LATITUDE<=(-26)~"SC",
                             LONGITUDE<115.5 & LATITUDE<=(-27)~"WC",
                             LONGITUDE<=114.834 & LATITUDE>(-27)~"Gascoyne",
                             LONGITUDE>=114.834 & LONGITUDE<=129 & LATITUDE>=(-27)~"NC",
                             TRUE ~ NA_character_),
         BIOREGION=ifelse(LATITUDE>(-34.5) & LATITUDE<(-30) &LONGITUDE< 118,"WC",BIOREGION))

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

#YR.selected=c(1995:1998,2002:2003)
YR.selected=table(DATA$YEAR)
YR.selected=as.numeric(names(YR.selected[YR.selected>1e3]))
DATA$Yr.dummy=with(DATA,ifelse(YEAR%in%YR.selected,YEAR,NA))
DATA$PER=as.character(with(DATA,paste(Yr.dummy,SEASON)))
DATA$PERIOD=with(DATA,ifelse(substr(PER,1,2)=="NA",NA,PER))
DATA=DATA[,-match(c("Yr.dummy","PER"),colnames(DATA))]

  #Add gillnet effort
DATA$EFFORT=DATA$NET_LENGTH*DATA$SOAK_TIME


  #Remove tropical species that belong to the Northern Territory and the "other sharks and scale fish" group
DATA=subset(DATA, FATE!="A" & SPECIES!="XX")

  #Add number caught
DATA$INDIVIDUALS=1

  #Fill in TL if FL available  
Len.cof=Len.cof%>%mutate(Species=capitalize(tolower(Species)),
                         Species=ifelse(Species=="Port jackson","Port Jackson",Species))
DATA=DATA%>%
      left_join(Len.cof%>%dplyr::select(Species,Intercept,Slope),
                by=c("COMMON_NAME"="Species"))%>%
      mutate(TL=ifelse(is.na(TL) & !is.na(FL),FL*Slope+Intercept,TL),
             TL=ifelse(COMMON_NAME%in%c("Southern eagle ray","Stingrays") & is.na(TL),
                       Disc.width,TL))%>%
      dplyr::select(-c(Intercept,Slope))

  #Fill in teleost FL with a constant
DATA$FL=with(DATA,{ifelse(NATURE=="T" & is.na(FL),TL*.85,FL)}) 

  #fix missing FL
Sp.len.dat=as.character(Len.cof$Species)
ID=with(DATA,which(COMMON_NAME%in%Sp.len.dat & is.na(FL) & !is.na(TL)))
if(length(ID)>0)
{
  DATA.fix=DATA[ID,]
  DATA=DATA[-ID,]
  DATA.fix=merge(DATA.fix,Len.cof[,-14],by.x="COMMON_NAME",by.y="Species",all.x=T)
  DATA.fix$FL=with(DATA.fix,(TL-Intercept)/Slope)
  DATA.fix=DATA.fix[,match(names(DATA),names(DATA.fix))]
  DATA=rbind(DATA,DATA.fix)
}

  #remove typo
DATA=DATA%>%filter(!COMMON_NAME=='Whale')

DATA$SKIPPER=with(DATA,ifelse(is.na(SKIPPER) & YEAR%in%c(2001:2002),"C. GULLOTTI",SKIPPER))
DATA=DATA%>%
  mutate(SKIPPER=tolower(SKIPPER),
         SKIPPER=case_when(grepl('lotti',SKIPPER)~"c. gullotti",
                           SKIPPER%in%c('carlo','carlo gulloti','c.gulloti')~"c. gullotti",
                           grepl('c00ke',SKIPPER)~"j.cooke",
                           grepl('parker',SKIPPER)~"r.parker",
                           grepl(paste(c('thorton','thornton','j.t'),collapse='|'),SKIPPER)~"j. thornton",
                           TRUE~SKIPPER),
         BOAT=case_when(is.na(BOAT) & SKIPPER=='c. gullotti'~'F505',
                        is.na(BOAT) & SKIPPER=='j.cooke'~'E35',
                        is.na(BOAT) & SKIPPER=='n.soulos'~'F505',
                        is.na(BOAT) & SKIPPER=='p.osborn'~'B22',
                        TRUE~BOAT))

  #Select relevant variables 
DATA=subset(DATA,select=c(SHEET_NO,YEAR,MONTH,BLOCK,BOAT,SKIPPER,SPECIES,FATE,
                          COMMON_NAME,NATURE,TL,SEX,TROPHIC_LEVEL,TL_SD,BOTDEPTH,
                          MESH_SIZE,NET_LENGTH,SOAK_TIME,EFFORT,LATITUDE,LONGITUDE,
                          ZONE,REGION,BIOREGION,AREA,SEASON,PERIOD,INDIVIDUALS,
                          PCS,SCIENTIFIC_NAME))

#remove records from shots with less than Min.shts per year-block
d=DATA[!duplicated(DATA$SHEET_NO),]
d$N=1
d1=aggregate(N~BLOCK+YEAR,d,sum)
N.base.case=subset(d1,N>Min.shts)
N.base.case$Yr.blk=with(N.base.case,paste(YEAR,BLOCK,sep="_"))
DATA$Yr.blk=with(DATA,paste(YEAR,BLOCK,sep="_"))
DATA=subset(DATA,Yr.blk%in%N.base.case$Yr.blk)

#Remove boats with less than Min.recs
AA=sort(table(DATA$BOAT))
AA=AA[AA>Min.recs]
DATA=subset(DATA,BOAT%in%names(AA))

#Remove years with less than Min.recs
AA=sort(table(DATA$YEAR))
AA=AA[AA>Min.recs]
DATA=subset(DATA,YEAR%in%names(AA))

#Show when species start appearing in data set
YEAR.min=min(DATA$YEAR)
a=DATA%>%
  mutate(Julian=julian(as.Date(paste(YEAR,MONTH,1,sep='-')),origin=as.Date(paste(YEAR.min,1,1,sep='-'))))%>%
  distinct(MONTH,YEAR,Julian,SCIENTIFIC_NAME)%>%
  mutate(N=1)
dummi=data.frame(SCIENTIFIC_NAME=names(sort(table(a$SCIENTIFIC_NAME))))%>%mutate(row_id=row_number())
LabL=dummi$SCIENTIFIC_NAME
names(LabL)=dummi$row_id
a=a%>%
  left_join(dummi,by="SCIENTIFIC_NAME")
a%>%
  ggplot(aes(Julian,row_id))+
  geom_point()+
  scale_y_continuous(labels=LabL,n.breaks=length(LabL),breaks=as.numeric(names(LabL)))+
  ylab('')+xlab('Julian day')+
  theme(axis.text.y = element_text(size = 7))
ggsave("Outputs/Exploratory/Species appearing in observer data.tiff",width = 6,height = 10,compression = "lzw")



  #3.3. Manipulated landings data (Monthly-aggregated records to 2021)
if(do.commercial)
{
  Data.monthly=Data.monthly%>%
            filter(METHOD=="GN" & zone%in%c('West','Zone1','Zone2') & LATITUDE<=(-26))%>%
            dplyr::select(-NETLEN.c)
  Data.daily=Data.daily%>%filter(METHOD=="GN" & zone%in%c('West','Zone1','Zone2') & LATITUDE<=(-26))
  if(combine.daily.monthly)
  {
    Data.daily.agg=Data.daily%>%
                    group_by(YEAR,MONTH,BOAT,METHOD,BLOCK,SPECIES,SNAME,YEAR.c,
                             LATITUDE,LONGITUDE,TYPE.DATA,zone,SKIPPER)%>%
                    summarise(LIVEWT.c=sum(LIVEWT.c))
    Effort.daily.agg=Data.daily%>%
                    distinct(SHEET_NO,EFFORT,YEAR,MONTH,BOAT,METHOD,BLOCK)%>%
                    group_by(YEAR,MONTH,BOAT,METHOD,BLOCK)%>%
                    summarise(EFFORT=sum(EFFORT))
    Data.monthly=rbind(Data.monthly,
                       Data.daily.agg%>%
                         left_join(Effort.daily.agg,by=c('YEAR','MONTH','BOAT','METHOD','BLOCK'))%>%
                         mutate(SHEET_NO=paste(YEAR,MONTH,BOAT,METHOD,BLOCK))%>%
                         relocate(names(Data.monthly)))
  }
  
  #remove records from shots with less than Min.shts per year-block
  d=Data.monthly[!duplicated(Data.monthly$SHEET_NO),]
  d$N=1
  d1=aggregate(N~BLOCK+YEAR,d,sum)
  N.base.case=subset(d1,N>Min.shts)
  N.base.case$Yr.blk=with(N.base.case,paste(YEAR,BLOCK,sep="_"))
  Data.monthly$Yr.blk=with(Data.monthly,paste(YEAR,BLOCK,sep="_"))
  Data.monthly=subset(Data.monthly,Yr.blk%in%N.base.case$Yr.blk)
  
  #Remove boats with less than Min.recs
  AA=sort(table(Data.monthly$BOAT))
  AA=AA[AA>Min.recs]
  Data.monthly=subset(Data.monthly,BOAT%in%names(AA))
  
  #Remove years with less than Min.recs
  AA=sort(table(Data.monthly$YEAR))
  AA=AA[AA>Min.recs]
  Data.monthly=subset(Data.monthly,YEAR%in%names(AA))
  
  #Fix duplicated SPECIES
  Data.monthly=Data.monthly%>%
    mutate(SNAME=ifelse(SPECIES==22999,'shark, other',
                 ifelse(SPECIES==18003,'Dusky Whaler',
                 ifelse(SPECIES==18014,'Blacktip Shark',
                        SNAME))))
  
  #add species trophic level,etc
  Data.Mon=merge(Data.monthly,SPECIES_PCS_FATE_Com,by="SPECIES",all.x=T,all.y=F)%>%  
                    rename(ZONE=zone)%>%
    mutate(BIOREGION=case_when(LONGITUDE>=115.5 & LONGITUDE<=129 & LATITUDE<=(-26)~"SC",
                               LONGITUDE<115.5 & LATITUDE<=(-27)~"WC",
                               LONGITUDE<=114.834 & LATITUDE>(-27)~"Gascoyne",
                               LONGITUDE>=114.834 & LONGITUDE<=129 & LATITUDE>=(-27)~"NC",
                               TRUE ~ NA_character_),
           BIOREGION=ifelse(LATITUDE>(-34.5) & LATITUDE<(-30) &LONGITUDE< 118,"WC",BIOREGION))
  Data.monthly=subset(Data.Mon,select=c(SHEET_NO,YEAR,MONTH,BLOCK,BOAT,TYPE.DATA,SPECIES,
                                        SCIENTIFIC_NAME,NATURE,LIVEWT.c,TROPHIC_LEVEL,TL_SD2,EFFORT,
                                        LATITUDE,LONGITUDE,ZONE,BIOREGION,
                                        Loo,K,tmax,tm,Lm))
  Data.monthly=subset(Data.monthly,!is.na(Data.monthly$TROPHIC_LEVEL))
  Data.monthly$SPECIES=as.factor(Data.monthly$SPECIES)
  Data.monthly$INDIVIDUALS=Data.monthly$LIVEWT.c
  Data.monthly$FINYEAR=Data.monthly$YEAR
  Data.monthly$YEAR=as.numeric(substr(Data.monthly$YEAR,1,4))
  
  if(get.traits)
  {
    library(FishLife)
    normal.scale=c('Temperature','rho','h','r','G') #Mean_pred parameters in normal scale
    fn.convert.normal.space=function(d,normal.scale)
    {
      id=match(normal.scale,names(d[[1]]$Mean_pred))
      Median=c(exp(d[[1]]$Mean_pred[-id]),d[[1]]$Mean_pred[id])
      Mean=c(exp(d[[1]]$Mean_pred[-id]+0.5*diag(d[[1]]$Cov_pred)[-id]),d[[1]]$Mean_pred[id]) #biased-corrected mean
      return(list(Median=Median,Mean=Mean))
    }
    
    AA=Data.monthly%>%distinct(SPECIES,SCIENTIFIC_NAME)
    
    sp.list=data.frame(SPECIES=AA$SPECIES,
                       Genus=gsub( " .*$", "", AA$SCIENTIFIC_NAME),
                       Species=sub("^\\S+\\s+", '', AA$SCIENTIFIC_NAME))
    Store.taxa.fishlife=vector('list',nrow(sp.list))
    names(Store.taxa.fishlife)=sp.list$Common.name
    for(s in 1:nrow(sp.list))
    {
      tryCatch({
        this=Search_species(Genus=sp.list$Genus[s],
                            Species=sp.list$Species[s],add_ancestors=FALSE)$match_taxonomy
        Store.taxa.fishlife[[s]]=Plot_taxa(this,mfrow=c(3,3))
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    Store.keypars.fishlife=Store.taxa.fishlife
    keypars=c('Loo','K','tmax','tm','M','Lm','ln_Fmsy_over_M','ln_Fmsy','h','r','G')
    for(s in 1:nrow(sp.list))
    {
      if(!is.null(Store.taxa.fishlife[[s]]))
      {
        Median.and.mean=fn.convert.normal.space(d=Store.taxa.fishlife[[s]],normal.scale)
        dummy=rbind(Median.and.mean$Median,Median.and.mean$Mean)%>%
          data.frame%>%
          dplyr::select(all_of(keypars))%>%
          summarise_all(list(mean))%>%
          rename(Fmsy_over_M=ln_Fmsy_over_M,
                 Fmsy=ln_Fmsy)%>%
          mutate(Species=sp.list$SPECIES[s])%>%
          relocate(Species)%>%
          mutate(across(where(is.numeric), round, 3))
        
        
      }
      if(is.null(Store.taxa.fishlife[[s]]))
      {
        dummy=data.frame(Species=sp.list$SPECIES[s],
                         Loo=NA,K=NA,tmax=NA,tm=NA,M=NA,Lm=NA,Fmsy_over_M=NA,Fmsy=NA, h=NA,r=NA,G=NA)
      }
      Store.keypars.fishlife[[s]]=dummy
      rm(dummy)
    }
    write.csv(do.call(rbind,Store.keypars.fishlife),"dummy.csv",row.names = F)
    
  }
  
  #Show when species start appearing in data set 
  YEAR.min=min(Data.monthly$YEAR)
  a=Data.monthly%>%
    mutate(Julian=julian(as.Date(paste(YEAR,MONTH,1,sep='-')),origin=as.Date(paste(YEAR.min,1,1,sep='-'))))%>%
    distinct(MONTH,YEAR,Julian,SCIENTIFIC_NAME)%>%
    mutate(N=1)
  dummi=data.frame(SCIENTIFIC_NAME=names(sort(table(a$SCIENTIFIC_NAME))))%>%mutate(row_id=row_number())
  LabL=dummi$SCIENTIFIC_NAME
  names(LabL)=dummi$row_id
  a=a%>%
    left_join(dummi,by="SCIENTIFIC_NAME")
  
  a%>%
    ggplot(aes(Julian,row_id))+
    geom_point()+
    scale_y_continuous(labels=LabL,n.breaks=length(LabL),breaks=as.numeric(names(LabL)))+
    ylab('')+xlab('Julian day')+
    theme(axis.text.y = element_text(size = 7))
  ggsave("Outputs/Exploratory/Species appearing in logbook.tiff",width = 6,height = 10,compression = "lzw")

  rm(Data.Mon,a)
}


# 4 Ecosystems indicators analysis-----------------------------------------------------------------------
Data.list=list(Observer=DATA)
if(do.commercial) Data.list$Logbook=Data.monthly
Store.data.list=vector('list',length(Data.list))
names(Store.data.list)=names(Data.list)
for(l in 1:length(Data.list))
{
  #1. Preliminary analyses
  if(do.exploratory=='YES')
  {
    #Number of sheet number per year
    surveys=data.frame(YEAR=as.factor(""), SURVEYS=as.numeric(0))#dummy
    for(i in unique(Data.list[[l]]$YEAR))
    {
      year=subset(Data.list[[l]], YEAR==i)
      n=length(unique(year$SHEET_NO))
      surveys=rbind(surveys, data.frame(YEAR=as.factor(i), SURVEYS=n))
    }
    surveys=surveys[-1,]#deleting dummy row
    surveys=surveys[order(as.numeric(as.character(surveys$YEAR))),]#sorting by year
    
    fn.fig(paste0("Outputs/Exploratory/Number of sheets by year_",names(Data.list)[l]),2000,2000)
    a=barplot(surveys$SURVEYS,ylim=c(0,700),main="No of sheet # per year",space=0.5)
    text(a[,1],surveys$SURVEYS+20,surveys$YEAR,srt=45,cex=0.75)
    dev.off()
    
    #Map the number of 'Sheet #' per block
    surveys=data.frame(BLOCK=as.factor(""), SURVEYS=as.numeric(0))#dummy
    for(i in unique(Data.list[[l]]$BLOCK))   {
      block=subset(Data.list[[l]], BLOCK==i)
      n=length(unique(block$SHEET_NO))
      surveys=rbind(surveys, data.frame(BLOCK=as.factor(i), SURVEYS=n))
    }
    surveys=surveys[-1,]#deleting dummy row
    surveys$LAT_BLOCK=-as.numeric(substr(surveys$BLOCK,1,2))
    surveys$LONG_BLOCK=as.numeric(substr(surveys$BLOCK,3,4))+100
    ocean.pal=colorRampPalette(brewer.pal(n=9,'YlOrRd'))(max(surveys$SURVEYS))#colors
    
    fn.fig(paste0("Outputs/Exploratory/Map of number of sheets per block_",names(Data.list)[l]),2000,2000) 
    par(las=1,mar=c(1.75,3.5,1.5,.1),oma=c(1.5,1,.1,.1),mgp=c(1,.8,0),cex.lab=1.25)
    
    maps::map("worldHires",xlim=c(112.95, 129.5),ylim=c(-36, -26.5),col="grey80",fill=F)
    for(i in seq(surveys$SURVEYS))
    {
      rect(surveys$LONG_BLOCK[i],surveys$LAT_BLOCK[i]-1,surveys$LONG_BLOCK[i]+1,surveys$LAT_BLOCK[i],border="blue",col=ocean.pal[surveys$SURVEYS[i]])
    }
    maps::map("worldHires",xlim=c(112.95, 129.5),ylim=c(-36, -26.5),col="grey80",fill=T,border="grey80",add=T)
    rect(surveys$LONG_BLOCK,surveys$LAT_BLOCK-1,surveys$LONG_BLOCK+1,surveys$LAT_BLOCK,border="blue")
    mtext("No of sheet # per block",line=1,cex=1.5)
    text(surveys$LONG_BLOCK+0.5,surveys$LAT_BLOCK-0.5,surveys$SURVEYS,cex=0.8)
    text(127.5,-27,paste("n =",sum(surveys[2])))
    dev.off()
    
    Data.list[[l]]%>% 
      distinct(SHEET_NO,YEAR,BLOCK,INDIVIDUALS)%>%
      group_by(YEAR,BLOCK)%>%
      tally()%>%
      mutate(LAT=-as.numeric(substr(BLOCK,1,2)),
             LON=100+as.numeric(substr(BLOCK,3,4)))%>%
      ggplot(aes(LON,LAT,size=n,color=n))+
      geom_point()+
      facet_wrap(~YEAR)+
      xlim(112.5,129)+ylim(-35.5,-27)
    ggsave(paste0("Outputs/Exploratory/Map of number of sheets per block per year_",names(Data.list)[l],'.tiff'),
           width = 8,height = 8,compression = "lzw")
    
    #No of shots per year-block-month
    a=Data.list[[l]][!duplicated(Data.list[[l]]$SHEET_NO),]
    a$N=1
    a=aggregate(N~YEAR+MONTH+BLOCK+LATITUDE+LONGITUDE,a,sum)
    
    yrs=sort(unique(a$YEAR))
    pdf(paste0('Outputs/Exploratory/Shots per year-block-month_',names(Data.list)[l],'.pdf'))
    for(i in 1:length(yrs)) fn.see(d=subset(a,YEAR==yrs[i]))
    dev.off()
    rm(a)
    
    #Export species by year
    if(names(Data.list)[l]=="Observer")dis.var='COMMON_NAME'
    if(names(Data.list)[l]=="Logbook")dis.var='SCIENTIFIC_NAME'
    TBLA=table(paste(Data.list[[l]][,dis.var],Data.list[[l]]$TROPHIC_LEVEL),Data.list[[l]]$YEAR)
    write.csv(TBLA,paste0('Outputs/Exploratory/TBLA_species.TL_year_',names(Data.list)[l],'.csv'),row.names=T)
    TBLA=table(Data.list[[l]][,dis.var],Data.list[[l]]$YEAR)
    write.csv(TBLA,paste0('Outputs/Exploratory/TBLA_species_year_',names(Data.list)[l],'.csv'),row.names=T)
    TBLA=table(Data.list[[l]]$TROPHIC_LEVEL,Data.list[[l]]$YEAR)
    colfunc <- colorRampPalette(c("red", "yellow"))
    CL=rep("grey40",nrow(TBLA))
    xx=as.numeric(colnames(TBLA))
    mltplr=5
    fn.fig(paste0('Outputs/Exploratory/TrophicLevel_by_yr_',names(Data.list)[l]),2000,2000) 
    par(las=1,mar=c(1.75,3.5,1.5,.1),oma=c(1.5,1,.1,.1),mgp=c(1,.8,0),cex.lab=1.25)
    plot(xx,rep(1,length(xx)),cex=(TBLA[1,]/max(TBLA))*mltplr,ylab="",pch=19,col=CL[1],xlab="",yaxt='n',ylim=c(0,nrow(TBLA)))
    for(n in 2:nrow(TBLA)) points(xx,rep(n,length(xx)),cex=(TBLA[n,]/max(TBLA))*mltplr,pch=19,col=CL[n])
    axis(2,1:nrow(TBLA),row.names(TBLA))
    mtext("Trophic level",2,2.5,cex=1.5,las=3)
    mtext("Year",1,2,cex=1.5)
    dev.off()
  }
  
  
  #2. Define factors and character variables
  Data.list[[l]]=Data.list[[l]]%>%
                  mutate(YEAR=as.character(YEAR),
                         MONTH=as.character(MONTH),
                         BLOCK=as.character(BLOCK),
                         SPECIES=as.factor(SPECIES),
                         Yr.blk=paste(YEAR,BLOCK,sep="_"))

  #3. Calculate ecological and functional indicators
  if(names(Data.list)[l]=="Observer")
  {
    idvarS=IDVAR
    resp.vars=subset(Ecol.Indicators,!Ecol.Indicators%in%c("MTL","FIB"))
  }
  if(names(Data.list)[l]=="Logbook")
  {
    idvarS=c("SHEET_NO","YEAR","BLOCK","BOAT","MONTH","LATITUDE","LONGITUDE")
    resp.vars=c(Ecol.Indicators%>%subset(Ecol.Indicators%in%c("Shannon","Pielou","Simpson","MTL")),
                traits)
  }
  n.rv=length(resp.vars)
  Main.title=resp.vars
  Res.var.in.log=rep("NO",length(resp.vars))  #fit response var in log space or not? 
    
  OUT.base.case=fn.calc.ecol.ind(DaTA=Data.list[[l]],
                                 dat.nm="BC",
                                 normalised="YES",
                                 Drop.yrs="NO",
                                 idvarS=idvarS,
                                 resp.vars=resp.vars)
  if(names(Data.list)[l]=="Observer")
  {
    OUT.base.case.retained=fn.calc.ecol.ind(DaTA=Data.list[[l]]%>%filter(FATE=='C'),
                                            dat.nm="BC",
                                            normalised="YES",
                                            Drop.yrs="NO",
                                            idvarS=idvarS,
                                            resp.vars=subset(resp.vars,!resp.vars=="Prop.Disc"))
    OUT.base.case.discarded=fn.calc.ecol.ind(DaTA=Data.list[[l]]%>%filter(FATE=='D'),
                                            dat.nm="BC",
                                            normalised="YES",
                                            Drop.yrs="NO",
                                            idvarS=idvarS,
                                            resp.vars=subset(resp.vars,!resp.vars=="Prop.Disc"))
  }
  
  #Exploratory analyses
  if(do.exploratory=="YES")
  {
    #Look at annual averages  
    fn.viol.box(d=OUT.base.case,byzone=TRUE,resp.vars)
    ggsave(paste0("Outputs/Exploratory/bxplt_year.zone_",names(Data.list)[l],'.tiff'),
           width = 8,height = 6,compression = "lzw")
    fn.viol.box(d=OUT.base.case,byzone=FALSE,filcol='lightsalmon2',resp.vars)
    ggsave(paste0("Outputs/Exploratory/bxplt_year_",names(Data.list)[l],'.tiff'),
           width = 8,height = 6,compression = "lzw")
    if(names(Data.list)[l]=="Observer")
    {
      fn.viol.box(d=OUT.base.case.retained,byzone=FALSE,filcol='lightsalmon2',resp.vars=subset(resp.vars,!resp.vars=="Prop.Disc"))
      ggsave(paste0("Outputs/Exploratory/bxplt_year_",names(Data.list)[l],'_retained.tiff'),
             width = 8,height = 6,compression = "lzw")
      fn.viol.box(d=OUT.base.case.discarded,byzone=FALSE,filcol='lightsalmon2',resp.vars=subset(resp.vars,!resp.vars=="Prop.Disc"))
      ggsave(paste0("Outputs/Exploratory/bxplt_year_",names(Data.list)[l],'_discarded.tiff'),
             width = 8,height = 6,compression = "lzw")
    }
    
     #Response var error structure
    fn.fig(paste0("Outputs/Exploratory/error_structure_Obs",names(Data.list)[l]),2000,2000) 
    par(mfrow=c(n.rv,2),las=1,mar=c(1,1.5,1,1.2),oma=c(1,1,.1,1),mgp=c(1,.5,0))
    for(i in 1:n.rv)
    {
      density.response(d=OUT.base.case,Var=resp.vars[i],Log="NO")
      density.response(d=OUT.base.case,Var=resp.vars[i],Log="YES")
    }
    dev.off()
    
    #Check if correlation between effort and indices
    fn.fig(paste0("Outputs/Exploratory/Correlation_effort_indices_",names(Data.list)[l]),2000,2000) 
    smart.par(n.plots=length(resp.vars),MAR=c(1.75,3.5,1.5,.1),OMA=c(3,1,.3,2),MGP=c(1,.8,0))
    for(i in 1:n.rv)
    {
      dd=OUT.base.case[,match(c("EFFORT",resp.vars[i]),names(OUT.base.case))]
      dd=dd[!is.na(dd[,2]),]
      plot(dd$EFFORT,dd[,2])
      mtext(paste(resp.vars[i]," cor=",round(cor(dd[,1],dd[,2]),2),sep=""),3,col=2)
      rm(dd)
    }
    dev.off()

  }
  
  #ACA apply models developed in 'fn.loop.over.obsrvr.data()', etc
  
  
  #Store results
  Store.data.list[[l]]=OUT.base.case
  rm(OUT.base.case)
}



  #--- 4.1 WA Fisheries observer data---

# 4.1.1 Preliminary analyses





#boats and skippers 
# A=table(as.character(DATA$BOAT),as.character(DATA$SKIPPER),useNA = "ifany")
# A[A>0]=1






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
