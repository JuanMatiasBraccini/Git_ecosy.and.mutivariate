
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

library(readxl)
library(maps)
library(mapdata)
library(RColorBrewer)
#library(gamlss)
library(reshape2)
#library(iNEXT) 
library(ggplot2)
library(caret)      #for all things data mining
library(plotrix)
#library(ReporteRs)
library(mvtnorm)      #for multivariate normal pdf
library(lme4) #mixed models
library(MuMIn)  #model selection and pseudoR2 mixed effect models
library(car)    #to get ANOVA from lmer model
library(data.table)
library(tidyverse)
library(Hmisc)
library(tictoc)
library(ggpubr)
require(grid)
library(ggpmisc)
library(ggstream)
library(ggtext)

#Define working directory
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

setwd(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch"))
User="Matias"


# 1 Data section-----------------------------------------------------------------------

# 1. Bring in WA shark observer data
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R"))
rm(DATA)
DATA=DATA.ecosystems%>%
      filter(!COMMON_NAME%in%c('Unidentified','','Other sharks'))%>%
  mutate(SCIENTIFIC_NAME=str_remove(SCIENTIFIC_NAME,paste(c("Families ","Family "),collapse='|')))



# 2. Bring in WA Species names-PCS-FATE
setwd(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch"))
SPECIES_PCS_FATE=read.csv("SPECIES+PCS+FATE.csv",stringsAsFactors=F)
SPECIES_PCS_FATE_Com=read.csv("SPECIES+PCS+FATE_Com.csv",stringsAsFactors=F)
Functional.traits=read.csv("Functional.traits.csv",stringsAsFactors=F)  #from FishBase get.fishbase.traits.R

SPECIES_PCS_FATE=left_join(SPECIES_PCS_FATE,
                               Functional.traits%>%
                                 dplyr::select(-SCIENTIFIC_NAME)%>%
                                 rename(SCIENTIFIC_NAME=SCIENTIFIC_NAME.original),
                               by=c('SCIENTIFIC_NAME'))

SPECIES_PCS_FATE_Com=left_join(SPECIES_PCS_FATE_Com,
                               Functional.traits%>%
                                 dplyr::select(-SCIENTIFIC_NAME)%>%
                                 rename(SCIENTIFIC_NAME=SCIENTIFIC_NAME.original),
                               by=c('SCIENTIFIC_NAME'))

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
Min.shts=5        #Use records with at least 5 shots per year-block
Min.shts.sens=c(Min.shts*2)
Min.recs=10        #select boats with at least this number of records
Min.individuals=5   #minimum number of individuals per shot to use
Min.shots.year=30  #select years with at least this number of shots
Min.years=3        #minimum number of years with Min.shots.year
  #vessel used as mixed effect
MixedEff="BOAT" 

# define if looking at landed and discarded separate
check.discards=FALSE

  #choose if using commercial data
do.commercial=TRUE
combine.daily.monthly=TRUE
use.NSF=FALSE

get.traits=FALSE  #get life history traits from FishLife
Min.YEAR=1988     # First year of analysis. before 1988 some species not reported in commercial logbooks; also some species become protected. Interpret results within this caveat
protected.species=c(8001,10003) 
names(protected.species)=c("greynurse shark","white shark")

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
Predictors=c("YEAR","BOAT","MONTH","LATITUDE","LONGITUDE") #"BLOCK"
Expl.varS=c("YEAR","BOAT","MONTH","BLOCK","BOTDEPTH","LATITUDE","LONGITUDE")
FactoRS=Expl.varS[-match(c("BOTDEPTH","LATITUDE","LONGITUDE","MONTH"),Expl.varS)]
OFFSETT=NA

response.var='cpue'  #use catch rates to calculate indicators

  #choose indicators 
Ecol.Indicators=c("Shannon","Pielou","Simpson","MTL","MML","Prop.Disc")
# FIB not applicable to Effort managed fishery (if TL is maintained by catches reduced due to Management
#                                               then FIB is <0)
traits=c() #"Age.mat.prop","K" "MaxAge" & "MaxLen" used in Functional diversity

Functional.diversity=c("FnRich_morph","FnRich_ecol")
traits_ecol=c('Trophic.level.fishbase','Habitat','Movement.scale','Feeding.group')  
traits_morph=c('Body.shape','Maximum.size','tmax','Body.depth','Head.length','Eye.diameter','Pre.orbital.length')


do.multivariate=TRUE  #focus only on ecological indicators

display.catch.effort=FALSE

#Predicted lats and longs for displaying year effect
Pred.lat=c(-29,-34,-34)
Pred.long=c(114,115,122)
names(Pred.lat)=names(Pred.long)=c('West','Zone 1','Zone 2')

# 3 Procedure section-----------------------------------------------------------------------

#3.1. Remove protected species from commercial catch
if(!is.null(protected.species))
{
  Data.daily=Data.daily%>%filter(!SPECIES%in%protected.species)
  Data.daily.LL=Data.daily.LL%>%filter(!SPECIES%in%protected.species)
  Data.daily.north=Data.daily.north%>%filter(!SPECIES%in%protected.species)
  Data.monthly=Data.monthly%>%filter(!SPECIES%in%protected.species)
  Data.monthly.north=Data.monthly.north%>%filter(!SPECIES%in%protected.species)
}


#3.2. Manipulate trophic levels
  #add SD from FishBase to Cortes'; use FishBase for Teleosts
SPECIES_PCS_FATE=SPECIES_PCS_FATE%>%
  mutate(TL_SD=ifelse(is.na(TL_SD)& NATURE%in%c("T"),TL_SD2,
               ifelse(is.na(TL_SD)& NATURE%in%c("S","R"),TROPHIC_LEVEL*TL_SD2/TROPHIC_LEVEL2,
                      TL_SD)),
         TROPHIC_LEVEL=ifelse(is.na(TROPHIC_LEVEL),TROPHIC_LEVEL2,TROPHIC_LEVEL),
         Trophic.level=ifelse(!is.na(TROPHIC_LEVEL),TROPHIC_LEVEL,Trophic.level))%>%
  rename(Trophic.level.fishbase=Trophic.level)

SPECIES_PCS_FATE_Com=SPECIES_PCS_FATE_Com%>%
        mutate(TL_SD=ifelse(is.na(TL_SD)& NATURE%in%c("T"),TL_SD2,
                            ifelse(is.na(TL_SD)& NATURE%in%c("S","R"),TROPHIC_LEVEL*TL_SD2/TROPHIC_LEVEL2,
                                   TL_SD)),
               TROPHIC_LEVEL=ifelse(is.na(TROPHIC_LEVEL),TROPHIC_LEVEL2,TROPHIC_LEVEL),
               Trophic.level=ifelse(!is.na(TROPHIC_LEVEL),TROPHIC_LEVEL,Trophic.level))%>%
        rename(Trophic.level.fishbase=Trophic.level)


#3.3. Manipulate WA shark observer data

  #Extract year and month
DATA$DATE=as.Date(DATA$date,format="%d/%m/%Y")
DATA$YEAR=as.numeric(strftime(DATA$DATE, format="%Y"))
DATA$MONTH=as.numeric(strftime(DATA$DATE, format="%m"))

  #Fix method and mesh size issues
DATA$Method=with(DATA,ifelse(is.na(Method) & BOAT%in%c("B67","F244","F517","E35"),"GN",Method))
DATA$MESH_SIZE=with(DATA,ifelse(is.na(MESH_SIZE) & BOAT=="B67" & YEAR>1997,"7",
                         ifelse(is.na(MESH_SIZE) & BOAT=="F517" & YEAR>2002,"7",MESH_SIZE)))

  #Select relevant gear 
used.gear='GN'
if(use.NSF) used.gear=c(used.gear,'LL')
DATA=subset(DATA,Method%in%used.gear) 
DATA=DATA%>%
       mutate(keep=ifelse(Method=='LL'| (Method=='GN' & MESH_SIZE%in%c("6","6.5","7")) ,'yes','no'))%>%
  filter(keep=='yes')%>%
  dplyr::select(-keep)
  
  
Research.vess=c("HAM","HOU","RV GANNET","RV BREAKSEA","NATT","NAT","FLIN","RV SNIPE 2")
DATA=subset(DATA,!BOAT%in%Research.vess)
ALL.yrs=sort(as.numeric(unique(DATA$YEAR)))

  #Get latitude and longitude from de first end of the net for location
DATA=DATA%>%
  mutate(LATITUDE=as.numeric(paste(END1LATD,".",
                                    ifelse(trunc(END1LATM*100/60)>=10,
                                               paste(substr(END1LATM*100/60,1,2),substr(END1LATM*100/60,4,5),sep=""),
                                               paste(0,substr(END1LATM*100/60,1,1),substr(END1LATM*100/60,3,4),sep="")),sep="")),
         LATITUDE=ifelse(is.na(LATITUDE),Mid.Lat,LATITUDE),
         LATITUDE=ifelse(SHEET_NO=='R00384',Mid.Lat,LATITUDE))

DATA=DATA%>%
  mutate(LONGITUDE=as.numeric(paste(END1LNGD,".",
                                     ifelse(trunc(END1LNGM*100/60)>=10,
                                             paste(substr(END1LNGM*100/60,1,2),substr(END1LNGM*100/60,4,5),sep=""),
                                            paste(0,substr(END1LNGM*100/60,1,1),substr(END1LNGM*100/60,3,4),sep="")),sep="")),
         LONGITUDE=ifelse(is.na(LONGITUDE),Mid.Long,LONGITUDE))

  #Data range
if(!use.NSF) DATA=subset(DATA, LATITUDE<(-26) | LATITUDE==0) #zero to include the dodgy sheet numbers F00001, F00002 and F00003 with zero position

  #Identify fishery
DATA=DATA%>%mutate(Fishery=ifelse(LATITUDE>(-26),'NSF','TDGDLF'),
                   Fishery=ifelse(LATITUDE==0,'TDGDLF',Fishery),
                   Keep=ifelse((Fishery=='NSF' & Method=='LL') | (Fishery=='TDGDLF' & Method=='GN'),'Yes','No'))%>%
            filter(Keep=='Yes')%>%
            dplyr::select(-Keep)

  #Fix sex
DATA$SEX=as.character(DATA$SEX)
DATA$SEX=with(DATA,ifelse(SEX%in%c("f","F"),"F",ifelse(SEX%in%c("m","M"),"M","U")))

  #Get TEPS interactions from Comments in boat.hdr  
Comments=subset(DATA,select=c(SHEET_NO,COMMENTS.hdr))
Comments=Comments[!is.na(Comments$COMMENTS.hdr),]
Comments=Comments[!duplicated(Comments$SHEET_NO),]

  #Add nature, fate, PCS and trophic level  
DATA=DATA%>%
  left_join(SPECIES_PCS_FATE%>%
            filter(SPECIES%in%unique(DATA$SPECIES))%>%
              dplyr::select(-c(COMMON_NAME,SCIENTIFIC_NAME)),
            by="SPECIES")

  #Special treatment for White shark (WP) and Grey nurse shark (GN) as they became protected throughout the period
if(is.null(protected.species))
{
  DATA$FATE=ifelse(DATA$SPECIES=="WP",ifelse(DATA$YEAR<1997,"C","D"), 
                   ifelse(DATA$SPECIES=="GN",ifelse(DATA$YEAR<2001,"C","D"),as.character(DATA$FATE)))
  
  DATA$NATURE=ifelse(DATA$SPECIES=="WP" | DATA$SPECIES=="GN",ifelse(DATA$FATE=="C","S","TEPS"),as.character(DATA$NATURE))           
  
}
                                               
  #remove unknown species codes
a=subset(DATA,is.na(FATE),selec=c(SPECIES,COMMON_NAME));a=a[!duplicated(a$SPECIES),]
NN.sp=a$SPECIES
if(length(NN.sp)>0) DATA=subset(DATA,!SPECIES%in%NN.sp)
DATA=DATA%>%filter(!SCIENTIFIC_NAME=='')

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
add.region=FALSE
if(add.region)
{
  DATA$REGION=as.character(with(DATA,ifelse(LONGITUDE>=124 & LONGITUDE<=129,"Region1",
                                            ifelse(LONGITUDE>=119 & LONGITUDE<124,"Region2",
                                                   ifelse(LONGITUDE>=116 & LONGITUDE<119,"Region3",
                                                          ifelse(LONGITUDE<116 & LATITUDE<=(-33),"Region4",
                                                                 ifelse(LATITUDE>(-33) & LATITUDE<=(-30),"Region5",
                                                                        ifelse(LATITUDE>(-30) & LATITUDE<=(-27),"Region6","Out.of.region"))))))))
  
}

  #Add area variable for temporal comparison
if(add.region)
{
  DATA$AREA=as.character(with(DATA,ifelse(BLOCK%in%c(2813,2814,2914),"Area1",
                                         ifelse(BLOCK%in%c(3114,3115,3214,3215),"Area2",
                                                ifelse(BLOCK%in%c(3314,3315,3414,3415),"Area3",
                                                       ifelse(BLOCK%in%c(3322,3323,3324),"Area4",NA))))))
  
}


  #Add period variable for spatial comparison
DATA$SEASON=as.character(with(DATA,ifelse(MONTH==12 | MONTH<3,"Summer",
         ifelse(MONTH>=3 & MONTH<6,"Autumn",
         ifelse(MONTH>=6 & MONTH<9,"Winter",
         ifelse(MONTH>=9 & MONTH<12,"Spring",NA))))))

#YR.selected=c(1995:1998,2002:2003) 
do.this=FALSE
if(do.this)  
{
  YR.selected=table(DATA$YEAR)
  YR.selected=as.numeric(names(YR.selected[YR.selected>1e3]))
  DATA$Yr.dummy=with(DATA,ifelse(YEAR%in%YR.selected,YEAR,NA))
  DATA$PER=as.character(with(DATA,paste(Yr.dummy,SEASON)))
  DATA$PERIOD=with(DATA,ifelse(substr(PER,1,2)=="NA",NA,PER))
  DATA=DATA[,-match(c("Yr.dummy","PER"),colnames(DATA))]
  
}

  #Add effort
DATA=DATA%>%
      mutate(EFFORT=ifelse(Method=='GN',NET_LENGTH*SOAK_TIME,
                    ifelse(Method=='LL',N.hooks*SOAK_TIME,
                    NA)))


  #Remove tropical species that belong to the Northern Territory and the "other sharks and scale fish" group
DATA=subset(DATA, FATE!="A" & SPECIES!="XX")

  #Add number caught
DATA$INDIVIDUALS=1
if(response.var=='cpue') DATA$INDIVIDUALS=DATA$INDIVIDUALS/DATA$EFFORT

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

  #Fill in teleost FL with a proportion of TL
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
                           grepl('d. rogers',SKIPPER)~"d.rogers", 
                           grepl('c. barnard',SKIPPER)~"c.barnard",
                           grepl('parker',SKIPPER)~"r.parker",
                           grepl('hayle',SKIPPER)~"p.hayler",
                           grepl('hoffha',SKIPPER)~"r. hoffhamer",
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
                          ZONE,BIOREGION,SEASON,INDIVIDUALS,
                          PCS,SCIENTIFIC_NAME,Fishery,Loo,K,tmax,tm,Lm,Body.shape,Habitat,Movement.scale,
                          Maximum.size,Body.depth,Head.length,Eye.diameter,Pre.orbital.length,Feeding.group,Trophic.level.fishbase))

#remove records with less than Min.shts per year-block
d=DATA[!duplicated(DATA$SHEET_NO),]
d$N=1
d1=aggregate(N~BLOCK+YEAR,d,sum)
N.base.case=subset(d1,N>=Min.shts)
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
fn.shw.appear=function(a)
{
  a=a%>%
    mutate(Julian=julian(as.Date(paste(YEAR,MONTH,1,sep='-')),origin=as.Date(paste(YEAR.min,1,1,sep='-'))))%>%
    distinct(MONTH,YEAR,Julian,SCIENTIFIC_NAME)%>%
    mutate(N=1)
  dummi=data.frame(SCIENTIFIC_NAME=names(sort(table(a$SCIENTIFIC_NAME))))%>%mutate(row_id=row_number())
  LabL=dummi$SCIENTIFIC_NAME
  names(LabL)=dummi$row_id
  a=a%>%
    left_join(dummi,by="SCIENTIFIC_NAME")
 p=a%>%
    ggplot(aes(Julian,row_id))+
    geom_point()+
    scale_y_continuous(labels=LabL,n.breaks=length(LabL),breaks=as.numeric(names(LabL)))+
    ylab('')+xlab('Julian day')+
    theme(axis.text.y = element_text(size = 7))
  print(p)
}
fn.shw.appear(a=DATA%>%filter(Fishery=='TDGDLF'))
ggsave("Outputs/Exploratory/Species appearing in observer data_TDGDLF.tiff",width = 6,height = 10,compression = "lzw")

if(use.NSF)
{
  fn.shw.appear(a=DATA%>%filter(Fishery=='NSF'))
  ggsave("Outputs/Exploratory/Species appearing in observer data_NSF.tiff",width = 6,height = 10,compression = "lzw")
}



  #3.4. Manipulated landings data (Monthly-aggregated records to 2021)
if(do.commercial)
{
  Data.monthly=Data.monthly%>%
            filter(METHOD=="GN" & zone%in%c('West','Zone1','Zone2') & LATITUDE<=(-26))%>%
            filter(NETLEN.c>=300 & EFFORT>0)%>%
            dplyr::select(-NETLEN.c)
  Data.daily=Data.daily%>%filter(METHOD=="GN" & zone%in%c('West','Zone1','Zone2') & LATITUDE<=(-26))%>%
    filter(NETLEN>=300 & EFFORT>0)
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
  
  #remove records with less than Min.shts per year-block
  d=Data.monthly[!duplicated(Data.monthly$SHEET_NO),]
  d$N=1
  d1=aggregate(N~BLOCK+YEAR,d,sum)
  N.base.case=subset(d1,N>=Min.shts)
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
  Data.monthly=subset(Data.Mon,!is.na(Loo),select=c(SHEET_NO,YEAR,MONTH,BLOCK,BOAT,TYPE.DATA,SPECIES,
                                        SCIENTIFIC_NAME,NATURE,LIVEWT.c,TROPHIC_LEVEL,TL_SD2,EFFORT,
                                        LATITUDE,LONGITUDE,ZONE,BIOREGION,
                                        Loo,K,tmax,tm,Lm,Movement.scale,
                                        Maximum.size,Body.depth,Head.length,Eye.diameter,Pre.orbital.length,Feeding.group,
                                        Trophic.level.fishbase,Body.shape,Habitat))
  
  Data.monthly=subset(Data.monthly,!is.na(Data.monthly$TROPHIC_LEVEL))
  Data.monthly$SPECIES=as.factor(Data.monthly$SPECIES)
  Data.monthly$INDIVIDUALS=Data.monthly$LIVEWT.c
  if(response.var=='cpue') Data.monthly$INDIVIDUALS=Data.monthly$INDIVIDUALS/Data.monthly$EFFORT
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
  ggsave("Outputs/Exploratory/Species appearing in logbook_TDGDLF.tiff",width = 6,height = 10,compression = "lzw")

  rm(Data.Mon,a)
  
  if(use.NSF)   
  {
    Data.monthly.north=Data.monthly.north%>%
      filter(METHOD=="LL" & zone%in%c('Joint','North') & LATITUDE>(-26))%>%
      dplyr::select(-NETLEN.c)
    Data.daily.north=Data.daily.north%>%filter(METHOD=="LL" & zone%in%c('Joint','North') & LATITUDE>(-26))
    if(combine.daily.monthly)
    {
      Data.daily.agg=Data.daily.north%>%
        group_by(YEAR,MONTH,BOAT,METHOD,BLOCK,SPECIES,SNAME,YEAR.c,
                 LATITUDE,LONGITUDE,TYPE.DATA,zone,SKIPPER)%>%
        summarise(LIVEWT.c=sum(LIVEWT.c))
      Effort.daily.agg=Data.daily.north%>%
        distinct(SHEET_NO,EFFORT,YEAR,MONTH,BOAT,METHOD,BLOCK)%>%
        group_by(YEAR,MONTH,BOAT,METHOD,BLOCK)%>%
        summarise(EFFORT=sum(EFFORT))
      Data.monthly.north=rbind(Data.monthly.north,
                         Data.daily.agg%>%
                           left_join(Effort.daily.agg,by=c('YEAR','MONTH','BOAT','METHOD','BLOCK'))%>%
                           mutate(SHEET_NO=paste(YEAR,MONTH,BOAT,METHOD,BLOCK))%>%
                           relocate(names(Data.monthly.north)))
    }
    
    #remove records with less than Min.shts per year-block
    d=Data.monthly.north[!duplicated(Data.monthly.north$SHEET_NO),]
    d$N=1
    d1=aggregate(N~BLOCK+YEAR,d,sum)
    N.base.case=subset(d1,N>=Min.shts)
    N.base.case$Yr.blk=with(N.base.case,paste(YEAR,BLOCK,sep="_"))
    Data.monthly.north$Yr.blk=with(Data.monthly.north,paste(YEAR,BLOCK,sep="_"))
    Data.monthly.north=subset(Data.monthly.north,Yr.blk%in%N.base.case$Yr.blk)
    
    #Remove boats with less than Min.recs
    AA=sort(table(Data.monthly.north$BOAT))
    AA=AA[AA>Min.recs]
    Data.monthly.north=subset(Data.monthly.north,BOAT%in%names(AA))
    
    #Remove years with less than Min.recs
    AA=sort(table(Data.monthly.north$YEAR))
    AA=AA[AA>Min.recs]
    Data.monthly.north=subset(Data.monthly.north,YEAR%in%names(AA))
    
    #Fix duplicated SPECIES
    Data.monthly.north=Data.monthly.north%>%
      mutate(SNAME=ifelse(SPECIES==22999,'shark, other',
                          ifelse(SPECIES==18003,'Dusky Whaler',
                                 ifelse(SPECIES==18014,'Blacktip Shark',
                                        SNAME))))
    
    #add species trophic level,etc 
    Data.Mon=merge(Data.monthly.north,SPECIES_PCS_FATE_Com,by="SPECIES",all.x=T,all.y=F)%>%  
      rename(ZONE=zone)%>%
      mutate(BIOREGION=case_when(LONGITUDE>=115.5 & LONGITUDE<=129 & LATITUDE<=(-26)~"SC",
                                 LONGITUDE<115.5 & LATITUDE<=(-27)~"WC",
                                 LONGITUDE<=114.834 & LATITUDE>(-27)~"Gascoyne",
                                 LONGITUDE>=114.834 & LONGITUDE<=129 & LATITUDE>=(-27)~"NC",
                                 TRUE ~ NA_character_),
             BIOREGION=ifelse(LATITUDE>(-34.5) & LATITUDE<(-30) &LONGITUDE< 118,"WC",BIOREGION))
    Data.monthly.north=subset(Data.Mon,!is.na(Loo),select=c(SHEET_NO,YEAR,MONTH,BLOCK,BOAT,TYPE.DATA,SPECIES,
                                                      SCIENTIFIC_NAME,NATURE,LIVEWT.c,TROPHIC_LEVEL,TL_SD2,EFFORT,
                                                      LATITUDE,LONGITUDE,ZONE,BIOREGION,
                                                      Loo,K,tmax,tm,Lm,Movement.scale,
                                                      Maximum.size,Body.depth,Head.length,Eye.diameter,Pre.orbital.length,Feeding.group,
                                                      Trophic.level.fishbase,Body.shape,Habitat))
    
    Data.monthly.north=subset(Data.monthly.north,!is.na(Data.monthly.north$TROPHIC_LEVEL))
    Data.monthly.north$SPECIES=as.factor(Data.monthly.north$SPECIES)
    Data.monthly.north$INDIVIDUALS=Data.monthly.north$LIVEWT.c
    if(response.var=='cpue') Data.monthly.north$INDIVIDUALS=Data.monthly.north$INDIVIDUALS/Data.monthly.north$EFFORT
    Data.monthly.north$FINYEAR=Data.monthly.north$YEAR
    Data.monthly.north$YEAR=as.numeric(substr(Data.monthly.north$YEAR,1,4))
    

    #Show when species start appearing in data set  
    YEAR.min=min(Data.monthly.north$YEAR)
    a=Data.monthly.north%>%
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
    ggsave("Outputs/Exploratory/Species appearing in logbook_NSF.tiff",width = 6,height = 10,compression = "lzw")
    
    rm(Data.Mon,a)
  }
}


  #3.5. Put data set in right format 
Data.list=list(Observer=DATA%>%
                          filter(Fishery=='TDGDLF' & SOAK_TIME>0.5 & NET_LENGTH>=0.3 &
                                        EFFORT>0 & MESH_SIZE%in%c(6.5,7))%>%
                          select(-Fishery)%>%
                          mutate(YEAR=ifelse(YEAR==2021,2020,YEAR))%>%
                          filter(!YEAR==1993))   #too few observations
if(do.commercial) Data.list$Logbook=Data.monthly
if(use.NSF)
{
  Data.list$Observer.north=DATA%>%filter(Fishery=='NSF')%>%select(-Fishery)
  if(do.commercial) Data.list$Logbook.north=Data.monthly.north
}

#drop dodgie effort
for(l in 1:length(Data.list))
{
  Data.list[[l]]=Data.list[[l]]%>%
    filter(!is.na(EFFORT))%>%
    filter(EFFORT>10 & EFFORT<2000)
  Data.list[[l]]%>%
    ggplot(aes(round(EFFORT)))+geom_bar()
}

#drop data sets if not meeting minimum number of record criteria  
Keep.dis=rep('no',length(Data.list))
for(l in 1:length(Data.list))
{
  Table.shots.year=with(Data.list[[l]]%>%distinct(SHEET_NO,YEAR),table(YEAR))
  iers=Table.shots.year[Table.shots.year>=Min.shots.year]
  Data.list[[l]]=Data.list[[l]]%>%filter(YEAR%in%as.numeric(names(iers)))
  if(length(iers)>=Min.years) Keep.dis[l]='yes'
}
Data.list=Data.list[which(Keep.dis=='yes')]

#starting year of analyses 
for(l in 1:length(Data.list))  Data.list[[l]]=Data.list[[l]]%>%filter(YEAR>=Min.YEAR) 

#fixing some dodgy scientific names
for(l in 1:length(Data.list))
  {
    Data.list[[l]]=Data.list[[l]]%>%
      mutate(SCIENTIFIC_NAME=case_when(SCIENTIFIC_NAME=='Pironace glauca'~'Prionace glauca',
                                       SCIENTIFIC_NAME=='Arripis truttaceus'~'Arripis truttacea',
                                       SCIENTIFIC_NAME=='Bothidae & Paralichthyidae'~'Bothidae & Pleuronectidae',
                                       SCIENTIFIC_NAME=='Centroberyx sp.'~'Centroberyx spp.',
                                       SCIENTIFIC_NAME=='Glaucosoma hebracium'~'Glaucosoma hebraicum',
                                       SCIENTIFIC_NAME=='Seriola hipos'~'Seriola hippos',
                                       SCIENTIFIC_NAME%in%c('Carcharhinus tilstoni','Carcharhinus limbatus')~'Carcharhinus limbatus/tilstoni',
                                       SCIENTIFIC_NAME=='Carcharodon Carcharias'~'Carcharodon carcharias',
                                       SCIENTIFIC_NAME=='Sphyrnidae'~'Sphyrna spp.',
                                       SCIENTIFIC_NAME=='Squatinidae'~'Squatina spp.',
                                       SCIENTIFIC_NAME=='Lethrinus laticaudis'~'Lethrinus spp.',
                                       SCIENTIFIC_NAME=='Carcharhiniformes'~'Carcharhinus spp.',
                                       TRUE~SCIENTIFIC_NAME),
             SCIENTIFIC_NAME=str_remove(SCIENTIFIC_NAME, '\\.'))
 }

#fix dodgy body shape
for(l in 1:length(Data.list))
{
  Data.list[[l]]=Data.list[[l]]%>%mutate(Body.shape=ifelse(Body.shape=="\tfusiform / normal","fusiform / normal",Body.shape))
}

#Get list of species by group
All.sp=vector('list',length(Data.list))
for(l in 1:length(All.sp))
{
  ss=Data.list[[l]]%>%distinct(SPECIES,SCIENTIFIC_NAME)
  if(any(grepl('.T',ss$SPECIES)))
  {
    ss=ss%>%
      mutate(Group=ifelse(grepl('.T',SPECIES),'Teleost','Elasmobranch'),
             Group=ifelse(SPECIES%in%c('BT','HT'),'Elasmobranch',Group))
  }else
  {
    ss=ss%>%
      mutate(SPECIES=as.numeric(as.character(SPECIES)),
             Group=ifelse(SPECIES>90030,'Teleost','Elasmobranch'))
  }

  

  
  All.sp[[l]]=ss
}
All.sp=do.call(rbind,All.sp)%>%
  distinct(SCIENTIFIC_NAME,Group)%>%arrange(Group,SCIENTIFIC_NAME)

#Export traits
a=Data.list$Logbook%>%distinct(SCIENTIFIC_NAME,
                               Body.shape,Maximum.size,tmax,Body.depth,Head.length,Eye.diameter,Pre.orbital.length,
                               Movement.scale,Feeding.group,TROPHIC_LEVEL,Habitat)%>%
  mutate(SCIENTIFIC_NAME=str_replace(SCIENTIFIC_NAME,'spp','spp.'),
         SCIENTIFIC_NAME=case_when(SCIENTIFIC_NAME=='Platycephalus laevigatus'~'Platycephalidae',
                                   SCIENTIFIC_NAME%in%c('Epinephelus daemelii','Epinephelus ergastularius')~'Epinephelus spp.',
                                   SCIENTIFIC_NAME=='Rhinobatidae'~'Rhinobatidae & Rhynchobatidae',
                                   TRUE~SCIENTIFIC_NAME))%>%
  distinct(SCIENTIFIC_NAME,.keep_all = T)
xx=read.csv(handl_OneDrive('Data/Species.code.csv'))%>%dplyr::select(COMMON_NAME,SCIENTIFIC_NAME)%>%
  mutate(SCIENTIFIC_NAME=str_remove(SCIENTIFIC_NAME,'Family '),
         SCIENTIFIC_NAME=str_remove(SCIENTIFIC_NAME,'Families '))%>%
  filter(SCIENTIFIC_NAME%in%a$SCIENTIFIC_NAME)%>%
  distinct(SCIENTIFIC_NAME,.keep_all = T)
a=a%>%left_join(xx,by='SCIENTIFIC_NAME')%>%
  relocate(COMMON_NAME,SCIENTIFIC_NAME)
write.csv(a,'Table functional diversity traits_logbook.csv',row.names = F)

  #3.6. Display raw species composition as proportions
setwd(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs"))
p.list=vector('list',length(Data.list))
for(l in 1:length(Data.list))
{
  NM=names(Data.list)[l]
  if(NM=="Logbook") NM="TDGDLF"
  if(NM=="Logbook.north") NM="NSF"
  #if(NM=="Observer") NM="Observer (TDGDLF)"
  p.list[[l]]=Catch.comp(ddd=Data.list[[l]]%>%
                           filter(!is.na(EFFORT))%>%
                           mutate(INDIVIDUALS=INDIVIDUALS*EFFORT)%>%
                           group_by(SCIENTIFIC_NAME,YEAR)%>%
                           summarise(INDIVIDUALS=sum(INDIVIDUALS,na.rm=T))%>%
                           ungroup()%>%
                           group_by(YEAR)%>%
                           mutate(Tot=sum(INDIVIDUALS),
                                  Prop=INDIVIDUALS/Tot)%>%
                           ungroup()%>%
                           mutate(Dataset=NM),
                         All.sp=All.sp,
                         Display='tile')
}
hei.vec=c(1, 0.8)
if(use.NSF) hei.vec=c(1, 0.8,0.6)
ggarrange(plotlist=p.list, common.legend = TRUE,ncol=1,heights = hei.vec)
ggsave(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Species composition.tiff"),
       width = 6.5,height = 10,compression = "lzw")

  #Streamgraph
do.this=FALSE  #very time consuming
if(do.this)
{
  p.list.stream=vector('list',length(Data.list))
  for(l in 1:length(Data.list))
  {
    NM=names(Data.list)[l]
    if(NM=="Logbook") NM="TDGDLF"
    if(NM=="Logbook.north") NM="NSF"
    
    p.list.stream[[l]]=Catch.comp.stream(ddd=Data.list[[l]]%>%
                                           filter(!is.na(EFFORT))%>%
                                           mutate(INDIVIDUALS=INDIVIDUALS*EFFORT)%>%
                                           group_by(SCIENTIFIC_NAME,YEAR)%>%
                                           summarise(INDIVIDUALS=sum(INDIVIDUALS,na.rm=T))%>%
                                           ungroup()%>%
                                           group_by(YEAR)%>%
                                           mutate(Tot=sum(INDIVIDUALS),
                                                  Prop=INDIVIDUALS/Tot)%>%
                                           ungroup()%>%
                                           mutate(Dataset=NM),
                                         All.sp=All.sp)
  }
  hei.vec=c(1, 0.8)
  if(use.NSF) hei.vec=c(1, 0.8,0.6)
  ggarrange(plotlist=p.list.stream, common.legend = TRUE,ncol=1,heights = hei.vec)
  ggsave(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Species composition_streamgraph.tiff"),
         width = 6.5,height = 10,compression = "lzw")
}


  #data sets combined
p.list=vector('list',length(Data.list))
names(p.list)=names(Data.list)
for(l in 1:length(Data.list))
{
  NM=names(Data.list)[l]
  if(NM=="Logbook") NM="TDGDLF"
  if(NM=="Logbook.north") NM="NSF"
  
  ddd=Data.list[[l]]%>%
    filter(!is.na(EFFORT))
  
  dumi=ddd%>%
    mutate(Ktch=INDIVIDUALS*EFFORT)%>%
    group_by(SCIENTIFIC_NAME)%>%
    summarise(Ktch=sum(Ktch))%>%
    ungroup()%>%
    arrange(-Ktch)%>%
    mutate(Cumktch=cumsum(Ktch),
           Tot=sum(Ktch),
           Cumktch2=Cumktch/Tot,
           SCIENTIFIC_NAME2=ifelse(Cumktch2<=0.96,SCIENTIFIC_NAME,'Other'))%>%
    left_join(All.sp,by='SCIENTIFIC_NAME')
  nelas=nrow(dumi%>%filter(SCIENTIFIC_NAME2=='Other' & Group== 'Elasmobranch'))
  ntel=nrow(dumi%>%filter(SCIENTIFIC_NAME2=='Other' & Group== 'Teleost'))
  dumi=dumi%>%
    mutate(SCIENTIFIC_NAME2=case_when(SCIENTIFIC_NAME2=='Other' & Group== 'Elasmobranch'~paste0('Other elasmobranchs (n=',nelas,' species)'),
                                      SCIENTIFIC_NAME2=='Other' & Group== 'Teleost'~paste0('Other teleosts (n=',ntel,' species)'),
                                      TRUE~SCIENTIFIC_NAME2))
  
  ddd=ddd%>%
    left_join(dumi%>%distinct(SCIENTIFIC_NAME,SCIENTIFIC_NAME2,Group),by='SCIENTIFIC_NAME')%>%
    dplyr::select(-SCIENTIFIC_NAME)%>%
    rename(SCIENTIFIC_NAME=SCIENTIFIC_NAME2)%>%
    mutate(Ktch=INDIVIDUALS*EFFORT)%>%
    group_by(SCIENTIFIC_NAME,Group)%>%
    summarise(Ktch=sum(Ktch))%>%
    ungroup()%>%
    mutate(color=ifelse(Group=='Teleost','dodgerblue4',
                        ifelse(Group=='Elasmobranch','firebrick3',NA)),
           color=ifelse(is.na(color),'forestgreen',color))
  
  LVLs=ddd%>%distinct(SCIENTIFIC_NAME,color,Group,Ktch)%>%arrange(desc(Group),Ktch)%>%pull(SCIENTIFIC_NAME)
  a=ddd%>%distinct(SCIENTIFIC_NAME,color)%>%arrange(factor(SCIENTIFIC_NAME,levels=LVLs))
  FILL=ddd%>%distinct(SCIENTIFIC_NAME,Group,Ktch)%>%arrange(Group,Ktch)
  colfunc <- colorRampPalette(c("lightpink", "firebrick4"))
  FILL.elas=colfunc(nrow(FILL%>%filter(Group=='Elasmobranch')))
  colfunc <- colorRampPalette(c("cadetblue2", "dodgerblue4"))
  FILL.tel=colfunc(nrow(FILL%>%filter(Group=='Teleost')))
  FILL=FILL%>%
        dplyr::select(-c(Group,Ktch))%>%
        mutate(FILL=c(FILL.elas,FILL.tel))
  fill.vec=FILL$FILL
  names(fill.vec)=FILL$SCIENTIFIC_NAME
  p.list[[l]]=ddd%>%
                  mutate(Tot=sum(Ktch),
                         Prop=100*Ktch/Tot,
                         SCIENTIFIC_NAME=factor(SCIENTIFIC_NAME,levels=LVLs))%>%
                  ggplot(aes(SCIENTIFIC_NAME,Prop))+
                  geom_bar(aes(fill=SCIENTIFIC_NAME),stat="identity", width=1) +
                  coord_flip()+
                  scale_fill_manual(values=fill.vec,drop=FALSE)+
                  ylab('Percentage')+xlab('')+
                  ggtitle(NM)+
                  theme_PA(axs.t.siz=10,Ttl.siz=15)+
                  theme(legend.position = 'none',
                        plot.title.position = "plot")+
                  scale_y_continuous(limits = c(0, NA))
  
}
hei.vec=c(1, 0.8)
if(use.NSF) hei.vec=c(1, 0.8,0.6)
ggarrange(plotlist=p.list, ncol=1,heights = hei.vec)
ggsave(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Univariate/Species composition_pie.tiff"),
       width = 6.5,height = 10,compression = "lzw")

    #by data set
for(l in 1:length(p.list))
{
  NM=names(Data.list)[l]
  if(NM=="Logbook") NM="TDGDLF"
  if(NM=="Logbook.north") NM="NSF"
  print(p.list[[l]]+labs(title = NULL) )+theme(axis.text = element_text(size = 11.5))
  ggsave(handl_OneDrive(paste0("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Univariate/Species composition_pie_",NM,".tiff")),
         width = 6,height = 8,compression = "lzw")
}

rm(DATA.bio,DATA.ecosystems,DATA,Boat_hdr)

# 4 Ecosystems indicators analyses-----------------------------------------------------------------------

  #4.1. calculate indicators

# Define which ecosystem indicator to consider
Resp.vars_observer=c(subset(Ecol.Indicators,!Ecol.Indicators%in%c("FIB")),traits,Functional.diversity) 
Resp.vars_logbook=c(Ecol.Indicators%>%subset(Ecol.Indicators%in%c("Shannon","Pielou","Simpson","MTL")),
                    traits,Functional.diversity)

Store.data.list=vector('list',length(Data.list))
names(Store.data.list)=names(Data.list)
Check.each.fun.rich=FALSE
for(l in 1:length(Data.list))  #takes 13 minutes for observer and logbook data sets
{
  print(paste("Ecosystems indicators calculation for -----------",names(Data.list)[l]))
  
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
    
    fn.fig(paste0("Exploratory/Number of sheets by year_",names(Data.list)[l]),2000,2000)
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
    
    fn.fig(paste0("Exploratory/Map of number of sheets per block_",names(Data.list)[l]),2000,2000) 
    par(las=1,mar=c(1.75,3.5,1.5,.1),oma=c(1.5,1,.1,.1),mgp=c(1,.8,0),cex.lab=1.25)
    MIN=max(-36,min(surveys$LAT_BLOCK))
    maps::map("worldHires",xlim=c(112.95, 129.5),ylim=c(MIN, max(surveys$LAT_BLOCK)),col="grey80",fill=F)
    for(i in seq(surveys$SURVEYS))
    {
      rect(surveys$LONG_BLOCK[i],surveys$LAT_BLOCK[i]-1,surveys$LONG_BLOCK[i]+1,surveys$LAT_BLOCK[i],border="blue",col=ocean.pal[surveys$SURVEYS[i]])
    }
    maps::map("worldHires",xlim=c(112.95, 129.5),ylim=c(MIN, max(surveys$LAT_BLOCK)),col="grey80",fill=T,border="grey80",add=T)
    rect(surveys$LONG_BLOCK,surveys$LAT_BLOCK-1,surveys$LONG_BLOCK+1,surveys$LAT_BLOCK,border="blue")
    mtext("Total sheet numbers per block",line=1,cex=1.5)
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
      xlim(112.5,129)+ylim(MIN,max(surveys$LAT_BLOCK))
    ggsave(paste0("Exploratory/Map of number of sheets per block per year_",names(Data.list)[l],'.tiff'),
           width = 8,height = 8,compression = "lzw")
    
    #Number of shots per year-block-month
    a=Data.list[[l]][!duplicated(Data.list[[l]]$SHEET_NO),]
    a$N=1
    a=aggregate(N~YEAR+MONTH+BLOCK+LATITUDE+LONGITUDE,a,sum)
    
    yrs=sort(unique(a$YEAR))
    pdf(paste0('Exploratory/Shots per year-block-month_',names(Data.list)[l],'.pdf'))
    for(i in 1:length(yrs)) fn.see(d=subset(a,YEAR==yrs[i]))
    dev.off()
    rm(a)
    
    #Export species by year
    if(grepl("Observer",names(Data.list)[l]))dis.var='COMMON_NAME'
    if( grepl("Logbook",names(Data.list)[l]))dis.var='SCIENTIFIC_NAME'
    TBLA=table(paste(Data.list[[l]][,dis.var],Data.list[[l]]$TROPHIC_LEVEL),Data.list[[l]]$YEAR)
    write.csv(TBLA,paste0('Exploratory/TBLA_species.TL_year_',names(Data.list)[l],'.csv'),row.names=T)
    TBLA=table(Data.list[[l]][,dis.var],Data.list[[l]]$YEAR)
    write.csv(TBLA,paste0('Exploratory/TBLA_species_year_',names(Data.list)[l],'.csv'),row.names=T)
    TBLA=table(Data.list[[l]]$TROPHIC_LEVEL,Data.list[[l]]$YEAR)
    colfunc <- colorRampPalette(c("red", "yellow"))
    CL=rep("grey40",nrow(TBLA))
    xx=as.numeric(colnames(TBLA))
    mltplr=5
    fn.fig(paste0('Exploratory/TrophicLevel_by_yr_',names(Data.list)[l]),2000,2000) 
    par(las=1,mar=c(1.75,3.5,1.5,.1),oma=c(1.5,1,.1,.1),mgp=c(1,.8,0),cex.lab=1.25)
    plot(xx,rep(1,length(xx)),cex=(TBLA[1,]/max(TBLA))*mltplr,ylab="",pch=19,col=CL[1],xlab="",yaxt='n',ylim=c(0,nrow(TBLA)))
    for(n in 2:nrow(TBLA)) points(xx,rep(n,length(xx)),cex=(TBLA[n,]/max(TBLA))*mltplr,pch=19,col=CL[n])
    axis(2,1:nrow(TBLA),row.names(TBLA))
    mtext("Trophic level",2,2.5,cex=1.5,las=3)
    mtext("Year",1,2,cex=1.5)
    dev.off()

   }
  

  #2. Calculate ecological and functional indicators
  idvarS=IDVAR[which(IDVAR%in%names(Data.list[[l]]))]
  if(grepl("Observer",names(Data.list)[l])) resp.vars=Resp.vars_observer
  if(grepl("Logbook",names(Data.list)[l]))  resp.vars=Resp.vars_logbook
  
  n.rv=length(resp.vars)
  Main.title=resp.vars
  Res.var.in.log=rep("NO",length(resp.vars))  #fit response var in log space or not? 
  
  ddd=Data.list[[l]]%>%
    mutate(YEAR=as.character(YEAR),
           MONTH=as.character(MONTH),
           BLOCK=as.character(BLOCK),
           SPECIES=as.factor(SPECIES),
           Yr.blk=paste(YEAR,BLOCK,sep="_"))
  
  OUT.base.case=fn.calc.ecol.ind(DaTA=ddd,  
                                 normalised="YES",
                                 Drop.yrs="NO",
                                 idvarS=idvarS,  
                                 resp.vars=resp.vars,
                                 check.each.fun.rich.var=Check.each.fun.rich)
  if(check.discards & grepl("Observer",names(Data.list)[l]))
  {
    OUT.base.case.retained=fn.calc.ecol.ind(DaTA=ddd%>%filter(FATE=='C'),
                                            dat.nm="BC",
                                            normalised="YES",
                                            Drop.yrs="NO",
                                            idvarS=idvarS,
                                            resp.vars=subset(resp.vars,!resp.vars=="Prop.Disc"),
                                            check.each.fun.rich.var=FALSE)
    OUT.base.case.discarded=fn.calc.ecol.ind(DaTA=ddd%>%filter(FATE=='D'),
                                            dat.nm="BC",
                                            normalised="YES",
                                            Drop.yrs="NO",
                                            idvarS=idvarS,
                                            resp.vars=subset(resp.vars,!resp.vars=="Prop.Disc"),
                                            check.each.fun.rich.var=FALSE)
  }
  
  #Exploratory analyses 
  if(do.exploratory=="YES")
  {
    #Show species composition and MTL  
      a=Data.list[[l]]
      
      SORT=a%>%distinct(SCIENTIFIC_NAME,.keep_all = T)%>%arrange(NATURE)%>%pull(SCIENTIFIC_NAME)
      a%>%
        mutate(SCIENTIFIC_NAME=factor(SCIENTIFIC_NAME,levels=SORT))%>%
        group_by(YEAR,SCIENTIFIC_NAME)%>%
        summarise(Tonnes=sum(LIVEWT.c,na.rm=T)/1000,
                  TL=mean(TROPHIC_LEVEL,na.rm=T))%>%
        ggplot(aes(YEAR,SCIENTIFIC_NAME,size=Tonnes,color=TL))+
        geom_point()+
        theme_PA(axs.t.siz=8)+ylab('Scientific name')
      ggsave(paste0("Univariate/catch comp by year_",names(Data.list)[l],'.tiff'),
             width = 6,height = 10.5,compression = "lzw")
      b=a%>%
        group_by(SCIENTIFIC_NAME)%>%
        mutate(Min.year=min(YEAR))
      write.csv(b%>%
                  select(all_of(c('SCIENTIFIC_NAME','Min.year',traits_ecol,traits_morph)))%>%
                  distinct(SCIENTIFIC_NAME,.keep_all = T)%>%arrange(Min.year),
                paste0("catch comp by year_traits_",names(Data.list)[l],'.csv'),row.names = F)
      
  
    
    #Traits thru time (annual averages)  
    fn.viol.box(d=OUT.base.case,byzone=FALSE,filcol='lightsalmon2',c(traits_ecol,traits_morph))
    ggsave(paste0("Exploratory/bxplt_year_traits_",names(Data.list)[l],'.tiff'),
           width = 10,height = 8,compression = "lzw")
    
    # Habitat         Rep
    # reef-associated   1
    # demersal          2
    # benthopelagic     3
    # pelagic-oceanic   4
    # pelagic-neritic   5
    
    # Movement.scale Rep
    #           0-100   1
    #             500   2
    #         100-500   3
            
  # Feeding.group                 Rep
  # hunting macrofauna (predator)   1
  
  # Body.shape Rep
  #    elongated   1
  #    fusiform / normal   2
  #  short and / or deep   3

    
    #Look at annual averages  
    fn.viol.box(d=OUT.base.case,byzone=FALSE,filcol='lightsalmon2',resp.vars)
    ggsave(paste0("Exploratory/bxplt_year_",names(Data.list)[l],'.tiff'),
           width = 10,height = 8,compression = "lzw")
    # fn.viol.box(d=OUT.base.case,byzone=TRUE,resp.vars=resp.vars)
    # ggsave(paste0("Outputs/Exploratory/bxplt_year.zone_",names(Data.list)[l],'.tiff'),
    #        width = 8,height = 6,compression = "lzw")
    if(check.discards & grepl("Observer",names(Data.list)[l]))
    {
      fn.viol.box(d=OUT.base.case.retained,byzone=FALSE,filcol='lightsalmon2',resp.vars=subset(resp.vars,!resp.vars=="Prop.Disc"))
      ggsave(paste0("Outputs/Exploratory/bxplt_year_",names(Data.list)[l],'_retained.tiff'),
             width = 8,height = 6,compression = "lzw")
      fn.viol.box(d=OUT.base.case.discarded,byzone=FALSE,filcol='lightsalmon2',resp.vars=subset(resp.vars,!resp.vars=="Prop.Disc"))
      ggsave(paste0("Outputs/Exploratory/bxplt_year_",names(Data.list)[l],'_discarded.tiff'),
             width = 8,height = 6,compression = "lzw")
    }
    
     #Response var error structure
    fn.fig(paste0("Exploratory/error_structure_Obs",names(Data.list)[l]),2000,2000) 
    par(mfrow=c(n.rv,2),las=1,mar=c(1,1.5,1,1.2),oma=c(1,1,.1,1),mgp=c(1,.5,0))
    for(i in 1:n.rv)
    {
      density.response(d=OUT.base.case,Var=resp.vars[i],Log="NO")
      density.response(d=OUT.base.case,Var=resp.vars[i],Log="YES")
    }
    dev.off()
    
    #Check if correlation between effort and indices
    fn.fig(paste0("Exploratory/Correlation_effort_indices_",names(Data.list)[l]),2000,2000) 
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
  
  #Store results
  Store.data.list[[l]]=OUT.base.case
  rm(OUT.base.case)
}

  # Display raw indicators as proportions 
p.list=vector('list',length(Store.data.list))
for(l in 1:length(Store.data.list))
{
  if(grepl("Observer",names(Store.data.list)[l])) Resp.vars=Resp.vars_observer
  if(grepl("Logbook",names(Store.data.list)[l]))  Resp.vars=Resp.vars_logbook
  
  NM=names(Store.data.list)[l]
  if(NM=="Logbook") NM="TDGDLF"
  if(NM=="Logbook.north") NM="NSF"
  #if(NM=="Observer") NM="Observer (TDGDLF)" max(ddd%>%filter(Indicator=='K')%>%pull(Prop))
  
  XLAB=''
  if(l==length(Store.data.list)) XLAB='Financial year'
  
  ddd=Store.data.list[[l]]%>%
    filter(!is.na(EFFORT))%>%
    dplyr::select(YEAR,all_of(Resp.vars))%>%
    gather(Indicator,Prop,-YEAR)%>%
    group_by(Indicator,YEAR)%>%
    summarise(Prop=mean(Prop,na.rm=T))%>%
    ungroup()%>%
    group_by(Indicator)%>%
    mutate(Prop.mean=mean(Prop,na.rm=T))%>%
    ungroup()%>%
    mutate(Rel.value=Prop/Prop.mean,
           Indicator=case_when(Indicator=="FnRich_morph"~"Functional richness (morph.)",
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
  
  p.list[[l]]=ddd%>%
    ggplot(aes(YEAR, Indicator , fill= Rel.value)) + 
    geom_tile()+
    scale_fill_gradient2(low="lightgoldenrodyellow",mid="darkgoldenrod2", high="brown4",
                         midpoint = mean(range(ddd$Rel.value)),name = "Relative value")+
    ylab('')+xlab(XLAB)+
    theme_PA(axs.t.siz=10,Ttl.siz=13)+
    theme(legend.position = 'top',
          plot.title.position = "plot",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(NM)
}
hei.vec=c(1, 1)
if(use.NSF) hei.vec=c(1, 1,0.6)
ggarrange(plotlist=p.list, common.legend = TRUE,ncol=1,heights = hei.vec)
ggsave(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Univariate/Indicators by year.tiff"),
       width = 6.5,height = 8,compression = "lzw")


  #4.2. Stats     
setwd(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Univariate"))
Store.out=vector('list',length(Data.list))
names(Store.out)=names(Data.list)
tic()    #takes 16 mins
for(l in 1:length(Store.out)) 
{
  print(paste("Ecosystems indicators stats for -----------",names(Data.list)[l]))
  
  IdvarS=IDVAR[which(IDVAR%in%names(Store.data.list[[l]]))]
  if(grepl("Observer",names(Store.data.list)[l])) Resp.vars=Resp.vars_observer
  if(grepl("Logbook",names(Store.data.list)[l]))  Resp.vars=Resp.vars_logbook
  
  inters="YES"
  if(names(Store.data.list)[l]=="Logbook.north") inters="NO"   #not enough degrees of freedom
  
  Store.out[[l]]=fn.apply.model(DaTA=Store.data.list[[l]],
                                dat.nm=names(Store.data.list)[l],
                                normalised="YES",
                                Drop.yrs="NO",
                                idvarS=IdvarS,
                                resp.vars=Resp.vars,
                                ADD.INTER=inters)
  
  Store.out[[l]]=Store.out[[l]]%>%
                    mutate(Data.set=names(Store.data.list)[l],
                           Indicator=case_when(Indicator=="FnRich_morph"~"Functional richness (morph.)",
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
}
toc()

#fix indicator label
for(l in 1:length(Store.out)) 
{
  Store.out[[l]]=Store.out[[l]]%>%mutate(Indicator=ifelse(Indicator=="Proportion mature",'Age mat. over Max age',Indicator))
}

#4.3. Display Year prediction  
p.list=vector('list',length(Data.list))
names(p.list)=names(Data.list)

  #Relative by data set
ADD.smoother=FALSE
for(l in 1:length(Store.out)) 
{
  p.list[[l]]=fn.plot.preds(d=Store.out[[l]]%>%
                              filter(Variable=='YEAR' & Relative=='YES')%>%
                              mutate(Value=as.numeric(Value),
                                     Zone=ifelse(LONGITUDE1<116 & LATITUDE1>(-33),'West',
                                                 ifelse(LONGITUDE1<116 & LATITUDE1<=(-33),'Zone 1',
                                                        'Zone 2'))),
                            add.smoother=ADD.smoother,
                            YLAB='Relative value',
                            add.zone=FALSE)
  
  print(p.list[[l]])
  ggsave(paste0("Year prediction_",names(Store.out)[l],"_relative.tiff"),width = 6,height = 6,compression = "lzw")
}

  #Relative data sets combined  
do.call(rbind,Store.out)%>%
  filter(Variable=='YEAR'& Relative=='YES')%>%
  mutate(Value=as.numeric(Value))%>%
  mutate(yr=case_when(Data.set=='Logbook'~Value+0.25,
                      Data.set=='Logbook.north'~Value-0.25,
                      TRUE~Value),
         Data.set=case_when(Data.set=='Logbook'~"TDGDLF",
                            Data.set=='Logbook.north'~'NSF',
                            TRUE~Data.set))%>%
  ggplot(aes(yr,MeAn,color=Data.set))+
# geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95,alpha=0.35)+
  geom_point(alpha=0.8,size=1.1)+
  geom_errorbar(aes(ymin = LowCI, ymax = UppCI))+
  geom_line(alpha=.6,linetype='dotted')+
  facet_wrap(~Indicator,scales='free_y',ncol=2)+
  #scale_y_continuous(limits = c(0, NA))+
  theme_PA(strx.siz=10)+ylab('Relative value')+xlab('Financial year')+
  theme(legend.position = 'top',
        legend.title=element_blank())
ggsave(paste0("Year prediction_Combined_relative.tiff"),width = 5,height = 8,compression = "lzw")


#Absolute by data set
for(l in 1:length(Store.out)) 
{
  fn.plot.preds(d=Store.out[[l]]%>%
                  filter(Variable=='YEAR' & Relative=='NO')%>%
                  mutate(Value=as.numeric(Value),
                         Zone=ifelse(LONGITUDE1<116 & LATITUDE1>(-33),'West',
                                     ifelse(LONGITUDE1<116 & LATITUDE1<=(-33),'Zone 1',
                                            'Zone 2'))),
                add.smoother=ADD.smoother,
                YLAB='Indicator value',
                add.zone=FALSE)
  ggsave(paste0("Year prediction_",names(Store.out)[l],".tiff"),width = 6,height = 6,compression = "lzw")
}

#Published observer
fn.plot.preds(d=Store.out$Observer%>%
                filter(!Indicator%in%c("Proportion of discards","Mean maximum length"))%>%
                filter(Variable=='YEAR' & Relative=='NO')%>%
                mutate(Value=as.numeric(Value),
                       Zone=ifelse(LONGITUDE1<116 & LATITUDE1>(-33),'West',
                                   ifelse(LONGITUDE1<116 & LATITUDE1<=(-33),'Zone 1',
                                          'Zone 2'))),
              add.smoother=ADD.smoother,
              YLAB='Indicator value',
              add.zone=FALSE)
ggsave(paste0("Year prediction_Observer_published.tiff"),width = 6,height = 6,compression = "lzw")

#Absolute data sets combined  
do.call(rbind,Store.out)%>%
  filter(Variable=='YEAR'& Relative=='NO')%>%
  mutate(Value=as.numeric(Value))%>%
  mutate(yr=case_when(Data.set=='Logbook'~Value+0.25,
                      Data.set=='Logbook.north'~Value-0.25,
                      TRUE~Value),
         Data.set=case_when(Data.set=='Logbook'~"TDGDLF",
                            Data.set=='Logbook.north'~'NSF',
                            TRUE~Data.set))%>%
  ggplot(aes(yr,MeAn,color=Data.set))+
# geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95,alpha=0.35)+
  geom_point(alpha=0.8,size=1.1)+
  geom_errorbar(aes(ymin = LowCI, ymax = UppCI))+
  geom_line(alpha=.6,linetype='dotted')+
  facet_wrap(~Indicator,scales='free_y',ncol=2)+
  #scale_y_continuous(limits = c(0, NA))+
  theme_PA(strx.siz=10)+ylab('Indicator value')+xlab('Financial year')+
  theme(legend.position = 'top',
        legend.title=element_blank())
ggsave(paste0("Year prediction_Combined.tiff"),width = 5,height = 8,compression = "lzw")


#4.4. Display Month prediction  
  #by data set
for(l in 1:length(Store.out)) 
{
  Store.out[[l]]%>%
    filter(Variable=='MONTH')%>%
    mutate(Value=as.numeric(Value))%>%
    ggplot(aes(Value,MeAn))+
    geom_point()+
    geom_errorbar(aes(ymin = LowCI, ymax = UppCI))+
    facet_wrap(~Indicator,scales='free_y')+
    scale_y_continuous(limits = c(0, NA))+
    theme_PA(strx.siz=8.5)+ylab('Relative value')+xlab('Month')+
    scale_x_continuous(breaks=seq(1,12,1))
  ggsave(paste0("Month prediction_",names(Store.out)[l],".tiff"),width = 8,height = 6,compression = "lzw")
}

  #data sets combined  
do.call(rbind,Store.out)%>%
  filter(Variable=='MONTH')%>%
  mutate(Value=as.numeric(Value))%>%
  mutate(yr=case_when(Data.set=='Logbook'~Value+0.25,
                      Data.set=='Logbook.north'~Value-0.25,
                      TRUE~Value),
         Data.set=case_when(Data.set=='Logbook'~"TDGDLF",
                            Data.set=='Logbook.north'~'NSF',
                            TRUE~Data.set))%>%
  ggplot(aes(yr,MeAn,color=Data.set))+
  geom_point(alpha=0.8,size=0.9)+
  geom_errorbar(aes(ymin = LowCI, ymax = UppCI))+
  #geom_line()+
  facet_wrap(~Indicator,scales='free_y',ncol=2)+
#  scale_y_continuous(limits = c(0, NA))+
  theme_PA(strx.siz=10)+ylab('Relative value')+xlab('Month')+
  theme(legend.position = 'top',
        legend.title=element_blank())+
  scale_x_continuous(breaks=seq(1,12,1))
ggsave(paste0("Month prediction_Combined.tiff"),width = 5,height = 8,compression = "lzw")


#4.5. Display Lat and Long prediction  
  #by data set
for(l in 1:length(Store.out)) 
{
  Store.out[[l]]%>%
    filter(Variable=='LATITUDE LONGITUDE')%>%
    mutate(Lat=as.numeric(gsub( " .*$", "", Value)),
           Long=as.numeric(gsub("^\\S+ ", "", Value)))%>%
    ggplot(aes(Long,Lat,size=MeAn,color=MeAn))+
    geom_point(alpha=0.6)+
    facet_wrap(~Indicator,scales='free_y')+
    theme_PA(strx.siz=8.5)+ylab('Latitude')+xlab('Longitude')+
    theme(legend.title = element_blank(),
          legend.position = 'top')
  ggsave(paste0("Lat long prediction_",names(Store.out)[l],".tiff"),width = 8,height = 6,compression = "lzw")
}

# 5 Multivariate analyses-------------------------------------------------------------------------
if(do.multivariate)
{
  model.type='lm'  #'gam' could not be used as anova.manyany summary.manyany not implemented in package
  if(model.type=='gam')
  {
    do.Gam=TRUE
    do.Lm=FALSE
  }
  if(model.type=='lm')
  {
    do.Gam=FALSE
    do.Lm=TRUE
  }
  
  Store.out.multi=vector('list',length(Data.list))
  names(Store.out.multi)=names(Data.list)
  Grouping.vars=c('sheet_no','YEAR','MONTH','LATITUDE','LONGITUDE','BOAT')
  if(do.Gam) TERM.form.g='year+s(month, k = 12, bs = "cc")+ s(boat, bs = "re")+s(latitude, longitude)'
  if(do.Lm) TERM.form.g=TERM.form.perm='year+latitude*longitude'
  Group.term.ord=labels(terms(formula(paste('y',TERM.form.g,sep='~'))))
  if(do.Gam) this=grep("s",Group.term.ord)
  if(do.Lm) this=grep(':',Group.term.ord)
  if(length(this)>0) Group.term.ord=Group.term.ord[-this]
  tic()
  for(l in 1:length(Store.out.multi))   
  {
    print(paste("Multivariate stats for -----------",names(Data.list)[l]))
    
    Store.out.multi[[l]]=multivariate.fn(d=Data.list[[l]]%>%
                                           mutate(YEAR=as.factor(YEAR),
                                                  MONTH=as.integer(MONTH),
                                                  BLOCK=as.character(BLOCK),
                                                  SPECIES=as.character(SPECIES)),
                                         Terms=tolower(Grouping.vars),
                                         Def.sp.term=c('species','year'),
                                         Transf='',
                                         Show.term='year',
                                         Group='95',
                                         hndl=handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Multivariate"),
                                         All.species.names=rbind(Data.list[[l]]%>%distinct(SPECIES,SCIENTIFIC_NAME),
                                                                 data.frame(SPECIES="Other",SCIENTIFIC_NAME="Other")),
                                         dat.name=names(Data.list)[l],
                                         do.MDS=TRUE,do.pcoa=FALSE, do.gam=FALSE, do.tweedie=TRUE, do.lm=FALSE,
                                         do.boral=FALSE, do.mvabund=TRUE, do.permanova=FALSE, do.simper=FALSE,
                                         group.ordination=TRUE,aggregate.monthly=TRUE,
                                         Group.term.ordination=Group.term.ord,
                                         term.form.model=TERM.form.g,
                                         term.form.permanova=TERM.form.perm,
                                         use.Other=TRUE,
                                         dis.lat= Pred.lat, dis.long=Pred.long,
                                         do.anova=FALSE)
    
  }
  
  
  #Display year effect 
    #data sets combined
  Min.yr=min(Data.list$Logbook$YEAR)
  Max.yr=max(Data.list$Logbook$YEAR)
  p.list=vector('list',length(Store.out.multi))
  names(p.list)=names(Store.out.multi)
  for(l in 1:length(p.list))
  {
    NM=names(Data.list)[l]
    if(NM=="Logbook") NM="TDGDLF"
    if(NM=="Logbook.north") NM="NSF"
    #if(NM=="Observer") NM="Observer (TDGDLF)"
    
    ddd=Store.out.multi[[l]]$Preds%>%
            filter(!species=="Other")%>%
            mutate(year=as.numeric(as.character(year)))%>%
      left_join(All.sp,by=c('scientific_name'='SCIENTIFIC_NAME'))%>%
      mutate(color=ifelse(Group=='Teleost','dodgerblue4',
                   ifelse(Group=='Elasmobranch','firebrick3',NA)),
             color=ifelse(is.na(color),'black',color))
    
    LVLs=ddd%>%distinct(scientific_name,color,Group)%>%arrange(desc(Group),scientific_name)%>%pull(scientific_name)
    a=ddd%>%distinct(scientific_name,color,Group)%>%arrange(desc(Group),factor(scientific_name,levels=LVLs))
    
    ddd=ddd%>%    
      rename(Proportion=Median)%>%
      mutate(Proportion=ifelse(Proportion<0,0,Proportion),
             scientific_name=factor(scientific_name,levels=LVLs))%>%
      mutate(Zone=ifelse(longitude<116 & latitude>(-33),'West',
                         ifelse(longitude<116 & latitude<=(-33),'Zone 1','Zone 2')))
    xlAb=''
    if(l==length(p.list)) xlAb='Financial year'
    Nzones=length(unique(ddd$Zone))
    p.list[[l]]=ddd%>%
                group_by(year,latitude,longitude)%>%
                mutate(Prop.total=sum(Proportion))%>%
                ungroup()%>%
                mutate(Proportion=Proportion/Prop.total)%>%
                ggplot(aes(year, scientific_name , fill= Proportion)) + 
                geom_tile()+
                facet_wrap(~Zone,ncol=Nzones)+
                 scale_fill_gradient(low="lightblue1", high="darkorange4",breaks = seq(0,0.5,length.out=3))+
                # scale_fill_gradient2(low="lightblue1",mid="darkgoldenrod2", high="brown4",midpoint = mean(range(ddd$Proportion)))+
                ylab('')+xlab(xlAb)+
                ggtitle(NM)+
                theme_PA(axs.t.siz=8,Ttl.siz=15)+
                theme(legend.position = 'top',
                      plot.title.position = "plot",
                      axis.text.y = element_text(colour = a%>%pull(color)))+
                xlim(Min.yr,Max.yr)

  }
  hei.vec=c(1, 0.8)
  if(use.NSF) hei.vec=c(1, 0.8,0.6)
  ggarrange(plotlist=p.list, common.legend = TRUE,ncol=1,heights = hei.vec)
  ggsave(handl_OneDrive(paste0("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Multivariate/","Year prediction_All.tiff")),
         width = 7,height = 6.5,compression = "lzw")
  
    #by data set  
  for(l in 1:length(p.list))
  {
    NM=names(Data.list)[l]
    if(NM=="Logbook") NM="TDGDLF"
    if(NM=="Logbook.north") NM="NSF"
  
    print(p.list[[l]]+ labs(title = NULL) )
    ggsave(handl_OneDrive(paste0("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Multivariate/Year prediction_",NM,".tiff")),
           width = 6,height = 4,compression = "lzw")
    
  }

  
  
  #Display lat and long 
    #data sets combined
  do.this=FALSE
  if(do.this)
  {
    p.list=vector('list',length(Store.out.multi))
    names(p.list)=names(Store.out.multi)
    for(l in 1:length(p.list))
    {
      NM=names(Data.list)[l]
      if(NM=="Logbook") NM="TDGDLF"
      if(NM=="Logbook.north") NM="NSF"
      #if(NM=="Observer") NM="Observer (TDGDLF)"
      
      ddd=Store.out.multi[[l]]$Preds.lat.lon%>%
        left_join(All.sp,by=c('scientific_name'='SCIENTIFIC_NAME'))%>%
        mutate(color=ifelse(Group=='Teleost','dodgerblue4',
                            ifelse(Group=='Elasmobranch','firebrick3',NA)),
               color=ifelse(is.na(color),'black',color))
      
      LVLs=ddd%>%distinct(scientific_name,color,Group)%>%arrange(desc(Group),scientific_name)%>%pull(scientific_name)
      ddd=ddd%>%
        filter(scientific_name%in%c('Carcharhinus obscurus','Heterodontus portusjacksoni',
                                    'Mustelus antarcticus','Carcharhinus plumbeus','Furgaleus macki',
                                    'Nemadactylus valenciennesi'))
      a=ddd%>%distinct(scientific_name,color,Group)%>%arrange(desc(Group),factor(scientific_name,levels=LVLs))
      
      ddd=ddd%>%    
        rename(Proportion=Median)%>%
        mutate(Proportion=ifelse(Proportion<0,0,Proportion),
               scientific_name=factor(scientific_name,levels=LVLs))
      ylAb=''
      if(l==length(p.list)) ylAb='Latitude'
      p.list[[l]]=ddd%>%
        ggplot(aes(longitude,latitude,size=Proportion))+
        geom_point()+
        facet_wrap(~scientific_name,ncol=1)+
        ylab('')+xlab('Longitude')+
        ggtitle(NM)+
        theme_PA(axs.t.siz=10,Ttl.siz=15)+
        theme(legend.position = 'top')
      
    }
    ggarrange(plotlist=p.list, common.legend = TRUE,ncol=2)
    ggsave(handl_OneDrive(paste0("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Multivariate/","Spatial prediction_All.tiff")),
           width = 6.5,height = 8,compression = "lzw")
    
  }
  toc()
}



# 6. Cluster analysis of Data.daily -----------------------------------------------------------------------
#note: variable is catch as proportion
Terms=c('SHEET_NO','BOAT','YEAR','NamE','MONTH','BLOCK','LATITUDE','LONGITUDE','zone')
Modelled.var='proportion' # 'proportion'  'cpue'
zonEs=sort(unique(Data.daily$zone))
Clus=vector('list',length(zonEs))
names(Clus)=zonEs
p1.list=p2.list=Clus

#Gillnet and Longline
daily.list=list(GN=Data.daily,LL=Data.daily.LL)
tic()
for(l in 1:length(daily.list))
{
  for(z in 1:length(zonEs))
  {
    print(paste('-----Cluster analysis for Method----',names(daily.list)[l], 'zone------',zonEs[z]))
    Clus[[z]]=Cluster.fn(d=daily.list[[l]]%>%
                           filter(zone==zonEs[z])%>%
                           mutate(NamE=SNAME)%>%
                           group_by_at(Terms)%>%
                           summarise(var=sum(LIVEWT.c,na.rm=T))%>%
                           #summarise(var=sum(LIVEWT.c/EFFORT,na.rm=T))%>%
                           ungroup(), 
                         Terms=tolower(Terms),
                         n.recent.years=5,
                         percent.ktch.explained=0.80,
                         var=Modelled.var)
    
    #display main species by cluster
    p1.list[[z]]=Clus[[z]]$dt%>%
      gather(key = species, value = var, -c(tolower(subset(Terms,!Terms=='NamE')),'cluster'))%>%
      mutate(cluster=as.character(cluster))%>%
      group_by(cluster,species)%>%
      summarise(Mean=round(mean(var),2))%>%
      mutate(species=ifelse(species=='Wobbegong','Wobbie',species))%>%
      ggplot(aes(species,Mean,fill=cluster))+
      geom_bar(stat='identity')+
      geom_text(aes(label=Mean),size = 2, position=position_dodge(width=0.9), vjust=0.5)+
      facet_wrap(~cluster,ncol=1)+ylab('Proportion')+xlab('')+
      scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
      theme(axis.text=element_text(size=8))
    
    #display boxplot species raw data by cluster
    p2.list[[z]]=Clus[[z]]$dt%>%
      gather(key = species, value = var, -c(tolower(subset(Terms,!Terms=='NamE')),'cluster'))%>%
      mutate(cluster=as.character(cluster))%>%
      ggplot(aes(species,var,fill=cluster))+
      geom_violin()+coord_flip()+
      geom_jitter(aes(color=cluster),height = 0, width = 0.1,alpha=0.2)+
      xlab('')+ylab('Proportion')
    
    #combine all figures
    ggarrange(plotlist=list(Clus[[z]]$num.clus,Clus[[z]]$cluster,p2.list[[z]],p1.list[[z]]),
              common.legend = TRUE,ncol=2,nrow=2)
    ggsave(handl_OneDrive(paste0("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Multivariate/Cluster_metier_",
                                 names(daily.list)[l],'_',zonEs[z],".tiff")),
           width = 7,height = 6,compression = "lzw")
  }
}
toc()
rm(daily.list)


# 7. Display.catch.effort -----------------------------------------------------------------------
if(display.catch.effort)
{
  library(rgdal)
  Dat.repository=handl_OneDrive('Analyses/Data_outs/')  #locations where all data are stored
  fn.in=function(NM) fread(paste(Dat.repository,NM,sep=""),data.table=FALSE)
  #TDGDLF
  Effort.monthly=fn.in(NM='Annual.total.eff.days.csv')          #all years
  Effort.monthly.zone=fn.in(NM='Annual.zone.eff.days.csv') 
  
  #NSF
  Effort.monthly.north=fn.in(NM='Annual.total.eff_NSF.csv') 
  
  Southern=fn.in(NM='Data.monthly.csv')%>%
    filter(!Shark.fishery=='non.shark.fishery' & SPECIES<99999)%>%
    filter(!is.na(BLOCKX))%>%
    group_by(FINYEAR,FisheryZone,BLOCKX,SPECIES)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
  Northern=fn.in(NM='Data.monthly.north.csv')%>%
    filter(!Shark.fishery=='non.shark.fishery' & SPECIES<99999)%>%
    filter(!is.na(BLOCKX))%>%
    group_by(FINYEAR,FisheryZone,BLOCKX,SPECIES)%>%
    summarise(LIVEWT.c=sum(LIVEWT.c,na.rm=T))
  if(use.NSF)
  {
    Spatio.temp.dat=rbind(Southern%>%
                            mutate(Fishery='Southern'),
                          Northern%>%
                            mutate(Fishery='Northern'))
  }else
  {
    Spatio.temp.dat=Southern%>%mutate(Fishery='Southern')
  }
  
  
  
  #spatial effort
  Relevant.yrs=use.dis.yrs=unique(Spatio.temp.dat$FINYEAR)
  Effort.monthly_blocks=fn.in("Effort.monthly.csv")
  Effort.daily_blocks=fn.in("Effort.daily.csv")
  Effort.monthly.north_blocks=fn.in("Effort.monthly.north.csv")
  Effort.daily.north_blocks=fn.in("Effort.daily.north.csv") 
  
  Effort.monthly_blocks=Effort.monthly_blocks%>%filter(FINYEAR%in%use.dis.yrs) 
  Effort.monthly.north_blocks=Effort.monthly.north_blocks%>%filter(FINYEAR%in%use.dis.yrs) 
  Effort.daily_blocks=Effort.daily_blocks%>%filter(finyear%in%use.dis.yrs) 
  Effort.daily.north_blocks=Effort.daily.north_blocks%>%filter(finyear%in%use.dis.yrs) 
  
  FINYrS=sort(unique(c(as.character(unique(Effort.monthly_blocks$FINYEAR))),sort(as.character(unique(Effort.daily_blocks$finyear)))))
  grouping=8
  FINYrS.gp=seq(1,length(FINYrS),by=grouping)
  FINYrS.gped=vector('list',length(FINYrS.gp))
  for(f in 1:length(FINYrS.gped))
  {
    if(f==length(FINYrS.gped))
    {
      FINYrS.gped[[f]]=FINYrS[FINYrS.gp[f]:length(FINYrS)]
      if(length(FINYrS.gped[[f]])==1)
      {
        names(FINYrS.gped)[f]=FINYrS.gped[[f]][1]
      }else
      {
        names(FINYrS.gped)[f]=paste(FINYrS.gped[[f]][1],"to",FINYrS.gped[[f]][length(FINYrS.gped[[f]])])
      }
      
    }else
    {
      FINYrS.gped[[f]]=FINYrS[FINYrS.gp[f]:(FINYrS.gp[f+1]-1)]
      names(FINYrS.gped)[f]=paste(FINYrS.gped[[f]][1],"to",FINYrS.gped[[f]][length(FINYrS.gped[[f]])])
    }
  }
  suppressWarnings({FINYrS.gped=do.call(cbind,FINYrS.gped)%>%
    data.frame%>%
    gather('Rango','FINYEAR')%>%
    mutate(Rango=substr(Rango,2,50),
           Rango=sub(".to.", " to ", Rango),
           Rango=str_replace_all(Rango, c("\\."), "-"))%>%
    distinct(Rango,FINYEAR)})
  daily.years=unique(Effort.daily_blocks$finyear)
  daily.years=subset(daily.years,!daily.years=="2005-06")
  
  Spatial.effort_monthly=Effort.monthly_blocks%>%
    filter(!Shark.fishery=='non.shark.fishery' &
             FINYEAR%in%Relevant.yrs&
             !FINYEAR%in%daily.years)%>%
    mutate(Effort=Km.Gillnet.Days.c,
           LAT1=-floor(abs(LAT)),
           LONG1=floor(LONG))%>%
    filter(!is.na(Effort) | !is.na(LAT))%>%
    filter(METHOD=='GN')%>%
    group_by(Same.return,zone,LAT1,LONG1,FINYEAR)%>%
    summarise(Effort=max(Effort))%>%
    ungroup()%>%
    group_by(LAT1,LONG1,zone,FINYEAR)%>%
    summarise(Effort=sum(Effort))
  Spatial.effort_daily=Effort.daily_blocks%>%
    filter(!Shark.fishery=='non.shark.fishery' &
             finyear%in%Relevant.yrs )%>%
    mutate(Effort=Km.Gillnet.Days.c,
           LAT1=-floor(abs(LAT)),
           LONG1=floor(LONG))%>%
    filter(!is.na(Effort) | !is.na(LAT))%>%
    filter(method=='GN')%>%
    group_by(Same.return.SNo,zone,LAT1,LONG1,finyear)%>%
    summarise(Effort=max(Effort))%>%
    ungroup()%>%
    rename(FINYEAR=finyear)%>%
    group_by(LAT1,LONG1,zone,FINYEAR)%>%
    summarise(Effort=sum(Effort))
  
  Spatial.effort_monthly.north=Effort.monthly.north_blocks%>%
    filter(FINYEAR%in%Relevant.yrs &
             !FINYEAR%in%unique(Effort.daily.north_blocks$finyear))%>%
    mutate(Effort=hook.days,
           LAT1=-floor(abs(as.numeric(substr(BLOCKX,1,2)))),
           LONG1=floor(100+as.numeric(substr(BLOCKX,3,4))),
           LAT1=ifelse(BLOCKX%in%c(96021),-25,LAT1),
           LONG1=ifelse(BLOCKX%in%c(96021),113,LONG1))%>%
    filter(!is.na(Effort) | !is.na(LAT1))%>%
    filter(METHOD=='LL')%>%
    group_by(Same.return,LAT1,LONG1,zone,FINYEAR)%>%
    summarise(Effort=max(Effort))%>%
    ungroup()%>%
    group_by(LAT1,LONG1,zone,FINYEAR)%>%
    summarise(Effort=sum(Effort))
  Spatial.effort_daily.north=Effort.daily.north_blocks%>%
    rename(METHOD=method,
           FINYEAR=finyear,
           BLOCKX=blockx)%>%
    filter(FINYEAR%in%Relevant.yrs)%>%
    mutate(Effort=hook.days,
           LAT1=-floor(abs(as.numeric(substr(BLOCKX,1,2)))),
           LONG1=floor(100+as.numeric(substr(BLOCKX,3,4))),
           LAT1=ifelse(BLOCKX%in%c(96021),-25,LAT1),
           LONG1=ifelse(BLOCKX%in%c(96021),113,LONG1))%>%
    filter(!is.na(Effort) | !is.na(LAT1))%>%
    filter(METHOD=='LL')%>%
    group_by(Same.return.SNo,LAT1,LONG1,zone,FINYEAR)%>%
    summarise(Effort=max(Effort))%>%
    ungroup()%>%
    group_by(LAT1,LONG1,zone,FINYEAR)%>%
    summarise(Effort=sum(Effort))
  
  if(use.NSF)
  {
    Spatial.effort=rbind(Spatial.effort_monthly%>%mutate(Fishery='Southern'),
                         Spatial.effort_daily%>%mutate(Fishery='Southern'),
                         Spatial.effort_monthly.north%>%mutate(Fishery='Northern'),
                         Spatial.effort_daily.north%>%mutate(Fishery='Northern'))%>%
      filter(!is.na(Effort))%>%
      left_join(FINYrS.gped,by='FINYEAR')
  }else
  {
    Spatial.effort=rbind(Spatial.effort_monthly%>%mutate(Fishery='Southern'),
                         Spatial.effort_daily%>%mutate(Fishery='Southern'))%>%
      filter(!is.na(Effort))%>%
      left_join(FINYrS.gped,by='FINYEAR')
  }
  
  WAcoast<-read.table(handl_OneDrive("Data/Mapping/WAcoastPointsNew.txt"), header=T)
  WAcoast=WAcoast%>%mutate(Latitude=abs(Latitude))
  Fishing.zones="Data/Mapping/Shark_shape_files/"
  SDGDLL_zone1=readOGR(handl_OneDrive(paste0(Fishing.zones,"SDGDLL_zone1.shp")), layer="SDGDLL_zone1") 
  SDGDLL_zone2=readOGR(handl_OneDrive(paste0(Fishing.zones,"SDGDLL_zone2.shp")), layer="SDGDLL_zone2") 
  WCDGDLL=readOGR(handl_OneDrive(paste0(Fishing.zones,"WCDGDLL.shp")), layer="WCDGDLL") 
  
  
  p1=Spatial.effort%>%
    group_by(LAT1,LONG1,Rango,Fishery)%>%
    summarise(Effort=sum(Effort))%>%
    group_by(Fishery)%>%
    mutate(Rel.effort=Effort/max(Effort))%>%
    ungroup()%>%
    mutate(LAT1=abs(LAT1))%>%
    ggplot(aes(LONG1,LAT1))+
    geom_polygon(data=WAcoast,aes(x=Longitude,y=Latitude),fill="grey90")+
    geom_raster(aes(fill = Rel.effort))+
    facet_wrap(~Rango,ncol=2)+
    scale_fill_gradientn(colours=c('ivory2','gold',"red2","darkred"),
                         guide = guide_legend(direction="vertical"))+
    theme_PA(leg.siz=7,axs.t.siz=10,axs.T.siz=14,Sbt.siz=14,strx.siz=11)+
    theme(
      legend.title=element_blank(),
      legend.key.size = unit(.65, 'cm'),
      legend.position = c(0.90, 0.9),
      plot.title.position = "plot")+
    ylab('')+xlab('')+
    scale_y_reverse()+xlim(113.5,130)+ggtitle('TDGDLF')
  
  #Observed shots
  Min.x=109
  zone.west=WCDGDLL@polygons[[1]]@Polygons[[1]]@coords%>%
    data.frame%>%
    mutate(X2=abs(X2))%>%
    filter(X1>Min.x)
  zone.zone1=SDGDLL_zone1@polygons[[1]]@Polygons[[1]]@coords%>%
    data.frame%>%
    mutate(X2=abs(X2))%>%
    filter(X1>Min.x)
  zone.zone2=SDGDLL_zone2@polygons[[1]]@Polygons[[1]]@coords%>%
    data.frame%>%
    mutate(X2=abs(X2))%>%
    filter(X1>Min.x)
  
  N.shots=Data.list$Observer%>%
    distinct(SHEET_NO,.keep_all = T)%>%
    left_join(FINYrS.gped%>%mutate(YEAR=as.numeric(substr(FINYEAR,1,4))),by='YEAR')%>%
    group_by(Rango,ZONE,LATITUDE,LONGITUDE)%>%
    tally()
  
  p2=N.shots%>% 
    mutate(Longitude=LONGITUDE,
           Latitude=abs(LATITUDE))%>%
    ggplot()+
    geom_polygon(data=WAcoast,aes(x=Longitude,y=Latitude),fill="grey90")+
    geom_polygon(data=zone.west,aes(x=X1,y=X2),fill="#F8766D")+
    geom_polygon(data=zone.zone1,aes(x=X1,y=X2),fill="#00BA38")+
    geom_polygon(data=zone.zone2,aes(x=X1,y=X2),fill="#619CFF")+
    geom_point(aes(Longitude,Latitude,size=n),shape=21, fill="grey30", color="black",alpha=0.4)+
    facet_wrap(~Rango,ncol=1)+
    theme_PA(leg.siz=7,axs.t.siz=10,axs.T.siz=14,Sbt.siz=14,strx.siz=11)+
    theme(legend.title=element_blank(),
          legend.key.size = unit(.45, 'cm'),
          legend.position = c(0.83, 0.9),
          #legend.position = c(0.23, 0.965),
          #legend.direction='horizontal',
          legend.key = element_rect(fill = "transparent"),
          plot.title.position = "plot")+
    ylab('')+xlab('')+
    scale_y_reverse()+xlim(Min.x,130)+ggtitle('Observer')
  
  #Temporal catch and effort
  Catch.agg=Spatio.temp.dat%>%
    mutate(Year=as.numeric(substr(FINYEAR,1,4)))%>%
    group_by(Year,FisheryZone)%>%
    summarise(Tons=sum(LIVEWT.c,na.rm=T)/1000)    
  
  Eff.agg=Effort.monthly.zone%>%
    gather(FisheryZone,Effort,-FINYEAR)%>%
    mutate(Year=as.numeric(substr(FINYEAR,1,4)))
  
  p3=Catch.agg%>%
    ggplot(aes(Year,Tons,fill=FisheryZone))+
    geom_bar(stat='identity')+
    #geom_point(size=2)+geom_line(linewidth=1.1,alpha=.4)+ #linetype ='dotted' 
    ylab('Total landings (tonnes)')+xlab('Financial year')+
    theme_PA(leg.siz=12,axs.t.siz=10,axs.T.siz=14,str.siz=10)+
    theme(legend.title = element_blank())
  
  CPUE=Catch.agg%>%
    left_join(Eff.agg,by=c('Year','FisheryZone'))%>%
    mutate(CPUE=Tons/Effort)
 
  p4=Eff.agg%>%
    ggplot(aes(Year,Effort,fill=FisheryZone))+
    geom_bar(stat='identity')+
    #geom_point(size=2)+geom_line(linewidth=1.1,alpha=.4)+
    ylab('Total effort (1000 km gn days)')+xlab('Financial year')+
    theme_PA(leg.siz=12,axs.t.siz=10,axs.T.siz=14,str.siz=10)+
    theme(legend.title = element_blank())

  pp1=ggarrange(plotlist=list(p2+ rremove("ylab") + rremove("xlab"),p1+ rremove("ylab") + rremove("xlab")),
                common.legend = FALSE,ncol=2,widths = c(.55,1))
  pp=ggarrange(plotlist=list(p4,p3), common.legend = TRUE,ncol=2)
  ggarrange(plotlist=list(annotate_figure(pp1, left = textGrob(expression('Latitude ('*~degree*S*')'), rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                                          bottom = textGrob(expression('Longitude ('*~degree*E*')'), gp = gpar(cex = 1.3))),
                          pp), ncol=1,heights = c(1.3,.7))
  ggsave(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Univariate/Catch and effort.tiff"),
         width = 7,height = 11.5,compression = "lzw")
  
  #TDGDLF only
  p6=N.shots%>% 
    mutate(Longitude=LONGITUDE,
           Latitude=abs(LATITUDE))%>%
    ggplot()+
    geom_polygon(data=WAcoast,aes(x=Longitude,y=Latitude),fill="grey90")+
    geom_polygon(data=zone.west,aes(x=X1,y=X2),fill="#F8766D")+
    geom_polygon(data=zone.zone1,aes(x=X1,y=X2),fill="#00BA38")+
    geom_polygon(data=zone.zone2,aes(x=X1,y=X2),fill="#619CFF")+
    theme_PA(leg.siz=7,axs.t.siz=10,axs.T.siz=14,Sbt.siz=14,strx.siz=11)+
    theme(plot.title.position = "plot")+
    ylab(expression('Latitude ('*~degree*S*')'))+
    xlab(expression('Longitude ('*~degree*E*')'))+
    scale_y_reverse()+xlim(Min.x,130)+
    ggtitle('TDGDLF')
  p3.published=p3+
    geom_line(data=CPUE,aes(Year,CPUE*10),linewidth=0.5) + #linetype='dotted'
    scale_y_continuous(sec.axis = sec_axis(~./10, name="CPUE"))+
    geom_point(data=CPUE,aes(Year,CPUE*10,fill=FisheryZone),shape=21,size=1,show.legend = FALSE)
  
  pp1=ggarrange(plotlist=list(p4+theme(axis.text = element_text(size = 9))+rremove("xlab"),
                              p3.published+theme(axis.text = element_text(size = 9))+rremove("xlab")),
                common.legend = TRUE,ncol=2,widths = c(0.75,1))
  pp1=annotate_figure(pp1, bottom = textGrob("Financial year", gp = gpar(cex = 1.3)))
  
  
  pp2=ggarrange(plotlist=list(p6,pp1),common.legend = TRUE,ncol=2,widths = c(0.4,1))
  
  ggarrange(plotlist=list(pp2,
                          p1+
                            facet_wrap(~Rango,ncol=3)+theme(legend.position = c(0.90, 0.825))+ggtitle('')+
                             ylab(expression('Latitude ('*~degree*S*')'))+
                             xlab(expression('Longitude ('*~degree*E*')'))),
                    ncol=1, heights = c(.7,1.3))
  
  ggsave(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Univariate/Catch and effort_logbook_with CPUE.tiff"),
         width = 7,height = 8,compression = "lzw")
  
  CPUE%>%
    ggplot(aes(Year,CPUE,color=FisheryZone))+
    geom_line(linewidth=1.1) +
    ylab('CPUE (tonnes/1000 km gn days)')+xlab('Financial year')+
    theme_PA(leg.siz=12,axs.t.siz=10,axs.T.siz=14,str.siz=10)+
    theme(legend.title = element_blank())
  ggsave(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/Univariate/CPUE_logbook.tiff"),
         width = 6,height = 6,compression = "lzw")
   
}


