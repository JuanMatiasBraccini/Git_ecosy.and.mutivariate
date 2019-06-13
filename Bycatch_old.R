
#SCRIPT FOR ANALYSING BYCATCH PATTERNS IN THE SHARK GILLNET FISHERIES OF AUSTRALIA


#INDEX
# ----- DATA SECTION
#         1. Bring in WA shark gillnet data
#         2. Bring in MAFFRI shark gillnet data
#         3. Bring in WA Species names

# ----- PROCEDURE SECTION
#         1. Manipulated WA shark gillnet data
#         2. Manipulated MAFFRI shark gillnet data
#         3. Start analysis by data set




require(lubridate)  #for dates manipulation
library(PBSmapping)
data(worldLLhigh)
library(plotrix)
library(lme4) #mixed models
library(vegan)



#DATA SECTION

setwd("C:/Matias/Analyses/Bycatch/Shark-bycatch")


# 1. Bring in WA shark gillnet data

DATA=read.csv("WA.csv",stringsAsFactors=FALSE)


# 2. Bring in MAFFRI shark gillnet data
DATA.MAFFRI=read.csv("MAFFRI_data.csv") 
DATA_TEPS.MAFFRI=read.csv("MAFFRI_TEPS_data.csv")




# 3. Bring in WA Species names
SPECIES.names=read.csv("Species.code.csv")


todas.las.especies=read.csv("dudas/todas.las.especies.csv")
NA.especies=subset(todas.las.especies,is.na(COMMON_NAME))

Net.len.missing=read.csv("dudas/NET_LENGTH_NA.csv")

#PROCEDURE SECTION


#1. Manipulated WA shark gillnet data

names(DATA)[match("NO HOOKS",names(DATA))]="N.hooks"
names(DATA)[match("SOAK TIME",names(DATA))]="SOAK_TIME"


#Select commercial vessels              
Research.vess=c("HAM","HOU","RV Gannet","RV GANNET","RV BREAKSEA",
                "NATT","NAT","FLIN","RV SNIPE 2")
DATA=subset(DATA,!BOAT%in%Research.vess)

#Fix lats
names(DATA)[match(c("MID.LAT","MID.LONG"),names(DATA))]=c('Mid.Lat','Mid.Long')
DATA$Mid.Lat=-DATA$Mid.Lat
DATA$END1LATD=-DATA$END1LATD
DATA$END2LATD=-DATA$END2LATD


#Extract month,year, day, hour
DATA$date=as.Date(DATA$DATE,format="%Y-%m-%d")
DATA$Day=mday(DATA$date)
DATA$Month=month(DATA$date)
DATA$year=year(DATA$date)

DATA$Set.time=strftime(DATA$START_SET, format='%H:%M')
DATA$Haul.time=strftime(DATA$START_HAUL, format='%H:%M')

#fix dodgy longitude                                      
DATA$END1LNGD=with(DATA,ifelse(SHEET_NO=="R00890",113.96,END1LNGD))
DATA$Mid.Long=with(DATA,ifelse(SHEET_NO=="R00890",113.9632,Mid.Long))
DATA$Mid.Long=with(DATA,ifelse(SHEET_NO=="N00099",113.4272,Mid.Long))
DATA$Mid.Long=with(DATA,ifelse(SHEET_NO=="N00597",113.2404,Mid.Long))
DATA$Mid.Long=with(DATA,ifelse(SHEET_NO=="N00558",113.2055,Mid.Long))

#fix dodgy latitude
DATA$Mid.Lat=with(DATA,ifelse(SHEET_NO=="N00401",-20.90,Mid.Lat))

#remove land blocks
DATA=subset(DATA,!BLOCK%in%c(2415,2315,2519))  

DATA$Lat.round = ceiling(DATA$Mid.Lat) -0.5
DATA$Long.round = floor(DATA$Mid.Long) +0.5

#fix method issues
DATA$Method=with(DATA,ifelse(is.na(Method) & MESH_SIZE%in%c("6","6.5","7","8"),"GN",
                             ifelse(is.na(Method)& N.hooks>25,"LL",
                                    ifelse(is.na(Method)& BOAT=="E35","GN",Method))))

DATA$MESH_SIZE=with(DATA,ifelse(is.na(MESH_SIZE) & BOAT=="E35" & year>=2011,"6.5",MESH_SIZE))
DATA$Method=with(DATA,ifelse(is.na(Method) & BOAT%in%c("B67","F244","F517"),"GN",Method))

DATA$MESH_SIZE=with(DATA,ifelse(is.na(MESH_SIZE) & BOAT=="B67" & year>1997,"7",
                                ifelse(is.na(MESH_SIZE) & BOAT=="F244" & year>2004,"6.5",            
                                       ifelse(is.na(MESH_SIZE) & BOAT=="F517" & year==2002,"6.5",
                                              ifelse(is.na(MESH_SIZE) & BOAT=="F517" & year>2002,"7",
                                                     MESH_SIZE)))))
#Get mid point of block
DATA$LAT=-as.numeric(substr(DATA$BLOCK,1,2))
DATA$LONG=100+as.numeric(substr(DATA$BLOCK,3,4))

#select commercial gillnet
DATA=subset(DATA,Method=="GN" & MESH_SIZE%in%c("6","6.5","605","6.587","7"))


#Select relevant variables
DATA=subset(DATA,select=c(SHEET_NO,SPECIES,TL,FL,PL,SEX,
            AVE.SET.TIME,START_SET,AVE.HAUL.TIME,END_SET,SOAK.TIME,
            Month,year,BOAT,BOTDEPTH,MESH_SIZE,NET_LENGTH,Mid.Lat,Mid.Long
  ))


#Add Number caught
DATA$Number=1
                            


#Fix sex
DATA$SEX=as.character(DATA$SEX)
DATA$SEX=with(DATA,ifelse(SEX=="f","F",ifelse(SEX=="m","M",
      ifelse(SEX%in%c(",","?","n","N","p","P","Y"),"U",SEX))))

DATA$FL=with(DATA,ifelse(FL<10,NA,FL))
DATA$FL=with(DATA,ifelse(SPECIES=="BW" & FL<60,60,ifelse(SPECIES=="WH" & FL<25,25,
             ifelse(SPECIES=="SD" & FL<20,20,FL))))
#Range analysis
DATA=subset(DATA,!(SPECIES=="GM" & Mid.Lat>(-26)))


#Add common name
names(SPECIES.names)[1]="SPECIES"
DATA=merge(DATA,SPECIES.names,by="SPECIES",all.x=T)
DATA=DATA[,-match("SPECIES",names(DATA))]

#2. Manipulated MAFFRI shark gillnet data
DATA.MAFFRI$SHEET_NO=with(DATA.MAFFRI,paste(Cruise,Station))

DATA.MAFFRI$Mid.Lat=-DATA.MAFFRI$StartLat/100
DATA.MAFFRI$Mid.Long=DATA.MAFFRI$StartLong/100

DATA_TEPS.MAFFRI$Mid.Lat=-DATA_TEPS.MAFFRI$StartLat/100
DATA_TEPS.MAFFRI$Mid.Long=DATA_TEPS.MAFFRI$StartLong/100

DATA.MAFFRI=DATA.MAFFRI[,-match(c("StartLat","StartLong","Dummy"),names(DATA.MAFFRI))]
DATA_TEPS.MAFFRI=DATA_TEPS.MAFFRI[,-match(c("StartLat","StartLong","Dummy"),names(DATA_TEPS.MAFFRI))]


# 3. Start analysis by data set

#3.1 WA Fisheries
DATA

#3.2 MAFFRI
DATA.MAFFRI
 

#TEPS analysis
DATA_TEPS.MAFFRI