#--- SCRIPT FOR EXTRACTING AND PUTTING DATA IN RIGHT FORMAT FOR ECOSYSTEM ANALYSIS OF CATCH ---#

#----DATA SECTION----


#Commercial catch and effort TDGDLF
Data.monthly=read.csv("C:/Matias/Analyses/Data_outs/Data.monthly.csv",stringsAsFactors=F)
Data.daily=read.csv("C:/Matias/Analyses/Data_outs/Data.daily.csv",stringsAsFactors=F)
Effort.monthly=read.csv("C:/Matias/Analyses/Data_outs/Effort.monthly.csv",stringsAsFactors=F)
Effort.daily=read.csv("C:/Matias/Analyses/Data_outs/Effort.daily.csv",stringsAsFactors=F)



#----PROCEDURE SECTION----

#Use only TDGDLF GN records
Data.monthly=subset(Data.monthly,METHOD=="GN" & Estuary=="NO" & NETLEN.c>=100 & LAT<=(-26),
                    select=c(Same.return,FINYEAR,MONTH,VESSEL,METHOD,
                             BLOCKX,SPECIES,SNAME,YEAR.c,LIVEWT.c,LAT,LONG,TYPE.DATA,
                             Bioregion,zone,NETLEN.c))

SHOTS_c=unique(Data.monthly$Same.return)

#Select corresponding effort (km gn days)

daily.yrs=unique(Effort.daily$finyear)

  #monthly
Effort.monthly=subset(Effort.monthly,!FINYEAR%in%daily.yrs)
Effort.monthly=subset(Effort.monthly,Same.return%in%SHOTS_c)
Effort.monthly=aggregate(Km.Gillnet.Days.c~Same.return,data=Effort.monthly,max,na.rm=T) #remove duplicates

#aggregate daily to monthly
Use.Date="YES"    #Rory's approach (aggregating by DATE)
#Use.Date="NO"     # aggregating by SNo and DSNo 
if(Use.Date=="NO")  Effort.daily=aggregate(Km.Gillnet.Days.c~ID+Same.return,data=Effort.daily,max,na.rm=T)
if(Use.Date=="YES") Effort.daily=aggregate(Km.Gillnet.Days.c~date+Same.return,data=Effort.daily,max,na.rm=T) 
Effort.daily=aggregate(Km.Gillnet.Days.c~Same.return,data=Effort.daily,sum,na.rm=T)
Effort.daily=subset(Effort.daily,Same.return%in%SHOTS_c)

Effort.monthly$Same.return=as.character(Effort.monthly$Same.return)
Effort.daily$Same.return=as.character(Effort.daily$Same.return)
Effort=rbind(Effort.monthly,Effort.daily)


#Aggregate catch for daily records
Daily.agg=subset(Data.monthly,FINYEAR%in%daily.yrs)
Data.monthly=subset(Data.monthly,!FINYEAR%in%daily.yrs)

Daily.agg$LAT=-as.numeric(substr(Daily.agg$BLOCKX,1,2))
Daily.agg$LONG=100+as.numeric(substr(Daily.agg$BLOCKX,3,4))
Daily.agg=aggregate(LIVEWT.c~Same.return+FINYEAR+MONTH+VESSEL+METHOD+BLOCKX+SPECIES+
                                SNAME+YEAR.c+LAT+LONG+TYPE.DATA,Daily.agg,sum)

Data.monthly=Data.monthly[,match(names(Daily.agg),names(Data.monthly))]
Data.monthly=rbind(Data.monthly,Daily.agg)



#Attach effort
Data.monthly$Same.return=as.character(Data.monthly$Same.return)
Data.monthly=merge(Data.monthly,Effort,by="Same.return",all.x=T)


#Change variable names to match observers data

Data.monthly$SHEET_NO=Data.monthly$Same.return
Data.monthly$YEAR=Data.monthly$FINYEAR
Data.monthly$BLOCK=Data.monthly$BLOCKX
Data.monthly$BOAT=Data.monthly$VESSEL
Data.monthly$SKIPPER=NA
  
Data.monthly$EFFORT=Data.monthly$Km.Gillnet.Days.c
Data.monthly$LATITUDE=Data.monthly$LAT
Data.monthly$LONGITUDE=Data.monthly$LONG


Data.monthly$BIOREGION=as.character(with(Data.monthly,ifelse(LONG>=115.5 & LONG<=129 & LAT<=(-26),"SC", 
               ifelse(LONG<115.5 & LAT<=(-27),"WC",
               ifelse(LONG<=114.834 & LAT>(-27),"Gascoyne",
               ifelse(LONG>114.834 & LAT>=(-27) & LONG<=129,"NC",NA))))))
    
Data.monthly$BIOREGION=with(Data.monthly,
               ifelse(BIOREGION=="SC"& LAT>(-34) & LONG <115.91 ,"WC",BIOREGION))
    
    
Data.monthly$ZONE=as.character(with(Data.monthly,ifelse(LONG>=116.5 & LAT<=(-26),"Zone2",
                ifelse(LONG<116.5 & LAT<=(-33),"Zone1",
                ifelse(LAT>(-33) & LAT<=(-26) & LONG<116.5,"West",
                ifelse(LAT>(-26) & LONG<114,"Closed",
                ifelse(LAT>(-23) & LONG>=114 & LONG<123.75,"North",
                ifelse(LAT>(-23) & LONG>=123.75,"Joint",NA))))))))
