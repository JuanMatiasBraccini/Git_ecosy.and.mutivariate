#--- SCRIPT FOR EXTRACTING AND PUTTING DATA IN RIGHT FORMAT FOR ECOSYSTEM ANALYSIS OF CATCH ---#
library(tidyverse)


#----DATA SECTION----
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

#Commercial catch and effort TDGDLF
Data.monthly=read.csv(handl_OneDrive("Analyses/Data_outs/Data.monthly.csv"),stringsAsFactors=F)
Data.daily=read.csv(handl_OneDrive("Analyses/Data_outs/Data.daily.csv"),stringsAsFactors=F)
Effort.monthly=read.csv(handl_OneDrive("Analyses/Data_outs/Effort.monthly.csv"),stringsAsFactors=F)
Effort.daily=read.csv(handl_OneDrive("Analyses/Data_outs/Effort.daily.csv"),stringsAsFactors=F)

#Commercial catch and effort NSF
Data.monthly.north=read.csv(handl_OneDrive("Analyses/Data_outs/Data.monthly.north.csv"),stringsAsFactors=F)
Data.daily.north=read.csv(handl_OneDrive("Analyses/Data_outs/Data.daily.north.csv"),stringsAsFactors=F)
Effort.monthly.north=read.csv(handl_OneDrive("Analyses/Data_outs/Effort.monthly.north.csv"),stringsAsFactors=F)
Effort.daily.north=read.csv(handl_OneDrive("Analyses/Data_outs/Effort.daily.north.csv"),stringsAsFactors=F)


#----NSF----
Data.monthly.north=Data.monthly.north%>%
                    filter(METHOD=="LL" & Estuary=="NO" & LAT>(-26) &
                           !FINYEAR%in%unique(Data.daily.north$FINYEAR))%>%
                    select(Same.return,FINYEAR,MONTH,VESSEL,METHOD,
                           BLOCKX,SPECIES,SNAME,YEAR.c,LIVEWT.c,LAT,LONG,TYPE.DATA,
                           zone)%>%mutate(NETLEN.c=NA)
Data.daily.north=Data.daily.north%>%
                    filter(METHOD=="LL" & Estuary=="NO" & LAT>(-26))%>%
                    select(Same.return.SNo,FINYEAR,MONTH,VESSEL,METHOD,
                           BLOCKX,SPECIES,SNAME,YEAR.c,LIVEWT.c,LAT,LONG,TYPE.DATA,
                           zone)%>%mutate(NETLEN=NA)

SHOTS_c=unique(Data.monthly.north$Same.return)

#Select corresponding effort (km gn days)
daily.yrs=unique(Effort.daily.north$finyear)

#monthly
Effort.monthly.north=Effort.monthly.north%>%
                        filter(!FINYEAR%in%daily.yrs)%>%
                        filter(Same.return%in%SHOTS_c)%>%
                        group_by(Same.return)%>%
                        summarise(hook.hours=max(hook.hours,na.rm=T))

#daily
Effort.daily.north=Effort.daily.north%>%
                      group_by(Same.return.SNo)%>%
                      filter(!is.na(hook.hours))%>%
                      summarise(hook.hours=max(hook.hours,na.rm=T))

##
#Attach effort
Data.monthly.north=Data.monthly.north%>%
                    left_join(Effort.monthly.north,by="Same.return")%>%
                    filter(!is.na(hook.hours))%>%
                    filter( LIVEWT.c>0)
Data.daily.north=Data.daily.north%>%
                    left_join(Effort.daily.north,by="Same.return.SNo")%>%
                    filter(!is.na(hook.hours))%>%
                    filter(LIVEWT.c>0)



#Change variable names to match observers data
Data.monthly.north=Data.monthly.north%>%
                      rename(SHEET_NO=Same.return,
                             YEAR=FINYEAR,
                             BLOCK=BLOCKX,
                             BOAT=VESSEL,
                             EFFORT=hook.hours,
                             LATITUDE=LAT,
                             LONGITUDE=LONG)%>%
                      mutate(SKIPPER=NA)

Data.daily.north=Data.daily.north%>%
                      rename(SHEET_NO=Same.return.SNo,
                             YEAR=FINYEAR,
                             BLOCK=BLOCKX,
                             BOAT=VESSEL,
                             EFFORT=hook.hours,
                             LATITUDE=LAT,
                             LONGITUDE=LONG)%>%
                      mutate(SKIPPER=NA)


#Remove nonsense catches (either too small or too large)
Data.monthly.north=Data.monthly.north%>%filter(LIVEWT.c>1 & LIVEWT.c<35000)    
Data.daily.north=Data.daily.north%>%filter(LIVEWT.c>.1 & LIVEWT.c<15000)  


##

#----TDGDLF----

#Use only TDGDLF GN records
Data.monthly=Data.monthly%>%filter(METHOD=="GN" & Estuary=="NO" & 
                                     NETLEN.c>=100 & LAT<=(-26) &
                                     !FINYEAR%in%unique(Data.daily$FINYEAR))%>%
                    select(Same.return,FINYEAR,MONTH,VESSEL,METHOD,
                             BLOCKX,SPECIES,SNAME,YEAR.c,LIVEWT.c,LAT,LONG,TYPE.DATA,
                             zone,NETLEN.c)
Data.daily.LL=Data.daily%>%
  filter(METHOD=="LL" & Estuary=="NO" & 
           HOOKS>=80 & LAT<=(-26))%>%
  select(Same.return.SNo,FINYEAR,MONTH,VESSEL,METHOD,
         BLOCKX,SPECIES,SNAME,YEAR.c,LIVEWT.c,LAT,LONG,TYPE.DATA,
         zone,HOOKS)
Data.daily=Data.daily%>%
              filter(METHOD=="GN" & Estuary=="NO" & 
                                     NETLEN>=100 & LAT<=(-26))%>%
              select(Same.return.SNo,FINYEAR,MONTH,VESSEL,METHOD,
                         BLOCKX,SPECIES,SNAME,YEAR.c,LIVEWT.c,LAT,LONG,TYPE.DATA,
                         zone,NETLEN)


SHOTS_c=unique(Data.monthly$Same.return)

#Select corresponding effort (km gn days)
daily.yrs=unique(Effort.daily$finyear)

  #monthly
Effort.monthly=Effort.monthly%>%
                  filter(!FINYEAR%in%daily.yrs)%>%
                  filter(Same.return%in%SHOTS_c)%>%
                  group_by(Same.return)%>%
                  summarise(Km.Gillnet.Hours.c=max(Km.Gillnet.Hours.c,na.rm=T))

  #daily
Effort.daily=Effort.daily%>%
              group_by(Same.return.SNo)%>%
              filter(!is.na(Km.Gillnet.Hours.c))%>%
              summarise(Km.Gillnet.Hours.c=max(Km.Gillnet.Hours.c,na.rm=T))


#Attach effort
Data.monthly=Data.monthly%>%
                  left_join(Effort.monthly,by="Same.return")%>%
                  filter(!is.na(Km.Gillnet.Hours.c))%>%
                  filter( LIVEWT.c>0)
Data.daily=Data.daily%>%
                  left_join(Effort.daily,by="Same.return.SNo")%>%
                  filter(!is.na(Km.Gillnet.Hours.c))%>%
                  filter(LIVEWT.c>0)



#Change variable names to match observers data
Data.monthly=Data.monthly%>%
                rename(SHEET_NO=Same.return,
                       YEAR=FINYEAR,
                       BLOCK=BLOCKX,
                       BOAT=VESSEL,
                       EFFORT=Km.Gillnet.Hours.c,
                       LATITUDE=LAT,
                       LONGITUDE=LONG)%>%
              mutate(SKIPPER=NA)

Data.daily=Data.daily%>%
            rename(SHEET_NO=Same.return.SNo,
                   YEAR=FINYEAR,
                   BLOCK=BLOCKX,
                   BOAT=VESSEL,
                   EFFORT=Km.Gillnet.Hours.c,
                   LATITUDE=LAT,
                   LONGITUDE=LONG)%>%
            mutate(SKIPPER=NA)
Data.daily.LL=Data.daily.LL%>%
            rename(SHEET_NO=Same.return.SNo,
                   YEAR=FINYEAR,
                   BLOCK=BLOCKX,
                   BOAT=VESSEL,
                   LATITUDE=LAT,
                   LONGITUDE=LONG)%>%
            mutate(SKIPPER=NA)
          

#Remove nonsense catches (either too small or too large)
Data.monthly=Data.monthly%>%filter(LIVEWT.c>1 & LIVEWT.c<25000)    
Data.daily=Data.daily%>%filter(LIVEWT.c>.1 & LIVEWT.c<8000)  
Data.daily.LL=Data.daily.LL%>%filter(LIVEWT.c>.1 & LIVEWT.c<8000)  

