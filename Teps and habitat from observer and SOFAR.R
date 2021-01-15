# SCRIPT FOR EXTRACTING PROTECTED SPECIES AND HABITAT OBSERVATIONS FROM SHARK BIO"

library(RODBC)
library(lunar)   #moon phases
library(lubridate)
library(tidyverse)
options(stringsAsFactors = FALSE)


#missing: Parks Australia and whatever I want to do with observer-habitat


# DATA SECTION -----------------------------------------------------------------------

#Sharks data base
setwd("U:/Shark")  # working directory
channel <- odbcConnectAccess2007("Sharks v20200323.mdb")  #new databased updated by Vero
Boat_hdr=sqlFetch(channel, "Boat_hdr", colnames = F)   
close(channel)


#McAuley & Simpfendorfer 2003
McAuley.Simp_rates=data.frame(
                    Species=c('marine.mammals','seabirds','turtles'),
                    number.caught.per.km.gn.hour=c(1/10000,4/100000,1/100000))
McAuley.Simp_Table.5=data.frame(
                    Species=c("Dolphin (common)","Dolphin (unspecified)","Dolphin (bottlenose)",
                              "Cormorants","Fairy penguin","Mutton bird","Sea Lion","Seal",
                              "Turtle (unspecified)"),
                    number.all.regions=c(4,2,2,rep(1,6)))
#Annual effort
Annual.total.eff.hours=read.csv('C:/Matias/Analyses/Data_outs/Annual.total.eff.hours.csv')
Annual.total.eff.hours$Total=1000*Annual.total.eff.hours$Total  #convert to 1000s of km gn h


#Sofar TEPS
Current.yr="2018_19"
SOFAR.TEPS=read.csv(paste("C:/Matias/Data/Catch and Effort",Current.yr,"TEPS_PROTECTEDSP.csv",sep='/'))


#Parks Australia


# TEPS_Get teps records from Shark survey boat header -----------------------------------------------------------------------
GN_hdr=Boat_hdr%>%
      filter(Method=='GN' & !BOAT%in%c("NAT"))%>%
      mutate(COMMENTS=tolower(COMMENTS))
search.teps=tolower(c("PRION","Mutton","shearwater",'shear','lion','seal','sealion',
                      'albatross','smooth','penguin','SM','manta','turtl',
                      'dead wp','dolphin','nurse'))
search.habitat=tolower(c('algae','kelp','weed','eklonia','coral','sargassum'))
Teps.in.GN.boat.header=GN_hdr%>%
      filter(grepl(paste(search.teps, collapse="|"), COMMENTS))


# TEPS_Scale up observed TEP interactions to whole fishery -----------------------------------------------------------------------
hndl.teps='C:/Matias/Analyses/Ecosystem indices and multivariate/Shark-bycatch/Outputs/TEPS/'

  #get total numbers by group
Scaled.obs=Annual.total.eff.hours%>%
            mutate(marine.mammals_n.per.km.gn.h=McAuley.Simp_rates$number.caught.per.km.gn.hour[match('marine.mammals',McAuley.Simp_rates$Species)],
                   seabirds_n.per.km.gn.h=McAuley.Simp_rates$number.caught.per.km.gn.hour[match('seabirds',McAuley.Simp_rates$Species)],
                   turtles_n.per.km.gn.h=McAuley.Simp_rates$number.caught.per.km.gn.hour[match('turtles',McAuley.Simp_rates$Species)])%>%
            mutate(total.number_marine.mammals=marine.mammals_n.per.km.gn.h*Total,
                   total.number_seabirds=seabirds_n.per.km.gn.h*Total,
                   total.number_turtles=turtles_n.per.km.gn.h*Total)%>%
  rename(Annual.effort_km.gn.h=Total)

  #reapportion group to species
marine.mammals=c('Dolphin (common)','Dolphin (unspecified)','Dolphin (bottlenose)','Sea Lion','Seal')
seabirds=c('Cormorants','Fairy penguin','Mutton bird')

Marine.mammals=data.frame(marine.mammals)%>%
                left_join(McAuley.Simp_Table.5,by=c('marine.mammals'='Species'))%>%
                mutate(prop=number.all.regions/sum(number.all.regions))%>%
                dplyr::select(-number.all.regions)%>%
                spread(marine.mammals,prop)
Marine.mammals=cbind(Marine.mammals,
                     Scaled.obs%>%
                       dplyr::select(FINYEAR,total.number_marine.mammals))
Marine.mammals=Marine.mammals%>%
  mutate_at(vars(-FINYEAR,-total.number_marine.mammals), function(x) x*Marine.mammals$total.number_marine.mammals)%>%
  dplyr::select(-total.number_marine.mammals)%>%
  mutate_at(vars(-FINYEAR),round)

Seabirds=data.frame(seabirds)%>%
  left_join(McAuley.Simp_Table.5,by=c('seabirds'='Species'))%>%
  mutate(prop=number.all.regions/sum(number.all.regions))%>%
  dplyr::select(-number.all.regions)%>%
  spread(seabirds,prop)
Seabirds=cbind(Seabirds,
                     Scaled.obs%>%
                       dplyr::select(FINYEAR,total.number_seabirds))
Seabirds=Seabirds%>%
  mutate_at(vars(-FINYEAR,-total.number_seabirds), function(x) x*Seabirds$total.number_seabirds)%>%
  dplyr::select(-total.number_seabirds)%>%
  mutate_at(vars(-FINYEAR),round)

Turtles=Scaled.obs%>%
  dplyr::select(FINYEAR,total.number_turtles)%>%
  mutate('Turtle (unspecified)'=round(total.number_turtles))%>%
  dplyr::select(-total.number_turtles)

Tab1=Scaled.obs%>%
        dplyr::select(FINYEAR,Annual.effort_km.gn.h)%>%
        left_join(Marine.mammals,by="FINYEAR")%>%
        left_join(Seabirds,by="FINYEAR")%>%
        left_join(Turtles,by="FINYEAR")%>%
        mutate(Annual.effort_km.gn.h=round(Annual.effort_km.gn.h))

write.csv(Tab1,paste(hndl.teps,"Table_scaled.teps.observer_McAuley.Simp_2003.csv",sep=''),row.names=F)



# TEPS_SOFAR records -----------------------------------------------------------------------
SOFAR.TEPS=SOFAR.TEPS%>%
  mutate(CommonName=case_when(SpeciesCode==40041050 ~'Mutton bird',
                              SpeciesCode==37010003 ~ "Shark, white",
                              SpeciesCode==37035000 ~ "Manta Rays",
                              TRUE ~ CommonName))%>%
  filter(!SpeciesCode%in%c(37020000,37018001,37018003,37018022,37038000,37990030))
sofar.tep.names=SOFAR.TEPS%>%distinct(SpeciesCode,.keep_all = T)%>%
  dplyr::select(SpeciesCode,CommonName,RSCommonName)
sofar.CategoryName=SOFAR.TEPS%>%distinct(SpeciesCode,.keep_all = T)%>%
  dplyr::select(CommonName,CategoryName)


Table.sofar.teps=SOFAR.TEPS%>%
          group_by(finyear,CommonName,Status)%>%
          summarise(Numbers=sum(Number))%>%
          mutate(Yr_status=paste(finyear,Status))%>%
          dplyr::select(-Status)

Table.sofar.teps=Table.sofar.teps[,-1]%>%
          spread(Yr_status,Numbers,fill='')

Table.sofar.teps=Table.sofar.teps%>%
  left_join(sofar.CategoryName,by='CommonName')%>%
  arrange(CategoryName)%>%
  dplyr::select(-CategoryName)
write.csv(Table.sofar.teps,paste(hndl.teps,"Table_SOFAR.csv",sep=''),row.names=F)


# Habitat_get habitat records from Shark survey boat header -----------------------------------------------------------------------
Habitat.in.GN.boat.header=GN_hdr%>%
  filter(grepl(paste(search.habitat, collapse="|"), COMMENTS))