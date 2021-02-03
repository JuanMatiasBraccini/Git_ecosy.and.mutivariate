# SCRIPT FOR EXTRACTING PROTECTED SPECIES AND HABITAT OBSERVATIONS FROM SHARK BIO"

library(RODBC)
library(lunar)   #moon phases
library(lubridate)
library(tidyverse)
library(reshape2)
library(Hmisc)
library(RColorBrewer)
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

#Plots
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}
fn.plt=function(d,mammals,seabirds,reptiles,elasmos,out.file,Year.max)
{
  d <- d%>%
    mutate(Year=as.numeric(substr(FINYEAR,1,4)))%>%
    dplyr::select(-Annual.effort_km.gn.h,-FINYEAR)%>%
    melt(id.vars="Year")%>%
    mutate(variable=as.character(variable),
           variable=case_when(variable=='Dolphin (bottlenose)'~'Bottlenose dolphin',
                              variable=='Dolphin (common)'~'Common dolphin',
                              variable=='Cormorants'~'Cormorant (unspecified)',
                              TRUE ~ variable))%>%
    mutate(group=case_when(variable%in%mammals~'Marine mammals',
                           variable%in%seabirds~'Seabirds',
                           variable%in%reptiles~'Reptiles',
                           variable%in%elasmos~'Elasmobranchs'),
           group=factor(group,levels=c('Elasmobranchs','Marine mammals','Reptiles','Seabirds')))%>%
    filter(Year<=Year.max)
  myColors <- brewer.pal(4, "Spectral")
  names(myColors) <- levels(d$group)
  
  d%>%
    ggplot(aes(Year,value, colour=group)) + 
    geom_point(size=3) + 
    facet_wrap(~variable,scales="free_y")+
    xlab('Financial year')+ylab('Number of interactions')+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          strip.text = element_text(size = 17),
          legend.title = element_blank(),
          legend.text = element_text(size = 18),
          legend.position="top",
          axis.text=element_text(size=15),
          axis.title=element_text(size=20),
          plot.margin=unit(c(.1,.5,.1,.1),"cm"))+
    scale_colour_manual(name = "group", values = myColors)+
    scale_y_continuous(breaks = integer_breaks())
  ggsave(paste(hndl.teps,out.file,sep=''),width = 10,height = 10,compression = "lzw")
  
}
fn.plt(d=Tab1,
       mammals=c("Bottlenose dolphin","Common dolphin","Dolphin (unspecified)","Sea Lion","Seal"),
       seabirds=c("Cormorant (unspecified)","Fairy penguin","Mutton bird"),
       reptiles=c("Turtle (unspecified)"),
       elasmos=NA,
       out.file="Plot_scaled.teps.observer_McAuley.Simp_2003.tiff",
       Year.max=2019)


  
# TEPS_SOFAR records -----------------------------------------------------------------------
SOFAR.TEPS=SOFAR.TEPS%>%
  mutate(CommonName=case_when(SpeciesCode==40041050 ~'Mutton bird',
                              SpeciesCode==37010003 ~ "Shark, white",
                              SpeciesCode==37035000 | CommonName%in%c('Manta Rays','Ray, Manta') ~ "Manta rays",
                              TRUE ~ CommonName))%>%
  filter(!SpeciesCode%in%c(37020000,37018001,37018003,37018022,37038000,37990030))%>%
  mutate(CommonName=capitalize(tolower(CommonName)),
         CommonName=case_when(CommonName=='New zealand fur-seal'~'NZ fur-seal',
                              CommonName=='Shark, grey nurse'~'Grey nurse shark',
                              CommonName=='Shark, white'~'White shark',
                              CommonName=='Snake, sea'~'Sea snakes',
                              TRUE ~ CommonName))

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

dd=Table.sofar.teps%>%
      melt(id.vars="CommonName")%>%
      mutate(value=as.numeric(value),
             FINYEAR=substr(variable,1,7))%>%
      group_by(CommonName,FINYEAR)%>%
      summarise(value=sum(value,na.rm=T))%>%
      spread(CommonName,value)%>%
      mutate(Annual.effort_km.gn.h=NA)


fn.plt(d=dd,
       mammals=c("Dolphins","NZ fur-seal","Sea lions","Seals","Whales"),
       seabirds=c("Mutton bird","Sea birds"),
       reptiles=c("Sea snakes","Turtles"),
       elasmos=c("Grey nurse shark","Manta rays","Sawfish","White shark"),
       out.file="Plot_SOFAR.tiff",
       Year.max=2019)

