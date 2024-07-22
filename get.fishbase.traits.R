library(rfishbase)
library(tidyverse)


a=read.csv(handl_OneDrive('Analyses/Ecosystem indices and multivariate/Shark-bycatch/SPECIES+PCS+FATE.csv'))
b=read.csv(handl_OneDrive('Analyses/Ecosystem indices and multivariate/Shark-bycatch/SPECIES+PCS+FATE_Com.csv'))

fish=sort(unique(c(a$SCIENTIFIC_NAME,b$SCIENTIFIC_NAME)))
fish=subset(fish,!fish=='')
fish=str_replace(fish, "\xa0", " ")
fish=str_replace(fish, "Family ", "")
fish=unique(fish)

Traits.col=c('Maximum.size','Body.depth','Head.length','Eye.diameter','Pre.orbital.length','Trophic.level',
                    'Habitat','Feeding.group')
Traits=as.data.frame(matrix(nrow=length(fish),ncol=length(Traits.col)+1))
colnames(Traits)=c('SCIENTIFIC_NAME',Traits.col)

Traits=Traits%>%
        mutate(SCIENTIFIC_NAME.original=fish,
               SCIENTIFIC_NAME=case_when(SCIENTIFIC_NAME.original=="Alopiidae"~"Alopias pelagicus",               
                                         SCIENTIFIC_NAME.original=="Balistidae"~"Xanthichthys lineopunctatus",
                                         SCIENTIFIC_NAME.original=="Bothidae & Paralichthyidae"~"Pseudorhombus argus",
                                         SCIENTIFIC_NAME.original=="Carangidae"~"Carangoides chrysophrys",
                                         SCIENTIFIC_NAME.original%in%c("Carcharhiniformes","Carcharhinus sp.")~"Carcharhinus plumbeus",
                                         SCIENTIFIC_NAME.original=="Centroberyx sp."~"Centroberyx australis",
                                         SCIENTIFIC_NAME.original=="Cheilodactylidae"~"Goniistius rubrolabiatus",
                                         SCIENTIFIC_NAME.original=="Dasyatidae & Urolophidae"~"Urolophus paucimaculatus",
                                         SCIENTIFIC_NAME.original=="Istiophoridae"~"Istiophorus platypterus",                                      
                                         SCIENTIFIC_NAME.original%in%c("Lethrinidae","Lethrinidae")~"Lethrinus nebulosus",
                                         SCIENTIFIC_NAME.original=="Mugilidae"~"Mugil cephalus",
                                         SCIENTIFIC_NAME.original=="Loliginidae"~"",
                                         SCIENTIFIC_NAME.original=="Monacanthidae"~"Aluterus monoceros",
                                         SCIENTIFIC_NAME.original=="Muraenidae"~"Gymnomuraena zebra",
                                         SCIENTIFIC_NAME.original=="Neosebastes sp."~"Neosebastes scorpaenoides",
                                         SCIENTIFIC_NAME.original=="Orectolobus spp."~"Orectolobus maculatus",
                                         SCIENTIFIC_NAME.original=="Ostraciidae & Anacaridae"~"Lactoria cornuta",
                                         SCIENTIFIC_NAME.original=="Pentacerotidae"~"Paristiopterus labiosus",
                                         SCIENTIFIC_NAME.original=="Platycephalidae"~"Platycephalus laevigatus",
                                         SCIENTIFIC_NAME.original=="Plotosidae"~"Euristhmus microceps",      
                                         SCIENTIFIC_NAME.original=="Pristiophoridae"~"Pristiophorus cirratus",
                                         SCIENTIFIC_NAME.original=="Rajidae"~"Dipturus acrobelus",
                                         SCIENTIFIC_NAME.original=="Rhinidae"~"Rhina ancylostomus",
                                         SCIENTIFIC_NAME.original%in%c("Rhinobatidae","Rhinobatidae & Rhynchobatidae")~"Rhynchobatus australiae",
                                         SCIENTIFIC_NAME.original=="Scaridae"~"Scarus frenatus",
                                         SCIENTIFIC_NAME.original%in%c("Scomberoides sp.","Scombridae","Scombridae")~"Scomberomorus commerson",
                                         SCIENTIFIC_NAME.original=="Scorpaenidae"~"Scorpaena sumptuosa",
                                         SCIENTIFIC_NAME.original%in%c("Scyliorhinidae ","Scyliorhinidae")~"Cephaloscyllium laticeps",
                                         SCIENTIFIC_NAME.original=="Sepiidae"~"",
                                         SCIENTIFIC_NAME.original=="Serranidae"~"Acanthistius serratus",
                                         SCIENTIFIC_NAME.original%in%c("Sphyraena sp.","Sphyraenidae")~"Sphyraena obtusata",
                                         SCIENTIFIC_NAME.original%in%c("Squalus sp.","Squalidae")~"Squalus megalops",
                                         SCIENTIFIC_NAME.original=="Squatinidae"~"Squatina australis",
                                         SCIENTIFIC_NAME.original=="Teleostei"~"Pagrus auratus",
                                         SCIENTIFIC_NAME.original=="Sphyrnidae"~"Sphyrna zygaena",
                                         SCIENTIFIC_NAME.original=="Trachinocephalus spp"~"Trachinocephalus trachinus",
                                         SCIENTIFIC_NAME.original=="Triglidae"~"Lepidotrigla papilio",
                                         SCIENTIFIC_NAME.original=="Carcharhinus limbatus/tilstoni"~"Carcharhinus limbatus",
                                         SCIENTIFIC_NAME.original=="Mustelus sp."~"Mustelus antarcticus",
                                         SCIENTIFIC_NAME.original=="Myliobatidae"~"Myliobatis tenuicaudatus",
                                         SCIENTIFIC_NAME.original=="Sphyrna sp."~"Sphyrna zygaena",
                                         SCIENTIFIC_NAME.original=="Uranoscopidae"~"Kathetostoma laeve",
                                         SCIENTIFIC_NAME.original=="Arripis georgianus"~"Arripis georgiana",
                                         SCIENTIFIC_NAME.original=="Aulopus purpurissatus"~"Latropiscis purpurissatus",
                                         SCIENTIFIC_NAME.original=="Carcharhinus amblyrhinchoides"~"Carcharhinus plumbeus",
                                         SCIENTIFIC_NAME.original=="Centrophorus acus"~"Centrophorus granulosus",
                                         SCIENTIFIC_NAME.original=="Cheilodactylus rubrolabiatus"~"Goniistius rubrolabiatus",
                                         SCIENTIFIC_NAME.original=="Chrysophrys auratus"~"Pagrus auratus",
                                         SCIENTIFIC_NAME.original=="Epinephelus ergastularius"~"Hyporthodus ergastularius",
                                         SCIENTIFIC_NAME.original=="Hypothalassia acerba"~"Lagocephalus sceleratus",
                                         SCIENTIFIC_NAME.original=="Moolgarda buchanani"~"Crenimugil buchanani",
                                         SCIENTIFIC_NAME.original=="Plectorhincus flavomaculatus"~"Plectorhinchus flavomaculatus",
                                         SCIENTIFIC_NAME.original=="Portunus armatus"~"",
                                         SCIENTIFIC_NAME.original=="Hemiscylliidae"~"Hemiscyllium ocellatum",
                                         SCIENTIFIC_NAME.original=="Rachycentron canadus"~"Rachycentron canadum",
                                         SCIENTIFIC_NAME.original=="Rhina ancylostoma"~"Rhina ancylostomus",
                                         SCIENTIFIC_NAME.original=="Schedophilus labyrinthica"~"Schedophilus velaini",
                                         SCIENTIFIC_NAME.original=="Scorpis georgianus"~"Scorpis georgiana",
                                         SCIENTIFIC_NAME.original=="Seriola hipos"~"Seriola hippos",
                                         SCIENTIFIC_NAME.original=="Stegastoma fasciatum"~"Stegostoma tigrinum",
                                         SCIENTIFIC_NAME.original=="Tilodon sexfasciatum"~"Tilodon sexfasciatus",
                                         SCIENTIFIC_NAME.original=="Lagocephalus scleratus"~"Lagocephalus sceleratus",
                                         TRUE~SCIENTIFIC_NAME.original))%>%
  filter(!SCIENTIFIC_NAME=='')
  
fishes=unique(Traits$SCIENTIFIC_NAME)

Fish.base.traits=fb_tbl("species") %>% 
                  mutate(SCIENTIFIC_NAME = paste(Genus, Species)) %>%
                  filter(SCIENTIFIC_NAME %in% fishes) %>% 
                  select(SCIENTIFIC_NAME, Length)%>%
                  rename(Maximum.size=Length)
#fishes[which(!fishes%in%Fish.base.traits$SCIENTIFIC_NAME)]


Traits=Traits%>%
  select(-c(Maximum.size))%>%
  left_join(Fish.base.traits,by='SCIENTIFIC_NAME')


Fish.base.morph=morphometrics(species_list=fishes)%>%
                    rename(SCIENTIFIC_NAME=Species)%>%
                    select(SCIENTIFIC_NAME,TL,BD,HL,ED,POL)%>%
                    mutate(HL=case_when(is.na(HL & SCIENTIFIC_NAME=='Pristiophorus nudipinnis')~213,
                                        is.na(HL & SCIENTIFIC_NAME=='Pristis pristis')~227.62,
                                        is.na(HL& SCIENTIFIC_NAME=='Pristis zijsron')~202.9,
                                        is.na(HL & SCIENTIFIC_NAME=='Rhina ancylostomus')~124,
                                        is.na(HL & SCIENTIFIC_NAME=='Squatina tergocellata')~107.4,
                                        is.na(HL & SCIENTIFIC_NAME=='Anoxypristis cuspidata')~207.86,
                                        is.na(HL & SCIENTIFIC_NAME=='Mustelus antarcticus')~110.968,
                                        TRUE~HL),
                           BD=100*BD/TL,
                           ED=100*ED/HL,
                           POL=100*POL/HL,
                           HL=100*HL/TL)%>%
                    group_by(SCIENTIFIC_NAME)%>%
                    summarise(TL=mean(TL,na.rm=T),
                              BD=mean(BD,na.rm=T),
                              HL=mean(HL,na.rm=T),
                              ED=mean(ED,na.rm=T),
                              POL=mean(POL,na.rm=T))%>%
                    ungroup()%>%
                    mutate(BD=case_when(is.na(BD & SCIENTIFIC_NAME=='Bathytoshia lata')~0,
                                        is.na(BD & SCIENTIFIC_NAME=='Kathetostoma laeve')~0,
                                        is.na(BD & SCIENTIFIC_NAME=='Pristiophorus cirratus')~9.8,
                                        is.na(BD & SCIENTIFIC_NAME=='Pristiophorus nudipinnis')~9.8,
                                        is.na(BD & SCIENTIFIC_NAME=='Pristis pristis')~9.8,
                                        is.na(BD & SCIENTIFIC_NAME=='Pristis zijsron')~9.8,
                                        is.na(BD & SCIENTIFIC_NAME=='Rhina ancylostomus')~9.8,
                                        is.na(BD & SCIENTIFIC_NAME=='Squatina tergocellata')~8,
                                        is.na(BD & SCIENTIFIC_NAME=='Anoxypristis cuspidata')~9.8,
                                        TRUE~BD))%>%
                      select(-TL)%>%
  rename(Body.depth=BD,
         Head.length=HL,
         Eye.diameter=ED,
         Pre.orbital.length=POL)

Traits=Traits%>%
  select(-c(Body.depth,Head.length,Eye.diameter,Pre.orbital.length))%>%
  left_join(Fish.base.morph,by='SCIENTIFIC_NAME')%>%
  mutate(Body.depth=case_when(is.na(Body.depth) & SCIENTIFIC_NAME=='Acanthistius serratus'~33.2,
                              is.na(Body.depth) & grepl(paste(c('Orectolobus','Aptychotrema','Trygonorrhina'),collapse = '|'),SCIENTIFIC_NAME)~9.8,
                              TRUE~Body.depth),
         Head.length=case_when(is.na(Head.length) & SCIENTIFIC_NAME=='Acanthistius serratus'~33,
                               is.na(Head.length) & grepl(paste(c('Orectolobus','Aptychotrema','Trygonorrhina'),collapse = '|'),SCIENTIFIC_NAME)~22.1,
                               TRUE~Head.length),
         Eye.diameter=case_when(is.na(Eye.diameter) & SCIENTIFIC_NAME=='Acanthistius serratus'~17.8,
                                is.na(Eye.diameter) & grepl(paste(c('Orectolobus','Aptychotrema','Trygonorrhina'),collapse = '|'),SCIENTIFIC_NAME)~11.7,
                                TRUE~Eye.diameter),
         Pre.orbital.length=case_when(is.na(Pre.orbital.length) & SCIENTIFIC_NAME=='Acanthistius serratus'~24.9,
                                      is.na(Pre.orbital.length) & grepl(paste(c('Orectolobus','Aptychotrema','Trygonorrhina'),collapse = '|'),SCIENTIFIC_NAME)~25.8,
                                      TRUE~Pre.orbital.length))


Fish.base.ecol=ecology(unique(Traits$SCIENTIFIC_NAME))%>%
                select(Species,DietTroph,FoodTroph,FeedingType)%>%
                mutate(Trophic.level=mean(c(DietTroph,FoodTroph),na.rm=T))%>%
                select(-c(DietTroph,FoodTroph))%>%
                rename(SCIENTIFIC_NAME=Species,
                       Feeding.group=FeedingType)

Traits=Traits%>%
  select(-c(Trophic.level,Feeding.group))%>%
  left_join(Fish.base.ecol,by='SCIENTIFIC_NAME')


Fish.base.morphol=morphology(species_list=fishes)%>%
  rename(SCIENTIFIC_NAME=Species,
         Body.shape=BodyShapeI)%>%
  select(SCIENTIFIC_NAME,Body.shape)

Traits=Traits%>%
  left_join(Fish.base.morphol,by='SCIENTIFIC_NAME')

Fish.base.habitat=species(fishes)%>%
  select(Species,DemersPelag)%>%
  rename(SCIENTIFIC_NAME=Species,
         Habitat=DemersPelag)

Traits=Traits%>%
  select(-Habitat)%>%
  left_join(Fish.base.habitat,by='SCIENTIFIC_NAME')

write.csv(Traits,handl_OneDrive('Analyses/Ecosystem indices and multivariate/Shark-bycatch/Functional.traits.csv'),row.names = F)
# Traits.col.pretty=c('Maximum size (cm)','Body depth (% of size)','Head length (% of size)',
#                     'Eye diameter (% of HL)','Pre-orbital length (% of HL)','Trophic level',
#                     'Habitat','Feeding group')
