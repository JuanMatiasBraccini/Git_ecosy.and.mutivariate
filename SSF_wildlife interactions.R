####### WILDLIFE INTERACTIONS IN THE AUSSIE SHARK FISHERIES #####################

#NOTE: script for analysing wildlife interactions in the Southern Shark Fishery
#      based on data from the 70s and 90-00s and the WA Shark Fishery
#       Terry has data for 1970s (162 sites  between Streaky bay and Vic/Tas, different nets)
#                        98-2001 commercial vessels, 153 sites, SA and Vic, 6 and 6.5 inch only


#MISSING:   change s$ for s4 in form 3; add other surveys


#Analysis to do:  (Interactions are anything other than commercial species catch)

# Figure 1. Map of area with management zones and fishery distribution (do densities of total efforts), and each shot, color-coded by study

# Figure 2. Heat Map for each interaction in which blocks the highest occurrences occur!!!! (use this to show spatial heterogeneity) leave blank unobserved blocks)

# Figure 3. Mesh size selectivity (all mesh sizes). Show 6/6.5 is best trade off
#           for economic and environmental considerations (economics=catch rates and size * average price for gummy, bronzy, etc)
#           environmental= number of species and number of individuals

# Table 1.  Study; Number of shots observed; Proportion of fishing effort observed; Area; Year; Interaction type: Discarding (# species, Freq. Occ),
#           Incidental capture of TEPS (# species, Freq. Occ), Food provisioning (# species, Freq. Occ), Seabird contacting vessel or gear (# species, Freq. Occ)
# Table 2.  Discards
# Table 3.  TEPS
# Table 4.  Seabird provisioning
# Table 5.  Contacts with vessel or gear

#Estimating total discard (factor PCS in) and total incidental capture of TEPS






library(RODBC)  			#include ODBC library for importing Acccess data
library(rmutil)	#generalised nonlinear regression (including beta distributions)
library(gnlm)
library(boot)	#for weighted correlation
library(MASS)				#for stepAIC model selection
library(relaimpo)			#relative importance of linear model terms
library(bootstrap)
 library(boot)
library(PBSmapping)			#needed to obtain maps

library(ade4)				#to add subplots
library(relaimpo)			#relative importance of linear model terms
library(plotrix)				#needed for graph legends
library(gam)

library(pscl)				#zero inflated binomial distribution
library(VGAM)		#zero inflated binomial distribution
  library(triangle)			#triangular pdf

  library(spam)				#packages for smoothing matrix
  library(fields)

library(gridBase)			#for inset map

library(plyr)				#for weighted mean aggregation

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

setwd(handl_OneDrive("Data/Shark survey/2007-2008"))
#setwd("F:/BackUps/Matias_1_12/Data/Shark survey/2007-2008")


#DATA SECTION

#XXXX survey data


#2007-2008 survey data
channel7 <- odbcConnectExcel("SharkSurveyData_30_09_2008") 
Data_08 <- sqlFetch(channel7, "F2_Sampling")
F1_Sampling_08<- sqlFetch(channel7,"F1_Sampling")
F1_SamplingTwo_08<- sqlFetch(channel7,"F1_SamplingTwo")
F3_Sampling_08<- sqlFetch(channel7,"F3_Sampling")
close(channel7) 

Sites=read.csv(handl_OneDrive("Analysis/ArcView/Shark survey 2007-2008/Sites.csv",sep=','))
#Sites=read.csv("F:/BackUps/Matias_1_12/Analysis/ArcView/Shark survey 2007-2008/Sites.csv",sep=',')


#Bathymetry data
Bathymetry_1=read.table(handl_OneDrive("Analysis/Post capture mortality/get_data_114_130.cgi"))
Bathymetry_2=read.table(handl_OneDrive("Analysis/Post capture mortality/get_data_130_140.cgi"))
Bathymetry_3=read.table(handl_OneDrive("Analysis/Post capture mortality/get_data_140_151.cgi"))

# Bathymetry_1=read.table("F:/BackUps/Matias_1_12/Analysis/Post capture mortality/get_data_114_130.cgi")
# Bathymetry_2=read.table("F:/BackUps/Matias_1_12/Analysis/Post capture mortality/get_data_130_140.cgi")
# Bathymetry_3=read.table("F:/BackUps/Matias_1_12/Analysis/Post capture mortality/get_data_140_151.cgi")

Bathymetry_1$pasted=with(Bathymetry_1,paste(V1, V2, V3))
Bathymetry_2$pasted=with(Bathymetry_2,paste(V1, V2, V3))
Bathymetry_3$pasted=with(Bathymetry_3,paste(V1, V2, V3))
Bathymetry=rbind(Bathymetry_1,Bathymetry_2,Bathymetry_3)
Bathymetry=Bathymetry[!duplicated(Bathymetry$pasted),]
Bathymetry=Bathymetry[,-4]




#PARAMETERS SECTION


#PROCEDURE SECTION

# 1. ---Mapping
F1_SamplingTwo$cruiseStation=paste(F1_SamplingTwo$Cruise,F1_SamplingTwo$Station)
latlong=F1_SamplingTwo[!duplicated(F1_SamplingTwo$cruiseStation),]
a=PCS_data[!(is.na(PCS_data$Stimuli)),]
a$cruiseStation=paste(a$Cruise,a$Station)
a=a[!duplicated(a$cruiseStation),]
b=a[,c(1:2,21)]
latlong=merge(latlong,b,"cruiseStation")
 latlong$StartLat=-(latlong$StartLat/100)
 latlong$StartLong=(latlong$StartLong/100)

latlong$StartLat=ifelse(latlong$Cruise.x%in%c(6,7),latlong$StartLat-0.31,latlong$StartLat)
latlong$StartLat=ifelse(latlong$Cruise.x==16 & latlong$StartLat>-35.779,35.785,latlong$StartLat)


extralong=c(131.1170, 131.2479, 131.3788, 131.5097, 131.6842, 131.8587, 132.7313, 131.9460,138.4901, 138.7519,
138.9264, 139.1445)  #add to move land sites to sea
extralat=c(-31.50102, -31.60555, -31.64039, -31.67523, -31.77976, -31.88429, -32.05850, -32.23271,
-35.64728, -35.64728, -35.78665, -35.89117)


insetOz <- function(){
opar <- par(mai = c(1.5,2,1,1))
on.exit(par(opar))
par(mar=rep(.1,4),xaxt="n",yaxt="n",plt=par("plt"))
plotMap(worldLLhigh, xlim=c(110,155), ylim=c(-44.5,-11),col="light grey", axes=F, xlab="", ylab="")
text(133,-25,("Australia"),col="black", cex=1.5)
points(153.2199,-11.53294,col="white",pch=21,bg="white")
  # add polygons
edgeX=c(127,151.5,151.5,127) 
edgeY=c(-30,-30,-44,-44)
polygon(x=edgeX,y=edgeY,lwd=2)
}

#Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
#Bathymetry=subset(Bathymetry,V3%in%c(-400:0))
#xbat=sort(unique(Bathymetry$V1))
#ybat=sort(unique(Bathymetry$V2))
#reshaped=as.matrix(reshape(Bathymetry,idvar="V1",  	#data frame transposed as a matrix (needed for image()		
#	timevar="V2",v.names="V3", direction="wide"))


# 2. ---Data manipulations

#2.1. combine data sets

#2007-08
  #add form 1 and 2
vars.08=c("Cruise","Station","Mesh1","Species","Csiro","Ret")
Data_08=Data_08[,match(vars.08,names(Data_08))]
vars.08=c("Cruise","Station","Year","Month","Day")
F1_Sampling_08=F1_Sampling_08[,match(vars.08,names(F1_Sampling_08))]
vars.08=c("Cruise","Station","Mesh","TimeSet","TimeUp",
          "Subs","Topo","Rock","Sponge","Coral","Other","Comments")
F1_SamplingTwo_08=F1_SamplingTwo_08[,match(vars.08,names(F1_SamplingTwo_08))]

DATA_08=merge(Data_08,F1_Sampling_08,by=c("Cruise","Station"),all.x=T)
DATA_08=merge(DATA_08,F1_SamplingTwo_08,by.x=c("Cruise","Station","Mesh1"),
              by.y=c("Cruise","Station","Mesh"),all.x=T)

  #add form 3
    #expand form 3 to include all combos of cruise, station, mesh
unique.cruise=unique(DATA_08$Cruise)
unique.cruise=unique.cruise[!(is.na(unique.cruise))]
unique.Station=unique(DATA_08$Station)
unique.Station=unique.Station[!(is.na(unique.Station))]
unique.Mesh1=unique(DATA_08$Mesh1)
unique.Mesh1=unique.Mesh1[!(is.na(unique.Mesh1))]
Dummy=expand.grid(unique.cruise,unique.Station)
#Dummy=expand.grid(unique.cruise,unique.Station,unique.Mesh1)
colnames(Dummy)=c("Cruise","Station")
#colnames(Dummy)=c("Cruise","Station","Mesh")
Dummy=Dummy[order(Dummy$Cruise),]
dummy.1=with(Dummy,paste(Cruise,Station))
#dummy.1=with(Dummy,paste(Cruise,Station,Mesh))
F3_Sampling_08$Mesh=ifelse(is.na(F3_Sampling_08$Mesh),"no mesh",F3_Sampling_08$Mesh)
dummy.2=with(F3_Sampling_08,paste(Cruise,Station))
dummy.2=unique(dummy.2)
remove.from.dummy=match(dummy.2,dummy.1)
remove.from.dummy=remove.from.dummy[!(is.na(remove.from.dummy))]
Dummy=Dummy[-remove.from.dummy,]
namesF3=names(F3_Sampling_08)[-(1:2)]
test=data.frame(t(namesF3))
names(test)=namesF3
test[1:nrow(Dummy),]=NA
Dummy=cbind(Dummy,test)
F3_Sampling_08=rbind(F3_Sampling_08,Dummy)
#ACA: what to do whith NA in mesh for interactions? how to merge????, not quite there yet!
colnames(F3_Sampling_08)[match("Mesh",colnames(F3_Sampling_08))]="Mesh1"
DATA_08=merge(DATA_08,F3_Sampling_08,by.x=c("Cruise","Station","Mesh1"),all.x=T)

# create bycatch as another interaction

#REPORT SECTION

#---FIGURE 1.
data(worldLLhigh)  					#high resolution version 

  #general map with points for recapture and releases and bathymetry




#X11(width=14,height=11)		#control width and height of graphic
#png(filename = "F:/BackUps/Matias_1_12/Analysis/Post capture mortality/Paper outputs/Figure1.png", width = 1000, height = 700)
#pdf(file="C:/Matias/Analysis/Post capture mortality/Paper outputs/Figure1.pdf")		#create pdf
#tiff(file="F:/BackUps/Matias_1_12/Analysis/Post capture mortality/Paper outputs/Figure1.tiff",width = 3200, height = 3200,units = "px", res = 300, compression = "lzw")    #create tiff
#tiff(file="C:/Matias/Analysis/Post capture mortality/Paper outputs/Figure1.tiff",width = 3200, height = 3200,units = "px", res = 300, compression = "lzw")    #create tiff
 
#par(mfcol=c(1,1),las=1,mai=c(.1, .1, .5, 0.5),omi=c(.1,.1,.5,0.1))
par(mfcol=c(1,1),las=1,mai=c(.2, .8, .5, 0.2),omi=c(.2,.9,.5,0.5))
plotMap(worldLLhigh, ylim=c(-44,-30), xlim=c(127,154.5),col="dark grey",tck=F,plt = c(0.01, 0.99, 0.05, 0.99), xlab="",ylab="", axes=F)
#contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=c(-44,-31), xlim=c(112,154.5), zlim=c(-1,-300),nlevels = 2,
#labcex=0.75,lty = 1,add=T)
par(new=T)
plotMap(worldLLhigh, ylim=c(-44,-30), xlim=c(127,154.5),col="dark grey",tck=F,plt = c(0.01, 0.99, 0.05, 0.99), xlab="",ylab="", axes=F)
points(latlong$StartLong,latlong$StartLat,col="black",pch=19,cex=1.25)
#points(extralong,extralat,col="black",pch=19)
points(131.2916,-35.9957,col="white",pch=19,cex=2)
#points(132.2033,-31.52928,col="dark grey",pch=19,cex=2.3)
#points(c(138.84,139.3126),c(-35.26307,-35.46168),col="dark grey",pch=19,cex=2.3)
arrows(153.4,-41.7,153.4,-40.7,length = 0.1, lwd=1.5)
text(153.4,-41.9,("N"),col="black", cex=1.25)

segments(153,-42.9,153.9,-42.9)
text(153.45,-42.625,("100 km"),col="black", cex=1.25)

text(131.3,-30.5,("Head of the"),col="black", cex=1.5)
text(131.6,-31.1,("Great Australian Bight"),col="black", cex=1.5)
text(145.7,-39.8,("Bass Strait"),col="black", cex=1.5)
mtext("Latitude (?S)",side=2,outer=T,line=2.5,font=1,las=0,cex=1.5)
mtext("Longitude (?E)",side=1,outer=T,line=-8,font=1,las=1,cex=1.5)
axis(2,at=seq(-42,-32,2),labels=rev(seq(32,42,2)),cex.axis=1.5)
axis(1,at=seq(128,154,2),labels=T,cex.axis=1.5)
box()

vp <- baseViewports()
pushViewport(vp$inner,vp$figure,vp$plot)
pushViewport(viewport(x=0.005,y=0.7,width=.4,height=.4,just=c("left","top")))
par(fig=gridFIG(),new=T)  
insetOz()
#dev.off()
#X11()  #reset with and height to defaut
