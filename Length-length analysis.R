
############################
#LENGTH-LENGTH RELATIONSHIP#
############################
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

#Define working directory
#Define user
User="Matias"
#User="Agustin"

if(User=="Matias")setwd(handl_OneDrive("Analyses/Ecosystem indices/Shark-bycatch"))
if(User=="Agustin")setwd("C:/Users/Agust?n/OneDrive/Research/2016/2. Length-length")

#Source functions
if(User=="Matias")source(handl_OneDrive("Analyses/Ecosystem indices/Shark-bycatch/Git_bycatch_TDGDLF/Length-length functions.R"))
if(User=="Agustin")source("Length-length functions.R")



#Data source

#Bring in WA shark gillnet data
DATA=read.csv("WA_Agustin.csv",stringsAsFactors=F)

#Bring in WA Species names
SPECIES.NAMES=read.csv("SPECIES NAMES.csv",stringsAsFactors=F)



#Get the final data frame

#Extract year and month
DATA$DATE=as.Date(DATA$date,format="%d/%m/%Y")
DATA$YEAR=as.numeric(strftime(DATA$DATE, format="%Y"))
DATA$MONTH=as.numeric(strftime(DATA$DATE, format="%m"))

#Fix method and mesh size issues
DATA$Method=with(DATA,ifelse(is.na(Method) & BOAT%in%c("B67","F244","F517","E35"),"GN",Method))
DATA$MESH_SIZE=with(DATA,ifelse(is.na(MESH_SIZE) & BOAT=="B67" & YEAR>1997,"7",
                                ifelse(is.na(MESH_SIZE) & BOAT=="F517" & YEAR>2002,"7",MESH_SIZE)))

#Fix sex
DATA$SEX=as.character(DATA$SEX)
DATA$SEX=with(DATA,ifelse(SEX%in%c("f","F"),"F",ifelse(SEX%in%c("m","M"),"M","U")))

#Add common and scientific names
DATA=merge(DATA,SPECIES.NAMES,by="SPECIES",all.x=T)

#Select relevant variables
#DATA=subset(DATA,select=c(SHEET_NO,YEAR,SPECIES,COMMON_NAME,NATURE,TL,FL,SEX))

#Remove really dodgy outliers that are clear typos
A=with(DATA,which(COMMON_NAME=="Graceful shark" & TL>160))         #R00711
B=with(DATA,which(COMMON_NAME=="Gummy Shark" & FL<45))             #C00099
C=with(DATA,which(COMMON_NAME=="Lemon Shark" & TL<150 & FL>150))   #R00386
D=with(DATA,which(COMMON_NAME=="Longnose Grey, Spinner Shark" & TL<150 & FL>150))#C00073
E=with(DATA,which(COMMON_NAME=="Milk Shark" & TL<75 & FL>70))      #N00694,N00816,N00712
F=with(DATA,which(COMMON_NAME=="Nervous Shark" & TL>150))          #R00951
G=with(DATA,which(COMMON_NAME=="Pencil Shark" & TL<20))            #T00006 (TL=100, no 10)
H=with(DATA,which(COMMON_NAME=="Port Jackson" & FL<10))            #S00014
I=with(DATA,which(COMMON_NAME=="Smooth Hammerhead" & TL>190 & FL<130))#C00005 
J=with(DATA,which(COMMON_NAME=="Spot tail Shark (Whaler)" & TL>800))#R00393 (TL=114, no 1014)         
K=with(DATA,which(COMMON_NAME=="Spur Dog" & TL>70 & FL<60))        #C00039
L=with(DATA,which(COMMON_NAME=="Thickskin Shark, Sandbar Shark" & TL<40))#13 from R (8) and S (5)
M=with(DATA,which(COMMON_NAME=="Thickskin Shark, Sandbar Shark" & FL<20))#4 from R (3) and S (1)
N=with(DATA,which(COMMON_NAME=="Tiger Shark" & TL<20))             #22 from R(21) and Q (1)
O=with(DATA,which(COMMON_NAME=="Tiger Shark" & TL>360 & FL<250))   #R00127
P=with(DATA,which(COMMON_NAME=="Whiskery Shark" & TL>155))         #X00061
Q=with(DATA,which(COMMON_NAME=="Whiskery Shark" & TL>140 & FL<115))#X00061
R=with(DATA,which(COMMON_NAME=="Gummy Shark" & TL>165 & FL<115))   #X00020
S=with(DATA,which(COMMON_NAME=="Milk Shark" & TL>85 & FL<63))      #N00272
T=with(DATA,which(COMMON_NAME=="Port Jackson" & TL>90 & FL<70))    #R00352
U=with(DATA,which(COMMON_NAME=="Western Wobbegong" & FL<10))       #S00030

DATA=DATA[-c(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U),] #n=57 removed

#Analysis

#List with robust linear regression models and r2 for species with both FL and TL measures  
shk.ray.sp=unique(subset(DATA,NATURE%in%c("R","S"))$COMMON_NAME)#only chondrichtyans
Morphom=subset(DATA,COMMON_NAME%in%shk.ray.sp & !is.na(FL) & !is.na(TL))#the ones with FL & TL measures at the same time
shk.ray.sp.Morph=sort(unique(Morphom$COMMON_NAME))#the species to apply function

OUT=vector('list',length(shk.ray.sp.Morph)) #dummie
names(OUT)=shk.ray.sp.Morph
for(s in 1:length(shk.ray.sp.Morph))
{
  OUT[[s]]=fn.get.TL(DAT=Morphom,SP=shk.ray.sp.Morph[s],Show.Text="NO",N=10,By.sex="NO")
}

#Summarize coefficients in a table
Tab=vector("list",length(shk.ray.sp.Morph))
for(s in 1:length(shk.ray.sp.Morph))
{
  if(length(OUT[[s]])>1)
    {
      a=as.data.frame(t(c(coef(OUT[[s]]$Model_FL.to.TL),OUT[[s]]$R2_FL.to.TL,OUT[[s]]$n)))
    }else
    {
      a=as.data.frame(t(rep(NA,4)))
    }
  colnames(a)=c("Intercept","slope","R2",'n')
  Tab[[s]]=a
}
Tab=do.call(rbind,Tab)
Tab=cbind(Species=shk.ray.sp.Morph,Tab)
Tab[,2:5]=round(Tab[,2:5],3)
Tab$Species=as.character(Tab$Species)

#Remove NA coefficients
Tab=Tab[!is.na(Tab$slope),]

#Remove wobbies and six gill
Wobbies=c("Wobbegong (General)","Cobbler Wobbegong","Banded Wobbegong",
          "Wobbegong (Floral Banded)","Spotted Wobbegong","Western Wobbegong")
Tab=subset(Tab,!Species%in%c(Wobbies,"Big Eye Six Gill Shark"))

#Add scientific name to table
Tab=merge(Tab,SPECIES.NAMES[,2:3],by.x="Species",by.y="COMMON_NAME")


#Output

#Export table
write.csv(Tab,"Coefficents table.csv",row.names=FALSE)



#Plotting

#Figure 1

tiff(file="Figure1.tiff",width=1600,height=2400,units="px",res=300,compression="lzw")

par(mfcol=c(2,1),mai=c(.4,.8,.1,.1),oma=c(1.2,.6,.1,.1),mgp=c(2.5,.5,0),las=1)

#A) Wobbies
Wobs=subset(Morphom,COMMON_NAME%in%Wobbies)
Wobs$CL=with(Wobs,ifelse(SPECIES=="WC","grey40",ifelse(SPECIES=="WD","black",
          ifelse(SPECIES=="WS","grey80",ifelse(SPECIES=="WW",'white',NA)))))
Wobs=Wobs[order(Wobs$COMMON_NAME),]
Legs=c("Banded wobbegong","Cobbler wobbegong","Spotted wobbegong","Western wobbegong")

plot(Wobs$TL,Wobs$FL,pch=21,col=1,bg=Wobs$CL,cex=1.5,ylab="",xlab="",cex.axis=1.5)
legend("topleft",Legs,col=rep(1,4),pt.bg=unique(Wobs$CL),bty='n',pch=21,pt.cex=1.15,cex=0.9)

#B) Six gill
test=subset(Morphom,COMMON_NAME=="Big Eye Six Gill Shark")
test.1=subset(test,YEAR==2003)
test.2=subset(test,YEAR==2004)
plot(test.1$TL,test.1$FL,cex=1.5,ylab="",xlab="",cex.axis=1.5)
points(test.2$TL,test.2$FL,pch=19,cex=1.5)
legend("topleft",c("Observer A","Observer B"),col="black",pt.bg=c("white","black"),bty='n',pch=21,pt.cex=1.15,cex=0.9)

mtext("Fork length (cm)",2,-1.5,outer=T,cex=1.75,las=3)
mtext("Total length (cm)",1,0,outer=T,cex=1.75)
dev.off()



