#Script for brining in MAFFRI shark survey and put in right format

library(RODBC)  			#include ODBC library for importing Acccess data

handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')

setwd(handl_OneDrive("Data/Shark survey/2007-2008"))
channel <- odbcConnectExcel("SharkSurveyData_30_09_2008")
Net_specs<- sqlFetch(channel,"Net_specification")
F1_Sampling<- sqlFetch(channel,"F1_Sampling")
F1_SamplingTwo<- sqlFetch(channel,"F1_SamplingTwo")
F2_Sampling<- sqlFetch(channel,"F2_Sampling")
TEPS<- sqlFetch(channel,"F3_Sampling")
close(channel)


#Combine data
names(Net_specs)[1]="Fisher"
names(Net_specs)[5]="VesselName"
names(Net_specs)[15]="Net_length"

Net_specs=subset(Net_specs,select=c(VesselName,Net_length))
Net_specs$VesselName=as.character(Net_specs$VesselName)
Net_specs$Net_length=as.numeric(as.character(Net_specs$Net_length))



F1_SamplingTwo$StartLat=with(F1_SamplingTwo,ifelse(is.na(StartLat),EndLat,StartLat))  
F1_SamplingTwo$StartLong=with(F1_SamplingTwo,ifelse(is.na(StartLong),EndLong,StartLong))     

Net_specs=Net_specs[!duplicated(Net_specs$VesselName),]
F1_Sampling$VesselName=as.character(F1_Sampling$VesselName)
F1_Sampling=subset(F1_Sampling,select=c(Cruise, Station, Year, Month, Day,VesselName))
F1_SamplingTwo=subset(F1_SamplingTwo,select=c(Cruise, Station,Mesh,StartLat, StartLong,MaxDepth,NetMaking))
F1_Sampling=merge(F1_Sampling,Net_specs,by="VesselName",all.x=T)
F1_Sampling=merge(F1_SamplingTwo,F1_Sampling,by=c("Cruise","Station"),all.x=T)
F2_Sampling=subset(F2_Sampling,select=c(Cruise, Station, Mesh1, Species, Sex,Ret,
                                        Length,LengthType))

F2_Sampling$Length=with(F2_Sampling,ifelse(Length==13690,1369,
        ifelse(Length==4250,425,ifelse(Length==42,420,ifelse(Length==168,1680,Length)))))
F2_Sampling$Length=F2_Sampling$Length/10

KeepMesh=c("S6","C6.00","C6.01","C6.02","C6.04","C6.51","C6.52","C6.53","C6.54")
F2_Sampling=subset(F2_Sampling,Mesh1%in%KeepMesh)
names(F2_Sampling)[3]="Mesh"

F1_Sampling$Dummy=with(F1_Sampling,paste(Cruise,Station,Mesh))
F2_Sampling$Dummy=with(F2_Sampling,paste(Cruise,Station,Mesh))
TEPS$Dummy=with(TEPS,paste(Cruise,Station,Mesh))
F1_Sampling=subset(F1_Sampling,Dummy%in%F2_Sampling$Dummy)
TEPS=subset(TEPS,Dummy%in%F2_Sampling$Dummy)

F2_Sampling$Mesh=as.character(F2_Sampling$Mesh)
F1_Sampling$Mesh=as.character(F1_Sampling$Mesh)
TEPS$Mesh=as.character(TEPS$Mesh)

F2_Sampling$Mesh=with(F2_Sampling,
            ifelse(Mesh%in%c("C6.00","C6.01","C6.02","C6.03","C6.04"),"6",
            ifelse(Mesh%in%c("C6.51","C6.52","C6.53","C6.54"),"6.5",Mesh)))
F1_Sampling$Mesh=with(F1_Sampling,
            ifelse(Mesh%in%c("C6.00","C6.01","C6.02","C6.03","C6.04"),"6",
           ifelse(Mesh%in%c("C6.51","C6.52","C6.53","C6.54"),"6.5",Mesh)))
TEPS$Mesh=with(TEPS,
          ifelse(Mesh%in%c("C6.00","C6.01","C6.02","C6.03","C6.04"),"6",
         ifelse(Mesh%in%c("C6.51","C6.52","C6.53","C6.54"),"6.5",Mesh)))

DATA_MAFFRI=merge(F2_Sampling,F1_Sampling,by=c("Dummy","Cruise","Station","Mesh"),all.x=T)
DATA_TEPS_MAFFRI=merge(TEPS,F1_Sampling,by=c("Dummy","Cruise","Station","Mesh"),all.x=T)


#export data to match WA format
setwd(handl_OneDrive("Analyses/Bycatch"))
write.csv(DATA_MAFFRI,"MAFFRI_data.csv",row.names=F)
write.csv(DATA_TEPS_MAFFRI,"MAFFRI_TEPS_data.csv",row.names=F)