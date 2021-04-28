

library(tidyverse)
library(dplyr)
library("readxl")


# Data section ------------------------------------------------------------

#SSF 2007-08
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')
setwd(handl_OneDrive('Data/SSF_survey_07_08'))
F1_Sampling <- read_excel("SharkSurveyData_30_09_2008.xls", sheet = "F1_Sampling")
F1_SamplingTwo <- read_excel("SharkSurveyData_30_09_2008.xls", sheet = "F1_SamplingTwo")
F3_Sampling <- read_excel("SharkSurveyData_30_09_2008.xls", sheet = "F3_Sampling")



#Parks Australia


# Procedure section ------------------------------------------------------------
SSF=F3_Sampling%>%left_join(F1_Sampling%>%
                              mutate(Station=as.numeric(Station))%>%
                              dplyr::select(VesselName,Cruise,Station,Year,Month,Day),
                            by=c('Cruise','Station'))%>%
  left_join(F1_SamplingTwo%>%
              distinct(Cruise,Station,.keep_all=T)%>%
              mutate(Station=as.numeric(Station))%>%
              dplyr::select(Cruise,Station,Mesh,TimeSet,StartLat,StartLong,EndLat,EndLong,TimeHaul,TimeUp,MaxDepth),
            by=c('Cruise','Station'))