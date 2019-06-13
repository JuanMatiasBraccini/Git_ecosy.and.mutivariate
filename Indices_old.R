#------------------- Diversity and Ecosystem-based indicators --------------------

#note: compare using just Numbers VS log(Numbers/effort), which standardises for using different efforts

library(vegan)


#1. Indices reported in Table 1 page 7-8 Hall & Wise 2011

communityI <- c(100, 100, 100, 100, 100,0,0,0)
communityII <- c(50, 60, 75, 90, 65,1,5,8)
communityIII <- c(30, 40, 35, 90, 60,10,20,50)

names(communityI)=paste("Sp.",1:length(communityI),sep='')
names(communityII)=names(communityIII)=names(communityI)
communities <- rbind(communityI,communityII,communityIII)


#Shannon-Wiener diversity index
SH_W=diversity(communities, index = "shannon")
Yrs=1:3
regression.result <- lm(SH_W ~ Yrs)

plot(Yrs,SH_W)  #annual trend
abline(regression.result, lty = 2)


#Margalef species richness
  #note: N = total number of individuals; S = the number of species 
Margalef=function(N,S) d=(S-1)/(log(N))

a=Margalef(N=sum(communityI),S=length(communityI))


#Evenness
  #note: S = number of species
Evenness=function(H_tilda,S) J=H_tilda/log(S)

a=Evenness(H_tilda=SH_W[1],S=length(communityI))


#Simpson
Simpson=diversity(communities, index = "simpson")


#Mean trophic level
Mean.TL=function(B,TL) sum(B*TL)/sum(B)



#Fishery in Balance



#Mean maximum length





#2. Other indices. rpackage Vegan


#Temporal MDS plots of species composition
#should apply transformation and calculate Bray curtis? see page 11 Hall & Wise 2011