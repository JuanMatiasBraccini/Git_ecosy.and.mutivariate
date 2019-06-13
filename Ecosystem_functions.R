
########### SCRIPT FOR COMPUTING ECOSYSTEM ANALYSIS OF FISHING IMPACTS ###################


#Reference: Hall & Wise (2011).

# "fun.plt.Indx"  --> for plotting indeces
# "fn.plt"        --> completing missing years for plots
# "Stats"         --> stats summary
# "Boot.eco.indx" --> to generate uncertainty in eco indeces through boostraping
# "Eco.ind.shot"  --> calculates Mean trophic level (MTL) and Mean maximum length (MML)
# "Div.ind.shot"  --> calculates Shannon's and Simpson's diversity, Pielou's evenness and Margalef's richness


require(vegan)
require(BEQI2)
require(asbio)


#--FUNCTIONS---

#Function for plotting
fun.plt.Indx = function(var.names, Indx, YLAB, XLAB, LABELS, YLIM, XLIM = NULL)       #function for plotting 
{
  plot(var.names, Indx, ylim = YLIM, xlim = XLIM, ylab = YLAB, xlab = XLAB, type = "o", xaxt = "n", pch = 20, cex = 2, lwd = 1.5, lty = 2, col = 1) #temporal trend
  axis(side = 1, at = var.names, labels = LABELS, las = 3)    
}

#Function for calculating Diversity indices
Div.ind.shot = function(data)
{
  dat.pos=subset(data,data>0)
  sp_names = names(dat.pos) 
  
  #Shannon-Wiener diversity index
  Shannon = unname(diversity(data, index = "shannon"))
  Shannon[Shannon == 0] = NA
      
  #Simpson's index
  Simpson = unname(diversity(data, index = "simpson"))
  Simpson[Simpson == 1] = NA
  
  #Margalef species richness
  Margalef= margalef(sp_names, dat.pos)
  Margalef[Margalef == 0] = NA
  
  #Pielou's evenness                        
  Evenness = function(H, S) {J = H / log(S)}  #where H = Shannon index for each community; S = Number of species
  Pielou = Evenness(H = Shannon, S = length(sp_names))
  Pielou[Pielou == 0] = NA
    
  #Summary
  return(list(Shannon = round(Shannon, 2), Margalef = round(Margalef, 2), Pielou = round(Pielou, 2), Simpson = round(Simpson, 2)))
}

#Function for calculating Ecosystem indicators
Eco.ind.shot = function(data, TROPHIC.LEVEL, MAX.BODY.LENGTH, TROPHIC.EFFICIENCY = 0.1)
{
  Trophic_level=Max_length=NA
  
  if(length(data)>1)
  {
    dat.pos=subset(data,data>0)
    sp_names = names(dat.pos) 
    
    ddd=data.frame(SPECIES=sp_names,Number=as.vector(data))
    TROPHIC.LEVEL=merge(TROPHIC.LEVEL,ddd,by="SPECIES")
   
    #Mean Trophic Level
    Trophic_level=with(TROPHIC.LEVEL,(sum(Number*TROPHIC_LEVEL)/sum(Number)))
    Trophic_level=round(Trophic_level, 3)
    
    #Mean Maximum Length
    if(!is.null(nrow(MAX.BODY.LENGTH)))
    {
      MAX.BODY.LENGTH=merge(MAX.BODY.LENGTH,ddd,by="SPECIES")
      Max_length=with(MAX.BODY.LENGTH,(sum(Number*FL)/sum(Number)))
      Max_length=round(Max_length, 3)
    }
  }
 
  return(list(MTL =Trophic_level , MML = Max_length))
}

#Bootstrapping for uncertainty in ecosystem indices
Boot.eco.indx=function(Boot.what,D)
{
  #loop over years 
  Yrs=sort(as.numeric(as.character(unique(D$YEAR))))
  
  dummy=D[!duplicated(D$SHEET_NO),]
  ShtN=sort(table(dummy$YEAR))
  if(Boot.what=="MIN")N.boot=ShtN[1]
  if(Boot.what=="MAX")N.boot=ShtN[length(ShtN)]
  
  boot.yr=vector('list',length(Yrs))
  names(boot.yr)=Yrs
  for(y in 1:length(Yrs))
  {
    d=subset(D,YEAR==Yrs[y])
    ShNy=unique(d$SHEET_NO)
    
    #boostrap
    bot.dat=vector('list',N.boots)
    names(bot.dat)=1:N.boots
    for(b in 1:N.boots)
    {
      ID=sample(ShNy, N.boot, replace = T)
      ID.list=vector('list',length(ID))
      for(x in 1:N.boot) ID.list[[x]]=subset(d,SHEET_NO==ID[x])
      bot.dat[[b]]=do.call(rbind,ID.list)
    }
    boot.yr[[y]]=bot.dat
  }
  return(boot.yr)
}

#Summary stats
Stats=function(Yr,DD,lw,up)
{
  MEAN=colMeans(DD,na.rm=T)
  LOW=apply(DD, 2, function(x) quantile(x, lw,na.rm=T))
  UP=apply(DD, 2, function(x) quantile(x, up,na.rm=T))
  names(MEAN)=names(LOW)=names(UP)=Yr
  return(list(MEAN=MEAN,LOW=LOW,UP=UP))
}

#Function for filling in missing year
fn.plt=function(d)
{
  #add missing year
  for(y in 1:length(d)) 
  {
    if(length(names(d[[y]]$MEAN))<length(YEARS))
    {
      id=which(!YEARS%in%names(d[[y]]$MEAN))
      addY=YEARS[id]
      dummy=rep(NA,length(addY))
      names(dummy)=addY
      mn=c(d[[y]]$MEAN,dummy)
      lw1=c(d[[y]]$LOW,dummy)
      up1=c(d[[y]]$UP,dummy)
      id=sort(names(mn))
      d[[y]]$MEAN=mn[id]
      d[[y]]$LOW=lw1[id]
      d[[y]]$UP=up1[id]
    }
  }
}








####################################################################
####################################################################
#Old stuff
# 
# #Function for calculating Diversity indices
# Div.ind = function(data)
# {
#   var_names = row.names(data)  
#   N = length(var_names)
#   sp_names = dimnames(data)[[2]]
#   
#   #Shannon-Wiener diversity index
#   Shannon = unname(diversity(data, index = "shannon"))
#   Shannon[Shannon == 0] = NA
#   
#   #Simpson's index
#   Simpson = unname(diversity(data, index = "simpson"))
#   Simpson[Simpson == 1] = NA
#   
#   #Margalef species richness
#   Margalef = rep(NA, N) 
#   for(i in 1:N) Margalef[i] = margalef(sp_names, data[i, ])
#   Margalef[Margalef == 0] = NA
#   
#   #Pielou's evenness                        
#   Evenness = function(H, S) {J = H / log(S)}  #where H = Shannon index for each community; S = Number of species
#   Pielou = rep(NA, N)
#   for(i in 1:N) Pielou[i] = Evenness(H = Shannon[[i]], S = length(sp_names))
#   Pielou[Pielou == 0] = NA
#   
#   #Summary
#   return(list(Shannon = round(Shannon, 2), Margalef = round(Margalef, 2), Pielou = round(Pielou, 2), Simpson = round(Simpson, 2), Var_names = var_names))
# }
# 
# 
# #Function for calculating Ecosystem indicators
# Eco.ind = function(data, TROPHIC.LEVEL, MAX.BODY.LENGTH, TROPHIC.EFFICIENCY = 0.1)
# {
#   var_names = row.names(data)  
#   N = length(var_names)
#   sp_names = dimnames(data)[[2]]
#   data[data == 0] = NA
#   
#   TROPHIC.LEVEL=subset(TROPHIC.LEVEL,SPECIES%in%colnames(data))
#   TROPHIC.LEVEL=TROPHIC.LEVEL[order(TROPHIC.LEVEL$SPECIES),]
#   
#   #Mean Trophic Level
#   Trophic_level = rep(NA,N)
#   for(i in 1:N)
#   {
#     a = data[i, ] * TROPHIC.LEVEL[, 2]
#     Trophic_level[i] =  sum(a, na.rm = T) / sum(data[i, ], na.rm = T)
#   }
#   
#   #Mean Maximum Length
#   Max_length = rep(NA,N)
#   a = data[,match(colnames(MAX.BODY.LENGTH),colnames(data))] * MAX.BODY.LENGTH
#   for(i in 1:N)
#   {
#     Max_length[i] = sum(a[i,], na.rm = T) / sum(data[i, ], na.rm = T)
#   }
#   
#   #Fishery in Balance
#   FIB = rep(NA,N)
#   for(i in 1:N) if(sum(data[i,], na.rm = T) > 0) {break}
#   
#   for(j in 1:N)
#   {
#     FIB[j] = log(sum(data[j,], na.rm = T) * (1 / TROPHIC.EFFICIENCY) ^ Trophic_level[j]) - log(sum(data[i,], na.rm = T) * (1 / TROPHIC.EFFICIENCY) ^ Trophic_level[i])
#   }
#   
#   #Summary
#   return(list(MTL = round(Trophic_level, 2), MML = round(Max_length, 2), FIB = round(FIB, 2)))
# }
# 
# 
