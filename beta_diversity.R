# source:      https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html
#vip resource!!
library(vegan)
load("C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Desktop/Argentina_Ecosystem effects paper/clean_data.Rdata")

sample_data$alpha <- diversity(obj$data$otu_rarefied[, sample_data$SampleID],
                               MARGIN = 2,
                               index = "invsimpson")
hist(sample_data$alpha)

# source: https://rdrr.io/rforge/vegan/man/betadiver.html

## Raw data and plotting
data(sipoo)
m <- betadiver(sipoo)
plot(m)
## The indices
betadiver(help=TRUE)
## The basic Whittaker index
d <- betadiver(sipoo, "w")
## This should be equal to Sorensen index (binary Bray-Curtis in
## vegan)
range(d - vegdist(sipoo, binary=TRUE))



#Source: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2012.00224.x
library(betapart)


data(ceram.s)
ceram.core.s<-betapart.core(ceram.s)
ceram.dist.jac<-beta.pair(ceram.core.s, index.family="jac")
ceram.dist.sor<-beta.pair(ceram.core.s, index.family="sor")
ceram.multi.jac<-beta.multi(ceram.core.s, index.family="jac")
ceram.multi.sor<-beta.multi(ceram.core.s, index.family="sor")




# source: https://www.youtube.com/watch?v=3ySrZFkNsU8
library(ade4)
library(adespatial)
data("doubs")
fish<-doubs$fish
fish<-fish[-8,]

#calculate beta diversity for jaccard
fish.bd.j<-beta.div.comp(fish,coef='J',quant=T)
fish.bd.j$part

#calculate beta diversity for sorensen
fish.bd.s<-beta.div.comp(fish,coef='S',quant=T)

#calculate local contriution to beta diversity
local.repl<-LCBD.comp(fish.bd.j$repl,sqrt.D = T)
