#pipeline for pairwise comparisons of methylation in the methylKit package and creation of figures
#Objects feeding into this pipeline were 

#install and load packages
BiocManager::install("genomation")
BiocManager::install("methylKit")
BiocManager::install("GenomicRanges")
BiocManager::install('bsseq')
install.packages("UpSetR")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
install.packages("cowplot")
install.packages("MASS")
install.packages("reshape2")
install.packages("dplyr")
install.packages("scales")
install.packages("matrixStats")
install.packages("ggfortify")
install.packages("ggplot2")
install.packages("cowplot")

library(genomation)
library(methylKit)
library(GenomicRanges)
library(bsseq)
library(UpSetR)
library(ComplexHeatmap)
library(MASS)
library(reshape2)
library(dplyr)
library(scales)
library(matrixStats)
library(ggfortify)
library(ggplot2)
library(cowplot)


#Load objects
Final.all.cov1 <- readRDS("Final.all.cov1")
Final.all.cov3 <- readRDS("Final.all.cov3")
Final.all.cov5 <- readRDS("Final.all.cov5")
Final.all.cov10 <- readRDS("Final.all.cov10")

#tiled by 1000bp
Final.all.cov3.tile1000 <-readRDS("Final.all.cov3.tile1000")
Final.all.cov5.tile1000 <-readRDS("Final.all.cov5.tile1000")
Final.all.cov10.tile1000 <-readRDS("Final.all.cov10.tile1000")

#filtered 
Final.all.cov3.tile1000.filt <-readRDS("Final.all.cov3.tile1000.filt")
Final.all.cov5.tile1000.filt <-readRDS("Final.all.cov5.tile1000.filt")
Final.all.cov10.tile1000.filt <-readRDS("Final.all.cov10.tile1000.filt")

#filtered and normalized
Final.all.cov3.tile1000.filt.nor <-readRDS("Final.all.cov3.tile1000.filt.nor")
Final.all.cov5.tile1000.filt.nor <-readRDS("Final.all.cov5.tile1000.filt.nor")
Final.all.cov10.tile1000.filt.nor <-readRDS("Final.all.cov10.tile1000.filt.nor")

#final united objects
Final.all.cov3.tile1000.filt.nor.unite <-readRDS("Final.all.cov3.tile1000.filt.nor.unite")
Final.all.cov5.tile1000.filt.nor.unite <-readRDS("Final.all.cov5.tile1000.filt.nor.unite")
Final.all.cov10.tile1000.filt.nor.unite <-readRDS("Final.all.cov10.tile1000.filt.nor.unite")

#model to assess sex effect that is filtered+normalized+united
Final.all.sexeffect.cov3.tile1000.filt.nor.unite <- readRDS("C:/Users/BOKVISTJ/Documents/UofC MSc/All Analysis December 2021/R Analysis All Sequencing Data/Final.all.sexeffect.cov3.tile1000.filt.nor.unite")

#model to assess family effect that is filtered+normalized+united
Final.all.3familyeffect.cov3.tile1000.filt.nor.unit<- readRDS("C:/Users/BOKVISTJ/Documents/UofC MSc/All Analysis December 2021/R Analysis All Sequencing Data/Final.all.3familyeffect.cov3.tile1000.filt.nor.unite")

#------------REATINED CPG SITES AND TILES--------------
nrows.Final.all.cov1 <- nrow(Final.all.cov1[[1]]) #4883950
nrows.Final.all.cov3 <-nrow(Final.all.cov3[[1]]) # 3768323
nrows.Final.all.cov5 <-nrow(Final.all.cov5[[1]]) #2901040
nrows.Final.all.cov10 <-nrow(Final.all.cov10[[1]]) #1463323
nrows.Final.all.cov3.tile1000 <-nrow(Final.all.cov3.tile1000[[1]]) #312723
nrows.Final.all.cov5.tile1000 <-nrow(Final.all.cov5.tile1000[[1]]) #258390
nrows.Final.all.cov10.tile1000 <-nrow(Final.all.cov10.tile1000[[1]]) #155783
nrows.Final.all.cov3.tile1000.filt <-nrow(Final.all.cov3.tile1000.filt[[1]]) #282814
nrows.Final.all.cov5.tile1000.filt <-nrow(Final.all.cov5.tile1000.filt[[1]]) #250487
nrows.Final.all.cov10.tile1000.filt <-nrow(Final.all.cov10.tile1000.filt[[1]]) #155622
nrows.Final.all.cov3.tile1000.filt.nor <-nrow(Final.all.cov3.tile1000.filt.nor[[1]]) # 282814
nrows.Final.all.cov5.tile1000.filt.nor <-nrow(Final.all.cov5.tile1000.filt.nor[[1]])#250487
nrows.Final.all.cov10.tile1000.filt.nor <-nrow(Final.all.cov10.tile1000.filt.nor[[1]]) #  155622


#tiles retained that are present in minimum 5 per treatment/lifestage coho (6), where there is at least 3 CpG per tile and 20x coverage overall
nrows.Final.all.cov3.tile1000.filt.nor.unite <-nrow(Final.all.cov3.tile1000.filt.nor.unite) #262268 
nrows.Final.all.cov5.tile1000.filt.nor.unite <-nrow(Final.all.cov5.tile1000.filt.nor.unite) #232366
nrows.Final.all.cov10.tile1000.filt.nor.unite <-nrow(Final.all.cov10.tile1000.filt.nor.unite) #148429


#---------------PAIRWISE CONTRASTS----------------------------
#subset the united file to subset adults and smolts
#compare at the smolt stage (ENR to WD, ENR to TD, TD to WD)
#compare at the adult stage (ENR to WD, ENR to TD, TD to WD)
#compare smotls to adults within a treatment 


#### Contrast smolts enriched vs wild - wild coded as 1
Smolts.ENRCvsWD.Final.all.cov3.tile1000.filt.nor.unite <- reorganize(Final.all.cov3.tile1000.filt.nor.unite, sample.ids=c("G3846",
                                                                                                                        "G3847",
                                                                                                                        "G3848",
                                                                                                                        "G3849",
                                                                                                                        "G3850",
                                                                                                                        "G3851",
                                                                                                                        "G3852",
                                                                                                                        "G3853",
                                                                                                                        "G3854",
                                                                                                                        "G3855",
                                                                                                                        "G3866",
                                                                                                                        "G3867",
                                                                                                                        "G3868",
                                                                                                                        "G3869",
                                                                                                                        "G3870",
                                                                                                                        "G3871",
                                                                                                                        "G3872",
                                                                                                                        "G3873",
                                                                                                                        "G3874",
                                                                                                                        "G3875",
                                                                                                                        "H656",
                                                                                                                        "H657",
                                                                                                                        "H658",
                                                                                                                        "H659",
                                                                                                                        "H660",
                                                                                                                        "H661",
                                                                                                                        "H662",
                                                                                                                        "H663",
                                                                                                                        "H664",
                                                                                                                        "H665",
                                                                                                                        "H676",
                                                                                                                        "H677",
                                                                                                                        "H678",
                                                                                                                        "H679",
                                                                                                                        "H680",
                                                                                                                        "H681",
                                                                                                                        "H682",
                                                                                                                        "H683",
                                                                                                                        "H684",
                                                                                                                        "H685")
                                                                   , treatment=c(2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1))

#covariate for this contrast
Smolts.ENRCvsWD.Covariate=data.frame(family=c(2,
                                              2,
                                              1,
                                              2,
                                              1,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              1,
                                              2,
                                              2,
                                              1,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2),
                                   sex=c(1,
                                         2,
                                         1,
                                         2,
                                         2,
                                         1,
                                         1,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         1,
                                         1,
                                         2,
                                         2,
                                         2,
                                         2,
                                         1,
                                         2,
                                         1,
                                         2,
                                         2,
                                         2,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         2,
                                         1,
                                         2))





#####contrast Smolts ENRCrich vs Traditional - Enriched coded as 1
Smolts.ENRCvsTD.Final.all.cov3.tile1000.filt.nor.unite <- reorganize(Final.all.cov3.tile1000.filt.nor.unite, sample.ids=c("G3846",
                                                                                                                        "G3847",
                                                                                                                        "G3848",
                                                                                                                        "G3849",
                                                                                                                        "G3850",
                                                                                                                        "G3851",
                                                                                                                        "G3852",
                                                                                                                        "G3853",
                                                                                                                        "G3854",
                                                                                                                        "G3855",
                                                                                                                        "G3856",
                                                                                                                        "G3858",
                                                                                                                        "G3859",
                                                                                                                        "G3860",
                                                                                                                        "G3861",
                                                                                                                        "G3862",
                                                                                                                        "G3863",
                                                                                                                        "G3864",
                                                                                                                        "G3865",
                                                                                                                        "H656",
                                                                                                                        "H657",
                                                                                                                        "H658",
                                                                                                                        "H659",
                                                                                                                        "H660",
                                                                                                                        "H661",
                                                                                                                        "H662",
                                                                                                                        "H663",
                                                                                                                        "H664",
                                                                                                                        "H665",
                                                                                                                        "H666",
                                                                                                                        "H667",
                                                                                                                        "H668",
                                                                                                                        "H669",
                                                                                                                        "H670",
                                                                                                                        "H671",
                                                                                                                        "H672",
                                                                                                                        "H673",
                                                                                                                        "H674",
                                                                                                                        "H675")
                                                                   , treatment=c(1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2))

#covariate for this contrast 
Smolts.ENRCvsTD.Covariate=data.frame(family=c(2,
                                              2,
                                              1,
                                              2,
                                              1,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              3,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              1,
                                              2,
                                              2,
                                              1,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2),
                                   sex=c(1,
                                         2,
                                         1,
                                         2,
                                         2,
                                         1,
                                         1,
                                         2,
                                         2,
                                         2,
                                         1,
                                         2,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         2,
                                         1,
                                         2,
                                         2,
                                         2,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         2,
                                         2,
                                         1,
                                         2,
                                         2,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2))


#####contrast Smolts wild vs Traditional - wild coded as 1
#no family covariate in this contrast as only family structure in Enriched smolts
Smolts.WDvsTD.Final.all.cov3.tile1000.filt.nor.unite <- reorganize(Final.all.cov3.tile1000.filt.nor.unite, sample.ids=c("G3856",
                                                                                                                        "G3858",
                                                                                                                        "G3859",
                                                                                                                        "G3860",
                                                                                                                        "G3861",
                                                                                                                        "G3862",
                                                                                                                        "G3863",
                                                                                                                        "G3864",
                                                                                                                        "G3865",
                                                                                                                        "G3866",
                                                                                                                        "G3867",
                                                                                                                        "G3868",
                                                                                                                        "G3869",
                                                                                                                        "G3870",
                                                                                                                        "G3871",
                                                                                                                        "G3872",
                                                                                                                        "G3873",
                                                                                                                        "G3874",
                                                                                                                        "G3875",
                                                                                                                        "H666",
                                                                                                                        "H667",
                                                                                                                        "H668",
                                                                                                                        "H669",
                                                                                                                        "H670",
                                                                                                                        "H671",
                                                                                                                        "H672",
                                                                                                                        "H673",
                                                                                                                        "H674",
                                                                                                                        "H675",
                                                                                                                        "H676",
                                                                                                                        "H677",
                                                                                                                        "H678",
                                                                                                                        "H679",
                                                                                                                        "H680",
                                                                                                                        "H681",
                                                                                                                        "H682",
                                                                                                                        "H683",
                                                                                                                        "H684",
                                                                                                                        "H685")
                                                                   , treatment=c(2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1))


# #covariate for this contrast
Smolts.WDvsTD.Covariate=data.frame(family=c(1,
                                            1,
                                            2,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1), sex=c(1,
                                         2,
                                         2,
                                         1,
                                         1,
                                         2,
                                         1,
                                         1,
                                         2,
                                         2,
                                         2,
                                         2,
                                         1,
                                         1,
                                         2,
                                         2,
                                         2,
                                         2,
                                         1,
                                         2,
                                         2,
                                         1,
                                         2,
                                         2,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         1,
                                         2,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         2,
                                         1,
                                         2))  






#CONTRASTING treatmentS AT THE ADULT STAGE

######contrast adults Enriched vs wild - wild coded as 1
Adults.ENRCvsWD.Final.all.cov3.tile1000.filt.nor.unite <- reorganize(Final.all.cov3.tile1000.filt.nor.unite, sample.ids=c("F159",
                                                                                                                        "F165",
                                                                                                                        "F166",
                                                                                                                        "F167",
                                                                                                                        "F33",
                                                                                                                        "F39",
                                                                                                                        "F43",
                                                                                                                        "F48",
                                                                                                                        "F73",
                                                                                                                        "M159",
                                                                                                                        "M161",
                                                                                                                        "M162",
                                                                                                                        "M164",
                                                                                                                        "M166",
                                                                                                                        "M167",
                                                                                                                        "M169",
                                                                                                                        "M34",
                                                                                                                        "M36",
                                                                                                                        "M38",
                                                                                                                        "M63",
                                                                                                                        "M71",
                                                                                                                        "M73")
                                                                   , treatment=c(2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 2,
                                                                                 1,
                                                                                 2,
                                                                                 2))

# #covariate for this contrast 
Adults.ENRCvsWD.Covariate=data.frame(family=c(1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              2,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1), sex=c(1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2)) 



######contrast adults Enriched vs traditional - Enriched coded as 1
Adults.ENRCvsTD.Final.all.cov3.tile1000.filt.nor.unite <- reorganize(Final.all.cov3.tile1000.filt.nor.unite, sample.ids=c("F159",
                                                                                                                        "F161",
                                                                                                                        "F163",
                                                                                                                        "F165",
                                                                                                                        "F166",
                                                                                                                        "F167",
                                                                                                                        "F33",
                                                                                                                        "F35",
                                                                                                                        "F39",
                                                                                                                        "F40",
                                                                                                                        "F41",
                                                                                                                        "F43",
                                                                                                                        "F45",
                                                                                                                        "F48",
                                                                                                                        "F49",
                                                                                                                        "F50",
                                                                                                                        "M159",
                                                                                                                        "M161",
                                                                                                                        "M162",
                                                                                                                        "M169",
                                                                                                                        "M34",
                                                                                                                        "M37",
                                                                                                                        "M38",
                                                                                                                        "M41",
                                                                                                                        "M71",
                                                                                                                        "M72",
                                                                                                                        "M73",
                                                                                                                        "M74",
                                                                                                                        "M77")
                                                                   , treatment=c(1,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 2,
                                                                                 1,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 2,
                                                                                 1,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 2,
                                                                                 1,
                                                                                 2,
                                                                                 1,
                                                                                 2,
                                                                                 1,
                                                                                 2,
                                                                                 2))

# #covariate for this contrast
Adults.ENRCvsTD.Covariate=data.frame(family=c(1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              2,
                                              2,
                                              2,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1), sex=c(1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2))


####contrast adults wild vs traditional - wild coded as 1
Adults.WDvsTD.Final.all.cov3.tile1000.filt.nor.unite <- reorganize(Final.all.cov3.tile1000.filt.nor.unite, sample.ids=c("F161",
                                                                                                                        "F163",
                                                                                                                        "F35",
                                                                                                                        "F40",
                                                                                                                        "F41",
                                                                                                                        "F45",
                                                                                                                        "F49",
                                                                                                                        "F50",
                                                                                                                        "F73",
                                                                                                                        "M164",
                                                                                                                        "M166",
                                                                                                                        "M167",
                                                                                                                        "M36",
                                                                                                                        "M37",
                                                                                                                        "M41",
                                                                                                                        "M63",
                                                                                                                        "M72",
                                                                                                                        "M74",
                                                                                                                        "M77")
                                                                   , treatment=c(2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 1,
                                                                                 2,
                                                                                 2,
                                                                                 1,
                                                                                 2,
                                                                                 2,
                                                                                 2))

# #covariate for this contrast
Adults.WDvsTD.Covariate=data.frame(family=c(1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            2,
                                            2,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1,
                                            1), sex=c(1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2))





#Contrasting treatments between smolts and adults

######contrast enriched smolts vs. adults - smolts coded as 1
Enriched.SMvsAD.Final.all.cov3.tile1000.filt.nor.unite <- reorganize(Final.all.cov3.tile1000.filt.nor.unite, sample.ids=c("F159",
                                                                                                                          "F165",
                                                                                                                          "F166",
                                                                                                                          "F167",
                                                                                                                          "F33",
                                                                                                                          "F39",
                                                                                                                          "F43",
                                                                                                                          "F48",
                                                                                                                          "G3846",
                                                                                                                          "G3847",
                                                                                                                          "G3848",
                                                                                                                          "G3849",
                                                                                                                          "G3850",
                                                                                                                          "G3851",
                                                                                                                          "G3852",
                                                                                                                          "G3853",
                                                                                                                          "G3854",
                                                                                                                          "G3855",
                                                                                                                          "H656",
                                                                                                                          "H657",
                                                                                                                          "H658",
                                                                                                                          "H659",
                                                                                                                          "H660",
                                                                                                                          "H661",
                                                                                                                          "H662",
                                                                                                                          "H663",
                                                                                                                          "H664",
                                                                                                                          "H665",
                                                                                                                          "M159",
                                                                                                                          "M161",
                                                                                                                          "M162",
                                                                                                                          "M169",
                                                                                                                          "M34",
                                                                                                                          "M38",
                                                                                                                          "M71",
                                                                                                                          "M73")
                                                                     , treatment=c(2,
                                                                                   2,
                                                                                   2,
                                                                                   2,
                                                                                   2,
                                                                                   2,
                                                                                   2,
                                                                                   2,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   2,
                                                                                   2,
                                                                                   2,
                                                                                   2,
                                                                                   2,
                                                                                   2,
                                                                                   2,
                                                                                   2))

#covariate for this contrast
Enriched.SMvsAD.Covariate=data.frame(family=c(2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              3,
                                              2,
                                              2,
                                              1,
                                              2,
                                              1,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              1,
                                              2,
                                              2,
                                              1,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2),
                                     sex=c(1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           1,
                                           2,
                                           1,
                                           2,
                                           2,
                                           1,
                                           1,
                                           2,
                                           2,
                                           2,
                                           2,
                                           1,
                                           2,
                                           2,
                                           2,
                                           2,
                                           1,
                                           1,
                                           1,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2,
                                           2))       


######contrast wild smolts vs. adults - smolts coded as 1
Wild.SMvsAD.Final.all.cov3.tile1000.filt.nor.unite <- reorganize(Final.all.cov3.tile1000.filt.nor.unite, sample.ids=c("F73",
                                                                                                                      "G3866",
                                                                                                                      "G3867",
                                                                                                                      "G3868",
                                                                                                                      "G3869",
                                                                                                                      "G3870",
                                                                                                                      "G3871",
                                                                                                                      "G3872",
                                                                                                                      "G3873",
                                                                                                                      "G3874",
                                                                                                                      "G3875",
                                                                                                                      "H676",
                                                                                                                      "H677",
                                                                                                                      "H678",
                                                                                                                      "H679",
                                                                                                                      "H680",
                                                                                                                      "H681",
                                                                                                                      "H682",
                                                                                                                      "H683",
                                                                                                                      "H684",
                                                                                                                      "H685",
                                                                                                                      "M164",
                                                                                                                      "M166",
                                                                                                                      "M167",
                                                                                                                      "M36",
                                                                                                                      "M63")
                                                                 , treatment=c(2,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               1,
                                                                               2,
                                                                               2,
                                                                               2,
                                                                               2,
                                                                               2))


# #covariate for this contrast 
Wild.SMvsAD.Covariate=data.frame(sex=c(1,
                                       2,
                                       2,
                                       2,
                                       1,
                                       1,
                                       2,
                                       2,
                                       2,
                                       2,
                                       1,
                                       1,
                                       2,
                                       1,
                                       1,
                                       1,
                                       1,
                                       2,
                                       2,
                                       1,
                                       2,
                                       2,
                                       2,
                                       2,
                                       2,
                                       2))  



######contrast traditional smolts vs. adults - smolts coded as 1
Traditional.SMvsAD.Final.all.cov3.tile1000.filt.nor.unite <- reorganize(Final.all.cov3.tile1000.filt.nor.unite, sample.ids=c("F161",
                                                                                                                             "F163",
                                                                                                                             "F35",
                                                                                                                             "F40",
                                                                                                                             "F41",
                                                                                                                             "F45",
                                                                                                                             "F49",
                                                                                                                             "F50",
                                                                                                                             "G3856",
                                                                                                                             "G3858",
                                                                                                                             "G3859",
                                                                                                                             "G3860",
                                                                                                                             "G3861",
                                                                                                                             "G3862",
                                                                                                                             "G3863",
                                                                                                                             "G3864",
                                                                                                                             "G3865",
                                                                                                                             "H666",
                                                                                                                             "H667",
                                                                                                                             "H668",
                                                                                                                             "H669",
                                                                                                                             "H670",
                                                                                                                             "H671",
                                                                                                                             "H672",
                                                                                                                             "H673",
                                                                                                                             "H674",
                                                                                                                             "H675",
                                                                                                                             "M37",
                                                                                                                             "M41",
                                                                                                                             "M72",
                                                                                                                             "M74",
                                                                                                                             "M77")
                                                                        , treatment=c(2,
                                                                                      2,
                                                                                      2,
                                                                                      2,
                                                                                      2,
                                                                                      2,
                                                                                      2,
                                                                                      2,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      1,
                                                                                      2,
                                                                                      2,
                                                                                      2,
                                                                                      2,
                                                                                      2))

# #covariate for this contrast
Traditional.SMvsAD.Covariate=data.frame(family=c(1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 2,
                                                 1,
                                                 1,
                                                 2,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1), sex=c(1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              2,
                                              2,
                                              1,
                                              1,
                                              2,
                                              1,
                                              1,
                                              2,
                                              2,
                                              2,
                                              1,
                                              2,
                                              2,
                                              2,
                                              1,
                                              1,
                                              1,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2))  






##CONTRAST TO ASSESS THE OVERALL EFFECT OF FAMILY
#the object Final.all.cov3.tile1000.filt.nor.unite is not used here, a different model is fitted completely to assess the effect of family

#covariate for this contrast
Family.Covariate=data.frame(treatment2=c(2,
                                         3,
                                         3,
                                         2,
                                         2,
                                         2,
                                         2,
                                         3,
                                         2,
                                         3,
                                         3,
                                         2,
                                         3,
                                         2,
                                         3,
                                         3,
                                         1,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         3,
                                         3,
                                         3,
                                         3,
                                         3,
                                         3,
                                         3,
                                         3,
                                         3,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         2,
                                         3,
                                         3,
                                         3,
                                         3,
                                         3,
                                         3,
                                         3,
                                         3,
                                         3,
                                         3,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         2,
                                         2,
                                         2,
                                         1,
                                         1,
                                         1,
                                         2,
                                         2,
                                         1,
                                         3,
                                         2,
                                         3,
                                         1,
                                         2,
                                         3,
                                         2,
                                         3,
                                         3),
                            sex=c(1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  2,
                                  1,
                                  2,
                                  2,
                                  1,
                                  1,
                                  2,
                                  2,
                                  2,
                                  1,
                                  2,
                                  2,
                                  1,
                                  1,
                                  2,
                                  1,
                                  1,
                                  2,
                                  2,
                                  2,
                                  2,
                                  1,
                                  1,
                                  2,
                                  2,
                                  2,
                                  2,
                                  1,
                                  2,
                                  1,
                                  2,
                                  2,
                                  2,
                                  2,
                                  1,
                                  1,
                                  1,
                                  2,
                                  2,
                                  2,
                                  1,
                                  2,
                                  2,
                                  2,
                                  1,
                                  1,
                                  1,
                                  2,
                                  1,
                                  2,
                                  1,
                                  1,
                                  1,
                                  1,
                                  2,
                                  2,
                                  1,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2,
                                  2))



##CONTRAST TO ASSESS THE OVERALL EFFECT OF SEX 
#the object Final.all.cov3.tile1000.filt.nor.unite is not used here, a different model is fitted completely to assess the effect of sex
#covariate for this contrast
Sex.Covariate=data.frame(treatment2=c(2,
                                      3,
                                      3,
                                      2,
                                      2,
                                      2,
                                      2,
                                      3,
                                      2,
                                      3,
                                      3,
                                      2,
                                      3,
                                      2,
                                      3,
                                      3,
                                      1,
                                      2,
                                      2,
                                      2,
                                      2,
                                      2,
                                      2,
                                      2,
                                      2,
                                      2,
                                      2,
                                      3,
                                      3,
                                      3,
                                      3,
                                      3,
                                      3,
                                      3,
                                      3,
                                      3,
                                      1,
                                      1,
                                      1,
                                      1,
                                      1,
                                      1,
                                      1,
                                      1,
                                      1,
                                      1,
                                      2,
                                      2,
                                      2,
                                      2,
                                      2,
                                      2,
                                      2,
                                      2,
                                      2,
                                      2,
                                      3,
                                      3,
                                      3,
                                      3,
                                      3,
                                      3,
                                      3,
                                      3,
                                      3,
                                      3,
                                      1,
                                      1,
                                      1,
                                      1,
                                      1,
                                      1,
                                      1,
                                      1,
                                      1,
                                      1,
                                      2,
                                      2,
                                      2,
                                      1,
                                      1,
                                      1,
                                      2,
                                      2,
                                      1,
                                      3,
                                      2,
                                      3,
                                      1,
                                      2,
                                      3,
                                      2,
                                      3,
                                      3),
                         family=c(1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  3,
                                  3,
                                  3,
                                  1,
                                  1,
                                  1,
                                  2,
                                  1,
                                  2,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  3,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  2,
                                  1,
                                  1,
                                  2,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1,
                                  1))




###---------------------CALLING SIGNIFICANT PARIWISE CONTRASTS---------------------#

######Smolts ENRC vs. WD
#control is treatment 1 WD (control)
smolts.ENRCvsWD.Tile1000.cov3<-calculateDiffMeth(Smolts.ENRCvsWD.Final.all.cov3.tile1000.filt.nor.unite,covariates=Smolts.ENRCvsWD.Covariate, overdispersion="MN", test="Chisq") 

#call significant regions
All.DMR.smolts.ENRCvsWD.Tile1000.cov3=getMethylDiff(smolts.ENRCvsWD.Tile1000.cov3, difference=15, qvalue=0.01, type="all")                                                              
head(All.DMR.smolts.ENRCvsWD.Tile1000.cov3)
nrow(All.DMR.smolts.ENRCvsWD.Tile1000.cov3) #15 

#hyper
Hyper.DMR.smolts.ENRCvsWD.Tile1000.cov3=getMethylDiff(smolts.ENRCvsWD.Tile1000.cov3, difference=15, qvalue=0.01, type="hyper")                                                              
head(Hyper.DMR.smolts.ENRCvsWD.Tile1000.cov3)
nrow(Hyper.DMR.smolts.ENRCvsWD.Tile1000.cov3) #2

#hypo
Hypo.DMR.smolts.ENRCvsWD.Tile1000.cov3=getMethylDiff(smolts.ENRCvsWD.Tile1000.cov3, difference=15, qvalue=0.01, type="hypo")                                                              
head(Hypo.DMR.smolts.ENRCvsWD.Tile1000.cov3)
nrow(Hypo.DMR.smolts.ENRCvsWD.Tile1000.cov3) #13



######Smolts ENRC vs. TD
smolts.ENRCvsTD.Tile1000.cov3<-calculateDiffMeth(Smolts.ENRCvsTD.Final.all.cov3.tile1000.filt.nor.unite, covariates=Smolts.ENRCvsTD.Covariate, overdispersion="MN", test="Chisq") 
nrow(smolts.ENRCvsTD.Tile1000.cov3) #262268

All.DMR.smolts.ENRCvsTD.Tile1000.cov3=getMethylDiff(smolts.ENRCvsTD.Tile1000.cov3, difference=15, qvalue=0.01, type="all")                                                              
head(All.DMR.smolts.ENRCvsTD.Tile1000.cov3)
nrow(All.DMR.smolts.ENRCvsTD.Tile1000.cov3) #8 
#hyper
Hyper.DMR.smolts.ENRCvsTD.Tile1000.cov3=getMethylDiff(smolts.ENRCvsTD.Tile1000.cov3, difference=15, qvalue=0.01, type="hyper")                                                              
head(Hyper.DMR.smolts.ENRCvsTD.Tile1000.cov3)
nrow(Hyper.DMR.smolts.ENRCvsTD.Tile1000.cov3) #5 

#hypo
Hypo.DMR.smolts.ENRCvsTD.Tile1000.cov3=getMethylDiff(smolts.ENRCvsTD.Tile1000.cov3, difference=15, qvalue=0.01, type="hypo")                                                              
head(Hypo.DMR.smolts.ENRCvsTD.Tile1000.cov3)
nrow(Hypo.DMR.smolts.ENRCvsTD.Tile1000.cov3) #3                                                              



#######Smolts WD vs. TD
#control is treatment 1 WD
smolts.WDvsTD.Tile1000.cov3<-calculateDiffMeth(Smolts.WDvsTD.Final.all.cov3.tile1000.filt.nor.unite, covariates=Smolts.WDvsTD.Covariate, overdispersion="MN", test="Chisq") 
head(smolts.WDvsTD.Tile1000.cov3) 
nrow(smolts.WDvsTD.Tile1000.cov3) #262268


All.DMR.smolts.WDvsTD.Tile1000.cov3=getMethylDiff(smolts.WDvsTD.Tile1000.cov3, difference=15, qvalue=0.01, type="all")                                                              
head(All.DMR.smolts.WDvsTD.Tile1000.cov3)
nrow(All.DMR.smolts.WDvsTD.Tile1000.cov3) #61

#hyper
Hyper.DMR.smolts.WDvsTD.Tile1000.cov3=getMethylDiff(smolts.WDvsTD.Tile1000.cov3, difference=15, qvalue=0.01, type="hyper")                                                              
head(Hyper.DMR.smolts.WDvsTD.Tile1000.cov3)
nrow(Hyper.DMR.smolts.WDvsTD.Tile1000.cov3) #24 

#hypo
Hypo.DMR.smolts.WDvsTD.Tile1000.cov3=getMethylDiff(smolts.WDvsTD.Tile1000.cov3, difference=15, qvalue=0.01, type="hypo")                                                              
head(Hypo.DMR.smolts.WDvsTD.Tile1000.cov3)
nrow(Hypo.DMR.smolts.WDvsTD.Tile1000.cov3) #38               



#########Adults ENRC vs. WD
#control is treatment 1 WD
adults.ENRCvsWD.Tile1000.cov3<-calculateDiffMeth(Adults.ENRCvsWD.Final.all.cov3.tile1000.filt.nor.unite, covariates=Adults.ENRCvsWD.Covariate, overdispersion="MN", test="Chisq") 
nrow(adults.ENRCvsWD.Tile1000.cov3) #262268

All.DMR.adults.ENRCvsWD.Tile1000.cov3=getMethylDiff(adults.ENRCvsWD.Tile1000.cov3, difference=15, qvalue=0.01, type="all")                                                              
head(All.DMR.adults.ENRCvsWD.Tile1000.cov3)
nrow(All.DMR.adults.ENRCvsWD.Tile1000.cov3) #14 

#hyper
Hyper.DMR.adults.ENRCvsWD.Tile1000.cov3=getMethylDiff(adults.ENRCvsWD.Tile1000.cov3, difference=15, qvalue=0.01, type="hyper")                                                              
head(Hyper.DMR.adults.ENRCvsWD.Tile1000.cov3)
nrow(Hyper.DMR.adults.ENRCvsWD.Tile1000.cov3) #8  

#hypo
Hypo.DMR.adults.ENRCvsWD.Tile1000.cov3=getMethylDiff(adults.ENRCvsWD.Tile1000.cov3, difference=15, qvalue=0.01, type="hypo")                                                              
head(Hypo.DMR.adults.ENRCvsWD.Tile1000.cov3)
nrow(Hypo.DMR.adults.ENRCvsWD.Tile1000.cov3) #6



##########Adults ENRC vs. TD3
adults.ENRCvsTD.Tile1000.cov3<-calculateDiffMeth(Adults.ENRCvsTD.Final.all.cov3.tile1000.filt.nor.unite, covariates=Adults.ENRCvsTD.Covariate, overdispersion="MN", test="Chisq") 
nrow(adults.ENRCvsTD.Tile1000.cov3) #262268


All.DMR.adults.ENRCvsTD.Tile1000.cov3=getMethylDiff(adults.ENRCvsTD.Tile1000.cov3, difference=15, qvalue=0.01, type="all")                                                              
head(All.DMR.adults.ENRCvsTD.Tile1000.cov3)
nrow(All.DMR.adults.ENRCvsTD.Tile1000.cov3) #8 


#hyper
Hyper.DMR.adults.ENRCvsTD.Tile1000.cov3=getMethylDiff(adults.ENRCvsTD.Tile1000.cov3, difference=15, qvalue=0.01, type="hyper")                                                              
head(Hyper.DMR.adults.ENRCvsTD.Tile1000.cov3)
nrow(Hyper.DMR.adults.ENRCvsTD.Tile1000.cov3) #3 

#hypo
Hypo.DMR.adults.ENRCvsTD.Tile1000.cov3=getMethylDiff(adults.ENRCvsTD.Tile1000.cov3, difference=15, qvalue=0.01, type="hypo")                                                              
head(Hypo.DMR.adults.ENRCvsTD.Tile1000.cov3)
nrow(Hypo.DMR.adults.ENRCvsTD.Tile1000.cov3)  #5 


#############Adults WD vs. TD
#control is treatment 1 WD
#this takes probably 10 mins
adults.WDvsTD.Tile1000.cov3<-calculateDiffMeth(Adults.WDvsTD.Final.all.cov3.tile1000.filt.nor.unite, covariates=Adults.WDvsTD.Covariate,overdispersion="MN", test="Chisq") 
head(adults.WDvsTD.Tile1000.cov3) 
nrow(adults.WDvsTD.Tile1000.cov3) # 262268

All.DMR.adults.WDvsTD.Tile1000.cov3=getMethylDiff(adults.WDvsTD.Tile1000.cov3, difference=15, qvalue=0.01, type="all")                                                              
head(All.DMR.adults.WDvsTD.Tile1000.cov3)
nrow(All.DMR.adults.WDvsTD.Tile1000.cov3) #470 


#hyper
Hyper.DMR.adults.WDvsTD.Tile1000.cov3=getMethylDiff(adults.WDvsTD.Tile1000.cov3, difference=15, qvalue=0.01, type="hyper")                                                              
head(Hyper.DMR.adults.WDvsTD.Tile1000.cov3)
nrow(Hyper.DMR.adults.WDvsTD.Tile1000.cov3)  #253

#hypo
Hypo.DMR.adults.WDvsTD.Tile1000.cov3=getMethylDiff(adults.WDvsTD.Tile1000.cov3, difference=15, qvalue=0.01, type="hypo")                                                              
head(Hypo.DMR.adults.WDvsTD.Tile1000.cov3)
nrow(Hypo.DMR.adults.WDvsTD.Tile1000.cov3) #217 



#################smolts to adults wild
#control is treatment 1 WD
#this takes probably 10 mins
DMR.WD.SMvsAD.Tile1000.cov3<-calculateDiffMeth(Wild.SMvsAD.Final.all.cov3.tile1000.filt.nor.unite, covariates=Wild.SMvsAD.Covariate, overdispersion="MN", test="Chisq") 
nrow(DMR.WD.SMvsAD.Tile1000.cov3) #] 262268

All.DMR.WD.SMvsAD.Tile1000.cov3=getMethylDiff(DMR.WD.SMvsAD.Tile1000.cov3, difference=15, qvalue=0.01, type="all")                                                              
nrow(All.DMR.WD.SMvsAD.Tile1000.cov3) #1552 


#hyper
Hyper.DMR.WD.SMvsAD.Tile1000.cov3=getMethylDiff(DMR.WD.SMvsAD.Tile1000.cov3, difference=15, qvalue=0.01, type="hyper")                                                              
nrow(Hyper.DMR.WD.SMvsAD.Tile1000.cov3) #618

#hypo
Hypo.DMR.WD.SMvsAD.Tile1000.cov3=getMethylDiff(DMR.WD.SMvsAD.Tile1000.cov3, difference=15, qvalue=0.01, type="hypo")                                                              
nrow(Hypo.DMR.WD.SMvsAD.Tile1000.cov3) #934




#####################smolts to adults Enriched
#this takes probably 10 mins
DMR.ENRC.SMvsAD.Tile1000.cov3<-calculateDiffMeth(Enriched.SMvsAD.Final.all.cov3.tile1000.filt.nor.unite, covariates=Enriched.SMvsAD.Covariate, overdispersion="MN", test="Chisq") 
nrow(DMR.ENRC.SMvsAD.Tile1000.cov3) #262268


All.DMR.ENRC.SMvsAD.Tile1000.cov3=getMethylDiff(DMR.ENRC.SMvsAD.Tile1000.cov3, difference=15, qvalue=0.01, type="all")                                                              
nrow(All.DMR.ENRC.SMvsAD.Tile1000.cov3) #2246

#hyper
Hyper.DMR.ENRC.SMvsAD.Tile1000.cov3=getMethylDiff(DMR.ENRC.SMvsAD.Tile1000.cov3, difference=15, qvalue=0.01, type="hyper")                                                              
nrow(Hyper.DMR.ENRC.SMvsAD.Tile1000.cov3) #1138

#hypo
Hypo.DMR.ENRC.SMvsAD.Tile1000.cov3=getMethylDiff(DMR.ENRC.SMvsAD.Tile1000.cov3, difference=15, qvalue=0.01, type="hypo")                                                              
nrow(Hypo.DMR.ENRC.SMvsAD.Tile1000.cov3) #1108



####################smolts to adults traditional
#this takes probably 10 mins
DMR.TD.SMvsAD.Tile1000.cov3<-calculateDiffMeth(Traditional.SMvsAD.Final.all.cov3.tile1000.filt.nor.unite, covariates=Traditional.SMvsAD.Covariate,overdispersion="MN", test="Chisq") 
nrow(DMR.TD.SMvsAD.Tile1000.cov3) #262268

All.DMR.TD.SMvsAD.Tile1000.cov3=getMethylDiff(DMR.TD.SMvsAD.Tile1000.cov3, difference=15, qvalue=0.01, type="all")                                                              
nrow(All.DMR.TD.SMvsAD.Tile1000.cov3) #2693

#hyper
Hyper.DMR.TD.SMvsAD.Tile1000.cov3=getMethylDiff(DMR.TD.SMvsAD.Tile1000.cov3, difference=15, qvalue=0.01, type="hyper")                                                              
nrow(Hyper.DMR.TD.SMvsAD.Tile1000.cov3) #1265

#hypo
Hypo.DMR.TD.SMvsAD.Tile1000.cov3=getMethylDiff(DMR.TD.SMvsAD.Tile1000.cov3, difference=15, qvalue=0.01, type="hypo")                                                              
nrow(Hypo.DMR.TD.SMvsAD.Tile1000.cov3) #1428 


##########Sex Effect
#model Final.all.sexeffect.cov3.tile1000.filt.nor.unite
Sex.Tile1000.cov3<-calculateDiffMeth(Final.all.sexeffect.cov3.tile1000.filt.nor.unite, covariates=Sex.Covariate, overdispersion="MN", test="Chisq") 
nrow(Sex.Tile1000.cov3) #342,938

All.DMR.Sex.Tile1000.cov3=getMethylDiff(Sex.Tile1000.cov3, difference=15, qvalue=0.01, type="all")                                                              
nrow(All.DMR.Sex.Tile1000.cov3) #17 

#hyper
Hyper.DMR.Sex.Tile1000.cov3=getMethylDiff(Sex.Tile1000.cov3, difference=15, qvalue=0.01, type="hyper")                                                              
nrow(Hyper.DMR.Sex.Tile1000.cov3) #11 

#hypo
Hypo.DMR.Sex.Tile1000.cov3=getMethylDiff(Sex.Tile1000.cov3, difference=15, qvalue=0.01, type="hypo")                                                              
nrow(Hypo.DMR.Sex.Tile1000.cov3) #6 


##############family effect
#model Final.all.3familyeffect.cov3.tile1000.filt.nor.unite
family.Tile1000.cov3<-calculateDiffMeth(Final.all.3familyeffect.cov3.tile1000.filt.nor.unite, covariates=Family.Covariate, overdispersion="MN", test="Chisq") 
nrow(family.Tile1000.cov3)#260251 tiles


All.DMR.family.Tile1000.cov3=getMethylDiff(family.Tile1000.cov3, difference=15, qvalue=0.01, type="all")                                                              
nrow(All.DMR.family.Tile1000.cov3) #94 

#hyper
Hyper.DMR.family.Tile1000.cov3=getMethylDiff(family.Tile1000.cov3, difference=15, qvalue=0.01, type="hyper")                                                              
nrow(Hyper.DMR.family.Tile1000.cov3) #94 



##----------------------------------Visual Plots of shared DMR's in pairwise comparisons-----------------------##

####UpSet plot to compare shared scaffolds between pairwise contrasts
All.DMR.smolts.ENRCvsWD.Tile1000.cov3.2 <-getData(All.DMR.smolts.ENRCvsWD.Tile1000.cov3)
All.DMR.smolts.ENRCvsTD.Tile1000.cov3.2<-getData(All.DMR.smolts.ENRCvsTD.Tile1000.cov3)
All.DMR.smolts.WDvsTD.Tile1000.cov3.2<-getData(All.DMR.smolts.WDvsTD.Tile1000.cov3)
All.DMR.adults.ENRCvsWD.Tile1000.cov3.2<-getData(All.DMR.adults.ENRCvsWD.Tile1000.cov3)
All.DMR.adults.ENRCvsTD.Tile1000.cov3.2<-getData(All.DMR.adults.ENRCvsTD.Tile1000.cov3)
All.DMR.adults.WDvsTD.Tile1000.cov3.2<-getData(All.DMR.adults.WDvsTD.Tile1000.cov3)
All.DMR.WD.SMvsAD.Tile1000.cov3.2<-getData(All.DMR.WD.SMvsAD.Tile1000.cov3)
All.DMR.ENRC.SMvsAD.Tile1000.cov3.2<-getData(All.DMR.ENRC.SMvsAD.Tile1000.cov3)
All.DMR.TD.SMvsAD.Tile1000.cov3.2<-getData(All.DMR.TD.SMvsAD.Tile1000.cov3)
All.DMR.Sex.Tile1000.cov3.2<-getData(All.DMR.Sex.Tile1000.cov3)
All.DMR.family.Tile1000.cov3.2<-getData(All.DMR.family.Tile1000.cov3)


#try with genomic regions
All.DMR.smolts.ENRCvsWD.Tile1000.cov3.genomic <- All.DMR.smolts.ENRCvsWD.Tile1000.cov3.2[ ,1]
All.DMR.smolts.ENRCvsTD.Tile1000.cov3.genomic <- All.DMR.smolts.ENRCvsTD.Tile1000.cov3.2[ ,1]
All.DMR.smolts.WDvsTD.Tile1000.cov3.genomic <- All.DMR.smolts.WDvsTD.Tile1000.cov3.2[ ,1]
All.DMR.adults.ENRCvsWD.Tile1000.cov3.genomic <- All.DMR.adults.ENRCvsWD.Tile1000.cov3.2[ ,1]
All.DMR.adults.ENRCvsTD.Tile1000.cov3.genomic <- All.DMR.adults.ENRCvsTD.Tile1000.cov3.2[ ,1]
All.DMR.adults.WDvsTD.Tile1000.cov3.genomic <- All.DMR.adults.WDvsTD.Tile1000.cov3.2[ ,1]
All.DMR.WD.SMvsAD.Tile1000.cov3.genomic <- All.DMR.WD.SMvsAD.Tile1000.cov3.2[ ,1]
All.DMR.ENRC.SMvsAD.Tile1000.cov3.genomic <- All.DMR.ENRC.SMvsAD.Tile1000.cov3.2[ ,1]
All.DMR.TD.SMvsAD.Tile1000.cov3.genomic <- All.DMR.TD.SMvsAD.Tile1000.cov3.2[ ,1]
All.DMR.Sex.Tile1000.cov3.genomic <- All.DMR.Sex.Tile1000.cov3.2[ ,1]
All.DMR.family.Tile1000.cov3.genomic <- All.DMR.family.Tile1000.cov3.2[ ,1]

list.genomic<-list(All.DMR.smolts.ENRCvsWD.Tile1000.cov3.genomic, All.DMR.smolts.ENRCvsTD.Tile1000.cov3.genomic, All.DMR.smolts.WDvsTD.Tile1000.cov3.genomic, All.DMR.adults.ENRCvsWD.Tile1000.cov3.genomic, 
                   All.DMR.adults.ENRCvsTD.Tile1000.cov3.genomic, All.DMR.adults.WDvsTD.Tile1000.cov3.genomic, All.DMR.WD.SMvsAD.Tile1000.cov3.genomic, All.DMR.ENRC.SMvsAD.Tile1000.cov3.genomic, All.DMR.TD.SMvsAD.Tile1000.cov3.genomic,All.DMR.Sex.Tile1000.cov3.genomic,All.DMR.family.Tile1000.cov3.genomic)
names(list.genomic) <-c("SM ERNCvsWD", "SM ENRCvsCV", "SM WDvsCV","AD ERNCvsWD", "AD ENRCvsCV", "AD WDvsCV", "WD SMvsAD","ENRC SMvsAD","TD SMvsAD", "SEX", "FAMILY")
matrix.genomic.chromosome = make_comb_mat(list.genomic)


ht = draw(UpSet(matrix.genomic.chromosome, right_annotation = upset_right_annotation(matrix.genomic.chromosome,
                                                                                     gp = gpar(fill = "#3182bd"))))
od = column_order(ht)
cs = comb_size(matrix.genomic.chromosome,)

decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(15, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})


####UpSet plot to compare shared DMRs between pairwise contrasts

All.DMR.smolts.ENRCvsWD.Tile1000.cov3.2 <-getData(All.DMR.smolts.ENRCvsWD.Tile1000.cov3)
All.DMR.smolts.ENRCvsWD.Tile1000.cov3.2 <- All.DMR.smolts.ENRCvsWD.Tile1000.cov3.2 %>%
  tidyr::unite('Merged', c("start","end"), sep="-", remove=FALSE)

All.DMR.smolts.ENRCvsTD.Tile1000.cov3.2<-getData(All.DMR.smolts.ENRCvsTD.Tile1000.cov3)
All.DMR.smolts.ENRCvsTD.Tile1000.cov3.2 <- All.DMR.smolts.ENRCvsTD.Tile1000.cov3.2 %>%
  tidyr::unite('Merged', c("start","end"), sep="-", remove=FALSE)

All.DMR.smolts.WDvsTD.Tile1000.cov3.2<-getData(All.DMR.smolts.WDvsTD.Tile1000.cov3)
All.DMR.smolts.WDvsTD.Tile1000.cov3.2 <- All.DMR.smolts.WDvsTD.Tile1000.cov3.2 %>%
  tidyr::unite('Merged', c("start","end"), sep="-", remove=FALSE)

All.DMR.adults.ENRCvsWD.Tile1000.cov3.2<-getData(All.DMR.adults.ENRCvsWD.Tile1000.cov3)
All.DMR.adults.ENRCvsWD.Tile1000.cov3.2 <- All.DMR.adults.ENRCvsWD.Tile1000.cov3.2 %>%
  tidyr::unite('Merged', c("start","end"), sep="-", remove=FALSE)

All.DMR.adults.ENRCvsTD.Tile1000.cov3.2<-getData(All.DMR.adults.ENRCvsTD.Tile1000.cov3)
All.DMR.adults.ENRCvsTD.Tile1000.cov3.2 <- All.DMR.adults.ENRCvsTD.Tile1000.cov3.2 %>%
  tidyr::unite('Merged', c("start","end"), sep="-", remove=FALSE)

All.DMR.adults.WDvsTD.Tile1000.cov3.2<-getData(All.DMR.adults.WDvsTD.Tile1000.cov3)
All.DMR.adults.WDvsTD.Tile1000.cov3.2 <- All.DMR.adults.WDvsTD.Tile1000.cov3.2 %>%
  tidyr::unite('Merged', c("start","end"), sep="-", remove=FALSE)

All.DMR.WD.SMvsAD.Tile1000.cov3.2<-getData(All.DMR.WD.SMvsAD.Tile1000.cov3)
All.DMR.WD.SMvsAD.Tile1000.cov3.2 <- All.DMR.WD.SMvsAD.Tile1000.cov3.2 %>%
  tidyr::unite('Merged', c("start","end"), sep="-", remove=FALSE)

All.DMR.ENRC.SMvsAD.Tile1000.cov3.2<-getData(All.DMR.ENRC.SMvsAD.Tile1000.cov3)
All.DMR.ENRC.SMvsAD.Tile1000.cov3.2 <- All.DMR.ENRC.SMvsAD.Tile1000.cov3.2 %>%
  tidyr::unite('Merged', c("start","end"), sep="-", remove=FALSE)

All.DMR.TD.SMvsAD.Tile1000.cov3.2<-getData(All.DMR.TD.SMvsAD.Tile1000.cov3)
All.DMR.TD.SMvsAD.Tile1000.cov3.2 <- All.DMR.TD.SMvsAD.Tile1000.cov3.2 %>%
  tidyr::unite('Merged', c("start","end"), sep="-", remove=FALSE)
nrow(All.DMR.TD.SMvsAD.Tile1000.cov3.2)

All.DMR.Sex.Tile1000.cov3.2<-getData(All.DMR.Sex.Tile1000.cov3)
All.DMR.Sex.Tile1000.cov3.2 <- All.DMR.Sex.Tile1000.cov3.2 %>%
  tidyr::unite('Merged', c("start","end"), sep="-", remove=FALSE)

All.DMR.family.Tile1000.cov3.2<-getData(All.DMR.family.Tile1000.cov3)
All.DMR.family.Tile1000.cov3.2 <- All.DMR.family.Tile1000.cov3.2 %>%
  tidyr::unite('Merged', c("start","end"), sep="-", remove=FALSE)


test_gr.1=All.DMR.smolts.ENRCvsWD.Tile1000.cov3.2$Merged
test_gr.2=All.DMR.smolts.ENRCvsTD.Tile1000.cov3.2$Merged
test_gr.3=All.DMR.smolts.WDvsTD.Tile1000.cov3.2$Merged
test_gr.4=All.DMR.adults.ENRCvsWD.Tile1000.cov3.2$Merged
test_gr.5=All.DMR.adults.ENRCvsTD.Tile1000.cov3.2$Merged
test_gr.6=All.DMR.adults.WDvsTD.Tile1000.cov3.2$Merged
test_gr.7=All.DMR.WD.SMvsAD.Tile1000.cov3.2$Merged
test_gr.8=All.DMR.ENRC.SMvsAD.Tile1000.cov3.2$Merged
test_gr.9=All.DMR.TD.SMvsAD.Tile1000.cov3.2$Merged
test_gr.10=All.DMR.Sex.Tile1000.cov3.2$Merged
test_gr.11=All.DMR.family.Tile1000.cov3.2$Merged
nrow(All.DMR.TD.SMvsAD.Tile1000.cov3.2$Merged)


list.genomic<-list(test_gr.1,test_gr.2,test_gr.3,test_gr.4,test_gr.5,test_gr.6,test_gr.7,test_gr.8,test_gr.9,test_gr.10,test_gr.11)
names(list.genomic) <-c("S SNvsWD", "S SNvsCV", "S WDvsCV","A SNvsWD", "A SNvsCV", "A WDvsCV", "WD SvsA","SN SvsA","CV SvsA","SEX","FAMILY")
matrix.genomic.DMR = make_comb_mat(list.genomic)


ht = draw(UpSet(matrix.genomic.DMR, right_annotation = upset_right_annotation(matrix.genomic.DMR,
                                                                              gp = gpar(fill = "#3182bd"))))          
od = column_order(ht)
cs = comb_size(matrix.genomic.DMR,)
#cs = formatC(cs, format = "e", digits = 0)

# rod = row_order(ht)
# ss = set_size(matrix.genomic.DMR)
# ss = formatC(ss, format = "e", digits = 0)

decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(4, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})

plot_grid(PCA7, PCA8, PCA9, align ="V", labels = "AUTO", ncol=1)



##----------------------------------Visual Plots of sequencing depth/coverage-----------------------##

#### Plot the raw reads comparing the nubmer of reads from HiSeq to Novaseq

Reads <- read.csv("C:/Users/BOKVISTJ/Documents/UofC MSc/Sequencing Read Stats Summarized.csv")
names(Reads)[names(Reads) == "NovaSeq.6000...2021"] <- "Novaseq 6000 - 2021"
names(Reads)[names(Reads) == "NovaSeq.6000...2020"] <- "Novaseq 6000 - 2020"
names(Reads)[names(Reads) == "HiSeq.4000..2019"] <- "HiSeq 4000 - 2019"

#gather data for subsequent processing
 Plot.data <- Reads %>%
   dplyr::select(Name,"Novaseq 6000 - 2021","Novaseq 6000 - 2020","HiSeq 4000 - 2019") %>%
   filter(Name!="M167") %>%
   filter(Name!="G3848") %>%
   mutate(Individual=rep(1:92))%>%
   dplyr::select(Individual, "Novaseq 6000 - 2021","Novaseq 6000 - 2020","HiSeq 4000 - 2019")
   
Plot.data$Individual <-as.factor(Plot.data$Individual)
Plot.data.melt<-melt.data.frame(Plot.data, id.vars="Individual") 

Plot.data.2 <- Reads %>%
  dplyr::select(Name,"Novaseq 6000 - 2021","Novaseq 6000 - 2020","HiSeq 4000 - 2019") %>%
 filter(Name=="M167" | Name=="G3848") %>%
  mutate(Individual=rep(93:94))%>%
  dplyr::select(Individual,"Novaseq 6000 - 2021","Novaseq 6000 - 2020","HiSeq 4000 - 2019")

Plot.data.2$Individual <-as.factor(Plot.data.2$Individual)
Plot.data.2.melt<-melt.data.frame(Plot.data.2, id.vars="Individual") 

#plot the samples
a<-ggplot(Plot.data.melt, aes(fill = variable,
                      y = value, x = Individual))+
  geom_bar(position = "stack", stat = "identity",colour="black")+
  scale_y_continuous(breaks = seq(0, 100000000, len = 8),labels = label_number(suffix = " M", scale = 1e-6))+
  labs(y="Number of reads (millions)", x="Individual")+
  scale_x_discrete(breaks = seq(1, 92, by = 2))+
 scale_fill_manual(values = c("#3182bd", "#9ecae1", "#deebf7"), name="Sequencing Run")+theme(
   legend.position = c(.95, .95),
   legend.justification = c("right", "top"),
   legend.box.just = "right",
   legend.margin = margin(6, 6, 6, 6),
   legend.key.size = unit(1, 'cm'))


#two individuals with over 100 million reads
b<-ggplot(Plot.data.2.melt, aes(fill = variable,
                           y = value, x = Individual))+
  geom_bar(position = "stack", stat = "identity",colour="black")+
  scale_y_continuous(breaks = seq(0, 325000000, len = 8),labels = label_number(suffix = " M", scale = 1e-6))+
  labs(y="Number of Rrads (millions)", x="Individual")+
  scale_fill_manual(values = c("#3182bd", "#9ecae1", "#deebf7"), name="Sequencing Run")+
  theme(legend.position="none")

plot_grid(a, b, labels = "AUTO", rel_widths = c(5, 1))


####coverage plot
data_all <- list.files(path = "C:/Users/BOKVISTJ/Documents/UofC MSc/All Analysis December 2021/R Analysis All Sequencing Data/Sequencing Depths",     
                       pattern = "*depth_genome", full.names = TRUE) 

myread<- function(fileName){
  data.frame(
    read.delim(fileName, header=F)
    , "Sample" = fileName
  )
}


data_all_files <- lapply(data_all, myread) %>%
  bind_rows %>%
  mutate(ID= case_when(
    grepl("F159",Sample)~"1",
    grepl("F161",Sample)~"2",
    grepl("F163",Sample)~"3",
    grepl("F165",Sample)~"4",
    grepl("F166",Sample)~"5",
    grepl("F167",Sample)~"6",
    grepl("F33",Sample)~"7",
    grepl("F35",Sample)~"8",
    grepl("F39",Sample)~"9",
    grepl("F40",Sample)~"10",
    grepl("F41",Sample)~"11",
    grepl("F43",Sample)~"12",
    grepl("F45",Sample)~"13",
    grepl("F48",Sample)~"14",
    grepl("F49",Sample)~"15",
    grepl("F50",Sample)~"16",
    grepl("F73",Sample)~"17",
    grepl("G3846",Sample)~"18",
    grepl("G3847",Sample)~"19",
    grepl("G3848",Sample)~"20",
    grepl("G3849",Sample)~"21",
    grepl("G3850",Sample)~"22",
    grepl("G3851",Sample)~"23",
    grepl("G3852",Sample)~"24",
    grepl("G3853",Sample)~"25",
    grepl("G3854",Sample)~"26",
    grepl("G3855",Sample)~"27",
    grepl("G3856",Sample)~"28",
    grepl("G3858",Sample)~"29",
    grepl("G3859",Sample)~"30",
    grepl("G3860",Sample)~"31",
    grepl("G3861",Sample)~"32",
    grepl("G3862",Sample)~"33",
    grepl("G3863",Sample)~"34",
    grepl("G3864",Sample)~"35",
    grepl("G3865",Sample)~"36",
    grepl("G3866",Sample)~"37",
    grepl("G3867",Sample)~"38",
    grepl("G3868",Sample)~"39",
    grepl("G3869",Sample)~"40",
    grepl("G3870",Sample)~"41",
    grepl("G3871",Sample)~"42",
    grepl("G3872",Sample)~"43",
    grepl("G3873",Sample)~"44",
    grepl("G3874",Sample)~"45",
    grepl("G3875",Sample)~"46",
    grepl("H656",Sample)~"47",
    grepl("H657",Sample)~"48",
    grepl("H658",Sample)~"49",
    grepl("H659",Sample)~"50",
    grepl("H660",Sample)~"51",
    grepl("H661",Sample)~"52",
    grepl("H662",Sample)~"53",
    grepl("H663",Sample)~"54",
    grepl("H664",Sample)~"55",
    grepl("H665",Sample)~"56",
    grepl("H666",Sample)~"57",
    grepl("H667",Sample)~"58",
    grepl("H668",Sample)~"59",
    grepl("H669",Sample)~"60",
    grepl("H670",Sample)~"61",
    grepl("H671",Sample)~"62",
    grepl("H672",Sample)~"63",
    grepl("H673",Sample)~"64",
    grepl("H674",Sample)~"65",
    grepl("H675",Sample)~"66",
    grepl("H676",Sample)~"67",
    grepl("H677",Sample)~"68",
    grepl("H678",Sample)~"69",
    grepl("H679",Sample)~"70",
    grepl("H680",Sample)~"71",
    grepl("H681",Sample)~"72",
    grepl("H682",Sample)~"73",
    grepl("H683",Sample)~"74",
    grepl("H684",Sample)~"75",
    grepl("H685",Sample)~"76",
    grepl("M159",Sample)~"77",
    grepl("M161",Sample)~"78",
    grepl("M162",Sample)~"79",
    grepl("M164",Sample)~"80",
    grepl("M166",Sample)~"81",
    grepl("M167",Sample)~"82",
    grepl("M169",Sample)~"83",
    grepl("M34",Sample)~"84",
    grepl("M36",Sample)~"85",
    grepl("M37",Sample)~"86",
    grepl("M38",Sample)~"87",
    grepl("M41",Sample)~"88",
    grepl("M63",Sample)~"89",
    grepl("M71",Sample)~"90",
    grepl("M72",Sample)~"91",
    grepl("M73",Sample)~"92",
    grepl("M74",Sample)~"93",
    grepl("M77",Sample)~"94"
  )) %>%
  filter(V2!=0)

 
#plot the coverage
line_plot <- data_all_files %>%
  filter(V2<101)

   ggplot(line_plot, aes(y =V3 , x = V2, fill="ID"))+
     geom_point()+
    scale_y_continuous(breaks = seq(0, 20000000, len = 8),labels = label_number(suffix = " M", scale = 1e-6))+
    labs(y="Number of bases (millions)", x="Coverage")+
    theme(legend.position="none")



##----------------------------------Making PCA's from MethylKit Data--------------------------------##

#create PCAs based on united file - Final.all.cov3.tile1000.filt.nor.unite

####PCA all samples colored by treatment
perc.mat<-percMethylation(Final.all.cov3.tile1000.filt.nor.unite)
perc.mat<-as.data.frame(perc.mat)
dim(perc.mat) #262268
perc.mat.2=perc.mat[ rowSums(is.na(perc.mat))==0, ]
dim(perc.mat.2) # 117185

sd.threshold <- 0.5 ## this takes only the position
sds=rowSds(as.matrix(perc.mat.2))
head(sds)
cutoff=quantile(sds,sd.threshold)
perc.mat.2=perc.mat.2[sds>cutoff,]
head(perc.mat.2)
dim(perc.mat.2) #58592    94

prcomp.perc.mat <-prcomp(t(perc.mat.2), scale. = TRUE, center = TRUE)
summary(prcomp.perc.mat)

myColors <- c("#fee0b6", "#f1a340", "#b35806", "#d8daeb","#998ec3","#542788")
names(myColors) <- c(1,2,3,4,5,6)
custom_colors <- scale_fill_manual(name = "Treatment", values = myColors, labels=c("WD-S","SN-S","CV-S","WD-A","SN-A","CV-A"), )

treatment=Final.all.cov3.tile1000.filt.nor.unite@treatment
treatment<-as.factor(treatment)
sample.ids=Final.all.cov3.tile1000.filt.nor.unite@sample.ids


PCA1<-ggplot(prcomp.perc.mat, aes(x = PC1, y =PC2, fill=treatment)) +
  geom_point(size=2, shape=21)+
  labs(x="PC1 10.04%", y="PC2 3.35%")+
  custom_colors


#what happens if I remove some PC's? Looking at the screen plot there to decide:

#looking at how many princirple compoenents to retain 
library(factoextra)

#nice screeplot 
#from the united object specifying treatment
scree1<-fviz_eig(prcomp.perc.mat, geom="bar", main = "" ,barfill = "#3182bd",
                 barcolor = "black", ncp=30, addlabels=T)

#note that screeplots for the united objects that specify treatement, sex or family are all the

#look at different combos of the PC's - first 3 are only important really

ggplot(prcomp.perc.mat, aes(x = PC2, y =PC3, fill=treatment)) +
  geom_point(size=2, shape=21)+
  labs(x="PC1 10.04%", y="PC2 3.35%")+
  custom_colors+
  geom_text(aes(label=sample.ids))


##### PCA all individuals colored by family
Final.all.3familyeffect.cov3.tile1000.filt.nor.unite <- readRDS("C:/Users/BOKVISTJ/Documents/UofC MSc/All Analysis December 2021/R Analysis All Sequencing Data/Final.all.3familyeffect.cov3.tile1000.filt.nor.unite")
nrow(Final.all.3familyeffect.cov3.tile1000.filt.nor.unite) #260251 tiles

All.PCA.family.perc.mat<-percMethylation(Final.all.3familyeffect.cov3.tile1000.filt.nor.unite)
All.PCA.family.perc.mat<-as.data.frame(All.PCA.family.perc.mat)
dim(All.PCA.family.perc.mat) #260251  
All.PCA.family.perc.mat.2=All.PCA.family.perc.mat[ rowSums(is.na(All.PCA.family.perc.mat))==0, ]
dim(All.PCA.family.perc.mat.2) # 117185     94

sd.threshold <- 0.5 ## this takes only the position
sds=rowSds(as.matrix(All.PCA.family.perc.mat.2))
head(sds)
cutoff=quantile(sds,sd.threshold)
All.PCA.family.perc.mat.2=All.PCA.family.perc.mat.2[sds>cutoff,]
head(All.PCA.family.perc.mat.2)
dim(All.PCA.family.perc.mat.2) # 58592    94

prcomp.All.PCA.family.perc.mat <-prcomp(t(All.PCA.family.perc.mat.2), scale. = TRUE, center = TRUE)
summary(prcomp.All.PCA.family.perc.mat)

myColors.2=c("black","yellow1","green2")  
names(myColors.2) <- c(1,2,3)
custom_colors.2 <- scale_fill_manual(name = "Family", values = myColors.2, labels=c("Unrelated","Family 1","Family 2"))
treatment.2=Final.all.3familyeffect.cov3.tile1000.filt.nor.unite@treatment
treatment.2<-as.factor(treatment.2)

PCA2<-ggplot(prcomp.All.PCA.family.perc.mat, aes(x = PC1, y =PC2, fill=treatment.2)) +
  geom_point(size=2, shape=21)+
  labs(x="PC1 10.04%", y="PC2 3.35%")+
  custom_colors.2


####PCA all individuals colored  by sex
Final.all.sexeffect.cov3.tile1000.filt.nor.unite <- readRDS("C:/Users/BOKVISTJ/Documents/UofC MSc/All Analysis December 2021/R Analysis All Sequencing Data/Final.all.sexeffect.cov3.tile1000.filt.nor.unite")
nrow(Final.all.sexeffect.cov3.tile1000.filt.nor.unite) # 342,938 larger number of regions than the original model had 

All.PCA.Sex.perc.mat<-percMethylation(Final.all.sexeffect.cov3.tile1000.filt.nor.unite)
All.PCA.Sex.perc.mat<-as.data.frame(All.PCA.Sex.perc.mat)
dim(All.PCA.Sex.perc.mat) #342938 
All.PCA.Sex.perc.mat.2=All.PCA.Sex.perc.mat[ rowSums(is.na(All.PCA.Sex.perc.mat))==0, ]
dim(All.PCA.Sex.perc.mat.2) # 117185     94

sd.threshold <- 0.5 ## this takes only the position
sds=rowSds(as.matrix(All.PCA.Sex.perc.mat.2))
head(sds)
cutoff=quantile(sds,sd.threshold)
All.PCA.Sex.perc.mat.2=All.PCA.Sex.perc.mat.2[sds>cutoff,]
head(All.PCA.Sex.perc.mat.2)
dim(All.PCA.Sex.perc.mat.2) # 58592    94
prcomp.All.PCA.Sex.perc.mat <-prcomp(t(All.PCA.Sex.perc.mat.2), scale. = TRUE, center = TRUE)
summary(prcomp.All.PCA.Sex.perc.mat)

myColors.3=c("white","black") 
names(myColors.3) <- c(1,2)
custom_colors.3 <- scale_fill_manual(name = "Sex", values = myColors.3, labels=c("Female","Male"))
treatment.3=Final.all.sexeffect.cov3.tile1000.filt.nor.unite@treatment
treatment.3<-as.factor(treatment.3)
sample.ids=Final.all.sexeffect.cov3.tile1000.filt.nor.unite@sample.ids


PCA3<-ggplot(prcomp.All.PCA.Sex.perc.mat, aes(x = PC1, y =PC2, fill=treatment.3)) +
  geom_point(size=2, shape=21)+
  labs(x="PC1 10.04%", y="PC2 3.35%")+
  custom_colors.3

plot_grid(PCA1, PCA2, PCA3, align ="V", labels = "AUTO", ncol=1)



####PCA of smolts colored by family

Final.all.3familyeffect.cov3.tile1000.filt.nor.unite <- readRDS("C:/Users/BOKVISTJ/Documents/UofC MSc/All Analysis December 2021/R Analysis All Sequencing Data/Final.all.3familyeffect.cov3.tile1000.filt.nor.unite")

Smolts.PCA.family <- reorganize(Final.all.3familyeffect.cov3.tile1000.filt.nor.unite, sample.ids=c("G3846",
                                                                              "G3847",
                                                                              "G3848",
                                                                              "G3849",
                                                                              "G3850",
                                                                              "G3851",
                                                                              "G3852",
                                                                              "G3853",
                                                                              "G3854",
                                                                              "G3855",
                                                                              "G3856",
                                                                              "G3858",
                                                                              "G3859",
                                                                              "G3860",
                                                                              "G3861",
                                                                              "G3862",
                                                                              "G3863",
                                                                              "G3864",
                                                                              "G3865",
                                                                              "G3866",
                                                                              "G3867",
                                                                              "G3868",
                                                                              "G3869",
                                                                              "G3870",
                                                                              "G3871",
                                                                              "G3872",
                                                                              "G3873",
                                                                              "G3874",
                                                                              "G3875",
                                                                              "H656",
                                                                              "H657",
                                                                              "H658",
                                                                              "H659",
                                                                              "H660",
                                                                              "H661",
                                                                              "H662",
                                                                              "H663",
                                                                              "H664",
                                                                              "H665",
                                                                              "H666",
                                                                              "H667",
                                                                              "H668",
                                                                              "H669",
                                                                              "H670",
                                                                              "H671",
                                                                              "H672",
                                                                              "H673",
                                                                              "H674",
                                                                              "H675",
                                                                              "H676",
                                                                              "H677",
                                                                              "H678",
                                                                              "H679",
                                                                              "H680",
                                                                              "H681",
                                                                              "H682",
                                                                              "H683",
                                                                              "H684",
                                                                              "H685")

                                                                     , treatment=c(1,
                                                                                   1,
                                                                                   2,
                                                                                   1,
                                                                                   2,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   3,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   2,
                                                                                   1,
                                                                                   1,
                                                                                   2,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1,
                                                                                   1))
Smolts.PCA.family.perc.mat<-percMethylation(Smolts.PCA.family)
Smolts.PCA.family.perc.mat<-as.data.frame(Smolts.PCA.family.perc.mat)
dim(Smolts.PCA.family.perc.mat) #260251
Smolts.PCA.family.perc.mat.2=Smolts.PCA.family.perc.mat[ rowSums(is.na(Smolts.PCA.family.perc.mat))==0, ]
dim(Smolts.PCA.family.perc.mat.2) # 120071      59

sd.threshold <- 0.5 ## this takes only the position
sds=rowSds(as.matrix(Smolts.PCA.family.perc.mat.2))
head(sds)
cutoff=quantile(sds,sd.threshold)
Smolts.PCA.family.perc.mat.2=Smolts.PCA.family.perc.mat.2[sds>cutoff,]
head(Smolts.PCA.family.perc.mat.2)
dim(Smolts.PCA.family.perc.mat.2) #60035    59

prcomp.Smolts.PCA.family.perc.mat <-prcomp(t(Smolts.PCA.family.perc.mat.2), scale. = TRUE, center = TRUE)
summary(prcomp.Smolts.PCA.family.perc.mat)

myColors.5=c("black","yellow1","green2")  
names(myColors.5) <- c(1,2,3)
custom_colors.5 <- scale_fill_manual(name = "Family", values = myColors.5, labels=c("Unrelated", "Family 1", "Family 2"))

screeplot(prcomp.Smolts.PCA.family.perc.mat, type="barplot", main=paste(context="CpG methylation PCA Screeplot"))
treatment.5=Smolts.PCA.family@treatment
treatment.5<-as.factor(treatment.5)
sample.ids=Smolts.PCA.family@sample.ids

PCA5<-ggplot(prcomp.Smolts.PCA.family.perc.mat, aes(x = PC1, y =PC2, fill=treatment.5)) +
  geom_point(size=2, shape=21)+
  labs(x="PC1 3.57%", y="PC2 3.35%")+
  custom_colors.5

#### PCA smolts colored by treatment
#have to re-organize the file and assign a new treatment coding
Smolts.PCA.treatment<- reorganize(Final.all.cov3.tile1000.filt.nor.unite, sample.ids=c("G3846",
                                                                                     "G3847",
                                                                                     "G3848",
                                                                                     "G3849",
                                                                                     "G3850",
                                                                                     "G3851",
                                                                                     "G3852",
                                                                                     "G3853",
                                                                                     "G3854",
                                                                                     "G3855",
                                                                                     "G3856",
                                                                                     "G3858",
                                                                                     "G3859",
                                                                                     "G3860",
                                                                                     "G3861",
                                                                                     "G3862",
                                                                                     "G3863",
                                                                                     "G3864",
                                                                                     "G3865",
                                                                                     "G3866",
                                                                                     "G3867",
                                                                                     "G3868",
                                                                                     "G3869",
                                                                                     "G3870",
                                                                                     "G3871",
                                                                                     "G3872",
                                                                                     "G3873",
                                                                                     "G3874",
                                                                                     "G3875",
                                                                                     "H656",
                                                                                     "H657",
                                                                                     "H658",
                                                                                     "H659",
                                                                                     "H660",
                                                                                     "H661",
                                                                                     "H662",
                                                                                     "H663",
                                                                                     "H664",
                                                                                     "H665",
                                                                                     "H666",
                                                                                     "H667",
                                                                                     "H668",
                                                                                     "H669",
                                                                                     "H670",
                                                                                     "H671",
                                                                                     "H672",
                                                                                     "H673",
                                                                                     "H674",
                                                                                     "H675",
                                                                                     "H676",
                                                                                     "H677",
                                                                                     "H678",
                                                                                     "H679",
                                                                                     "H680",
                                                                                     "H681",
                                                                                     "H682",
                                                                                     "H683",
                                                                                     "H684",
                                                                                     "H685")
                                
                                , treatment=c(2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              2,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              3,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1))

Smolts.PCA.treatment.perc.mat<-percMethylation(Smolts.PCA.treatment)
Smolts.PCA.treatment.perc.mat<-as.data.frame(Smolts.PCA.treatment.perc.mat)
dim(Smolts.PCA.treatment.perc.mat) #262268
Smolts.PCA.treatment.perc.mat.2=Smolts.PCA.treatment.perc.mat[ rowSums(is.na(Smolts.PCA.treatment.perc.mat))==0, ]
dim(Smolts.PCA.treatment.perc.mat.2) #  120071     59
head(Smolts.PCA.treatment.perc.mat.2)

sd.threshold <- 0.5 ## this takes only the position
sds=rowSds(as.matrix(Smolts.PCA.treatment.perc.mat.2))
cutoff=quantile(sds,sd.threshold)
Smolts.PCA.treatment.perc.mat.2=Smolts.PCA.treatment.perc.mat.2[sds>cutoff,]
head(Smolts.PCA.treatment.perc.mat.2)
dim(Smolts.PCA.treatment.perc.mat.2) # 60035   59

prcomp.Smolts.PCA.treatment.perc.mat <-prcomp(t(Smolts.PCA.treatment.perc.mat.2), scale. = TRUE, center = TRUE)
summary(prcomp.Smolts.PCA.treatment.perc.mat)

myColors.4=c("#fee0b6", "#f1a340", "#b35806") #orange shades
names(myColors.4) <- c(1,2,3)
custom_colors.4 <- scale_fill_manual(name = "Treatment", values = myColors.4, labels=c("WD-S","SN-S","CV-S") )

treatment.4=Smolts.PCA.treatment@treatment
treatment.4<-as.factor(treatment.4)
sample.ids=Smolts.PCA.treatment@sample.ids

PCA4<-ggplot(prcomp.Smolts.PCA.treatment.perc.mat, aes(x = PC1, y =PC2, fill=treatment.4)) +
  geom_point(size=2, shape=21)+
  labs(x="PC1 3.57%", y="PC2 3.35%")+
  custom_colors.4


#### PCA smolts colored by sex
#have to re-organize the file and assign a new treatment coding
Final.all.sexeffect.cov3.tile1000.filt.nor.unite <- readRDS("C:/Users/BOKVISTJ/Documents/UofC MSc/All Analysis December 2021/R Analysis All Sequencing Data/Final.all.sexeffect.cov3.tile1000.filt.nor.unite")
nrow(Final.all.sexeffect.cov3.tile1000.filt.nor.unite) # 342,938 larger number of regions than the original model had 

Smolts.PCA.sex<- reorganize(Final.all.sexeffect.cov3.tile1000.filt.nor.unite, sample.ids=c("G3846",
                                                                                       "G3847",
                                                                                       "G3848",
                                                                                       "G3849",
                                                                                       "G3850",
                                                                                       "G3851",
                                                                                       "G3852",
                                                                                       "G3853",
                                                                                       "G3854",
                                                                                       "G3855",
                                                                                       "G3856",
                                                                                       "G3858",
                                                                                       "G3859",
                                                                                       "G3860",
                                                                                       "G3861",
                                                                                       "G3862",
                                                                                       "G3863",
                                                                                       "G3864",
                                                                                       "G3865",
                                                                                       "G3866",
                                                                                       "G3867",
                                                                                       "G3868",
                                                                                       "G3869",
                                                                                       "G3870",
                                                                                       "G3871",
                                                                                       "G3872",
                                                                                       "G3873",
                                                                                       "G3874",
                                                                                       "G3875",
                                                                                       "H656",
                                                                                       "H657",
                                                                                       "H658",
                                                                                       "H659",
                                                                                       "H660",
                                                                                       "H661",
                                                                                       "H662",
                                                                                       "H663",
                                                                                       "H664",
                                                                                       "H665",
                                                                                       "H666",
                                                                                       "H667",
                                                                                       "H668",
                                                                                       "H669",
                                                                                       "H670",
                                                                                       "H671",
                                                                                       "H672",
                                                                                       "H673",
                                                                                       "H674",
                                                                                       "H675",
                                                                                       "H676",
                                                                                       "H677",
                                                                                       "H678",
                                                                                       "H679",
                                                                                       "H680",
                                                                                       "H681",
                                                                                       "H682",
                                                                                       "H683",
                                                                                       "H684",
                                                                                       "H685")
                                  
                                  , treatment=c(1,
                                                2,
                                                1,
                                                2,
                                                2,
                                                1,
                                                1,
                                                2,
                                                2,
                                                2,
                                                1,
                                                2,
                                                2,
                                                1,
                                                1,
                                                2,
                                                1,
                                                1,
                                                2,
                                                2,
                                                2,
                                                2,
                                                1,
                                                1,
                                                2,
                                                2,
                                                2,
                                                2,
                                                1,
                                                2,
                                                1,
                                                2,
                                                2,
                                                2,
                                                2,
                                                1,
                                                1,
                                                1,
                                                2,
                                                2,
                                                2,
                                                1,
                                                2,
                                                2,
                                                2,
                                                1,
                                                1,
                                                1,
                                                2,
                                                1,
                                                2,
                                                1,
                                                1,
                                                1,
                                                1,
                                                2,
                                                2,
                                                1,
                                                2))
Smolts.PCA.sex.perc.mat<-percMethylation(Smolts.PCA.sex)
Smolts.PCA.sex.perc.mat<-as.data.frame(Smolts.PCA.sex.perc.mat)
dim(Smolts.PCA.sex.perc.mat) #342938 
Smolts.PCA.sex.perc.mat.2=Smolts.PCA.sex.perc.mat[ rowSums(is.na(Smolts.PCA.sex.perc.mat))==0, ]
dim(Smolts.PCA.sex.perc.mat.2) #  120071     59

sd.threshold <- 0.5 ## this takes only the position
sds=rowSds(as.matrix(Smolts.PCA.sex.perc.mat.2))
head(sds)
cutoff=quantile(sds,sd.threshold)
Smolts.PCA.sex.perc.mat.2=Smolts.PCA.sex.perc.mat.2[sds>cutoff,]
head(Smolts.PCA.sex.perc.mat.2)
dim(Smolts.PCA.sex.perc.mat.2) # 60035    59

prcomp.Smolts.PCA.sex.perc.mat <-prcomp(t(Smolts.PCA.sex.perc.mat.2), scale. = TRUE, center = TRUE)
summary(prcomp.Smolts.PCA.sex.perc.mat)

myColors.6=c("white","black")
names(myColors.6) <- c(1,2)
custom_colors.6 <- scale_fill_manual(name = "Sex", values = myColors.6, labels=c("Female","Male") )

treatment.6=Smolts.PCA.sex@treatment
treatment.6<-as.factor(treatment.6)
sample.ids=Smolts.PCA.sex@sample.ids

PCA6<-ggplot(prcomp.Smolts.PCA.sex.perc.mat, aes(x = PC1, y =PC2, fill=treatment.6)) +
  geom_point(size=2, shape=21)+
  labs(x="PC1 3.57%", y="PC2 3.35%")+
  custom_colors.6

plot_grid(PCA4, PCA5, PCA6, align ="V", labels = "AUTO", ncol=1)


#### PCA adults colored by family
#have to re-organize the file and assign a new treatment coding
Final.all.3familyeffect.cov3.tile1000.filt.nor.unite <- readRDS("C:/Users/BOKVISTJ/Documents/UofC MSc/All Analysis December 2021/R Analysis All Sequencing Data/Final.all.3familyeffect.cov3.tile1000.filt.nor.unite")

Adults.PCA.family <- reorganize(Final.all.3familyeffect.cov3.tile1000.filt.nor.unite, sample.ids=c("F159",
                                                                                     "F161",
                                                                                     "F163",
                                                                                     "F165",
                                                                                     "F166",
                                                                                     "F167",
                                                                                     "F33",
                                                                                     "F35",
                                                                                     "F39",
                                                                                     "F40",
                                                                                     "F41",
                                                                                     "F43",
                                                                                     "F45",
                                                                                     "F48",
                                                                                     "F49",
                                                                                     "F50",
                                                                                     "F73",
                                                                                     "M159",
                                                                                     "M161",
                                                                                     "M162",
                                                                                     "M164",
                                                                                     "M166",
                                                                                     "M167",
                                                                                     "M169",
                                                                                     "M34",
                                                                                     "M36",
                                                                                     "M37",
                                                                                     "M38",
                                                                                     "M41",
                                                                                     "M63",
                                                                                     "M71",
                                                                                     "M72",
                                                                                     "M73",
                                                                                     "M74",
                                                                                     "M77")
                                
                                , treatment=c(1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              2,
                                              2,
                                              2,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1,
                                              1))

Adults.PCA.family.perc.mat<-percMethylation(Adults.PCA.family)
Adults.PCA.family.perc.mat<-as.data.frame(Adults.PCA.family.perc.mat)
dim(Adults.PCA.family.perc.mat) #262268
Adults.PCA.family.perc.mat.2=Adults.PCA.family.perc.mat[ rowSums(is.na(Adults.PCA.family.perc.mat))==0, ]
dim(Adults.PCA.family.perc.mat.2) #174790      35


sd.threshold <- 0.5 ## this takes only the position
sds=rowSds(as.matrix(Adults.PCA.family.perc.mat.2))
head(sds)
cutoff=quantile(sds,sd.threshold)
Adults.PCA.family.perc.mat.2=Adults.PCA.family.perc.mat.2[sds>cutoff,]
head(Adults.PCA.family.perc.mat.2)
dim(Adults.PCA.family.perc.mat.2) # 87395    35

prcomp.Adults.PCA.family.perc.mat <-prcomp(t(Adults.PCA.family.perc.mat.2), scale. = TRUE, center = TRUE)
summary(prcomp.Adults.PCA.family.perc.mat)

myColors.8=c("black","green2")  #other family "#fee08b"
names(myColors.8) <- c(1,2)
custom_colors.8 <- scale_fill_manual(name = "Family", values = myColors.8, labels=c("Unrelated","Family 2"))

treatment.8=Adults.PCA.family@treatment
treatment.8<-as.factor(treatment.8)
sample.ids=Adults.PCA.family@sample.ids

PCA8<-ggplot(prcomp.Adults.PCA.family.perc.mat, aes(x = PC1, y =PC2, fill=treatment.8)) +
  geom_point(size=2, shape=21)+
  labs(x="PC1 11.54%", y="PC2 5.06%")+
  custom_colors.8


#### PCA Adults colored by treatment
#have to re-organize the file and assign a new treatment coding
Adults.PCA.treatment<- reorganize(Final.all.cov3.tile1000.filt.nor.unite, sample.ids=c("F159",
                                                                                       "F161",
                                                                                       "F163",
                                                                                       "F165",
                                                                                       "F166",
                                                                                       "F167",
                                                                                       "F33",
                                                                                       "F35",
                                                                                       "F39",
                                                                                       "F40",
                                                                                       "F41",
                                                                                       "F43",
                                                                                       "F45",
                                                                                       "F48",
                                                                                       "F49",
                                                                                       "F50",
                                                                                       "F73",
                                                                                       "M159",
                                                                                       "M161",
                                                                                       "M162",
                                                                                       "M164",
                                                                                       "M166",
                                                                                       "M167",
                                                                                       "M169",
                                                                                       "M34",
                                                                                       "M36",
                                                                                       "M37",
                                                                                       "M38",
                                                                                       "M41",
                                                                                       "M63",
                                                                                       "M71",
                                                                                       "M72",
                                                                                       "M73",
                                                                                       "M74",
                                                                                       "M77")
                                  
                                  , treatment=c(2,
                                                3,
                                                3,
                                                2,
                                                2,
                                                2,
                                                2,
                                                3,
                                                2,
                                                3,
                                                3,
                                                2,
                                                3,
                                                2,
                                                3,
                                                3,
                                                1,
                                                2,
                                                2,
                                                2,
                                                1,
                                                1,
                                                1,
                                                2,
                                                2,
                                                1,
                                                3,
                                                2,
                                                3,
                                                1,
                                                2,
                                                3,
                                                2,
                                                3,
                                                3))
Adults.PCA.treatment.perc.mat<-percMethylation(Adults.PCA.treatment)
Adults.PCA.treatment.perc.mat<-as.data.frame(Adults.PCA.treatment.perc.mat)
dim(Adults.PCA.treatment.perc.mat) #262268
Adults.PCA.treatment.perc.mat.2=Adults.PCA.treatment.perc.mat[ rowSums(is.na(Adults.PCA.treatment.perc.mat))==0, ]
dim(Adults.PCA.treatment.perc.mat.2) # 174790    35

sd.threshold <- 0.5 ## this takes only the position
sds=rowSds(as.matrix(Adults.PCA.treatment.perc.mat.2))
head(sds)
cutoff=quantile(sds,sd.threshold)
Adults.PCA.treatment.perc.mat.2=Adults.PCA.treatment.perc.mat.2[sds>cutoff,]
head(Adults.PCA.treatment.perc.mat.2)
dim(Adults.PCA.treatment.perc.mat.2) #  87395    35

prcomp.Adults.PCA.treatment.perc.mat <-prcomp(t(Adults.PCA.treatment.perc.mat.2), scale. = TRUE, center = TRUE)
summary(prcomp.Adults.PCA.treatment.perc.mat)

names(myColors.7) <- c(1,2,3)
custom_colors.7 <- scale_fill_manual(name = "Treatment", values = myColors.7, labels=c("WD-A","SN-A","CV-A"))

treatment.7=Adults.PCA.treatment@treatment
treatment.7<-as.factor(treatment.7)
sample.ids=Adults.PCA.treatment@sample.ids

PCA7<-ggplot(prcomp.Adults.PCA.treatment.perc.mat, aes(x = PC1, y =PC2, fill=treatment.7)) +
  geom_point(size=2, shape=21)+
  labs(x="PC1 11.54%", y="PC2 5.06%")+
  custom_colors.7




#### PCA adults colored by sex
#have to re-organize the file and assign a new treatment coding
Final.all.sexeffect.cov3.tile1000.filt.nor.unite <- readRDS("C:/Users/BOKVISTJ/Documents/UofC MSc/All Analysis December 2021/R Analysis All Sequencing Data/Final.all.sexeffect.cov3.tile1000.filt.nor.unite")
nrow(Final.all.sexeffect.cov3.tile1000.filt.nor.unite) # 342,938 

Adults.PCA.sex<- reorganize(Final.all.sexeffect.cov3.tile1000.filt.nor.unite, sample.ids=c("F159",
                                                                                 "F161",
                                                                                 "F163",
                                                                                 "F165",
                                                                                 "F166",
                                                                                 "F167",
                                                                                 "F33",
                                                                                 "F35",
                                                                                 "F39",
                                                                                 "F40",
                                                                                 "F41",
                                                                                 "F43",
                                                                                 "F45",
                                                                                 "F48",
                                                                                 "F49",
                                                                                 "F50",
                                                                                 "F73",
                                                                                 "M159",
                                                                                 "M161",
                                                                                 "M162",
                                                                                 "M164",
                                                                                 "M166",
                                                                                 "M167",
                                                                                 "M169",
                                                                                 "M34",
                                                                                 "M36",
                                                                                 "M37",
                                                                                 "M38",
                                                                                 "M41",
                                                                                 "M63",
                                                                                 "M71",
                                                                                 "M72",
                                                                                 "M73",
                                                                                 "M74",
                                                                                 "M77")
                            
                            , treatment=c(1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          1,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2,
                                          2))

Adults.PCA.sex.perc.mat<-percMethylation(Adults.PCA.sex)
Adults.PCA.sex.perc.mat<-as.data.frame(Adults.PCA.sex.perc.mat)
dim(Adults.PCA.sex.perc.mat) #262268
Adults.PCA.sex.perc.mat.2=Adults.PCA.sex.perc.mat[ rowSums(is.na(Adults.PCA.sex.perc.mat))==0, ]
dim(Adults.PCA.sex.perc.mat.2) #  174790     35

sd.threshold <- 0.5 ## this takes only the position
sds=rowSds(as.matrix(Adults.PCA.sex.perc.mat.2))
head(sds)
cutoff=quantile(sds,sd.threshold)
Adults.PCA.sex.perc.mat.2=Adults.PCA.sex.perc.mat.2[sds>cutoff,]
head(Adults.PCA.sex.perc.mat.2)
dim(Adults.PCA.sex.perc.mat.2) # 87395    35

prcomp.Adults.PCA.sex.perc.mat <-prcomp(t(Adults.PCA.sex.perc.mat.2), scale. = TRUE, center = TRUE)
summary(prcomp.Adults.PCA.sex.perc.mat)
myColors.9=c("white","black")
names(myColors.9) <- c(1,2)
custom_colors.9 <- scale_fill_manual(name = "Sex", values = myColors.9, labels=c("Female","Male"))

treatment.9=Adults.PCA.sex@treatment
treatment.9<-as.factor(treatment.9)
sample.ids=Adults.PCA.sex@sample.ids

PCA9<-ggplot(prcomp.Adults.PCA.sex.perc.mat, aes(x = PC1, y =PC2, fill=treatment.9)) +
  geom_point(size=2, shape=21)+
  labs(x="PC1 11.54%", y="PC2 5.06%")+
  custom_colors.9

plot_grid(PCA7, PCA8, PCA9, align ="V", labels = "AUTO", ncol=1)



####Scree plots for PCAs
scree1<-fviz_eig(prcomp.perc.mat, geom="bar", main = "" ,barfill = "#3182bd",
                 barcolor = "black", ncp=30, addlabels=T,xlab = "Principle components")+theme_gray()  #for the one with all adults

scree2<-fviz_eig(prcomp.Smolts.PCA.treatment.perc.mat, geom="bar", main = "" ,barfill = "#3182bd",
                 barcolor = "black", ncp=30, addlabels=T,xlab = "Principle components")+theme_gray()  #for the one with smolts only

scree3<-fviz_eig(prcomp.Adults.PCA.treatment.perc.mat, geom="bar", main = "" ,barfill = "#3182bd",
                 barcolor = "black", ncp=30, addlabels=T,xlab = "Principle components")+theme_gray()  #for the one with adults only

plot_grid(scree1, scree2, scree3, align ="V", labels = "AUTO", ncol=1)

