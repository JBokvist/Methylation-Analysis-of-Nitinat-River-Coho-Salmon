#pipeline for pairwise comparisons of methylation in the DSS package and creation of figures

install.packages("UpSetR")
library(UpSetR)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

package.version(UpSetR)
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

library("DSS")
library("bsseq")
library(GenomicRanges)
library(ComplexUpset)
library(tidyr)

#load objects
DSSmodel.all.3family<- readRDS("DSSmodel.all.globus.3family")
#(Intercept) TreatmentADTRAD TreatmentADWILD TreatmentSMENRC TreatmentSMTRAD TreatmentSMWILD FamilyF2 FamilyF3 SexM
DSSmodel.all.3family$formula

#contrast matrices
#compare Smolts ENRC vs. WD
contrast.smolts.ENRCvsWD <- matrix(c(0,0,0,1,0,-1,0,0,0), ncol=1)

#Smolts ENRC vs. TD
contrast.smolts.ENRCvsTD <- matrix(c(0,0,0,1,-1,0,0,0,0), ncol=1)

#Smolts WD vs. TD
contrast.smolts.WDvsTD <- matrix(c(0,0,0,0,1,-1,0,0,0), ncol=1)


#COMPARISON OF ADULTS
#compare adults ENRC vs. WD
contrast.adults.ENRCvsWD <- matrix(c(0,0,1,0,0,0,0,0,0), ncol=1)

#Smolts ENRC vs. TD
contrast.adults.ENRCvsTD <- matrix(c(0,1,0,0,0,0,0,0,0), ncol=1)

#Smolts WD vs. TD
contrast.adults.WDvsTD <- matrix(c(0,1,-1,0,0,0,0,0,0), ncol=1)
#COMPARISON OF ADULTS
#compare adults ENRC vs. WD
contrast.adults.ENRCvsWD <- matrix(c(0,0,1,0,0,0,0,0,0), ncol=1)

#Smolts ENRC vs. TD
contrast.adults.ENRCvsTD <- matrix(c(0,1,0,0,0,0,0,0,0), ncol=1)

#Smolts WD vs. TD
contrast.adults.WDvsTD <- matrix(c(0,1,-1,0,0,0,0,0,0), ncol=1)


#COMPARISON OF SMOLTS TO ADULTS
#compare smolt to adult WD
contrast.WD.SMvsAD <- matrix(c(0,0,1,0,0,-1,0,0,0), ncol=1)

#compare smolt to adult ENRC
contrast.ENRC.SMvsAD <- matrix(c(0,0,0,1,0,0,0,0,0), ncol=1)

#compare smolt to adult TD
contrast.TD.SMvsAD <- matrix(c(0,1,0,0,-1,0,0,0,0), ncol=1)

#EFFECT OF BEING MALE VS FEMALE
contrast.SEXM <- matrix(c(0,0,0,0,0,0,0,0,1), ncol=1)


#Get DML's and add adjusted p-value columns
smolts.ENRCvsWD<- DMLtest.multiFactor(DSSmodel.all.3family, Contrast=contrast.smolts.ENRCvsWD) #10366649 obs. of 5 variables
smolts.ENRCvsTD<- DMLtest.multiFactor(DSSmodel.all.3family, Contrast=contrast.smolts.ENRCvsTD) 
smolts.WDvsTD<- DMLtest.multiFactor(DSSmodel.all.3family, Contrast=contrast.smolts.WDvsTD) 
adults.ENRCvsWD<- DMLtest.multiFactor(DSSmodel.all.3family, Contrast=contrast.adults.ENRCvsWD)
adults.ENRCvsTD<- DMLtest.multiFactor(DSSmodel.all.3family, Contrast=contrast.adults.ENRCvsTD)
adults.WDvsTD<- DMLtest.multiFactor(DSSmodel.all.3family, Contrast=contrast.adults.WDvsTD)
WD.SMvsAD<- DMLtest.multiFactor(DSSmodel.all.3family, Contrast=contrast.WD.SMvsAD)
ENRC.SMvsAD<- DMLtest.multiFactor(DSSmodel.all.3family, Contrast=contrast.ENRC.SMvsAD)
TD.SMvsAD<- DMLtest.multiFactor(DSSmodel.all.3family, Contrast=contrast.TD.SMvsAD)
SEXM<- DMLtest.multiFactor(DSSmodel.all.3family, Contrast=contrast.SEXM)
FAMILY<- DMLtest.multiFactor(DSSmodel.all.3family, term="Family")


#CALL SIGNIFICANT DMRS

DMR50.adults.ENRCvsWD <- callDMR(DMLresult=adults.ENRCvsWD,minlen=50,p.threshold=0.00001,minCG=3, dis.merge=50, pct.sig=0.5) #12 
DMR50.adults.ENRCvsTD <- callDMR(DMLresult=adults.ENRCvsTD,minlen=50,p.threshold=0.00001,minCG=3, dis.merge=50, pct.sig=0.5) #2
DMR50.adults.WDvsTD <- callDMR(DMLresult=adults.WDvsTD,minlen=50,p.threshold=0.00001,minCG=3, dis.merge=50, pct.sig=0.5) #191 
DMR50.WD.SMvsAD <- callDMR(DMLresult=WD.SMvsAD,minlen=50,p.threshold=0.00001,minCG=3, dis.merge=50, pct.sig=0.5) #564 
DMR50.ENRC.SMvsAD <- callDMR(DMLresult=ENRC.SMvsAD,minlen=50,p.threshold=0.00001,minCG=3, dis.merge=50, pct.sig=0.5)#648
DMR50.TD.SMvsAD <- callDMR(DMLresult=TD.SMvsAD,minlen=50,p.threshold=0.00001,minCG=3, dis.merge=50, pct.sig=0.5) #695 
DMR50.SEXM <- callDMR(DMLresult=SEXM,minlen=50,p.threshold=0.00001,minCG=3, dis.merge=50, pct.sig=0.5) #31 
DMR50.FAMILY <- callDMR(DMLresult=FAMILY,minlen=50,p.threshold=0.00001,minCG=3, dis.merge=50, pct.sig=0.5) #26 


##-----------------------------------------UpSetR Plots from DSS Objects----------------------------------------##

#want to compare between contrasts if the same scaffolds were compared

#with genomic regions
DMR50.smolts.ENRCvsWD.genomic <- DMR50.smolts.ENRCvsWD[ ,1]
DMR50.smolts.ENRCvsTD.genomic <- DMR50.smolts.ENRCvsTD[ ,1]
DMR50.smolts.WDvsTD.genomic <- DMR50.smolts.WDvsTD[ ,1]
DMR50.adults.ENRCvsWD.genomic <- DMR50.adults.ENRCvsWD[ ,1]
DMR50.adults.ENRCvsTD.genomic <- DMR50.adults.ENRCvsTD[ ,1]
DMR50.adults.WDvsTD.genomic <- DMR50.adults.WDvsTD[ ,1]
DMR50.WD.SMvsAD.genomic <- DMR50.WD.SMvsAD[ ,1]
DMR50.ENRC.SMvsAD.genomic <- DMR50.ENRC.SMvsAD[ ,1]
DMR50.TD.SMvsAD.genomic <- DMR50.TD.SMvsAD[ ,1]
DMR50.SEX.genomic <- DMR50.SEXM[ ,1]
DMR50.FAMILY.genomic <- DMR50.FAMILY[ ,1]

list.genomic<-list(DMR50.smolts.ENRCvsWD.genomic, DMR50.smolts.ENRCvsTD.genomic, DMR50.smolts.WDvsTD.genomic, DMR50.adults.ENRCvsWD.genomic, DMR50.adults.ENRCvsTD.genomic, DMR50.adults.WDvsTD.genomic,
                   DMR50.WD.SMvsAD.genomic, DMR50.ENRC.SMvsAD.genomic, DMR50.TD.SMvsAD.genomic, DMR50.SEX.genomic, DMR50.FAMILY.genomic)
names(list.genomic) <-c("S SNvsWD", "S SNvsCV", "S WDvsCV","A SNvsWD", "A SNvsCV", "A WDvsCV", "WD SvsA","SN SvsA","CV SvsA", "SEX", "FAMILY")

matrix.genomic.chromosome = make_comb_mat(list.genomic)

                                                                      
ht = draw(UpSet(matrix.genomic.chromosome, right_annotation = upset_right_annotation(matrix.genomic.chromosome,
                                                                                     gp = gpar(fill = "#3182bd"))))          
od = column_order(ht)
cs = comb_size(matrix.genomic.chromosome,)
#cs = formatC(cs, format = "e", digits = 0)

rod = row_order(ht)
ss = set_size(matrix.genomic.chromosome)
ss = formatC(ss, format = "e", digits = 0)

decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(4, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})


#now make the same plot for the shared DMRs
DMR50.smolts.ENRCvsWD.genomic <- DMR50.smolts.ENRCvsWD
DMR50.smolts.ENRCvsWD.genomic <- DMR50.smolts.ENRCvsWD.genomic %>%
  unite('Merged', c("start","end"), sep="-", remove=FALSE)

DMR50.smolts.ENRCvsTD.genomic <- DMR50.smolts.ENRCvsTD
DMR50.smolts.ENRCvsTD.genomic <- DMR50.smolts.ENRCvsTD.genomic %>%
  unite('Merged', c("start","end"), sep="-", remove=FALSE)

DMR50.smolts.WDvsTD.genomic <- DMR50.smolts.WDvsTD
DMR50.smolts.WDvsTD.genomic <- DMR50.smolts.WDvsTD.genomic %>%
  unite('Merged', c("start","end"), sep="-", remove=FALSE)

DMR50.adults.ENRCvsWD.genomic <- DMR50.adults.ENRCvsWD
DMR50.adults.ENRCvsWD.genomic <- DMR50.adults.ENRCvsWD.genomic %>%
  unite('Merged', c("start","end"), sep="-", remove=FALSE)

DMR50.adults.ENRCvsTD.genomic <- DMR50.adults.ENRCvsTD
DMR50.adults.ENRCvsTD.genomic <- DMR50.adults.ENRCvsTD.genomic %>%
  unite('Merged', c("start","end"), sep="-", remove=FALSE)

DMR50.adults.WDvsTD.genomic <- DMR50.adults.WDvsTD
DMR50.adults.WDvsTD.genomic <- DMR50.adults.WDvsTD.genomic %>%
  unite('Merged', c("start","end"), sep="-", remove=FALSE)

DMR50.WD.SMvsAD.genomic <- DMR50.WD.SMvsAD
DMR50.WD.SMvsAD.genomic <- DMR50.WD.SMvsAD.genomic %>%
  unite('Merged', c("start","end"), sep="-", remove=FALSE)

DMR50.ENRC.SMvsAD.genomic <- DMR50.ENRC.SMvsAD
DMR50.ENRC.SMvsAD.genomic <- DMR50.ENRC.SMvsAD.genomic %>%
  unite('Merged', c("start","end"), sep="-", remove=FALSE)

DMR50.TD.SMvsAD.genomic <- DMR50.TD.SMvsAD
DMR50.TD.SMvsAD.genomic <- DMR50.TD.SMvsAD.genomic %>%
  unite('Merged', c("start","end"), sep="-", remove=FALSE)

DMR50.SEX.genomic <- DMR50.SEXM
DMR50.SEX.genomic <- DMR50.SEX.genomic %>%
  unite('Merged', c("start","end"), sep="-", remove=FALSE)

DMR50.FAMILY.genomic <- DMR50.FAMILY
DMR50.FAMILY.genomic <- DMR50.FAMILY.genomic %>%
  unite('Merged', c("start","end"), sep="-", remove=FALSE)


test_gr.1=DMR50.smolts.ENRCvsWD.genomic$Merged
test_gr.2=DMR50.smolts.ENRCvsTD.genomic$Merged
test_gr.3=DMR50.smolts.WDvsTD.genomic$Merged
test_gr.4=DMR50.adults.ENRCvsWD.genomic$Merged
test_gr.5=DMR50.adults.ENRCvsTD.genomic$Merged
test_gr.6=DMR50.adults.WDvsTD.genomic$Merged
test_gr.7=DMR50.WD.SMvsAD.genomic$Merged
test_gr.8=DMR50.ENRC.SMvsAD.genomic$Merged
test_gr.9=DMR50.TD.SMvsAD.genomic$Merged
test_gr.10=DMR50.SEX.genomic$Merged
test_gr.11=DMR50.FAMILY.genomic$Merged

list.genomic<-list(test_gr.1,test_gr.2,test_gr.3,test_gr.4,test_gr.5,test_gr.6,test_gr.7,test_gr.8,test_gr.9,test_gr.10,test_gr.11)
names(list.genomic) <-c("S SNvsWD", "S SNvsCV", "S WDvsCV","A SNvsWD", "A SNvsCV", "A WDvsCV", "WD SvsA","SN SvsA","CV SvsA", "SEX", "FAMILY")

matrix.genomic.DMR = make_comb_mat(list.genomic)

ht = draw(UpSet(matrix.genomic.DMR, right_annotation = upset_right_annotation(matrix.genomic.DMR,
                                                                              gp = gpar(fill = "#3182bd"))))          
od = column_order(ht)
cs = comb_size(matrix.genomic.DMR,)
#cs = formatC(cs, format = "e", digits = 0)

# rod = row_order(ht)
# ss = set_size(matrix.genomic.DMR)
# #ss = formatC(ss, format = "e", digits = 0)

decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(4, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})





