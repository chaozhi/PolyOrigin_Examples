# Installing and using the PedigreeSimR to simulate the cross and GBS data
# devtools::install_github("rramadeu/PedigreeSimR")
library(PedigreeSimR)
packageVersion("PedigreeSimR")

library(rstudioapi)    
workdir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
getwd()

simmaphaplo <- function(nsnp,np,ploidy=4,chrlen=100,seed=1234){
  haplo = fake_haplo(n=np*ploidy,m=nsnp,seed=seed)
  cm0=rexp(nsnp-1,rate=1)
  cm0=cm0/(sum(cm0)/chrlen)
  cm0=round(cm0,digits=2)
  cm=c(0,cumsum(cm0))
  list(map=cm,hap=haplo)
}

ploidy=4
nsnp=120
np=5
popsize=200
sim=simmaphaplo(nsnp,np,ploidy=ploidy)
# pedigree = star_pedigree(parents=np,popsize=popsize)
# pedigree =linear_pedigree(parents=np,popsize=popsize)
# pedigree = round_pedigree(parents=np,popsize=popsize)
pedigree = diallel_pedigree(parents=np,selfs=0,popsize=popsize)
pedigreesimR(sim$map,sim$hap,
             ploidy=ploidy,
             sampleHap = FALSE,filename="",
             pedigree,workingfolder = getwd(),
             missingFreq=c(0.1,0.1),
             GBS=TRUE,GBSavgdepth = 20,GBSseq = 0.001
)
# delete all pedsim*.* files. 
sapply(list.files(pattern = "pedsim*"), unlink)
sapply(list.files(pattern = "snparray*"), unlink)
unlink("*truevalue.csv")
unlink("*geno.csv")


file.rename(list.files(pattern="*geno_GBS.csv")[1], "tetraploid_simgbs_geno.csv")
file.rename(list.files(pattern="*truevalue_ancestral.csv")[1], "tetraploid_simgbs_true.csv")
file.rename(list.files(pattern="*pedigree.csv")[1], "tetraploid_simgbs_ped.csv")

