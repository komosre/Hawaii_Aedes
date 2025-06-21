########################
# Base Tree Aanlysis on HovenWeeP
# 6/20/25
# Adam E. Vorsino
########################

dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER")) # add to the path


rm(list = ls())

# Fackages to use
pckg <- c('rlang', 'glue','R.utils','stringr', 'stringr', 'parallel', # 'rgdal','rgeos',
          'vcfR', 'adegenet', 'poppr', 'pinfsc50',  'ggplot2', 'reshape', 'tidyr', 'phangorn')#,'ggpubr'` 'devtools'  'easystats',)

for(i in 1:length(pckg)) {

  if (!(pckg[i] %in% installed.packages())) {
    install.packages(pckg[i],repos = "http://cran.us.r-project.org", dependencies = T)
  }
  print(pckg[i])
  do.call("library", list(pckg[i]))
}
# IterToUse <- 100 # making it so the numcores makes the model 500

### BASE
Loc <- "" # Location of the work director with the .vcf.gz file
ReRun = T


procfile <- sort(list.files(path = Loc, pattern = "VARSELnumpro", full.names = T))[1]
NumCores<- read.table(procfile)[[1]]# # reads in the number of cores set in slurm code
options(mc.cores = NumCores)
options(mc.cores = Sys.getenv("SLURM_CPUS_PER_TASK"))
cat(paste0('\n CPUs:', Sys.getenv("SLURM_CPUS_PER_TASK")))

### READING IN WITH VCFR
ReRun <- F
if(file.exists(paste0(Loc, "VCFandBase.RData")) == F || ReRun == T){ 
  #Reading in vcf 
  vcf <- read.vcfR(paste0(Loc,"AalbF5_June2024_DP8_mm0.05.vcf.gz")) # AalbF3nuc_2023-02_DP8_mm0.05.vcf.gz"))
  save.image(file = paste0(Loc, "VCFandBase.RData"))
}else{
  load(paste0(Loc, "VCFandBase.RData")) # .RData preserves file names
  
}

###Getting sample names from VCF for latter use

SampleNames <- colnames(vcf@gt) # get sample names

# Loc <- "C:/Users/avorsino/OneDrive - DOI/Desktop/RCode/GitHub/AlboGenData/"
# SampleNames <- read.csv(paste0(Loc, 'SampleNames.csv'))[,2]
MetaData <- read.csv(paste0(Loc, "Albo_Pacific_Samples_11_2024_Info.csv")) # Aalb_popgen_sample metadata_AalbF3.csv"))

AllSampleMetaData <- c()
for(SN in 2:length(SampleNames)){ # The first is "FORMAT"
  # SN  <- 2
  SampNM <- SampleNames[SN]
  
  SampleMetaData <- MetaData[which(MetaData$Sample.ID == SampNM),]
  if(nrow(SampleMetaData)>1){
    cat(pasteo("\n there are more then one rows for sample name: ", SampleNames[SN]))
    break()
  }
  if(nrow(SampleMetaData)<1){
    cat(paste0("\n there are no rows with sample name: ", SampleNames[SN]))
  }
  AllSampleMetaData <- rbind(SampleMetaData, AllSampleMetaData)
}

SampleLoc <- AllSampleMetaData$Collection_Location

######  CONVERTING FROM VCF TO GENLIGHT for analysis in adegenet
ReRun = F
if(file.exists(paste0(Loc, "Albogenlight.RData")) == F || ReRun == T){ 
  AlboGenLight <- vcfR2genlight(vcf, n.cores = NumCores)
  pop(AlboGenLight) <- as.factor(SampleLoc)
  cat(paste0('\n', popNames(AlboGenLight), ' \n change start'))
  popNames(AlboGenLight) <- iconv(popNames(AlboGenLight), "UTF-8", "UTF-8",sub='')
  cat(paste0('\n', popNames(AlboGenLight), ' \n change done'))
  save(AlboGenLight, file = paste0(Loc, "Albogenlight.RData"))
}else{
  load(file = paste0(Loc, "Albogenlight.RData"))# .RData preserves file names
}

ploidy(AlboGenLight) <- 2

# Distance tree
ReRun <- T

if(file.exists(paste0(Loc, "Tree_PopGenPhyloOutput.RData")) == F || ReRun == T){ 
  cat('\n Running NJ trees using bitwise distance \n')
  treeNJ_gDNA <- aboot(AlboGenLight, tree = "nj", distance = bitwise.dist, sample = 500,
                       showtree = F, cutoff = 50, quiet = F, threads = NumCores)
  save(SampleLoc, treeNJ_gDNA, file =paste0(Loc, 'Tree_PopGenPhyloOutput.RData'))
}

if(file.exists(paste0(Loc, "NJ_UPGMA_Tree_PopGenPhyloOutput.RData")) == F || ReRun == T){ 
  cat('\n Running UPGMA trees using bitwise distance \n')
  load(paste0(Loc, 'Tree_PopGenPhyloOutput.RData'))
  treeUPGMA_gDNA <- aboot(AlboGenLight, tree = "upgma", distance = bitwise.dist, sample = 500, 
                          showtree = F, cutoff = 50, quiet = F, threads = NumCores)
  
  save(SampleLoc, treeNJ_gDNA, treeUPGMA_gDNA, file =paste0(Loc, 'NJ_UPGMA_Tree_PopGenPhyloOutput.RData'))
}


