library(stringr)


############### INPUT THE FOLLOWING PRESETS AND DIRECTORIES ###############


#specify directories to mixcr and vdjtools (optional)
path_to_mixcr <- "/software/mixcr/mixcr.jar"
path_to_vdj <- "software/vdjtools/vdjtools.jar"

#sample batch file (.tsv, .txt) of the following structure: column1 - sampleID, column2 - SMART (6 letters), column3 - R1 fastq file name (without path), column 4 - R2 fastq file name (without path)
samples <- read.csv("pipeline_trial/mixcr_preset.tsv", sep = "\t", header = F)

#set directories for input and output files
path_to_fastq <- "shared-with-me/Px_Lupus1-SpA-JIA_NextSeq-IBCH_20042024/"
path_to_output_dir <- "/pipeline_trial/"

#set preset and setup for analyze (automatic effective threshold during refineTagsAndSort vs set recordsPerConsensus during assemble)
preset_analyze <- "milab-human-rna-tcr-umi-race --species hsa -f --assemble-clonotypes-by cdr3"
preset_analyze <- "generic-amplicon --floating-left-alignment-boundary V --floating-right-alignment-boundary V
  --floating-left-alignment-boundary J --floating-right-alignment-boundary J
  --rigid-left-alignment-boundary C --rigid-right-alignment-boundary C --rna --species hsa -f --assemble-clonotypes-by cdr3"

# "auto" (str) - with effective threshold
# n (int) - with recordsPerConsensus=n and postFilter=null
setup <- "auto"

#set export preset
preset_export <- "--chains TRB --drop-default-fields -uniqueTagCount Molecule -uniqueTagFraction Molecule -nFeature CDR3 -aaFeature CDR3 -vHit -dHit -jHit -positionOf VEnd -positionOf DBegin -positionOf DEnd -positionOf JBegin"

#The default tag pattern in this script is ^N{19}(R1:*)\^tggtatcaacgcagag{SMART}(UMI:TNNNNTNNNNTNNNN)TCT(R2:*). Tag pattern can be changed on line 161

############### NO NEED TO CHANGE ANYTHING PAST THIS POINT, ONLY EXECUTE COMMANDS ###############

if (setup == "auto") {
  mode <- " "
  md <- ""
} else {
  mode <- str_glue(" -M refineTagsAndSort.parameters.postFilter=null -M assemble.consensusAssemblerParameters.assembler.minRecordsPerConsensus={setup} ")
  md <- str_glue("_NULL_T{setup}")
}


#setup analyze report
#WARNING! the column names in the exported table are offset by one column to the left (didn't come up with a way to fix that yet)
run_stats = function(folderpath){
  al_report_features <- c("Total sequencing reads",         
                          "Successfully aligned reads",
                          "Alignment failed: no hits (not TCR/IG?):", 
                          "Alignment failed: absence of V hits:", 
                          "Alignment failed: absence of J hits:", 
                          "Alignment failed: no target with both V and J alignments:", 
                          "Alignment failed: absent barcode:", 
                          "Matched reads:")
  ref_report_features <- c("UMI input diversity", 
                           "UMI output diversity", 
                           "UMI input reads:", 
                           "UMI output reads:", 
                           "Number of groups", 
                           "Number of groups accepted", 
                           "UMI mean reads per tag:",
                           "Effective threshold")
  asmbl_report_features<- c("Final clonotype count",
                            "Number of input groups",
                            "Reads used in clonotypes, percent of total",
                            "Average number of reads per clonotype",
                            "Reads dropped due to the lack of a clone sequence",
                            "Reads dropped due to low quality",
                            "Reads dropped due to failed mapping",
                            "Reads dropped with low quality clones",
                            "Reads used in clonotypes before clustering, percent of total",
                            "Number of reads used as a core, percent of used",
                            "Mapped low quality reads, percent of used",
                            "Reads clustered in PCR error correction, percent of used",
                            "Reads pre-clustered due to the similar VJC-lists, percent of used",
                            "Clonotypes eliminated by PCR error correction",
                            "Clonotypes dropped as low quality",
                            "Clonotypes pre-clustered due to the similar VJC-lists")
  
  filelist_al = list.files(path = folderpath, pattern = "*align.report.txt", full.names = T)
  align_stats = lapply(filelist_al, function(x) read.table(file = x, header=F, sep = '\t', fill = T))
  names(align_stats) = list.files(path = folderpath, pattern = "*align.report.txt", full.names = F)
  filelist_ref = list.files(path =paste(folderpath, "", sep=''), pattern = "*refine.report.txt", full.names = T)
  ref_stats = lapply(filelist_ref, function(x) read.table(file = x, header=F, sep = '\t', fill = T))
  names(ref_stats) = list.files(path = folderpath, pattern = "*refine.report.txt", full.names = F)
  filelist_as = list.files(path = folderpath, pattern = "*assemble.report.txt", full.names = T)
  asmbl_stats = lapply(filelist_as, function(x) read.table(file = x, header=F, sep = '\t', fill = T))
  names(asmbl_stats) = list.files(path = folderpath, pattern = "*assemble.report.txt", full.names = F)
  
  al = list()
  for (i in 1:length(align_stats)){
    al[[i]] = as.data.frame(t(sapply(al_report_features, 
                                     function(x) {
                                       if (length(grep(align_stats[[i]][,1],pattern = x,value = T,fixed=T))==0) {
                                         return(0)
                                       } else {
                                         list = strsplit(grep(align_stats[[i]][,1],
                                                              pattern = x,
                                                              value = T,
                                                              fixed=T),
                                                         split = ": ")
                                         return(list[[1]][length(list[[1]])])
                                       }
                                     })))
  }
  names(al) = list.files(path =paste(folderpath, "",sep = ''), pattern = "*align.report.txt", full.names = F)
  al = do.call(rbind, al)
  
  re = list()
  for (i in 1:length(ref_stats)){
    re[[i]] = as.data.frame(t(sapply(ref_report_features, 
                                     function(x) {
                                       if (length(grep(ref_stats[[i]][,1],pattern = x,value = T,fixed=T))==0) {
                                         return(0)
                                       } else {
                                         list = strsplit(grep(ref_stats[[i]][,1],
                                                              pattern = x,
                                                              value = T,
                                                              fixed=T),
                                                         split = ": ")
                                         return(list[[1]][length(list[[1]])])
                                       }
                                     })))
  }
  names(re) = list.files(path =paste(folderpath, "", sep=''), pattern = "*refine.report.txt", full.names = F)
  re = do.call(rbind, re)
  
  as = list()
  for (i in 1:length(asmbl_stats)){
    as[[i]]=as.data.frame(t(sapply(asmbl_report_features, 
                                   function(x) {
                                     if (length(grep(asmbl_stats[[i]][,1],pattern = x,value = T,fixed=T))==0) {
                                       return(0)
                                     } else {
                                       list = strsplit(grep(asmbl_stats[[i]][,1],
                                                            pattern = x,
                                                            value = T,
                                                            fixed=T),
                                                       split = ": ")
                                       return(list[[1]][length(list[[1]])])
                                     }
                                   }
    )
    )
    )
  }
  names(as) = list.files(path = folderpath, pattern = "*assemble.report.txt", full.names = F)
  as = do.call(rbind, as)
  
  all <- cbind(al, re, as)
  write.table(t(all), file = paste(folderpath,'mixcr_stats.tsv', sep='/'), quote = F, sep = '\t', row.names = T, col.names = T)
  
}


############### PIPELINE STARTS HERE ###############


#preform analyze (align, refine and assemble)
for (i in 1:nrow(samples)) {
  tagPattern <- paste(r'(--tag-pattern "^N{19}(R1:*)\^)', str_glue('tggtatcaacgcagag{samples$V2[i]}(UMI:TNNNNTNNNNTNNNN)TCT(R2:*)\"'), sep ="")
  R1 <- str_glue("{path_to_fastq}{samples$V3[i]}")
  R2 <- str_glue("{path_to_fastq}{samples$V4[i]}")
  output <- str_glue("{path_to_output_dir}{samples$V1[i]}{md}")
  system(str_glue("java -jar {path_to_mixcr} analyze {preset_analyze}{mode}{tagPattern} {R1} {R2} {output}"))
}


#generate analyze report from align, refine and assemble reports
#WARNING! the column names in the exported table are offset by one column to the left (didn't come up with a way to fix that yet)
run_stats(str_glue("{path_to_output_dir}"))

#export .clns into _export.tsv 
for (i in 1:nrow(samples)) {
  input <- str_glue("{path_to_output_dir}{samples$V1[i]}{md}.clns")
  output <- str_glue("{path_to_output_dir}{samples$V1[i]}{md}_export.tsv")
  system(str_glue("java -jar {path_to_mixcr} exportClones {preset_export} {input} {output}"))
}

#make .tsv VDJtools friendly
for (i in 1:nrow(samples)) {
  rm(ctable)
  ctable <- read.csv(paste(str_glue("{path_to_output_dir}{samples$V1[i]}{md}_export"), ".tsv", sep = ""), sep = '\t')
  colnames(ctable) <- c("count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J", "Vend", "Dstart", "Dend", "Jstart")
  write.table(ctable, file = paste(str_glue("{path_to_output_dir}{samples$V1[i]}{md}_export"), ".txt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
}

#compute basic statistics and create a joint file for all samples 
for (i in 1:nrow(samples)) {
  system(str_glue("java -jar {path_to_vdj} CalcBasicStats {path_to_output_dir}{samples$V1[i]}{md}_export.txt {path_to_output_dir}{samples$V1[i]}{md}"))
}
rm(stats)
for (i in 1:nrow(samples)) {
  if (i == 1) {
    stats <- read.table(str_glue("{path_to_output_dir}{samples$V1[i]}{md}.basicstats.txt"), sep = '\t', header = TRUE)
  } else {
    stats <- rbind(read.table(str_glue("{path_to_output_dir}{samples$V1[i]}{md}.basicstats.txt"), sep = '\t', header = TRUE), stats)
  }
}
stats$metadata_blank <- NULL
write.table(stats, file = str_glue("{path_to_output_dir}basic_stats.tsv"), quote = F, sep = '\t', row.names = F, col.names = T)


############### END OF PIPELINE ###############
