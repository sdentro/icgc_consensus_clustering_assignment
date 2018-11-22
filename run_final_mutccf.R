#' Produce per mutation observed CCF values - before clustering
########################################################################
# Command line options
########################################################################
library(optparse)
option_list = list(
  make_option(c("-l", "--libpath"), type="character", default=NULL, help="Path to pipeline installation directory", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL, help="Directory where input (the output of the main pipeline) is stored", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default=NULL, help="Directory where output will be written", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

libpath = opt$libpath
indir = opt$input
outputdir = opt$outputdir

library(MutationTimeR)
source(file.path(libpath, "util.R"))

########################################################################
# Functions
########################################################################
get_ccf = function(vcf, mt_output, purity) {
  num_muts = nrow(mt_output)
  output = data.frame(chromosome=as.character(seqnames(vcf)), position=as.numeric(start(vcf)), type=rep(NA, num_muts), 
                      ccf=rep(NA, num_muts), major_cn=mt_output$MajCN, minor_cn=mt_output$MinCN, mcn=rep(NA, num_muts), mult=mt_output$MutCN,
                      chromosome2=rep(NA, num_muts), position2=rep(NA, num_muts))  
  output$mcn = mutationBurdenToMutationCopyNumber(burden=mt_output$altCount / (mt_output$altCount + mt_output$wtCount), 
                                                  cellularity=purity, 
                                                  normalCopyNumber=rep(2, num_muts), 
                                                  totalCopyNumber=output$major_cn + output$minor_cn)
  output$ccf = output$mcn / output$mult
  return(output)
}

########################################################################
# Create the output
########################################################################
for (infile in list.files(indir, pattern="_assignment.RData", full.names=T)) {
  load(infile)
  if (file.exists(file.path(outdir, paste0(samplename, "_mutation_ccf.txt")))) { next }
  print(samplename)
  
  output_snv = get_ccf(vcf_snv, MCN$D, purity)
  output_snv$type = "SNV"
  if (nrow(MCN_indel$D) > 0) {
    output_indel = get_ccf(vcf_indel, MCN_indel$D, purity)
    output_indel$type = "indel"
  } else {
    output_indel = NULL
  }
  if (!is.null(vcf_sv) && !all(is.na(MCN_sv$D$MutCN))) {
    output_sv = rbind(get_ccf(vcf_sv, MCN_sv$D, purity), get_ccf(vcf_sv_alt, MCN_sv_alt$D, purity))
    # sync chromosome and position
    output_sv[, c("chromosome", "position", "chromosome2", "position2")] = sv_output[,c("chromosome", "position", "chromosome2", "position2")]
    output_sv$type = "SV"
  } else {
    output_sv = NULL
  }
  
  output = rbind(output_snv, output_indel, output_sv)
  write.table(output, file=file.path(outputdir, paste0(samplename, "_mutation_ccf.txt")), row.names=F, sep="\t", quote=F)
}
