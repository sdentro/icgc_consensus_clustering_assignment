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
  if ("altCount" %in% colnames(mt_output)) {
    output$mcn = mutationBurdenToMutationCopyNumber(burden=mt_output$altCount / (mt_output$altCount + mt_output$wtCount), 
                                                    cellularity=purity, 
                                                    normalCopyNumber=rep(2, num_muts), 
                                                    totalCopyNumber=output$major_cn + output$minor_cn)
  } else {
    output$mcn = NA
  }
  output$ccf = output$mcn / output$mult
  return(output)
}

########################################################################
# Create the output
########################################################################
for (infile in list.files(indir, pattern="_assignment.RData", full.names=T)) {
  # clear anything from previous sample
  try(rm(MCN_indel, MCN, MCN_sv, MCN_sv_alt, vcf_snv, vcf_indel, vcf_sv, vcf_sv_alt), silent=T)
  
  load(infile)
  if (file.exists(file.path(outdir, paste0(samplename, "_mutation_ccf.txt")))) { next }
  print(samplename)
  
  output_snv = get_ccf(vcf_snv, MCN$D, purity)
  output_snv$type = "SNV"
  output_snv$svid = NA
  if (!is.null(vcf_indel) && nrow(MCN_indel$D) > 0) {
    output_indel = get_ccf(vcf_indel, MCN_indel$D, purity)
    output_indel$type = "indel"
    output_indel$svid = NA
  } else {
    output_indel = NULL
  }
  if (!is.null(vcf_sv) && !all(is.na(MCN_sv$D$MutCN))) {
    output_sv = rbind(get_ccf(vcf_sv, MCN_sv$D, purity), get_ccf(vcf_sv_alt, MCN_sv_alt$D, purity))
    # sync chromosome and position
    output_sv[, c("chromosome", "position", "chromosome2", "position2")] = sv_output[,c("chromosome", "position", "chromosome2", "position2")]
    output_sv$type = "SV"
    output_sv$svid = sv_output$svid
    # add missing SV ids that did not make it into the SVclone output
    output_sv = add_missing_entries(sv_vcf_file, "GRCh37", output_sv)
  } else {
    output_sv = NULL
  }
  
  output = rbind(output_snv, output_indel, output_sv)
  write.table(output, file=file.path(outputdir, paste0(samplename, "_mutation_ccf.txt")), row.names=F, sep="\t", quote=F)
}
