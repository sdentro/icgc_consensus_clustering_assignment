#' Produce per mutation observed CCF values - before clustering

source("~/repo/moritz_mut_assignment/MutationTime.R")
source("~/repo/moritz_mut_assignment/util.R")
source("~/repo/dpclust3p/R/interconvertMutationBurdens.R")

# produce the final output
indir = "output_wm"
outdir = "final_output"

#' Map SVs back onto their input location and merge with the existing assignments.
#' @param consensus_vcf_file
#' @param svid_map_file
#' @param output_sv
#' @return A list with two data.frames, one with hard assignments and one with probabilities. Every consensus SV is reported with their consensus location
remap_svs_ccf = function(consensus_vcf_file, svid_map_file, output_sv) {
  cons_sv = readVcf(consensus_vcf_file, "GRCh37")
  svmap = read.table(svid_map_file, header=T, stringsAsFactors=F)
  
  #' Helper function to split a string by multiple characters
  strsplits <- function(x, splits, ...) {
    for (split in splits) {
      x <- unlist(strsplit(x, split, ...))
    }
    return(x[!x == "" & sapply(x, function(y) !any(sapply(c("A", "C", "G", "T", "N"), grepl, y)))]) # Remove empty values
  }
  
  r = lapply(fixed(cons_sv)$ALT, function(x) strsplits(x, c("\\[", "\\]", ":", fixed=T)))
  m = matrix(unlist(r), ncol=2, byrow=T)
  all_sv_data = data.frame(chr=as.character(seqnames(cons_sv)), pos=start(cons_sv), chr2=m[,1], pos2=as.numeric(m[,2]), id=rownames(as.data.frame(rowRanges(cons_sv))),
                           type=rep(NA, nrow(m)), ccf=rep(NA, nrow(m)), mcn=rep(NA, nrow(m)), mult=rep(NA, nrow(m)), major_cn=rep(NA, nrow(m)), minor_cn=rep(NA, nrow(m)))
  
  # iterate over all svs and replace
  for (i in 1:nrow(output_sv)) {
    if (any(output_sv$position[i] == svmap$pos1)) {
      hit = which(output_sv$position[i] == svmap$pos1)
    } else {
      hit = which(output_sv$position[i] == svmap$pos2)
    }
    
    svid = unlist(strsplit(svmap$original_id[hit], "_", fixed=T))[1]
    all_sv_data_row = which(grepl(svid, all_sv_data$id))
    # now have mapped i onto all_sv_data_row, which contains both end points of the SV
    # save the assignments into the all_data tables
    
    all_sv_data$ccf[all_sv_data_row] = as.character(output_sv$ccf[i])
    all_sv_data$mcn[all_sv_data_row] = as.character(output_sv$mcn[i])
    all_sv_data$mult[all_sv_data_row] = as.character(output_sv$mult[i])
    all_sv_data$major_cn[all_sv_data_row] = as.character(output_sv$major_cn[i])
    all_sv_data$minor_cn[all_sv_data_row] = as.character(output_sv$minor_cn[i])
  }
  return(list(output_sv=all_sv_data))
}


for (infile in list.files(indir, pattern="_assignment.RData", full.names=T)) {
  load(infile)
  if (file.exists(file.path(outdir, paste0(samplename, "_mutation_ccf.txt")))) { next }
  print(samplename)
  
  get_ccf = function(vcf, mt_output, purity) {
    num_muts = nrow(mt_output)
    output = data.frame(chromosome=as.character(seqnames(vcf)), position=as.numeric(start(vcf)), type=rep(NA, num_muts), 
                        ccf=rep(NA, num_muts), major_cn=mt_output$MajCN, minor_cn=mt_output$MinCN, mcn=rep(NA, num_muts), mult=mt_output$MutCN,
                        chromosome2=rep(NA, num_muts), position2=rep(NA, num_muts))  
    output$mcn = mutationBurdenToMutationCopyNumber(burden=mt_output$altCount / (mt_output$altCount + mt_output$wtCount), cellularity=purity, normalCopyNumber=rep(2, num_muts), totalCopyNumber=output$major_cn + output$minor_cn)
    output$ccf = output$mcn / output$mult
    return(output)
  }
  
  output_snv = get_ccf(vcf_snv, MCN$D, purity)
  output_snv$type = "SNV"
  if (nrow(MCN_indel$D) > 0) {
    output_indel = get_ccf(vcf_indel, MCN_indel$D, purity)
    output_indel$type = "indel"
  } else {
    output_indel = NULL
  }
  if (!is.null(vcf_sv) && !all(is.na(MCN_sv$D$MutCN))) {
    output_sv = get_ccf(vcf_sv, MCN_sv$D, purity)
    output_sv = remap_svs_ccf(sv_vcf_file, svid_map_file, output_sv)$output_sv
    output_sv = output_sv[, c("chr", "pos", "type", "ccf", "major_cn", "minor_cn", "mcn", "mult", "chr2", "pos2")]
    colnames(output_sv) = colnames(output_snv)
    output_sv$type = "SV"
  } else {
    output_sv = NULL
  }
  
  output = rbind(output_snv, output_indel, output_sv)
  write.table(output, file=file.path(outdir, paste0(samplename, "_mutation_ccf.txt")), row.names=F, sep="\t", quote=F)
}