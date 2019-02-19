#' Produce per mutation probabilities of being gained
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
# Build the timing
########################################################################
timing_disagree = data.frame()
for (infile in list.files(indir, pattern="_assignment.RData", full.names=T)) {
  # clear anything from previous sample
  try(rm(MCN_indel, MCN, MCN_sv, MCN_sv_alt, vcf_snv, vcf_indel, vcf_sv, vcf_sv_alt), silent=T)
  
  count_timing_disagree = 0
  
  load(infile)
  print(samplename)
  snv_timing = data.frame(chromosome=as.character(seqnames(vcf_snv)),
                          position=as.numeric(start(vcf_snv)),
                          mut_type=rep("SNV", nrow(MCN$D)),
                          timing=MutationTimeR:::classifyMutations(MCN$D),
                          chromosome2=rep(NA, nrow(MCN$D)),
                          position2=rep(NA, nrow(MCN$D)),
                          svid=rep(NA, nrow(MCN$D)),
                          prob_gained=MCN$D$pGain,
                          prob_not_gained=MCN$D$pSingle,
                          prob_subclonal=MCN$D$pSub,
                          stringsAsFactors=F)
  
  if (!is.null(vcf_indel)) {
    indel_timing = data.frame(chromosome=as.character(seqnames(vcf_indel)),
                              position=as.numeric(start(vcf_indel)),
                              mut_type=rep("indel", nrow(MCN_indel$D)),
                              timing=MutationTimeR:::classifyMutations(MCN_indel$D),
                              chromosome2=rep(NA, nrow(MCN_indel$D)),
                              position2=rep(NA, nrow(MCN_indel$D)),
                              svid=rep(NA, nrow(MCN_indel$D)),
                              prob_gained=MCN_indel$D$pGain,
                              prob_not_gained=MCN_indel$D$pSingle,
                              prob_subclonal=MCN_indel$D$pSub,
                              stringsAsFactors=F)
  }
  
  if (!is.null(vcf_sv)) {
    # use one SV breakpoint for assignment
    sv_timing = data.frame(chromosome=as.character(seqnames(vcf_sv)),
                           position=start(vcf_sv),
                           mut_type=rep("SV", nrow(MCN_sv$D)),
                           timing=MutationTimeR:::classifyMutations(MCN_sv$D),
                           chromosome2=info(vcf_sv)$chr2,
                           position2=info(vcf_sv)$pos2,
                           svid=info(vcf_sv)$id,
                           prob_gained=MCN_sv$D$pGain,
                           prob_not_gained=MCN_sv$D$pSingle,
                           prob_subclonal=MCN_sv$D$pSub,
                           stringsAsFactors=F)
    
    if (!is.null(vcf_sv_alt)) {
      # here we use the other SV breakpoint for a separate assignment
      sv_alt_timing = data.frame(chromosome=as.character(seqnames(vcf_sv_alt)),
                                 position=start(vcf_sv_alt),
                                 mut_type=rep("SV", nrow(MCN_sv_alt$D)),
                                 timing=MutationTimeR:::classifyMutations(MCN_sv_alt$D),
                                 chromosome2=info(vcf_sv_alt)$chr2,
                                 position2=info(vcf_sv_alt)$pos2,
                                 svid=info(vcf_sv_alt)$id,
                                 prob_gained=MCN_sv_alt$D$pGain,
                                 prob_not_gained=MCN_sv_alt$D$pSingle,
                                 prob_subclonal=MCN_sv_alt$D$pSub,
                                 stringsAsFactors=F)

      sv_timing = rbind(sv_timing, sv_alt_timing)
      
    } else {
      # Remap SVs into their correct position
      res = remap_svs(sv_vcf_file, svclone_file, final_pcawg11_output$sv_assignments, final_pcawg11_output$sv_assignments_prob, sv_timing)
      final_pcawg11_output$sv_assignments = res$sv_assignments
      final_pcawg11_output$sv_assignments_prob = res$sv_assignments_prob
      sv_timing = res$sv_timing
    }
    
    
    sv_output = data.frame(chromosome=as.character(final_pcawg11_output$sv_assignments$chr),
                           position=final_pcawg11_output$sv_assignments$pos,
                           mut_type=rep("SV", nrow(final_pcawg11_output$sv_assignments_prob)),
                           final_pcawg11_output$sv_assignments_prob[, grepl("cluster", colnames(final_pcawg11_output$sv_assignments_prob)), drop=F],
                           chromosome2=final_pcawg11_output$sv_assignments$chr2,
                           position2=final_pcawg11_output$sv_assignments$pos2,
                           svid=final_pcawg11_output$sv_assignments$id,
                           stringsAsFactors=F)
    
    # If the sample is not clonal, then we look to merge the SV breakpoint assignment probabilities (i.e. each SV up until now
    # has two probabilities, one for each breakpoint, here they are merged into one)
    if (nrow(final_pcawg11_output$final_clusters) > 1) {
      before = colSums(sv_output[, grepl("cluster_", colnames(sv_output))], na.rm=T) / 2
      masked = c()
      # problematic = data.frame()
      # assignment probabilties for some sv breakpoint pairs differ a lot. here those are masked, because we don't really know where these belong
      for (svid in unlist(lapply(sv_output$svid, function(x) unlist(strsplit(x, "_"))[1]))) {
        probs_bp_a = sv_output[sv_output$svid==paste0(svid, "_1"), grepl("cluster_", colnames(sv_output))]
        probs_bp_b = sv_output[sv_output$svid==paste0(svid, "_2"), grepl("cluster_", colnames(sv_output))]
        prob_diff = mean(as.numeric(abs(probs_bp_a-probs_bp_b)))
        
        if (prob_diff > max_allowed_sv_prob_diff | is.na(prob_diff)) {
          # mask probabilities of breakpoint pairs that differ by more than the established threshold
          sv_output[grepl(svid, sv_output$svid), grepl("cluster_", colnames(sv_output))] = NA
          sv_timing[sv_timing$svid==paste0(svid, "_1"), "timing"] = NA
          sv_timing[sv_timing$svid==paste0(svid, "_2"), "timing"] = NA
          sv_timing[sv_timing$svid==paste0(svid, "_1"), grepl("prob_", colnames(sv_timing))] = NA
          sv_timing[sv_timing$svid==paste0(svid, "_2"), grepl("prob_", colnames(sv_timing))] = NA
          masked = c(masked, svid)
        } else {
          # Assignment probabilities - if there is a good concordance, combine the two separate estimates into one
          probs_combined = sapply(1:length(probs_bp_a), function(i) mean(c(as.numeric(probs_bp_a[i]), as.numeric(probs_bp_b[i])), na.rm=T))
          if (sum(probs_combined) >= (1+.Machine$double.eps)) {
            print(paste0(svid, " probs not equal 1, diff is ", 1-sum(probs_combined)))
          }
          sv_output[sv_output$svid==paste0(svid, "_1"), grepl("cluster_", colnames(sv_output))] = probs_combined
          sv_output[sv_output$svid==paste0(svid, "_2"), grepl("cluster_", colnames(sv_output))] = probs_combined
          
          # Synchronise the timing information - probabilities first
          timing_a = sv_timing[sv_timing$svid==paste0(svid, "_1"), grepl("prob_", colnames(sv_timing))]
          timing_b = sv_timing[sv_timing$svid==paste0(svid, "_2"), grepl("prob_", colnames(sv_timing))]
          timing_combined = sapply(1:length(timing_a), function(i) mean(c(as.numeric(timing_a[i]), as.numeric(timing_b[i])), na.rm=T))
          timing_combined[!is.finite(timing_combined)] = NA
          sv_timing[sv_timing$svid==paste0(svid, "_1"), grepl("prob_", colnames(sv_timing))] = timing_combined
          sv_timing[sv_timing$svid==paste0(svid, "_2"), grepl("prob_", colnames(sv_timing))] = timing_combined
          
          # If the original classification did not agree, here we set the timing verdict to NA
          # this is because timing also depends on the copy number, which remains unchanged
          # it is therefore not possible to get a joint timing estimate that describes both breakpoints
          if (!is.na(sv_timing[sv_timing$svid==paste0(svid, "_1"), "timing"]) & sv_timing[sv_timing$svid==paste0(svid, "_1"), "timing"] != sv_timing[sv_timing$svid==paste0(svid, "_2"), "timing"]) {
            sv_timing[sv_timing$svid==paste0(svid, "_1"), "timing"] = NA
            sv_timing[sv_timing$svid==paste0(svid, "_2"), "timing"] = NA
            sv_timing[sv_timing$svid==paste0(svid, "_1"), grepl("prob_", colnames(sv_timing))] = NA
            sv_timing[sv_timing$svid==paste0(svid, "_2"), grepl("prob_", colnames(sv_timing))] = NA
            
            count_timing_disagree = count_timing_disagree + 1
            # problematic = c(problematic, svid)
            
            res = data.frame(svid=svid,
                             timing[which(grepl(paste(svid, "_", sep=""), timing$svid)),c("timing", "prob_gained", "prob_not_gained", "prob_subclonal")],
                             output[which(grepl(paste(svid, "_", sep=""), timing$svid)),c("ccf", "major_cn", "minor_cn", "mcn", "mult")],
                             assign_probs[which(grepl(paste(svid, "_", sep=""), timing$svid)),grepl("cluster_", colnames(assign_probs))])
            # problematic = rbind(problematic, res)
          }
        }
      }
      
      # re-establish SV cluster sizes
      final_pcawg11_output$final_clusters[,c("cluster", "proportion", "ccf", "n_snvs", "n_indels", "n_svs")]
      # if there are no SVs with probabilities anymore, then set the cluster sizes to NA
      if (all(is.na(sv_output[, grepl("cluster_", colnames(sv_output))]))) {
        final_pcawg11_output$final_clusters[,c("n_svs")] = NA
      } else {
        final_pcawg11_output$final_clusters[,c("n_svs")] = colSums(sv_output[, grepl("cluster_", colnames(sv_output))], na.rm=T)
      }
    }
    
    # add missing SV ids that did not make it into the SVclone output
    sv_timing = add_missing_entries(sv_vcf_file, "GRCh37", sv_timing)
    sv_output = add_missing_entries(sv_vcf_file, "GRCh37", sv_output)
    
    timing = do.call(rbind, list(snv_timing, indel_timing, sv_timing))
    assign_probs = do.call(rbind, list(snv_output, indel_output, sv_output))
    
  } else {
    sv_output = NULL
    timing = rbind(snv_timing, indel_timing)
    assign_probs = rbind(snv_output, indel_output)
  }
  timing[, c("chromosome", "position", "chromosome2", "position2")] = assign_probs[, c("chromosome", "position", "chromosome2", "position2")]
  write.table(timing, file=file.path(outputdir, paste0(samplename, "_prob_gained.txt")), quote=F, row.names=F, sep="\t")
  
  timing_disagree = rbind(timing_disagree, data.frame(samplename=samplename, count=count_timing_disagree, stringsAsFactors=F))
}
write.table(timing_disagree, file=file.path(outputdir, paste0(samplename, "_timing_disagree.txt")), quote=F, row.names=F, sep="\t")



svid = "SVMERGE727"



