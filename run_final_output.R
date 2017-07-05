#' Produce PCAWG-wide output files

source("~/repo/moritz_mut_assignment/MutationTime.R")
source("~/repo/moritz_mut_assignment/util.R")
source("~/repo/dpclust3p/R/interconvertMutationBurdens.R")

# produce the final output
indir = "output_wm"
outdir = "final_output"

for (infile in list.files(indir, pattern="_assignment.RData", full.names=T)) {
  load(infile)
  print(samplename)
  subcl_struct = final_pcawg11_output$final_clusters[,c("cluster", "proportion", "ccf", "n_snvs", "n_indels", "n_svs")]
  colnames(subcl_struct)[2] = "fraction_total_cells"
  colnames(subcl_struct)[3] = "fraction_cancer_cells"
  
  if (!is.null(vcf_sv)) {
    # Bug in pipeline, fixed post-hoc
    sv_output = data.frame(chromosome=final_pcawg11_output$sv_assignments_prob$chr,
                           position=final_pcawg11_output$sv_assignments_prob$pos,
                           mut_type=rep("SV", nrow(final_pcawg11_output$sv_assignments_prob)),
                           final_pcawg11_output$sv_assignments_prob[, grepl("cluster", colnames(final_pcawg11_output$sv_assignments_prob)), drop=F],
                           chromosome2=final_pcawg11_output$sv_assignments_prob$chr2,
                           position2=final_pcawg11_output$sv_assignments_prob$pos2,
                           stringsAsFactors=F)
    assign_probs = do.call(rbind, list(snv_output, indel_output, sv_output))
  }
  
  for (i in which(is.na(timing$timing))) {
    assign_probs[i, grepl("cluster", colnames(assign_probs))] = NA
  }
  
  write.table(subcl_struct, file=file.path(outdir, paste0(samplename, "_subclonal_structure.txt")), quote=F, row.names=F, sep="\t")
  write.table(assign_probs, file=file.path(outdir, paste0(samplename, "_cluster_assignments.txt")), quote=F, row.names=F, sep="\t")
  write.table(timing, file=file.path(outdir, paste0(samplename, "_mutation_timing.txt")), quote=F, row.names=F, sep="\t")

  png(file.path(paste0(samplename, "_qq.png")), height=500, width=1800)
  if (!is.null(vcf_sv) && !all(is.na(MCN_sv$D$pMutCNTail))) {
  	grid.arrange(qqnorm(qnorm(MCN$D$pMutCNTail[!is.na(MCN$D$pMutCNTail)])), qqnorm(qnorm(MCN_indel$D$pMutCNTail[!is.na(MCN_indel$D$pMutCNTail)])), qqnorm(qnorm(MCN_sv$D$pMutCNTail[!is.na(MCN_sv$D$pMutCNTail)])), nrow=1, top=samplename)
  } else {
	grid.arrange(qqnorm(qnorm(MCN$D$pMutCNTail[!is.na(MCN$D$pMutCNTail)])), qqnorm(qnorm(MCN_indel$D$pMutCNTail[!is.na(MCN_indel$D$pMutCNTail)])), make_dummy_figure(), nrow=1, top=samplename)
  }
  dev.off()
}
