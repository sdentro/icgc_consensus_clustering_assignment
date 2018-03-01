#' Produce PCAWG-wide output files


args = commandArgs(T)
libpath = args[1]
indir = args[2]
outputdir = args[3]
has_indel = as.logical(args[4])
has_sv = as.logical(args[5])
make_plot = as.logical(args[6])
project = args[7]
infile = args[8]

source(file.path(libpath, "MutationTime.R"))
source(file.path(libpath, "util.R"))
source("~/repo/dpclust3p/R/interconvertMutationBurdens.R")
library(ggplot2)
library(gridExtra)

write_output_summary_table = function(structure_df, outdir, samplename, project, purity) {
  #' Determine which cluster is clonal
  is_clonal = structure_df$fraction_cancer_cells > 0.9
  if (any(is_clonal)) {
    if (sum(is_clonal, na.rm=T) > 1) {
      clone_id = which.min(abs(structure_df$fraction_cancer_cells[is_clonal]-1))
      is_clonal = rep(F, length(is_clonal))
      is_clonal[clone_id] = T
    }
    clonal_cluster = structure_df$cluster[is_clonal]
  } else {
    clonal_cluster = NA
  }
  
  if (!is.na(clonal_cluster)) {
    num_clonal = structure_df$n_snvs[structure_df$cluster==clonal_cluster]
    num_subclonal = ifelse(nrow(structure_df) > 1, structure_df$n_snvs[structure_df$cluster!=clonal_cluster], 0)
  } else {
    num_clonal = 0
    num_subclonal = sum(structure_df$n_snvs)
  }
  
  
  #' Write out the summary table entry
  summary_table = data.frame(cancer_type=project, 
                             samplename=samplename, 
                             num_subclones=nrow(structure_df)-1, 
                             purity=purity, 
                             ploidy=NA,
                             num_clonal=num_clonal,
                             num_subclonal=num_subclonal,
                             frac_clonal=num_clonal / (num_clonal+num_subclonal),
                             noCNA=NA,
                             clonal=NA,
                             subclonal=NA)
  write.table(summary_table, file=file.path(outdir, paste0(samplename, "_summary_table.txt")), row.names=F, sep="\t", quote=F)
}

#for (infile in list.files(indir, pattern="_assignment.RData", full.names=T)) {
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

  write.table(subcl_struct, file=file.path(outputdir, paste0(samplename, "_subclonal_structure.txt")), quote=F, row.names=F, sep="\t")
  write.table(assign_probs, file=file.path(outputdir, paste0(samplename, "_cluster_assignments.txt")), quote=F, row.names=F, sep="\t")
  write.table(timing, file=file.path(outputdir, paste0(samplename, "_mutation_timing.txt")), quote=F, row.names=F, sep="\t")
  
  
   write_output_summary_table(subcl_struct, outputdir, samplename, project, purity) 
  
  
  if (make_plot) {
  png(file.path(paste0(samplename, "_qq.png")), height=500, width=1800)
  #if (!is.null(vcf_sv) && !all(is.na(MCN_sv$D$pMutCNTail))) {
  if (has_indel & has_sv) {
  	grid.arrange(qqnorm(qnorm(MCN$D$pMutCNTail[!is.na(MCN$D$pMutCNTail)])), qqnorm(qnorm(MCN_indel$D$pMutCNTail[!is.na(MCN_indel$D$pMutCNTail)])), qqnorm(qnorm(MCN_sv$D$pMutCNTail[!is.na(MCN_sv$D$pMutCNTail)])), nrow=1, top=samplename)
  } else if (!has_indel & !has_sv) {
	  grid.arrange(qqnorm(qnorm(MCN$D$pMutCNTail[!is.na(MCN$D$pMutCNTail)])), make_dummy_figure(), make_dummy_figure(), nrow=1, top=samplename)
  } else {
	grid.arrange(qqnorm(qnorm(MCN$D$pMutCNTail[!is.na(MCN$D$pMutCNTail)])), qqnorm(qnorm(MCN_indel$D$pMutCNTail[!is.na(MCN_indel$D$pMutCNTail)])), make_dummy_figure(), nrow=1, top=samplename)
  }
  dev.off()
  }
#}
