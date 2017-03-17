args = commandArgs(T)

samplename = args[1]
outdir = args[2]
snv_vcf_file = args[3]
indel_vcf_file = args[4]
sv_vcf_file = args[5]
bb_file = args[6]
clust_file = args[7]
purity_file = args[8]
summary_table = args[9]
svclone_file = args[10]

merge_clusters = F
filter_small_clusters = F

vcf_template = "~/repo/moritz_mut_assignment/template_icgc_consensus.vcf"

# samplename = "55e520f4-0e4b-41a2-9951-c4e9f323100b"
# cons_method = "wm"
# clust_postfix = ".txt"
# outdir = paste0("output_", cons_method, "/")
# snv_vcf_file = paste0("../..//processed_data/consensusCalls/consensusCalls_2016_10_12/filtered/", samplename, ".consensus.20160830.somatic.snv_mnv.vcf.gz")
# indel_vcf_file = paste0("../..//processed_data/consensusCalls/consensusCalls_2016_10_12/indel/", samplename, ".consensus.20161006.somatic.indel.vcf.gz")
# sv_vcf_file = paste0("../../processed_data/consensusSVs/pcawg_consensus_1.5.160912/", samplename, ".pcawg_consensus_1.5.160912.somatic.sv.vcf.gz")
# bb_file = paste0("../../processed_data/copynumberConsensus/final_consensus_20170119/working_group_output_full_profile/", samplename, "_segments.txt.gz")
# clust_file = paste0("consensus_clusters_", cons_method, "/", samplename, "_subclonal_structure", clust_postfix)
# purity_file = "consensus.20170119.purity.ploidy.txt.gz"
# summary_table = "summary_table_2_20170303.txt"
# svclone_file = paste0("../../processed_data/sv_vafs_geoff/vafs/sv_vafs/", samplename, "_filtered_svs.tsv")
# 
# if (!file.exists(sv_vcf_file)) {
#   sv_vcf_file = "NA"
# }


# mult_file = args[9]
# dpc_assign_file = args[10]
# merge_clusters = F


#vcf_file = "final/final_consensus_12oct_passonly/snv_mnv/0040b1b6-b07a-4b6e-90ef-133523eaf412.consensus.20160830.somatic.snv_mnv.vcf.gz"
#bb_file = "dp/20161213_vanloo_wedge_consSNV_prelimConsCNAallStar/4_copynumber/0040b1b6-b07a-4b6e-90ef-133523eaf412_segments.txt.gz"
#clust_file = "dp/20161213_vanloo_wedge_consSNV_prelimConsCNAallStar/2_subclones/0040b1b6-b07a-4b6e-90ef-133523eaf412_subclonal_structure.txt.gz"
#purity_file = "dp/20161213_vanloo_wedge_consSNV_prelimConsCNAallStar/1_purity_ploidy/purity_ploidy.txt"

source("~/repo/moritz_mut_assignment/MutationTime.R")
source("~/repo/moritz_mut_assignment/util.R")
source("~/repo/dpclust3p/R/interconvertMutationBurdens.R")
library(ggplot2)
library(gridExtra)
library(grid)

########################################################################
# Parse the input
########################################################################
bb <- loadBB(bb_file)
clusters = read.table(clust_file, header=TRUE, sep="\t")

vcf_snv <- readVcf(snv_vcf_file, genome="GRCh37")
vcf_indel <- readVcf(indel_vcf_file, genome="GRCh37")
if (svclone_file=="NA") {
  vcf_sv = NULL
} else {
  vcf_sv <- prepare_svclone_output(svclone_file, vcf_template, genome="GRCh37")
}

purityPloidy <- read.table(purity_file, header=TRUE, sep="\t")
purity = purityPloidy$purity[purityPloidy$samplename==samplename]
ploidy = purityPloidy$ploidy[purityPloidy$samplename==samplename]

summary_table = read.table(summary_table, header=T, stringsAsFactors=F)
sex = summary_table$inferred_sex[summary_table$samplename==samplename]
is_wgd = purityPloidy$wgd_status[purityPloidy$samplename==samplename]=="wgd"

#' If not all clusters are there we need to renumber them
if (max(clusters$cluster) > nrow(clusters)) {
	# dpc_assign_input = dpc_assign
	for (i in 1:nrow(clusters)) {
		clusterid = clusters$cluster[i]
		# dpc_assign$cluster[dpc_assign_input$cluster==clusterid] = i
		clusters$cluster[i] = i
	}
}	

#' Merge clusters if requested
if (merge_clusters) { clusters = mergeClusters(clusters) }

#' Calculate CCF for each cluster
clusters$ccf = clusters$proportion/purity

########################################################################
# Assignments
########################################################################
#' Assign using Moritz' approach
MCN <- computeMutCn(vcf_snv, bb, clusters, purity, gender=sex, isWgd=is_wgd)
MCN_indel <- computeMutCn(vcf_indel, bb, clusters, purity, gender=sex, isWgd=is_wgd)
if (!is.null(vcf_sv)) {
  MCN_sv <- computeMutCn(vcf_sv, bb, clusters, purity, gender=sex, isWgd=is_wgd)
}

########################################################################
# Create the assignment - binom probability
########################################################################
snv_binom = assign_binom_ll(MCN, clusters, purity)
indel_binom = assign_binom_ll(MCN_indel, clusters, purity)
if (!is.null(vcf_sv)) {
  sv_binom = assign_binom_ll(MCN_sv, clusters, purity)
}

########################################################################
# Create the assignment - Kaixians approach
########################################################################
snv_moritz = assign_moritz(MCN, clusters, purity)
indel_moritz = assign_moritz(MCN_indel, clusters, purity)
if (!is.null(vcf_sv)) {
  sv_moritz = assign_moritz(MCN_sv, clusters, purity)
}

########################################################################
# Summary table entry
########################################################################
if (!is.null(vcf_sv)) {
  sv_assignment_table = sv_moritz$plot_data
} else {
  sv_assignment_table = NULL
}

sample_entry = get_summary_table_entry(samplename=samplename, 
                                       summary_table=summary_table, 
                                       cluster_info=snv_moritz$clusters_new, 
                                       snv_assignment_table=snv_moritz$plot_data, 
                                       indel_assignment_table=indel_moritz$plot_data, 
                                       sv_assignment_table=sv_assignment_table,
                                       do_filter=filter_small_clusters)

write.table(sample_entry, file.path(outdir, paste0(samplename, "_summary_table_entry.txt")), row.names=F, sep="\t", quote=F)

########################################################################
# Output to share with PCAWG
########################################################################
snv_timing = data.frame(chromosome=as.character(seqnames(vcf_snv)),
                        position=as.numeric(start(vcf_snv)),
                        mut_type=rep("SNV", nrow(MCN$D)),
                        timing=classifyMutations(MCN$D),
                        chromosome2=rep(NA, nrow(MCN$D)),
                        position2=rep(NA, nrow(MCN$D)),
                        stringsAsFactors=F)

snv_output = data.frame(chromosome=as.character(seqnames(vcf_snv)),
                        position=as.numeric(start(vcf_snv)),
                        mut_type=rep("SNV", nrow(MCN$D)),
                        ccf=snv_moritz$clusters$ccf[match(snv_moritz$plot_data$cluster, snv_moritz$clusters$cluster)],
                        chromosome2=rep(NA, nrow(MCN$D)),
                        position2=rep(NA, nrow(MCN$D)),
                        stringsAsFactors=F)

indel_timing = data.frame(chromosome=as.character(seqnames(vcf_indel)),
                          position=as.numeric(start(vcf_indel)),
                          mut_type=rep("indel", nrow(MCN_indel$D)),
                          timing=classifyMutations(MCN_indel$D),
                          chromosome2=rep(NA, nrow(MCN_indel$D)),
                          position2=rep(NA, nrow(MCN_indel$D)),
                          stringsAsFactors=F)

indel_output = data.frame(chromosome=as.character(seqnames(vcf_indel)),
                          position=as.numeric(start(vcf_indel)),
                          mut_type=rep("indel", nrow(MCN_indel$D)),
                          ccf=indel_moritz$clusters$ccf[match(indel_moritz$plot_data$cluster, indel_moritz$clusters$cluster)],
                          chromosome2=rep(NA, nrow(MCN_indel$D)),
                          position2=rep(NA, nrow(MCN_indel$D)),
                          stringsAsFactors=F)

if (!is.null(vcf_sv)) {
  #' Some magic required to map back to chr/pos
  sv_timing = data.frame(chromosome=info(vcf_sv)$chr1,
                         position=info(vcf_sv)$pos1,
                         mut_type=rep("SV", nrow(MCN_sv$D)),
                         timing=classifyMutations(MCN_sv$D),
                         chromosome2=info(vcf_sv)$chr2,
                         position2=info(vcf_sv)$pos2,
                         stringsAsFactors=F)
  
  sv_output = data.frame(chromosome=info(vcf_sv)$chr1,
                          position=info(vcf_sv)$pos1,
                          mut_type=rep("SV", nrow(MCN_sv$D)),
                          ccf=sv_moritz$clusters$ccf[match(sv_moritz$plot_data$cluster, sv_moritz$clusters$cluster)],
                          chromosome2=info(vcf_sv)$chr2,
                          position2=info(vcf_sv)$pos2,
                         stringsAsFactors=F)
}

if (!is.null(vcf_sv)) {
  timing = do.call(rbind, list(snv_timing, indel_timing, sv_timing))
  ccfs = do.call(rbind, list(snv_output, indel_output, sv_output))
} else {
  timing = rbind(snv_timing, indel_timing)
  ccfs = rbind(snv_output, indel_output)
}

ccfs$ccf = round(ccfs$ccf, 4)

# TODO disabled for now
# write.table(timing, file.path(outdir, paste0(samplename, "_timing_snv_indel_sv.txt")), row.names=F, sep="\t", quote=F)
# write.table(ccfs, file.path(outdir, paste0(samplename, "_ccfs_snv_indel_sv.txt")), row.names=F, sep="\t", quote=F)
final_pcawg11_output = pcawg11_output(snv_moritz, indel_moritz, sv_moritz, MCN, MCN_indel, MCN_sv, vcf_sv)
save(final_pcawg11_output, timing, ccfs, file=file.path(outdir, paste0(samplename, "_pcawg11_output.RData")))

########################################################################
# Plot
########################################################################
#' Make assignment plot for both assignment strategies
p = base_plot(snv_binom$plot_data, "ccf", "Consensus binomial assignment (control)") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf)) + xlab("ccf - snv")
p = p + scale_fill_hue(labels = rev(paste0(" ", 
                                           snv_binom$clusters$cluster, " : ", 
                                           round(snv_binom$clusters$ccf, 2), "  ", 
                                           snv_binom$clusters$n_ssms, "  ")))

p3 = base_plot(snv_moritz$plot_data, "ccf", "Consensus closest cluster assignment") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf)) + xlab("ccf - snv")
p3 = p3 + scale_fill_hue(labels = rev(paste0(" ", 
                                           snv_binom$clusters$cluster, " : ", 
                                           round(snv_moritz$clusters$ccf, 2), "  ", 
                                           snv_moritz$clusters$n_ssms, "  ")))

if (any(indel_moritz$plot_data$ccf < 1.5)) {
  p4 = base_plot(indel_moritz$plot_data, "ccf", "Consensus closest cluster assignment") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf)) + xlab("ccf - indel")
  p4 = p4 + scale_fill_hue(labels = rev(paste0(" ", 
                                             indel_binom$clusters$cluster, " : ", 
                                             round(indel_moritz$clusters$ccf, 2), "  ", 
                                             indel_moritz$clusters$n_ssms, "  ")))
} else {
  p4 = make_dummy_figure()
}

if (!is.null(vcf_sv)) {
  p5 = base_plot(sv_moritz$plot_data, "ccf", "Consensus closest cluster assignment") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf)) + xlab("ccf - sv")
  p5 = p5 + scale_fill_hue(labels = rev(paste0(" ", 
                                               sv_binom$clusters$cluster, " : ", 
                                               round(sv_moritz$clusters$ccf, 2), "  ", 
                                               sv_moritz$clusters$n_ssms, "  ")))
} else {
  p5 = make_dummy_figure()
}

if (!is.null(vcf_sv) && any(sv_moritz$plot_data$ccf < 1.5)) {
  all_data = do.call(rbind, list(snv_binom$plot_data, indel_binom$plot_data, sv_binom$plot_data))
  all_data$type = factor(c(rep("SNV", nrow(snv_binom$plot_data)), rep("indel", nrow(indel_binom$plot_data)), rep("sv", nrow(sv_binom$plot_data))), levels=c("SNV", "indel", "sv"))
} else {
  all_data = do.call(rbind, list(snv_binom$plot_data, indel_binom$plot_data))
  all_data$type = factor(c(rep("SNV", nrow(snv_binom$plot_data)), rep("indel", nrow(indel_binom$plot_data))), levels=c("SNV", "indel", "sv"))
}
p6 = base_plot(all_data, "ccf", "All data", fill="type") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf)) + xlab("ccf")

# my_legend = g_legend(p)

if (samplename %in% summary_table$samplename) {
  power = summary_table$nrpcc[summary_table$samplename==samplename]
  histology = summary_table$histology_abbreviation[summary_table$samplename==samplename]
} else {
  power = NA
  histology = NA
}
title = paste0(samplename, " - ",
               histology, " - ",
               "SNVs ", nrow(snv_binom$plot_data), " - ",
               "Purity ", purity, " - ",
               "Ploidy ", ploidy, " - ",
               "Power ", power, " - ",
               "Clust size input ", paste(rev(clusters$n_ssms), collapse=", "))

png(file.path(outdir, paste0(samplename, "_final_assignment.png")), height=400, width=2000)
grid.arrange(p6, p, p3, p4, p5, nrow=1, top=title)
dev.off()

save.image(file.path(outdir, paste0(samplename, "_assignment.RData")))

# No longer used
# #' DPClust output - sync the data frames
# assign_chr_pos = paste0(dpc_assign$chr, "_", dpc_assign$pos)
# mult_chr_pos = paste0(mult$chr, "_", mult$pos)
# inboth = intersect(mult_chr_pos, assign_chr_pos)
# dpc_assign = dpc_assign[assign_chr_pos %in% inboth,]
# mult = mult[mult_chr_pos %in% inboth,]
# 
# 
# ccf = mult$mutation.copy.number / mult$multiplicity
# plot_data_4 = data.frame(mcn=mult$mutation.copy.number, ccf=ccf, cluster=factor(dpc_assign$cluster, levels=rev(unique(sort(clusters$cluster)))))
# p4 = base_plot(plot_data_4, "ccf", "DPClust assignment") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf))
# 
# 
# png(paste0(samplename, "_final_assignment_all.png"), height=800, width=1000)
# grid.arrange(arrangeGrob(p + theme(legend.position="none"), 
# 			 p3 + theme(legend.position="none"),
# 			 p4 + theme(legend.position="none"),
# 			 p2 + theme(legend.position="none"), ncol=2), 
# 	     arrangeGrob(my_legend), nrow=2, heights=c(18,1), top=title)
# dev.off()

