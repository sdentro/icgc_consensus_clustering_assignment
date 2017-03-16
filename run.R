args = commandArgs(T)

samplename = args[1]
outdir = args[2]
snv_vcf_file = args[3]
indel_vcf_file = args[4]
bb_file = args[5]
clust_file = args[6]
purity_file = args[7]
summary_table = args[8]
mult_file = args[9]
dpc_assign_file = args[10]
merge_clusters = F


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
vcf_snv <- readVcf(snv_vcf_file, genome="GRCh37")
vcf_indel <- readVcf(indel_vcf_file, genome="GRCh37")

# TODO SV
# data_sv = parse_sv_data("test_data/1e27cc8a-5394-4958-9af6-5ece1fe24516_most_likely_copynumbers.txt", "test_data/1e27cc8a-5394-4958-9af6-5ece1fe24516_vaf_ccf.txt")

bb <- loadBB(bb_file)
clusters = read.table(clust_file, header=TRUE, sep="\t")

purityPloidy <- read.table(purity_file, header=TRUE, sep="\t")
purity = purityPloidy$purity[purityPloidy$samplename==samplename]
ploidy = purityPloidy$ploidy[purityPloidy$samplename==samplename]

summary_table = read.table(summary_table, header=T, stringsAsFactors=F)
# Commented out for now
# mult = read.table(mult_file, header=T, stringsAsFactors=F)
# dpc_assign = read.table(dpc_assign_file, header=T, stringsAsFactors=F)
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


#' Convert copy number clonal frequency to CCF
bb$clonal_frequency = bb$clonal_frequency / purity

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

#' TODO: add SV

########################################################################
# Create the assignment - binom probability
########################################################################
snv_binom = assign_binom_ll(MCN, clusters, purity)
# plot_data = res$plot_data
# clusters_new = res$clusters_new

indel_binom = assign_binom_ll(MCN_indel, clusters, purity)

#' TODO: add SV

########################################################################
# Create the assignment - Kaixians approach
########################################################################
snv_moritz = assign_moritz(MCN, clusters, purity)
# plot_data_2 = res$plot_data
# clusters_new_2 = res$clusters_new

indel_moritz = assign_moritz(MCN_indel, clusters, purity)

#' TODO: add SV

########################################################################
# Summary table entry
########################################################################
sample_entry = get_summary_table_entry(samplename=samplename, 
                                       summary_table=summary_table, 
                                       cluster_info=snv_moritz$clusters_new, 
                                       snv_assignment_table=snv_moritz$plot_data, 
                                       indel_assignment_table=indel_moritz$plot_data, 
                                       sv_assignment_table=NULL)

readr::write_tsv(sample_entry, file.path(outdir, paste0(samplename, "_summary_table_entry.txt")))

########################################################################
# Output to share with PCAWG
########################################################################
snv_timing = data.frame(chromosome=as.character(seqnames(vcf_snv)),
                        position=as.numeric(start(vcf_snv)),
                        mut_type=rep("SNV", nrow(MCD$D)),
                        timing=classifyMutations(MCN$D))

snv_output = data.frame(chromosome=as.character(seqnames(vcf_snv)),
                        position=as.numeric(start(vcf_snv)),
                        ccf=snv_moritz$clusters$ccf[match(snv_moritz$plot_data$cluster, snv_moritz$clusters$cluster)])

indel_timing = data.frame(chromosome=as.character(seqnames(vcf_indel)),
                          position=as.numeric(start(vcf_indel)),
                          mut_type=rep("indel", nrow(MCN_indel$D)),
                          timing=classifyMutations(MCN_indel$D))

indel_output = data.frame(chromosome=as.character(seqnames(vcf_indel)),
                          position=as.numeric(start(vcf_indel)),
                          ccf=clusters_new_2$ccf[match(indel_moritz$plot_data$cluster, indel_moritz$clusters$cluster)])

timing = rbind(snv_timing, indel_timing)
ccfs = rbind(snv_output, indel_output)

readr::write_tsv(timing, file.path(outdir, paste0(samplename, "_timing_snv_indel.txt")))
readr::write_tsv(ccfs, file.path(outdir, paste0(samplename, "_ccfs_snv_indel.txt")))

#' TODO: add SV

########################################################################
# Plot
########################################################################
base_plot = function(plot_data, x_variable, title=NA) {
  p = ggplot(plot_data) + aes_string(x=x_variable, y="..count..", fill="cluster") + geom_histogram(binwidth=0.02, colour="black", position="stack") + ylab("Count") + theme(legend.position="bottom") + scale_fill_discrete(drop = FALSE)
  if (!is.na(title)) {
	  p = p + ggtitle(title)
  }
  return(p)
}

g_legend <- function(a.gplot) { 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

#' Make assignment plot for both assignment strategies
p = base_plot(snv_binom$plot_data, "ccf", "Consensus binomial assignment") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf))
# p2 = base_plot(plot_data, "mcn", "Consensus binomial assignment") + xlim(0, 4)
p3 = base_plot(snv_moritz$plot_data, "ccf", "Consensus closest cluster assignment") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf))
p = p + scale_fill_hue(labels = rev(paste0(" ", 
                                           snv_binom$clusters$cluster, " : ", 
                                           round(snv_binom$clusters$ccf, 2), "  ", 
                                           snv_binom$clusters$n_ssms, "  ")))
my_legend = g_legend(p)

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
               "Power ", power)


png(file.path(outdir, paste0(samplename, "_final_assignment.png")), height=400, width=1000)
grid.arrange(arrangeGrob(p + theme(legend.position="none"), 
                         p3 + theme(legend.position="none"), ncol=2), 
             arrangeGrob(my_legend), nrow=2, heights=c(9,1), top=title)
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

