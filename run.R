args = commandArgs(T)

samplename = args[1]
vcf_file = args[2]
bb_file = args[3]
clust_file = args[4]
purity_file = args[5]
summary_table = args[6]
mult_file = args[7]
dpc_assign_file = args[8]
merge_clusters = F


#vcf_file = "final/final_consensus_12oct_passonly/snv_mnv/0040b1b6-b07a-4b6e-90ef-133523eaf412.consensus.20160830.somatic.snv_mnv.vcf.gz"
#bb_file = "dp/20161213_vanloo_wedge_consSNV_prelimConsCNAallStar/4_copynumber/0040b1b6-b07a-4b6e-90ef-133523eaf412_segments.txt.gz"
#clust_file = "dp/20161213_vanloo_wedge_consSNV_prelimConsCNAallStar/2_subclones/0040b1b6-b07a-4b6e-90ef-133523eaf412_subclonal_structure.txt.gz"
#purity_file = "dp/20161213_vanloo_wedge_consSNV_prelimConsCNAallStar/1_purity_ploidy/purity_ploidy.txt"

source("~/repo/moritz_mut_assignment/MutationTime.R")
source("~/repo/dpclust3p/R/interconvertMutationBurdens.R")
library(ggplot2)
library(gridExtra)
library(grid)

vcf <- readVcf(vcf_file, genome="GRCh37")
bb <- loadBB(bb_file)
clusters = read.table(clust_file, header=TRUE, sep="\t")
purityPloidy <- read.table(purity_file, header=TRUE, sep="\t")
purity = purityPloidy$purity[purityPloidy$samplename==samplename]
ploidy = purityPloidy$ploidy[purityPloidy$samplename==samplename]
summary_table = read.table(summary_table, header=T, stringsAsFactors=F)
mult = read.table(mult_file, header=T, stringsAsFactors=F)
dpc_assign = read.table(dpc_assign_file, header=T, stringsAsFactors=F)
sex = summary_table$inferred_sex[summary_table$samplename==samplename]
is_wgd = purityPloidy$wgd_status[purityPloidy$samplename==samplename]=="wgd"

#' If not all clusters are tehre we need to renumber them
if (max(clusters$cluster) > nrow(clusters)) {
	dpc_assign_input = dpc_assign
	for (i in 1:nrow(clusters)) {
		clusterid = clusters$cluster[i]
		dpc_assign$cluster[dpc_assign_input$cluster==clusterid] = i
		clusters$cluster[i] = i
	}
}	


#' Convert copy number clonal frequency to CCF
bb$clonal_frequency = bb$clonal_frequency / purity

if (merge_clusters) {
	  clusters = mergeClusters(clusters)
}
clusters$ccf = clusters$proportion/purity

MCN <- computeMutCn(vcf, bb, clusters, purity, gender=sex, isWgd=is_wgd)

#' Calculate bionomial probability of each mutation belonging in each cluster
getClustLL = function(data, cluster_locations, purity) {
	assignment_ll = sapply(1:length(cluster_locations), function(c) {
			        mutBurdens = mutationCopyNumberToMutationBurden(cluster_locations[c] * data$MutCN, data$MajCN+data$MinCN, purity, rep(2, nrow(data)))
				data$altCount*log(mutBurdens) + data$wtCount*log(1-mutBurdens)
				})
	return(assignment_ll)
}

clust_assign_ll = getClustLL(MCN$D, clusters$proportion/purity, purity)
best_cluster = unlist(apply(clust_assign_ll, 1, function(x) if (all(is.na(x))) NA else which.max(x)))

cluster_counts = table(best_cluster)
clusters_new = data.frame(clusters$cluster, sapply(clusters$cluster, function(x) cluster_counts[[as.character(x)]]), clusters$proportion, clusters$ccf)
colnames(clusters_new) = colnames(clusters)

#' Calculate MCN and CCF for plotting
mcn = mutationBurdenToMutationCopyNumber(burden=MCN$D$altCount / (MCN$D$altCount + MCN$D$wtCount), cellularity=purity, normalCopyNumber=rep(2, nrow(MCN$D)), totalCopyNumber=MCN$D$MajCN + MCN$D$MinCN)
ccf = mcn / MCN$D$MutCN

base_plot = function(plot_data, x_variable, title=NA) {
  p = ggplot(plot_data) + aes_string(x=x_variable, y="..count..", fill="cluster") + geom_histogram(binwidth=0.02, colour="black", position="stack") + ylab("Count") + theme(legend.position="bottom") + scale_fill_discrete(drop = FALSE)
  if (!is.na(title)) {
	  p = p + ggtitle(title)
  }
  return(p)
}

plot_data = data.frame(mcn=mcn, ccf=ccf, cluster=factor(best_cluster, levels=rev(unique(sort(clusters_new$cluster)))))
p = base_plot(plot_data, "ccf", "Consensus binomial assignment") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf))
p2 = base_plot(plot_data, "mcn", "Consensus binomial assignment") + xlim(0, 4)

p = p + scale_fill_hue(labels = rev(paste0(" ", clusters_new$cluster, " : ", round(clusters_new$ccf, 2), "  ", clusters_new$n_ssms, "  ")))

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

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
               "SNVs ", nrow(plot_data), " - ",
               "Purity ", purity, " - ",
               "Ploidy ", ploidy, " - ",
               "Power ", power)

#' Kaixians approach
best_cluster = sapply(MCN$D$CNF, function(x) if (is.na(x)) NA else which.min(abs(x-clusters_new$proportion)))
cluster_counts = table(best_cluster)
clusters_new_2 = data.frame(clusters$cluster, sapply(clusters$cluster, function(x) cluster_counts[[as.character(x)]]), clusters$proportion, clusters$ccf)
colnames(clusters_new_2) = colnames(clusters)

plot_data_2 = data.frame(mcn=mcn, ccf=MCN$D$CNF, cluster=factor(best_cluster, levels=rev(unique(sort(clusters_new$cluster)))))
plot_data_2$ccf = plot_data_2$ccf / purity
p3 = base_plot(plot_data_2, "ccf", "Consensus closest cluster assignment") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf))

#' DPClust output - sync the data frames
assign_chr_pos = paste0(dpc_assign$chr, "_", dpc_assign$pos)
mult_chr_pos = paste0(mult$chr, "_", mult$pos)
inboth = intersect(mult_chr_pos, assign_chr_pos)
dpc_assign = dpc_assign[assign_chr_pos %in% inboth,]
mult = mult[mult_chr_pos %in% inboth,]


ccf = mult$mutation.copy.number / mult$multiplicity
plot_data_4 = data.frame(mcn=mult$mutation.copy.number, ccf=ccf, cluster=factor(dpc_assign$cluster, levels=rev(unique(sort(clusters$cluster)))))
p4 = base_plot(plot_data_4, "ccf", "DPClust assignment") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf))

png(paste0(samplename, "_final_assignment.png"), height=400, width=1000)
grid.arrange(arrangeGrob(p + theme(legend.position="none"), 
			 p3 + theme(legend.position="none"), ncol=2), 
	     arrangeGrob(my_legend), nrow=2, heights=c(9,1), top=title)
dev.off()

png(paste0(samplename, "_final_assignment_all.png"), height=800, width=1000)
grid.arrange(arrangeGrob(p + theme(legend.position="none"), 
			 p3 + theme(legend.position="none"),
			 p4 + theme(legend.position="none"),
			 p2 + theme(legend.position="none"), ncol=2), 
	     arrangeGrob(my_legend), nrow=2, heights=c(18,1), top=title)
dev.off()

save.image(paste0(samplename, "_assignment.RData"))
#' # Check output
#' # Mutation Annotation
#' head(MCN$D)
# Classify as basic clonal states
# table(classifyMutations(MCN$D))
# Timing parameters
# MCN$P[[1]]
# Extract timing of segments
# bb$timing_param <- MCN$P
# bbToTime(bb)
