#
# This is the main function to run the ICGC PCAWG-11 consensus assignment procedure. It takes
# VCF files with mutations (at least SNVs required, others can be NA) copy number, cluster locations/sizes,
# a purity value, a PCAWG-11 summary table (with sex and WGD denomiation per patient at least)
# and two SVclone files with SV mappings (can be NA as well)
#
# Out of the box it cannot handle subclonal copy number without introducing additional cluster locations
# so there are switches: To round subclonal copy number or to remove it.
#

########################################################################
# Command line options
########################################################################
library(optparse)
option_list = list(
  make_option(c("-l", "--libpath"), type="character", default=NULL, help="Path to pipeline installation directory", metavar="character"),
  make_option(c("--sam"), type="character", default=NULL, help="Samplename", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default=NULL, help="Directory where output will be written", metavar="character"),
  make_option(c("--snv"), type="character", default=NULL, help="SNV VCF file", metavar="character"),
  make_option(c("--ind"), type="character", default=NULL, help="Indel VCF file", metavar="character"),
  make_option(c("--sv"), type="character", default=NULL, help="SV VCF file", metavar="character"),
  make_option(c("--cna"), type="character", default=NULL, help="CNA segments file", metavar="character"),
  make_option(c("--struct"), type="character", default=NULL, help="Subclonal structure file", metavar="character"),
  make_option(c("--pur"), type="numeric", default=NULL, help="Sample purity", metavar="character"),
  make_option(c("--ploi"), type="numeric", default=NULL, help="Sample ploidy", metavar="character"),
  make_option(c("--iswgd"), type="logical", action="store_true", default=FALSE, help="Provide when the sample has undergone a whole genome doubling", metavar="character"),
  make_option(c("--sv_vaf"), type="character", default=NULL, help="SVclone VAF file", metavar="character"),
  make_option(c("--round_subclonal_cna"), default=FALSE, type="logical", help="Round subclonal CNAs", metavar="logical", action="store_true"),
  make_option(c("--remove_subclonal_cna"), default=FALSE, type="logical", help="Remove subclonal CNAs", metavar="logical", action="store_true"),
  make_option(c("--sex"), type="character", default=NULL, help="Sex of the donor", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

libpath = opt$libpath
samplename = opt$sam
outdir = opt$outputdir
snv_vcf_file = opt$snv
bb_file = opt$cna
clust_file = opt$struct
purity = opt$pur
ploidy = opt$ploi
is_wgd = opt$iswgd
round_subclonal_cna = opt$round_subclonal_cna
remove_subclones = opt$remove_subclonal_cna
sex = opt$sex

# Possible not supplied arguments
indel_vcf_file = NULL
if (!is.null(opt$ind)) {
	if (opt$ind!="NA") {
		indel_vcf_file = opt$ind
	}
}

sv_vcf_file = NULL
if (!is.null(opt$sv)) {
	if (opt$sv!="NA") {
		sv_vcf_file = opt$sv
	}
}

svclone_file = NULL
if (!is.null(opt$sv_vaf)) {
	if (opt$sv_vaf!="NA") {
		svclone_file = opt$sv_vaf
	}
}

# check everything is there
if (is.null(libpath) | is.null(samplename) | is.null(outdir) | is.null(snv_vcf_file) | is.null(bb_file) | is.null(clust_file) | is.null(purity) | is.null(ploidy) | is.null(sex)) {
  print_help(opt_parser)
  stop("Missing required parameters\n", call.=FALSE)
}

if (round_subclonal_cna & remove_subclones) {
  stop("round_subclonal_cna and remove_subclonal_cna cannot both be set to TRUE, please chose one\n", call.=FALSE)
}


# parameters
merge_clusters = T
filter_small_clusters = F # no longer used
deltaFreq <- 0.00 # merge clusters withing deltaFreq
min_read_diff = 2 # merge clusters within this number of mutant reads
xmin = 0 # overdispersion parameter for beta-binomial distribution
min_ccf_clone = 0.95 # minimum CCF of the clonal cluster
max_allowed_sv_prob_diff = 0.5 # if the pair of breakpoints that underpin a SV differ by more than this, we mask the SV as its VAFs are deemed unreliable

vcf_template = file.path(libpath, "template_icgc_consensus.vcf")

library(MutationTimeR)
source(file.path(libpath, "util.R")) # this must occur after loading mtimer as it overloads the loadBB function
library(ggplot2)
library(gridExtra)
library(grid)

# Overdispersion parameter
rho_snv = 0.01
rho_indel = 0.01
rho_sv = 0.05
q = 0.05

########################################################################
# Parse the input
########################################################################
bb = loadBB(bb_file, round_subclones=round_subclonal_cna, remove_subclones=remove_subclones)
clusters = read.table(clust_file, header=TRUE, sep="\t")

if (!any(colnames(clusters)=="proportion")) {
	if (!any(colnames(clusters)=="ccf")) { print("Structure file requires at least proportion or ccf column"); q(save="no") }
	clusters$proportion = clusters$ccf * purity
}

if (!any(colnames(clusters)=="ccf")) {
  clusters$ccf = clusters$proportion / purity
}

# not al pipeliones adhere to the prescribed standard
if (any(colnames(clusters)=="n_ssm")) {
  colnames(clusters)[colnames(clusters)=="n_ssm"] = "n_ssms"
}

# bring clusters in the expected order of columns
clusters = clusters[,c("cluster", "proportion", "ccf", "n_ssms")]

# sort the clusters and renumber
clusters = clusters[with(clusters, order(proportion, decreasing=T)),]
clusters$cluster = 1:nrow(clusters)

# Read in VCF files
vcf_snv = readVcf(snv_vcf_file, genome="GRCh37")
if (is.null(indel_vcf_file)) {
	vcf_indel = NULL
} else {
	vcf_indel = readVcf(indel_vcf_file, genome="GRCh37")
	if (nrow(vcf_indel)==0) {
		vcf_indel = NULL
	}
}
if (is.null(svclone_file) | is.null(svclone_file)) {
  vcf_sv = NULL
  vcf_sv_alt = NULL
} else {
  vcf_sv = prepare_svclone_output(svclone_file, vcf_template, genome="GRCh37", sv_vcf_file=sv_vcf_file)
  vcf_sv_alt = prepare_svclone_output(svclone_file, vcf_template, genome="GRCh37", sv_vcf_file=sv_vcf_file, take_preferred_breakpoint=F)
  if (nrow(vcf_sv)==0) {
	  vcf_sv = NULL
	  vcf_sv_alt = NULL
  }
}

# Annotate copy number cell frequencies
if ("ccf" %in% colnames(elementMetadata(bb))) {
	bb$clonal_frequency = purity*bb$ccf
} else {
	bb$clonal_frequency = purity
}

# If not all clusters are there we need to renumber them
if (max(clusters$cluster) > nrow(clusters)) {
	for (i in 1:nrow(clusters)) {
		clusterid = clusters$cluster[i]
		clusters$cluster[i] = i
	}
}	

# Alt merge clusters
if (merge_clusters & nrow(clusters) > 1) { clusters = mergeClustersByMutreadDiff(clusters, purity, ploidy, vcf_snv, min_read_diff) }

########################################################################
# Assignments - MutationTimeR
########################################################################
	MCN <- computeMutCn(vcf_snv, bb, clusters, purity, gender=sex, isWgd=is_wgd, rho=rho_snv, n.boot=0, xmin=xmin, deltaFreq=deltaFreq)

	if (!is.null(vcf_indel)) {
	  temp = clusters
	  temp$n_ssms = estimate_cluster_size(clusters$ccf, vcf_indel, bb, purity, sex, is_wgd, rho_snv, xmin, deltaFreq)
		MCN_indel <- computeMutCn(vcf_indel, bb, temp, purity, gender=sex, isWgd=is_wgd, rho=rho_indel, n.boot=0, xmin=xmin, deltaFreq=deltaFreq)
	}

if (!is.null(vcf_sv)) {
  # SVclone has mapped SVs to particular copy number segments. This may not directly correspond to the
  # segment at the exact base of the SV breakpoint. Here we therefore create a custom copy number profile
  # where the copy number corresponds at the exact base of the SV.
  temp_bb = copynumber_at_sv_locations(bb, vcf_sv)
  temp_clusters = clusters
  temp_clusters$n_ssms = estimate_cluster_size(clusters$ccf, vcf_sv, temp_bb, purity, sex, is_wgd, rho_sv, xmin, deltaFreq)
  MCN_sv <- computeMutCn(vcf_sv, temp_bb, temp_clusters, purity, gender=sex, isWgd=is_wgd, rho=rho_sv, n.boot=0, xmin=xmin, deltaFreq=deltaFreq)
  
  # now do the same for the other SV allele, the non-preferred one by SVclone
  temp_bb = copynumber_at_sv_locations(bb, vcf_sv_alt)
  #temp_clusters_alt = clusters
  #temp_clusters_alt$n_ssms = estimate_cluster_size(clusters$ccf, vcf_sv_alt, temp_bb, purity, sex, is_wgd, rho_sv, xmin, deltaFreq)
  MCN_sv_alt <- computeMutCn(vcf_sv_alt, temp_bb, temp_clusters, purity, gender=sex, isWgd=is_wgd, rho=rho_sv, n.boot=0, xmin=xmin, deltaFreq=deltaFreq)
}

snv_mtimer = assign_mtimer(MCN, clusters, purity)
if (!is.null(vcf_indel)) {
  indel_mtimer = assign_mtimer(MCN_indel, clusters, purity)
} else {
  indel_mtimer = NULL
}

if (!is.null(vcf_sv)) {
  sv_mtimer = assign_mtimer(MCN_sv, clusters, purity)
}

if (!is.null(vcf_sv_alt)) {
  sv_alt_mtimer = assign_mtimer(MCN_sv_alt, clusters, purity)
}

########################################################################
# Assignments - binom probability (used as comparison for plotting)
########################################################################
snv_binom = assign_binom_ll(MCN, clusters, purity)
if (!is.null(vcf_indel)) {
  indel_binom = assign_binom_ll(MCN_indel, clusters, purity)
}
if (!is.null(vcf_sv)) {
  sv_binom = assign_binom_ll(MCN_sv, clusters, purity)
}
if (!is.null(vcf_sv_alt)) {
  sv_alt_binom = assign_binom_ll(MCN_sv_alt, clusters, purity)
}

########################################################################
# Obtain final output
########################################################################
final_pcawg11_output = pcawg11_output(snv_mtimer, indel_mtimer, sv_mtimer, MCN, MCN_indel, MCN_sv, vcf_snv, vcf_indel, vcf_sv, MCN_sv_alt, vcf_sv_alt, sv_alt_mtimer) #sv_vcf_file

########################################################################
# Format the output - including merging of SV breakpoint probs
########################################################################
snv_timing = data.frame(chromosome=as.character(seqnames(vcf_snv)),
                        position=as.numeric(start(vcf_snv)),
                        mut_type=rep("SNV", nrow(MCN$D)),
                        timing=MutationTimeR:::classifyMutations(MCN$D),
                        chromosome2=rep(NA, nrow(MCN$D)),
                        position2=rep(NA, nrow(MCN$D)),
                        svid=rep(NA, nrow(MCN$D)),
                        prob_clonal_early=MCN$D$pGain,
                        prob_clonal_late=MCN$D$pSingle,
                        prob_subclonal=MCN$D$pSub,
                        stringsAsFactors=F)

snv_output = data.frame(chromosome=final_pcawg11_output$snv_assignments_prob$chr,
                        position=final_pcawg11_output$snv_assignments_prob$pos,
                        mut_type=rep("SNV", nrow(MCN$D)),
                        final_pcawg11_output$snv_assignments_prob[, grepl("cluster", colnames(final_pcawg11_output$snv_assignments_prob)), drop=F],
                        chromosome2=rep(NA, nrow(MCN$D)),
                        position2=rep(NA, nrow(MCN$D)),
                        svid=rep(NA, nrow(MCN$D)),
                        stringsAsFactors=F)

if (!is.null(vcf_indel)) {
  indel_timing = data.frame(chromosome=as.character(seqnames(vcf_indel)),
                            position=as.numeric(start(vcf_indel)),
                            mut_type=rep("indel", nrow(MCN_indel$D)),
                            timing=MutationTimeR:::classifyMutations(MCN_indel$D),
                            chromosome2=rep(NA, nrow(MCN_indel$D)),
                            position2=rep(NA, nrow(MCN_indel$D)),
                            svid=rep(NA, nrow(MCN_indel$D)),
                            prob_clonal_early=MCN_indel$D$pGain,
                            prob_clonal_late=MCN_indel$D$pSingle,
                            prob_subclonal=MCN_indel$D$pSub,
                            stringsAsFactors=F)
  
  indel_output = data.frame(chromosome=final_pcawg11_output$indel_assignments_prob$chr,
                            position=final_pcawg11_output$indel_assignments_prob$pos,
                            mut_type=rep("indel", nrow(MCN_indel$D)),
                            final_pcawg11_output$indel_assignments_prob[, grepl("cluster", colnames(final_pcawg11_output$indel_assignments_prob)), drop=F],
                            chromosome2=rep(NA, nrow(MCN_indel$D)),
                            position2=rep(NA, nrow(MCN_indel$D)),
                            svid=rep(NA, nrow(MCN_indel$D)),
                            stringsAsFactors=F)
} else {
  indel_timing = NULL
  indel_output = NULL
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
                         prob_clonal_early=MCN_sv$D$pGain,
                         prob_clonal_late=MCN_sv$D$pSingle,
                         prob_subclonal=MCN_sv$D$pSub,
                         stringsAsFactors=F)
  # reset probabilities if one cannot distinguish between clonal early and late
  sv_timing$prob_clonal_early[as.character(sv_timing$timing)=="clonal [NA]"] = NA
  sv_timing$prob_clonal_late[as.character(sv_timing$timing)=="clonal [NA]"] = NA

  if (!is.null(vcf_sv_alt)) {
    # here we use the other SV breakpoint for a separate assignment
    sv_alt_timing = data.frame(chromosome=as.character(seqnames(vcf_sv_alt)),
                           position=start(vcf_sv_alt),
                           mut_type=rep("SV", nrow(MCN_sv_alt$D)),
                           timing=MutationTimeR:::classifyMutations(MCN_sv_alt$D),
                           chromosome2=info(vcf_sv_alt)$chr2,
                           position2=info(vcf_sv_alt)$pos2,
                           svid=info(vcf_sv_alt)$id,
                           prob_clonal_early=MCN_sv_alt$D$pGain,
                           prob_clonal_late=MCN_sv_alt$D$pSingle,
                           prob_subclonal=MCN_sv_alt$D$pSub,
                           stringsAsFactors=F)
    # reset probabilities if one cannot distinguish between clonal early and late
    sv_alt_timing$prob_clonal_early[as.character(sv_alt_timing$timing)=="clonal [NA]"] = NA
    sv_alt_timing$prob_clonal_late[as.character(sv_alt_timing$timing)=="clonal [NA]"] = NA
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
    # assignment probabilties for some sv breakpoint pairs differ a lot. here those are masked, because we don't really know where these belong
    for (svid in unlist(lapply(sv_output$svid, function(x) unlist(strsplit(x, "_"))[1]))) {
      probs_bp_a = sv_output[sv_output$svid==paste0(svid, "_1"), grepl("cluster_", colnames(sv_output))]
      probs_bp_b = sv_output[sv_output$svid==paste0(svid, "_2"), grepl("cluster_", colnames(sv_output))]
      prob_diff = mean(as.numeric(abs(probs_bp_a-probs_bp_b)))
      
      if (prob_diff > max_allowed_sv_prob_diff | is.na(prob_diff)) {
        # mask probabilities of breakpoint pairs that differ by more than the established threshold
        sv_output[grepl(svid, sv_output$svid), grepl("cluster_", colnames(sv_output))] = NA
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

  timing = do.call(rbind, list(snv_timing, indel_timing, sv_timing))
  assign_probs = do.call(rbind, list(snv_output, indel_output, sv_output))
  
} else {
  sv_output = NULL
  timing = rbind(snv_timing, indel_timing)
  assign_probs = rbind(snv_output, indel_output)
}

# set timing probs to NA for clonal mutations that could not be classified as early or late
# if this isn't done then prob_clonal_late is set to 1-prob_subclonal, which would be misleading
timing[!is.na(as.character(timing$timing)) & as.character(timing$timing)=="clonal [NA]", c("prob_clonal_early", "prob_clonal_late")] = NA

########################################################################
# Checking
########################################################################
# check for negative probabilities
# These can occur naturally due to rounding errors. here we check whether they exist
# if a probability is within twice the machine precision we round it of
# across pcawg, this occurs for only a single mutation
all_probs = assign_probs[,grepl("cluster_", colnames(assign_probs)), drop=F]
if (any(assign_probs[,grepl("cluster_", colnames(assign_probs))] < 0, na.rm=T)) {
  
  if (all(all_probs[!is.na(all_probs) & all_probs < 0] > -2*.Machine$double.eps)) {
    print("Encountered negative probabilities due to rounding")
    all_probs[!is.na(all_probs) & all_probs < 0] = round(all_probs[!is.na(all_probs) & all_probs < 0])
    assign_probs[,grepl("cluster_", colnames(assign_probs))] = all_probs
  } else {
    stop("Encountered major negative probabilities")
  }
}

all_timing_probs = timing[,grepl("prob_", colnames(timing)), drop=F]
if (any(all_timing_probs < 0, na.rm=T)) {
  if (all(all_timing_probs[!is.na(all_timing_probs) & all_timing_probs < 0] > -2*.Machine$double.eps)) {
    print("Encountered negative timing probabilities due to rounding")
    all_timing_probs[!is.na(all_timing_probs) & all_timing_probs < 0] = round(all_timing_probs[!is.na(all_timing_probs) & all_timing_probs < 0])
    timing[,grepl("cluster_", colnames(timing))] = all_timing_probs
  } else {
    stop("Encountered major negative timing probabilities")
  }
}

# check that all mutations have been accounted for
if (length(vcf_snv) != sum(assign_probs$mut_type=="SNV")) { print("Did not assign all SNVs") }
if (length(vcf_indel) != sum(assign_probs$mut_type=="indel")) { print("Did not assign all indels") }
if ((2*length(vcf_sv)) != sum(assign_probs$mut_type=="SV")) { print("Did not assign all SVs") }

assign_chrpos = paste(assign_probs$chromosome, "_", assign_probs$position, sep="")
timing_chrpos = paste(timing$chromosome, "_", timing$position, sep="")

orig_vcf_sv = readVcf(sv_vcf_file, "hg19")
orig_chrpos = paste(as.character(seqnames(orig_vcf_sv)), "_", start(orig_vcf_sv), sep="")
if (any(assign_chrpos!=orig_chrpos) | any(timing_chrpos!=orig_chrpos)) { 
  print("SV position discrepancy detected") 
  
  print(which(!assign_chrpos==orig_chrpos))
  print(which(!timing_chrpos==orig_chrpos))
  
  
}

########################################################################
# Summary table entry
########################################################################

# summarise tail probabilities for the beta-binomials fit to each mutation type
qq_snv <- mean(MCN$D$pMutCNTail < q/2 | MCN$D$pMutCNTail > 1-q/2, na.rm=T)
if (!is.null(vcf_indel)) {
  qq_indel <- mean(MCN_indel$D$pMutCNTail < q/2 | MCN_indel$D$pMutCNTail > 1-q/2, na.rm=T)
} else {
  qq_indel <- NA
}
if (!is.null(vcf_sv)) {
  qq_sv <- mean(MCN_sv$D$pMutCNTail < q/2 | MCN_sv$D$pMutCNTail > 1-q/2, na.rm=T)
} else {
  qq_sv = NA
}

p_snv = pbinom(sum(MCN$D$pMutCNTail < q/2 | MCN$D$pMutCNTail > 1-q/2, na.rm=T), nrow(MCN$D), 0.05, lower.tail=TRUE)
if (!is.null(vcf_indel)) {
  p_indel = pbinom(sum(MCN_indel$D$pMutCNTail < q/2 | MCN_indel$D$pMutCNTail > 1-q/2, na.rm=T), nrow(MCN_indel$D), 0.05, lower.tail=TRUE)
} else {
  p_indel = NA
}
if (!is.null(vcf_sv)) {
  p_sv = pbinom(sum(MCN_sv$D$pMutCNTail < q/2 | MCN_sv$D$pMutCNTail > 1-q/2, na.rm=T), nrow(MCN_sv$D), 0.05, lower.tail=TRUE)
} else {
  p_sv = NA
}
posthoc_stats = data.frame(samplename, qq_snv=qq_snv, qq_indel=qq_indel, qq_sv=qq_sv, p_snv=p_snv, p_indel=p_indel, p_sv=p_sv)

if (!is.null(vcf_sv)) {
  sv_assignment_table = sv_mtimer$plot_data
} else {
  sv_assignment_table = NULL
}
sample_entry = get_summary_table_entry(samplename=samplename, 
                                       purity=purity, 
                                       ploidy=ploidy, 
                                       sex=sex, 
                                       is_wgd=is_wgd,
                                       cluster_info=snv_mtimer$clusters_new, 
                                       # snv_assignment_table=snv_mtimer$plot_data,
                                       snv_assignment_table=snv_output, 
                                       min_clonal_ccf=min_ccf_clone,
                                       # indel_assignment_table=indel_mtimer$plot_data, 
                                       indel_assignment_table=indel_output,
                                       # sv_assignment_table=sv_assignment_table,
                                       sv_assignment_table=sv_output,
                                       do_filter=filter_small_clusters)
sample_entry = data.frame(sample_entry, posthoc_stats, stringsAsFactors=F)

write.table(sample_entry, file.path(outdir, paste0(samplename, "_summary_table_entry.txt")), row.names=F, sep="\t", quote=F)

########################################################################
# produce pcawg wide output files
########################################################################
subcl_struct = final_pcawg11_output$final_clusters[,c("cluster", "proportion", "ccf", "n_snvs", "n_indels", "n_svs")]
colnames(subcl_struct)[2] = "fraction_total_cells"
colnames(subcl_struct)[3] = "fraction_cancer_cells"

# check assignment probabilities and timing are consistently available
for (i in 1:nrow(all_probs)) {
  # if any prob is NA, then all should be NA
  if ((any(is.na(all_probs[i,])) & !all(is.na(all_probs[i,])))) { print(paste0("PROBS INCONSISTENT:", i)) }
  # if the assignment prob is NA, then timing information should be NA aswel
  if (all(is.na(all_probs[i,])) & !is.na(timing$timing[i])) { print(paste0("TIMING INCONSISTENT:", i)) }
}

write.table(subcl_struct, file=file.path(outdir, paste0(samplename, "_subclonal_structure.txt")), quote=F, row.names=F, sep="\t")
write.table(assign_probs, file=file.path(outdir, paste0(samplename, "_cluster_assignments.txt")), quote=F, row.names=F, sep="\t")
write.table(timing, file=file.path(outdir, paste0(samplename, "_mutation_timing.txt")), quote=F, row.names=F, sep="\t")
save.image(file.path(outdir, paste0(samplename, "_assignment.RData")))

# Save the PCAWG data - does not require loading of any libraries
save(samplename, final_pcawg11_output, timing, assign_probs, posthoc_stats, file=file.path(outdir, paste0(samplename, "_pcawg11_output.RData")))

########################################################################
# Plot
########################################################################
#' Make assignment plot for both assignment strategies
p = base_plot(snv_binom$plot_data, "ccf", "Consensus binomial assignment (control)") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf)) + xlab("ccf - snv")
p = p + scale_fill_hue(labels = rev(paste0(" ", 
                                           snv_binom$clusters$cluster, " : ", 
                                           round(snv_binom$clusters$ccf, 2), "  ", 
                                           snv_binom$clusters$n_ssms, "  ")))

p3 = base_plot(snv_mtimer$plot_data, "ccf", "Consensus closest cluster assignment") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf)) + xlab("ccf - snv")
p3 = p3 + scale_fill_hue(labels = rev(paste0(" ", 
                                           snv_binom$clusters$cluster, " : ", 
                                           round(snv_mtimer$clusters$ccf, 2), "  ", 
                                           snv_mtimer$clusters$n_ssms, "  ")))

if (!is.null(vcf_indel) && any(indel_mtimer$plot_data$ccf < 1.5)) {
  p4 = base_plot(indel_mtimer$plot_data, "ccf", "Consensus closest cluster assignment") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf)) + xlab("ccf - indel")
  p4 = p4 + scale_fill_hue(labels = rev(paste0(" ", 
                                             indel_binom$clusters$cluster, " : ", 
                                             round(indel_mtimer$clusters$ccf, 2), "  ", 
                                             indel_mtimer$clusters$n_ssms, "  ")))
} else {
  p4 = make_dummy_figure()
}

if (!is.null(vcf_sv) && any(!is.na(sv_mtimer$plot_data$ccf))) {
  p5 = base_plot(sv_mtimer$plot_data, "ccf", "Consensus closest cluster assignment") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf)) + xlab("ccf - sv")
  p5 = p5 + scale_fill_hue(labels = rev(paste0(" ", 
                                               sv_binom$clusters$cluster, " : ", 
                                               round(sv_mtimer$clusters$ccf, 2), "  ", 
                                               sv_mtimer$clusters$n_ssms, "  ")))
} else {
  p5 = make_dummy_figure()
}

if (!is.null(vcf_sv) && any(!is.na(sv_mtimer$plot_data$ccf)) && any(sv_mtimer$plot_data$ccf < 1.5)) {
  all_data = do.call(rbind, list(snv_binom$plot_data, indel_binom$plot_data, sv_binom$plot_data))
  all_data$type = factor(c(rep("SNV", nrow(snv_binom$plot_data)), rep("indel", nrow(indel_binom$plot_data)), rep("sv", nrow(sv_binom$plot_data))), levels=c("SNV", "indel", "sv"))
} else if(!is.null(vcf_indel)) {
  all_data = do.call(rbind, list(snv_binom$plot_data, indel_binom$plot_data))
  all_data$type = factor(c(rep("SNV", nrow(snv_binom$plot_data)), rep("indel", nrow(indel_binom$plot_data))), levels=c("SNV", "indel", "sv"))
} else {
  all_data = snv_binom$plot_data
  all_data$type = factor(c(rep("SNV", nrow(snv_binom$plot_data))), levels=c("SNV", "indel", "sv"))
}
p6 = base_plot(all_data, "ccf", "All data", fill="type") + xlim(0, 1.5) + geom_vline(data=clusters, mapping=aes(xintercept=ccf)) + xlab("ccf")

# if (samplename %in% summary_table$samplename) {
#   power = summary_table$nrpcc[summary_table$samplename==samplename]
#   histology = summary_table$histology_abbreviation[summary_table$samplename==samplename]
# } else {
  power = NA
  histology = NA
# }
title = paste0(samplename, " - ",
               histology, " - ",
               "SNVs ", nrow(snv_binom$plot_data), " - ",
               "Purity ", round(purity, 2), " - ",
               "Ploidy ", ploidy, " - ",
               "Power ", power, " - ",
               "Clust size input ", paste(rev(clusters$n_ssms), collapse=", "))

png(file.path(outdir, paste0(samplename, "_final_assignment.png")), height=400, width=2000)
grid.arrange(p6, p, p3, p4, p5, nrow=1, top=title)
dev.off()
