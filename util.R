########################################################################
# Conversion functions borrowed from dpclust3p
########################################################################

#' Mutation burden to mutation copy number
#' 
#' Function to convert mutation burdens into mutation copy number
#' @param burden A vector containing mutation burdens
#' @param totalCopyNumber A vector with total tumour copynumber
#' @param cellularity Float with the cellularity of this sample
#' @param normalCopyNumber A vector with the total copy number of the normal
#' @return A vector with mutation copy number
#' @author dw9
#' @export
mutationBurdenToMutationCopyNumber = function(burden, totalCopyNumber, cellularity, normalCopyNumber=rep(2, length(burden))) {
  mutCopyNumber = burden/cellularity*(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
  mutCopyNumber[is.nan(mutCopyNumber)] = 0
  return(mutCopyNumber)
}

#' Mutation copy number to mutation burden
#' 
#' Function to convert mutation copy number to mutation burden
#' @param copyNumber A vector containing mutation copy number
#' @param totalCopyNumber A vector with total tumour copynumber
#' @param cellularity Float with the cellularity of this sample
#' @param normalCopyNumber A vector with the total copy number of the normal
#' @return A vector with mutation burdens
#' @author dw9
#' @export
mutationCopyNumberToMutationBurden = function(copyNumber, totalCopyNumber, cellularity, normalCopyNumber=rep(2, length(copyNumber))){
  burden = copyNumber*cellularity/(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
  burden[is.nan(burden) | (burden<0.000001)] = 0.000001
  burden[burden>0.999999] = 0.999999
  return(burden)	
}


########################################################################
# Parsers
########################################################################

#' Obtain for every SV altCount, tumDepth and majorCN
#' 
#' @param cn_assignment_file An SVclone SV to copy number assignment file
#' @param vafs_file An SVclone file reporting VAFs per SV
#' @return A data.frame with altCount, tumDepth and majorCN columns
#' 
#' @author sd11
#' @export
parse_sv_data = function(cn_assignment_file, vafs_file) {
  cn_assignment = read.table(cn_assignment_file, header=T, stringsAsFactors=F, sep="\t")
  vafs = read.table(vafs_file, header=T, stringsAsFactors=F, sep="\t")
  
  parse_entry = function(cn_assignment, vafs, i) {
    vaf = vafs$adjusted_VAF[i]
    pos1_cn = cn_assignment$pos1_bb_CN[i]
    pos2_cn = cn_assignment$pos2_bb_CN[i]
    total_cn = as.numeric(unlist(stringr::str_split(ifelse(pos1_cn!="", pos1_cn, pos2_cn), ","))[1])
    depth = cn_assignment$depth[i]
    # sv_data = rbind(sv_data, data.frame(altCount=vaf*depth, tumDepth=depth, majorCN=total_cn))
    return(data.frame(altCount=vaf*depth, tumDepth=depth, majorCN=total_cn))
  }
  
  res = lapply(1:nrow(vafs), parse_entry, cn_assignment=cn_assignment, vafs=vafs)
  return(as.data.frame(do.call(rbind, res)))
}

#' Transform the SVclone raw data into a vcf object
#' 
#' @param svclone_file SVclone vaf file with copy number mapping
#' @param vcf_template Path to a vcf template file
#' @param genome The reference genome to pass to readVcf
#' 
#' @return A basic vcf file with chromosome/position according to the SVclone
#' preferred mapping (take the end that corresponds to the one mapped to copy
#' number) and t_alt_count and t_ref_count filled in according to the given
#' adjusted support and depth. For convenience the original chr1/pos1 and chr2/pos2
#' columns are also provided.
#' @author sd11
prepare_svclone_output = function(svclone_file, vcf_template, genome) {
  dat = read.table(svclone_file, header=T, stringsAsFactors=F, sep="\t")
  
  mutCount = dat$adjusted_support
  WTCount = dat$adjusted_depth-dat$adjusted_support
  
  #' Select the preferred SV end from SVclone
  sv_chrom_pos = data.frame()
  for (i in 1:nrow(dat)) {
    # Preferred copy number
    if (dat$preferred_side[i]==0) {
      sv_chrom_pos = rbind(sv_chrom_pos, data.frame(chrom=as.character(dat$chr1[i]), pos=dat$pos1[i]))
    } else {
      sv_chrom_pos = rbind(sv_chrom_pos, data.frame(chrom=as.character(dat$chr2[i]), pos=dat$pos2[i]))
    }
  }
  
  # Now push this into a VCF format with just alt and ref counts
  v <- readVcf(snv_vcf_file, genome=genome)
  d = data.frame(chromosome=sv_chrom_pos$chrom, position=sv_chrom_pos$pos)
  d.gr = makeGRangesFromDataFrame(d, start.field="position", end.field="position")
  d.info = DataFrame(t_alt_count=mutCount, t_ref_count=WTCount, chr1=dat$chr1, pos1=dat$pos1, chr2=dat$chr2, pos2=dat$pos2)
  d.v = VCF(rowRanges=d.gr, exptData=metadata(v), geno=geno(v), fixed=rep(fixed(v)[1,], nrow(d)), colData=colData(v), info=d.info)
  
  return(d.v)
}

########################################################################
# Mutation Assignment
########################################################################
#' Calculate binomial probability of each mutation belonging in each cluster
#' 
#' @param data
#' @param cluster_locations
#' @param purity The sample purity
#' 
getClustLL = function(data, cluster_locations, purity) {
  assignment_ll = sapply(1:length(cluster_locations), function(c) {
    mutBurdens = mutationCopyNumberToMutationBurden(cluster_locations[c] * data$MutCN, data$MajCN+data$MinCN, purity, rep(2, nrow(data)))
    data$altCount*log(mutBurdens) + data$wtCount*log(1-mutBurdens)
  })
  return(matrix(assignment_ll, ncol=length(cluster_locations)))
}

#' Assign mutations to the cluster with highest binomial likelihood
#' 
#' @param MCN Output from computeMutCn
#' @param clusters data.frame with cluster number, size, proportion and ccf
#' @param purity The sample purity
#' 
#' @return A list containing two data.frames: (1) The mutations with mcn, ccf and assigned cluster and (2) Cluster number, size, proportion and ccf
#' @author sd11
assign_binom_ll = function(MCN, clusters, purity) {
  clust_assign_ll = getClustLL(MCN$D, clusters$proportion/purity, purity)
  best_cluster = unlist(apply(clust_assign_ll, 1, function(x) if (all(is.na(x))) NA else which.max(x)))
  
  cluster_counts = table(factor(best_cluster, levels=clusters$cluster))
  clusters_new = data.frame(clusters$cluster, sapply(clusters$cluster, function(x) cluster_counts[[as.character(x)]]), clusters$proportion, clusters$ccf)
  colnames(clusters_new) = colnames(clusters)
  
  #' Calculate MCN and CCF for plotting
  mcn = mutationBurdenToMutationCopyNumber(burden=MCN$D$altCount / (MCN$D$altCount + MCN$D$wtCount), cellularity=purity, normalCopyNumber=rep(2, nrow(MCN$D)), totalCopyNumber=MCN$D$MajCN + MCN$D$MinCN)
  ccf = mcn / MCN$D$MutCN
  plot_data = data.frame(mcn=mcn, ccf=ccf, cluster=factor(best_cluster, levels=rev(unique(sort(clusters_new$cluster)))))
  return(list(plot_data=plot_data, clusters_new=clusters_new))
}

#' Assign mutations to the best cluster identified by Moritz method
#' 
#' @param MCN Output from computeMutCn
#' @param clusters data.frame with cluster number, size, proportion and ccf
#' @param purity The sample purity
#' 
#' @return A list containing two data.frames: (1) The mutations with mcn, ccf and assigned cluster and (2) Cluster number, size, proportion and ccf
#' @author sd11
assign_mtimer = function(MCN, clusters, purity) {
  best_cluster = sapply(MCN$D$CNF, function(x) if (is.na(x)) NA else which.min(abs(x-clusters$proportion)))
  cluster_counts = table(factor(best_cluster, levels=clusters$cluster))
  clusters_new_2 = data.frame(clusters$cluster, sapply(clusters$cluster, function(x) cluster_counts[[as.character(x)]]), clusters$proportion, clusters$ccf)
  colnames(clusters_new_2) = colnames(clusters)
  
  mcn = mutationBurdenToMutationCopyNumber(burden=MCN$D$altCount / (MCN$D$altCount + MCN$D$wtCount), cellularity=purity, normalCopyNumber=rep(2, nrow(MCN$D)), totalCopyNumber=MCN$D$MajCN + MCN$D$MinCN)
  if (all(is.na(best_cluster))) {
    mcn = rep(NA, length(best_cluster))
  }
  ccf = MCN$D$CNF / purity
  plot_data_2 = data.frame(mcn=mcn, ccf=ccf, cluster=factor(best_cluster, levels=rev(unique(sort(clusters$cluster)))))
  return(list(plot_data=plot_data_2, clusters_new=clusters_new_2))
}

########################################################################
# PCAWG11 Summary table
########################################################################

get_clusters_entry = function(clusters, assignments_table, indel_assignments=NULL, sv_assignments=NULL, min_clonal_ccf=0.9, max_clonal_ccf=1.1, do_filter=T) {
  assignments = table(assignments_table$cluster)
  # if (!is.null(indel_assignments)) { indel_assignments = table(indel_assignments$cluster) }
  # if (!is.null(sv_assignments)) { sv_assignments = table(sv_assignments$cluster) }
  
  if (do_filter) {
    total_muts = sum(assignments)
    if (total_muts < 100) {
      threshold = total_muts*FRAC_SNVS_CLUSTER
    } else {
      threshold = 30
    }
    
    # kept_clusters = names(assignments)[assignments > (total_muts*FRAC_SNVS_CLUSTER) & assignments > MIN_NUM_SNVS_CLUSTER]
    kept_clusters = names(assignments)[assignments > threshold]
    superclones_to_merge = names(assignments)[clusters$location > max_clonal_ccf & clusters$no.of.mutations < total_muts*FRAC_SNVS_CLUSTER & rep(length(kept_clusters) > 1, length(kept_clusters))]
  } else {
    kept_clusters = names(assignments)
    superclones_to_merge = c()
  }  
  
  # Count the number of subclones and SNVs assigned
  num_clonal = 0
  num_subclonal = 0
  num_subclones = 0
  num_superclones = 0
  num_superclonal = 0
  indel_clonal = 0
  indel_subclonal = 0
  indel_superclonal = 0
  sv_clonal = 0
  sv_subclonal = 0
  sv_superclonal = 0
  cluster_locations = c()
  cluster_sizes = c()
  for (cluster in rev(kept_clusters)) {
    if (clusters[clusters$cluster.no==cluster,]$location > min_clonal_ccf) {
      # Clonal
      num_clonal = num_clonal + assignments[cluster]
      # if (!is.null(indel_assignments)) { indel_clonal = indel_clonal + indel_assignments[cluster] }
      if (!is.null(indel_assignments)) { indel_clonal = indel_clonal + sum(indel_assignments$cluster==cluster, na.rm=T) }
      # if (!is.null(sv_assignments)) { sv_clonal = sv_clonal + sv_assignments[cluster] }
      if (!is.null(sv_assignments)) { sv_clonal = sv_clonal + sum(sv_assignments$cluster==cluster, na.rm=T) }
      
      if (clusters[clusters$cluster.no==cluster,]$location > max_clonal_ccf & ! cluster %in% superclones_to_merge) {
        # Superclonal
        num_superclones = num_superclones + 1
        num_superclonal = num_superclonal + assignments[cluster]
        # if (!is.null(indel_assignments)) { indel_superclonal = indel_superclonal + indel_assignments[cluster] }
        if (!is.null(indel_assignments)) { indel_superclonal = indel_superclonal + sum(indel_assignments$cluster==cluster, na.rm=T) }
        if (!is.null(sv_assignments)) { sv_superclonal = sv_superclonal + sum(sv_assignments$cluster==cluster, na.rm=T) }
      }
    } else {
      # Subclonal
      num_subclonal = num_subclonal + assignments[cluster]
      num_subclones = num_subclones + 1
      # if (!is.null(indel_assignments)) { indel_subclonal = indel_subclonal + indel_assignments[cluster] }
      if (!is.null(indel_assignments)) { indel_subclonal = indel_subclonal + sum(indel_assignments$cluster==cluster, na.rm=T) }
      # if (!is.null(sv_assignments)) { sv_subclonal = sv_subclonal + sv_assignments[cluster] }
      if (!is.null(sv_assignments)) { sv_subclonal = sv_subclonal + sum(sv_assignments$cluster==cluster, na.rm=T) }
    }
    
    if (! cluster %in% superclones_to_merge) {
      cluster_locations = c(cluster_locations, clusters[clusters$cluster.no==cluster,]$location)
      
      temp_cluster_size = assignments[cluster]
      cluster_sizes = c(cluster_sizes, temp_cluster_size)
    }
  }
  return(list(clust_stats=list(num_subclones, 
                               num_clonal, num_subclonal, num_superclones, num_superclonal, 
                               indel_clonal, indel_subclonal, indel_superclonal, 
                               sv_clonal, sv_subclonal, sv_superclonal), 
              clust_details=list(cluster_locations, cluster_sizes)))
}


FRAC_SNVS_CLUSTER = 0.01

#' Creates a summary table entry
#' 
#' @param samplename The samplename
#' @param summary_table An existing summary table with an entry for this sample
#' @param clusters_new_2
#' @param plot_data_2
#' 
#' @return 
#' @author sd11
get_summary_table_entry = function(samplename, summary_table, cluster_info, snv_assignment_table, indel_assignment_table=NULL, sv_assignment_table=NULL, do_filter=T) {
  sample_entry = summary_table[summary_table$samplename==samplename, ,drop=F]
  # temp because using old summary table
  sample_entry$num_subclones = 0
  sample_entry$num_clonal = 0
  sample_entry$num_subclonal = 0
  sample_entry$num_superclones = 0
  sample_entry$num_superclonal = 0
  sample_entry$indel_clonal = 0
  sample_entry$indel_subclonal = 0
  sample_entry$indel_superclonal = 0
  sample_entry$sv_clonal = 0
  sample_entry$sv_subclonal = 0
  sample_entry$sv_superclonal = 0
  colnames(cluster_info) = c("cluster.no", "no.of.muts", "proportion", "location")
  
  # Get the values for various SNV cluster related columns
  res = get_clusters_entry(cluster_info, 
                           snv_assignment_table, 
                           indel_assignments=indel_assignment_table, 
                           sv_assignments=sv_assignment_table)
  clust_stats = lapply(list(res), function(x) x$clust_stats)
  
  # Saving these for later in the script to be appended to the table
  clust_details = lapply(list(res), function(x) x$clust_details)
  
  # Combine into single row matrix and then copy the updated values into the sample_entry
  res = as.data.frame(matrix(unlist(clust_stats), ncol=11, byrow=T))
  colnames(res) = c("num_subclones", "num_clonal", "num_subclonal", "num_superclones", "num_superclonal", "indel_clonal", "indel_subclonal", "indel_superclonal", "sv_clonal", "sv_subclonal", "sv_superclonal")
  sample_entry[, colnames(sample_entry) %in% colnames(res)] = res
  sample_entry$frac_clonal = round(res$num_clonal / (res$num_subclonal+res$num_clonal), 3)
  
  # Append cluster details
  res = lapply(clust_details, function(x) data.frame(paste(round(x[[1]], 4), collapse=";"), paste(x[[2]], collapse=";")))
  res = do.call(rbind, res)
  colnames(res) = c("cluster_locations", "cluster_sizes")
  sample_entry[, colnames(sample_entry) %in% colnames(res)] = res
  return(sample_entry)
}

########################################################################
# PCAWG11 Calibration format
########################################################################

pcawg11_output = function(snv_mtimer, indel_mtimer, sv_mtimer, MCN, MCN_indel, MCN_sv, vcf_sv, consensus_vcf_file, svid_map_file) {
  # Cluster locations
  final_clusters = snv_mtimer$clusters
  
  # Assignments
  snv_assignments = data.frame(chr=as.character(seqnames(vcf_snv)), pos=as.numeric(start(vcf_snv)), cluster=snv_mtimer$plot_data$cluster)
  if (!is.null(indel_mtimer)) {
    indel_assignments = data.frame(chr=as.character(seqnames(vcf_indel)), pos=as.numeric(start(vcf_indel)), cluster=indel_mtimer$plot_data$cluster)
  } else {
    indel_assignments = NULL
  }
  if (!is.null(vcf_sv)) {
    sv_assignments = data.frame(chr=info(vcf_sv)$chr1, pos=info(vcf_sv)$pos1, chr2=info(vcf_sv)$chr2, pos2=info(vcf_sv)$pos2, cluster=sv_mtimer$plot_data$cluster)
  } else {
    sv_assignments = NULL
  }
  
  # Multiplicities
  snv_mult = data.frame(chr=snv_assignments$chr, 
                        pos=snv_assignments$pos, 
                        tumour_copynumber=MCN$D$MajCN+MCN$D$MinCN, 
                        multiplicity=MCN$D$MutCN, multiplicity_options=NA, probabilities=NA)
  
  if (!is.null(indel_mtimer) && nrow(indel_assignments) > 0) {
    indel_mult = data.frame(chr=indel_assignments$chr, 
                          pos=indel_assignments$pos, 
                          tumour_copynumber=MCN_indel$D$MajCN+MCN_indel$D$MinCN, 
                          multiplicity=MCN_indel$D$MutCN, multiplicity_options=NA, probabilities=NA)
  } else {
    indel_mult = NULL
  }
  
  if (!is.null(vcf_sv)) {
    sv_mult = data.frame(chr=sv_assignments$chr, 
                            pos=sv_assignments$pos, 
                            chr2=sv_assignments$chr2,
                            pos2=sv_assignments$pos2,
                            tumour_copynumber=MCN_sv$D$MajCN+MCN_sv$D$MinCN, 
                            multiplicity=MCN_sv$D$MutCN, multiplicity_options=NA, probabilities=NA)
  } else {
    sv_mult = NULL
  }
  
  get_probs = function(final_clusters, MCN, vcf_snv) {
    
    n_subclones = nrow(final_clusters)-1
    if (n_subclones==0) {
      r = t(t(sapply(MCN$D$pAllSubclones, function(x) 0)))
    } else if (n_subclones==1) {
      #r = t(t(sapply(MCN$D$pAllSubclones, function(x) if(length(x)!=0) x else rep(NA, n_subclones))))
      r = matrix(unlist(sapply(MCN$D$pAllSubclones, function(x) if(length(x)!=0) x else rep(NA, n_subclones))))
    } else {
      # r = t(sapply(MCN$D$pAllSubclones, function(x) if(length(x)!=0) x else rep(1, n_subclones)))
      # r = matrix(unlist(sapply(MCN$D$pAllSubclones, function(x) if(length(x)!=0) x else rep(1, n_subclones))), ncol=n_subclones, byrow=T)
      r = matrix(unlist(lapply(MCN$D$pAllSubclones, function(x) if (is.null(x)) { rep(NA, n_subclones) } else { x })), ncol=n_subclones, byrow=T)
    }
    snv_assignments_prob = data.frame(chr=as.character(seqnames(vcf_snv)), 
                                      pos=as.numeric(start(vcf_snv)), 
                                      clone=1-rowSums(r), 
                                      r)
    
    if (n_subclones==0) {
      snv_assignments_prob = snv_assignments_prob[,1:3]
    }
    
    # set cluster number in the header
    colnames(snv_assignments_prob) = c("chr", "pos", paste0("cluster_", final_clusters$cluster))
    return(snv_assignments_prob)
  }
  
  # Obtain probabilities of assignments - snv
  snv_assignments_prob = get_probs(final_clusters, MCN, vcf_snv)
  if (!is.null(indel_mtimer) && nrow(indel_assignments) > 0) {
    # Obtain probabilities - indel
    indel_assignments_prob = get_probs(final_clusters, MCN_indel, vcf_indel)
  } else {
    indel_assignments_prob = NULL
  }
  
  if (!is.null(vcf_sv)) {
    # Obtain probabilities - SV
    sv_assignments_prob = get_probs(final_clusters, MCN_sv, vcf_sv)
    sv_assignments_prob$chr2 = sv_assignments$chr2
    sv_assignments_prob$pos2 = sv_assignments$pos2
  } else {
    sv_assignments_prob = NULL
  }
  
  # Recalculate the size of the clusters
  final_clusters$n_snvs = colSums(snv_assignments_prob[, grepl("cluster", colnames(snv_assignments_prob)), drop=F], na.rm=T)
  if (!is.null(indel_mtimer) && nrow(indel_assignments) > 0) {
    final_clusters$n_indels = colSums(indel_assignments_prob[, grepl("cluster", colnames(indel_assignments_prob)), drop=F], na.rm=T)
  } else {
    final_clusters$n_indels = NA
  }
  if (!is.null(vcf_sv)) {
    final_clusters$n_svs = sv_mtimer$clusters$n_ssms
    if (nrow(final_clusters)==1) {
	final_clusters$n_svs = sum(sv_assignments_prob[, grepl("cluster", colnames(sv_assignments_prob))], na.rm=T)
    } else {
    	final_clusters$n_svs = colSums(sv_assignments_prob[, grepl("cluster", colnames(sv_assignments_prob)), drop=F], na.rm=T)
    }
  } else {
    final_clusters$n_svs = NA 
  }
   
  return(list(final_clusters=final_clusters, 
              snv_assignments=snv_assignments, 
              indel_assignments=indel_assignments,
              sv_assignments=sv_assignments,
              snv_mult=snv_mult,
              indel_mult=indel_mult,
              sv_mult=sv_mult,
              snv_assignments_prob=snv_assignments_prob,
              indel_assignments_prob=indel_assignments_prob,
              sv_assignments_prob=sv_assignments_prob))
}

#' Map SVs back onto their input location and merge with the existing assignments.
#' @param consensus_vcf_file
#' @param svid_map_file
#' @param sv_assignments
#' @param sv_assignments_prob
#' @return A list with two data.frames, one with hard assignments and one with probabilities. Every consensus SV is reported with their consensus location
remap_svs = function(consensus_vcf_file, svid_map_file, sv_assignments, sv_assignments_prob, sv_timing) {
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
  all_sv_data = data.frame(chr=as.character(seqnames(cons_sv)), pos=start(cons_sv), chr2=m[,1], pos2=as.numeric(m[,2]), id=rownames(as.data.frame(rowRanges(cons_sv))))
  all_sv_data_probs = all_sv_data
  all_sv_timing = data.frame(chromosome=all_sv_data$chr, 
                         position=all_sv_data$pos,
                         mut_type=rep("SV", nrow(all_sv_data)),
                         timing=NA, 
                         chromosome2=all_sv_data$chr2,
                         position2=all_sv_data$pos2, 
                         stringsAsFactors=F)
  
  # add extra columns for annotations
  all_sv_data$cluster = NA
  for (i in which(grepl("cluster", colnames(sv_assignments_prob)))) {
    all_sv_data_probs[,colnames(sv_assignments_prob)[i]] = NA
  }
  
  # iterate over all svs and replace
  for (i in 1:nrow(sv_assignments)) {
    if (any(sv_assignments$pos[i] == svmap$pos1)) {
      hit = which(sv_assignments$pos[i] == svmap$pos1)
    } else {
      hit = which(sv_assignments$pos[i] == svmap$pos2)
    }
    
    svid = unlist(strsplit(svmap$original_id[hit], "_", fixed=T))[1]
    all_sv_data_row = which(grepl(svid, all_sv_data$id))
    # now have mapped i onto all_sv_data_row, which contains both end points of the SV
    # save the assignments into the all_data tables
    
    all_sv_data$cluster[all_sv_data_row] = as.character(sv_assignments$cluster[i])
    all_sv_data_probs[all_sv_data_row, grepl("cluster", colnames(all_sv_data_probs))] = sv_assignments_prob[i, grepl("cluster", colnames(sv_assignments_prob))]
    all_sv_timing$timing[all_sv_data_row] = as.character(sv_timing$timing[i])
  }
  return(list(sv_assignments=all_sv_data, sv_assignments_prob=all_sv_data_probs, sv_timing=all_sv_timing))
}

########################################################################
# Plotting
########################################################################
make_dummy_figure = function() { return(ggplot(dat=data.frame(x=rnorm(10), y=rnorm(10))) + aes(x=x, y=y) + geom_point(size=0.00001) + theme_bw()) }

base_plot = function(plot_data, x_variable, title=NA, fill="cluster") {
  p = ggplot(plot_data) + aes_string(x=x_variable, y="..count..", fill=fill) + 
    geom_histogram(binwidth=0.05, colour="black", position="stack") + 
    ylab("Count") + theme(legend.position="bottom") + scale_fill_discrete(drop = FALSE)
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

calc_exp_mutreads_ccf = function(ccf, purity, ploidy, mean_depth) {
  return(ccf / (purity*ploidy + (1-purity)*2) * mean_depth)
}

mergeClustersByMutreadDiff = function(clusters, purity, ploidy, vcf_snv, min_read_diff) {
	print(clusters)
  clusters_new = clusters
  exp_reads = sapply(clusters$ccf, calc_exp_mutreads_ccf, purity=purity, ploidy=ploidy, mean_depth=mean(getTumorDepth(vcf_snv), na.rm=T))
  ccf_diff = exp_reads[1:(length(exp_reads)-1)] - exp_reads[2:length(exp_reads)]
  print("MERGING")
  print(exp_reads)
  print(ccf_diff)





  if (any(ccf_diff < min_read_diff)) {
    
    #' Iteratively merge a pair of clusters untill no more pairs within distance can be found
    merged = T
    while(merged) {
      merged = F
      
      exp_reads = sapply(clusters_new$ccf, calc_exp_mutreads_ccf, purity=purity, ploidy=ploidy, mean_depth=mean(getTumorDepth(vcf_snv), na.rm=T))
      ccf_diff = exp_reads[1:(length(exp_reads)-1)] - exp_reads[2:length(exp_reads)]
      to_merge = which(ccf_diff < min_read_diff)
      
      if (length(to_merge)==0) {
        merged = F
        break
      } else {
        i = to_merge[1]
        print(paste0("merging ", i, " and ", i+1))
        clusters_new$ccf[i] = sum(clusters_new$ccf[c(i, i+1)]*clusters_new$n_ssms[c(i, i+1)]) / sum(clusters_new$n_ssms[c(i, i+1)])
        clusters_new$n_ssms[i] = sum(clusters_new$n_ssms[c(i, i+1)])
        clusters_new = clusters_new[-(i+1),]
        merged = T
      }
      
      if (nrow(clusters_new)==1) {
        merged = F
        break
      }
    }
  }
  clusters_new$proportion = clusters_new$ccf * purity
  print(clusters_new)
  # sorting and renumbering clusters
  clusters_new = clusters_new[with(clusters_new, order(proportion, decreasing=T)),]
  clusters_new$cluster = 1:nrow(clusters_new)
  return(clusters_new)
}




