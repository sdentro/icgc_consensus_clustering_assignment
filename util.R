########################################################################
# Overloaded function from MutationTimeR
########################################################################
loadBB <- function(file, round_subclones=F, remove_subclones=F) {
  if (round_subclones & remove_subclones) {
    print("When supplying both rounding and removing to loadBB subclones are removed")
  }


        tab <- read.table(file, header=TRUE, sep='\t')
        r = GRanges(tab$chromosome, IRanges(tab$start, tab$end), strand="*", tab[-3:-1])

        if (remove_subclones) {
          o = findOverlaps(r, r)
          c = countSubjectHits(o)
          subclonal_segments = which(c > 1)
          r = r[-subclonal_segments,]

        } else if (round_subclones) {

          # Check for the ccf column
          if (!"ccf" %in% colnames(tab)) {
            print("No CCF column in segments supplied, stopping")
            q(save="no", status=1)
          }

          # Round subclonal copy number by taking the maximum CCF state
          o = findOverlaps(r, r)
          c = countSubjectHits(o)

          merged_subclonal = data.frame()
          if (any(c > 1)) {
	    subclonal_segments = which(c > 1)[seq(1,length(which(c > 1)), 2)]
            for (i in subclonal_segments) {
              tab_curr = tab[subjectHits(o)[queryHits(o)==i],]
              tab_select = tab_curr[which.max(tab_curr$ccf),]
              tab_select$ccf = 1
              merged_subclonal = rbind(merged_subclonal, tab_select)
            }
          }
          subclonal_segments = subjectHits(o)[which(c > 1)]
          tab_merged = rbind(tab[c==1,], merged_subclonal)
          r = sort(GRanges(tab_merged$chromosome, IRanges(tab_merged$start, tab_merged$end), strand="*", tab_merged[-3:-1]))
        }

        if (length(r)==0) {
          print("No copy number left after filtering, exit now")
          q(save="no", status=0)
        }


        return(r)
}

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

#' Calculate binomial probability of each mutation belonging in each cluster
#' 
#' @param data
#' @param cluster_locations
#' @param purity The sample purity
#' 
getClustLL2 = function(vaf, multiplicity, total_cn, mutcount, wtcount, cluster_locations, purity) {
  assignment_ll = sapply(1:length(cluster_locations), function(c) {
    mutBurdens = mutationCopyNumberToMutationBurden(cluster_locations[c] * multiplicity, total_cn, purity, rep(2, length(mutcount)))
    mutcount*log(mutBurdens) + wtcount*log(1-mutBurdens)
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

#' Assign mutations to the best cluster identified by mtimer
#' 
#' @param MCN Output from computeMutCn
#' @param clusters data.frame with cluster number, size, proportion and ccf
#' @param purity The sample purity
#' 
#' @return A list containing two data.frames: (1) The mutations with mcn, ccf and assigned cluster and (2) Cluster number, size, proportion and ccf
#' @author sd11
assign_mtimer = function(MCN, clusters, purity) {
  print("ASSIGN START")
  print(clusters)
  best_cluster = sapply(MCN$D$CNF, function(x) if (is.na(x)) NA else which.min(abs(x-clusters$proportion)))
  cluster_counts = table(factor(best_cluster, levels=clusters$cluster))
  clusters_new_2 = data.frame(clusters$cluster, clusters$proportion, clusters$ccf, sapply(clusters$cluster, function(x) cluster_counts[[as.character(x)]]))
  colnames(clusters_new_2) = colnames(clusters)
  print("ASSIGN END")
  print(clusters_new_2)
  
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
  print(clusters)
  for (cluster in rev(kept_clusters)) {
    print(paste0("considering cluster ", cluster))
    if (clusters[clusters$cluster.no==cluster,]$location > min_clonal_ccf) {
      print("is clonal")
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
      print("is subclonal")
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
  return(list(clust_stats=data.frame(num_subclones=num_subclones, 
                               num_clonal=num_clonal, num_subclonal=num_subclonal, num_superclones=num_superclones, num_superclonal=num_superclonal, 
                               indel_clonal=indel_clonal, indel_subclonal=indel_subclonal, indel_superclonal=indel_superclonal, 
                               sv_clonal=sv_clonal, sv_subclonal=sv_subclonal, sv_superclonal=sv_superclonal), 
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
get_summary_table_entry = function(samplename, cluster_info, snv_assignment_table, purity, ploidy, sex, is_wgd, indel_assignment_table=NULL, sv_assignment_table=NULL, do_filter=T) {
  sample_entry = data.frame(samplename = samplename, purity = purity, ploidy = ploidy, sex = sex, wgd_status = ifelse(is_wgd, "wgd", "no_wgd"), stringsAsFactors=F)
  # set default values
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
  sample_entry$cluster_locations = NA
  sample_entry$cluster_sizes = NA
  colnames(cluster_info) = c("cluster.no", "no.of.muts", "proportion", "location")
  
  # Get the values for various SNV cluster related columns
  res = get_clusters_entry(cluster_info, 
                           snv_assignment_table, 
                           indel_assignments=indel_assignment_table, 
                           sv_assignments=sv_assignment_table,
                           do_filter=do_filter)
  
  # Saving these for later in the script to be appended to the table
  clust_details = lapply(list(res), function(x) x$clust_details)
  
  # Combine into single row matrix and then copy the updated values into the sample_entry
  res = res$clust_stats
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
pcawg11_output = function(snv_mtimer, indel_mtimer, sv_mtimer, MCN, MCN_indel, MCN_sv, vcf_snv, vcf_indel, vcf_sv, MCN_sv_alt=NULL, vcf_sv_alt=NULL, sv_alt_mtimer=NULL) {
  # Cluster locations
  final_clusters = snv_mtimer$clusters
  
  ###################################################
  # Assignments
  ###################################################
  snv_assignments = data.frame(chr=as.character(seqnames(vcf_snv)), pos=as.numeric(start(vcf_snv)), cluster=snv_mtimer$plot_data$cluster, stringsAsFactors=F)
  if (!is.null(indel_mtimer)) {
    indel_assignments = data.frame(chr=as.character(seqnames(vcf_indel)), pos=as.numeric(start(vcf_indel)), cluster=indel_mtimer$plot_data$cluster, stringsAsFactors=F)
  } else {
    indel_assignments = NULL
  }
  if (!is.null(vcf_sv)) {
    sv_assignments = data.frame(chr=info(vcf_sv)$chr1, pos=info(vcf_sv)$pos1, chr2=info(vcf_sv)$chr2, pos2=info(vcf_sv)$pos2, cluster=sv_mtimer$plot_data$cluster, id=info(vcf_sv)$id, stringsAsFactors=F)
  } else {
    sv_assignments = NULL
  }
  
  if (!is.null(vcf_sv_alt)) {
    sv_alt_assignments = data.frame(chr=info(vcf_sv_alt)$chr2, pos=info(vcf_sv_alt)$pos2, chr2=info(vcf_sv_alt)$chr1, pos2=info(vcf_sv_alt)$pos1, cluster=sv_alt_mtimer$plot_data$cluster, id=info(vcf_sv_alt)$id, stringsAsFactors=F)
  } else {
    sv_alt_assignments = NULL
  }
  
  ###################################################
  # Multiplicities
  ###################################################
  snv_mult = data.frame(chr=snv_assignments$chr, 
                        pos=snv_assignments$pos, 
                        tumour_copynumber=MCN$D$MajCN+MCN$D$MinCN, 
                        multiplicity=MCN$D$MutCN, multiplicity_options=NA, probabilities=NA, stringsAsFactors=F)
  
  if (!is.null(indel_mtimer) && nrow(indel_assignments) > 0) {
    indel_mult = data.frame(chr=indel_assignments$chr, 
                          pos=indel_assignments$pos, 
                          tumour_copynumber=MCN_indel$D$MajCN+MCN_indel$D$MinCN, 
                          multiplicity=MCN_indel$D$MutCN, multiplicity_options=NA, probabilities=NA, stringsAsFactors=F)
  } else {
    indel_mult = NULL
  }
  
  if (!is.null(vcf_sv)) {
    sv_mult = data.frame(chr=sv_assignments$chr, 
                            pos=sv_assignments$pos, 
                            chr2=sv_assignments$chr2,
                            pos2=sv_assignments$pos2,
                            tumour_copynumber=MCN_sv$D$MajCN+MCN_sv$D$MinCN, 
                            multiplicity=MCN_sv$D$MutCN, multiplicity_options=NA, probabilities=NA, stringsAsFactors=F)
  } else {
    sv_mult = NULL
  }
  
  if (!is.null(vcf_sv_alt)) {
    sv_alt_mult = data.frame(chr=sv_alt_assignments$chr, 
                         pos=sv_alt_assignments$pos, 
                         chr2=sv_alt_assignments$chr2,
                         pos2=sv_alt_assignments$pos2,
                         tumour_copynumber=MCN_sv_alt$D$MajCN+MCN_sv_alt$D$MinCN, 
                         multiplicity=MCN_sv_alt$D$MutCN, multiplicity_options=NA, probabilities=NA, stringsAsFactors=F)
  } else {
    sv_alt_mult = NULL
  }
  
  ###################################################
  # Assignment probabilities
  ###################################################
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
                                      r, stringsAsFactors=F)
    
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
    same_end_selected = sv_assignments$chr==sv_assignments_prob$chr & sv_assignments$pos==sv_assignments_prob$pos
    sv_assignments_prob$chr2 = sv_assignments$chr2
    sv_assignments_prob$pos2 = sv_assignments$pos2
    sv_assignments_prob$chr[!same_end_selected] = sv_assignments$chr[!same_end_selected]
    sv_assignments_prob$pos[!same_end_selected] = sv_assignments$pos[!same_end_selected]
  } else {
    sv_assignments_prob = NULL
  }
  
  if (!is.null(vcf_sv_alt)) {
    # Obtain probabilities - SV, non preferred allele
    sv_alt_assignments_prob = get_probs(final_clusters, MCN_sv_alt, vcf_sv_alt)
    same_end_selected = sv_alt_assignments$chr==sv_alt_assignments_prob$chr & sv_alt_assignments$pos==sv_alt_assignments_prob$pos
    sv_alt_assignments_prob$chr2 = sv_alt_assignments$chr2
    sv_alt_assignments_prob$pos2 = sv_alt_assignments$pos2
    sv_alt_assignments_prob$chr[!same_end_selected] = sv_alt_assignments$chr[!same_end_selected]
    sv_alt_assignments_prob$pos[!same_end_selected] = sv_alt_assignments$pos[!same_end_selected]
  } else {
    sv_alt_assignments_prob = NULL
  }
  
  ###################################################
  # Recalculate the size of the clusters
  ###################################################
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
  
  ###################################################
  # Combine SV breakpoint data
  ###################################################
  
  # If both SV breakpoints have been supplied, then merge the output here
  # cluster sizes have been updated above already by using just a single
  # breakpoint, as each SV is represented by multiple
  if (!is.null(vcf_sv) & !is.null(vcf_sv_alt)) {
    sv_assignments_prob = rbind(sv_assignments_prob, sv_alt_assignments_prob)
    sv_assignments = rbind(sv_assignments, sv_alt_assignments)
    sv_mult = rbind(sv_mult, sv_alt_mult)
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


########################################################################
# Structural variants
########################################################################
#' Get copy number at exactly the SV breakpoint location
#' 
#' This creates a dummy copy number profile where the SV position corresponds
#' to its mapped copy number state, which by original coordinates may not
#' necesarily map. This function exists to circumvent this
#' @param bb A copy number profile loaded by loadBB
#' @param vcf_sv A VCF format object as output by \code{prepare_svclone_output}
#' @return A GRanges object with copynumber at exactly the SV bases
copynumber_at_sv_locations = function(bb, vcf_sv) {
  temp_bb = rep(bb[1,c("total_cn", "major_cn", "minor_cn", "clonal_frequency")], length(vcf_sv))
  for (i in 1:length(vcf_sv)) {

    # Catch rare case where there is no copy number for Y, but there is an SV on Y
    if (!as.character(seqnames(vcf_sv)[i]) %in% seqlevels(temp_bb)) {
    	seqlevels(temp_bb) = c(seqlevels(temp_bb), "Y")
    }
    seqnames(temp_bb)[i] = seqnames(vcf_sv)[i]
    # assign seg boundaries in right order as to not create negative length segments
    if (end(temp_bb[i]) < start(vcf_sv)[i]) {
      start(temp_bb)[i] = end(temp_bb[i]) = start(vcf_sv)[i]
    } else {
      end(temp_bb[i]) = start(temp_bb)[i] = start(vcf_sv)[i]
    }
    temp_bb$major_cn[i] = info(vcf_sv)$major_cn[i]
    temp_bb$minor_cn[i] = info(vcf_sv)$minor_cn[i]
    temp_bb$total_cn[i] = temp_bb$major_cn[i] + temp_bb$minor_cn[i]
  }
  return(temp_bb)
}

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
#' @param take_preferred_breakpoint Supply TRUE when to use the preferred breakpoint, set to FALSE when the not preferred breakpoint is to be used (Default: TRUE)
#' 
#' @return A basic vcf file with chromosome/position according to the SVclone
#' preferred mapping (take the end that corresponds to the one mapped to copy
#' number) and t_alt_count and t_ref_count filled in according to the given
#' adjusted support and depth. For convenience the original chr1/pos1 and chr2/pos2
#' columns are also provided.
#' @author sd11
prepare_svclone_output = function(svclone_file, vcf_template, genome, take_preferred_breakpoint=T) {
  dat = read.table(svclone_file, header=T, stringsAsFactors=F, sep="\t")
  
  mutCount = dat$adjusted_support
  # WTCount = dat$adjusted_depth-dat$adjusted_support
  WTCount = array(NA, length(mutCount))
  major_cn = array(NA, length(mutCount))
  minor_cn = array(NA, length(mutCount))
  ids = array(NA, length(mutCount))
  original_pos = array(NA, length(mutCount))
  
  #' Select the preferred SV end from SVclone
  sv_chrom_pos = data.frame()
  for (i in 1:nrow(dat)) {
    # Preferred copy number
    if ((dat$preferred_side[i]==0 & take_preferred_breakpoint) | (dat$preferred_side[i]==1 & !take_preferred_breakpoint)) {
      # sv_chrom_pos = rbind(sv_chrom_pos, data.frame(chrom=as.character(dat$chr1[i]), pos=dat$pos1[i]))
      
      copy_number = as.numeric(unlist(stringr::str_split(dat$gtype1[i], ","))[1:2])
      major_cn[i] = copy_number[1]
      minor_cn[i] = copy_number[2]

      # save the original position, and renumber the position with the row index. it cannot be guaranteed that 
      # pos1, pos2, original_pos1 or original_pos2 are unique. Hence the renumbering      
      original_pos[i] = dat$original_pos1[i]
      sv_chrom_pos = rbind(sv_chrom_pos, data.frame(chrom=as.character(dat$chr1[i]), pos=i, stringsAsFactors=F))
      WTCount[i] = dat$adjusted_norm1[i]
      ids[i] = dat$original_ID[i]
    } else if ((dat$preferred_side[i]==1 & take_preferred_breakpoint) | (dat$preferred_side[i]==0 & !take_preferred_breakpoint)) {
      # sv_chrom_pos = rbind(sv_chrom_pos, data.frame(chrom=as.character(dat$chr2[i]), pos=dat$pos2[i]))
      
      copy_number = as.numeric(unlist(stringr::str_split(dat$gtype2[i], ","))[1:2])
      major_cn[i] = copy_number[1]
      minor_cn[i] = copy_number[2]
      
      # save the original position, and renumber the position with the row index. it cannot be guaranteed that
      # pos1, pos2, original_pos1 or original_pos2 are unique. Hence the renumbering
      original_pos[i] = dat$original_pos2[i]
      sv_chrom_pos = rbind(sv_chrom_pos, data.frame(chrom=as.character(dat$chr2[i]), pos=i, stringsAsFactors=F))
      WTCount[i] = dat$adjusted_norm2[i]
      # replace XYZ1_1 with XYZ1_2 or XYZ1_2 with XYZ1_1
      ids[i] = ifelse(grepl("_1", dat$original_ID[i]), gsub("_1", "_2", dat$original_ID[i]), gsub("_2", "_1", dat$original_ID[i]))
    } 
  }
  
  # Now push this into a VCF format with just alt and ref counts
  v <- readVcf(snv_vcf_file, genome=genome)
  d = data.frame(chromosome=sv_chrom_pos$chrom, position=sv_chrom_pos$pos, stringsAsFactors=F)
  d.gr = makeGRangesFromDataFrame(d, start.field="position", end.field="position")
  d.info = DataFrame(t_alt_count=mutCount, t_ref_count=WTCount, chr1=dat$chr1, pos1=dat$original_pos1, chr2=dat$chr2, pos2=dat$original_pos2, id=ids, major_cn=major_cn, minor_cn=minor_cn, original_pos=original_pos)
  d.v = VCF(rowRanges=d.gr, exptData=metadata(v), geno=geno(v), fixed=rep(fixed(v)[1,], nrow(d)), colData=colData(v), info=d.info)
  return(d.v)
}

#' Map SVs back onto their input location and merge with the existing assignments.
#' @param consensus_vcf_file
#' @param svid_map_file
#' @param sv_assignments
#' @param sv_assignments_prob
#' @return A list with two data.frames, one with hard assignments and one with probabilities. Every consensus SV is reported with their consensus location
remap_svs = function(consensus_vcf_file, svclone_file, sv_assignments, sv_assignments_prob, sv_timing) {
  cons_sv = readVcf(consensus_vcf_file, "GRCh37")
  svclone_dat = read.table(svclone_file, header=T, stringsAsFactors=F, sep="\t")
  
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
                         svid=all_sv_data$id,
                         prob_clonal_early=NA,
                         prob_clonal_late=NA,
                         prob_subclonal=NA,
                         stringsAsFactors=F)
  
  # add extra columns for annotations
  all_sv_data$cluster = NA
  for (i in which(grepl("cluster", colnames(sv_assignments_prob)))) {
    all_sv_data_probs[,colnames(sv_assignments_prob)[i]] = NA
  }
  
  # iterate over all svs and replace
  #for (i in 1:nrow(sv_assignments)) {
  #  if (any(sv_assignments$pos[i] == svmap$pos1)) {
  #    hit = which(sv_assignments$pos[i] == svmap$pos1)
  #  } else {
  #    hit = which(sv_assignments$pos[i] == svmap$pos2)
  #  }
  #  
  #  svid = unlist(strsplit(svmap$original_id[hit], "_", fixed=T))[1]
  #  all_sv_data_row = which(grepl(svid, all_sv_data$id))
  #  # now have mapped i onto all_sv_data_row, which contains both end points of the SV
  #  # save the assignments into the all_data tables
  #  
  #  all_sv_data$cluster[all_sv_data_row] = as.character(sv_assignments$cluster[i])
  #  all_sv_data_probs[all_sv_data_row, grepl("cluster", colnames(all_sv_data_probs))] = sv_assignments_prob[i, grepl("cluster", colnames(sv_assignments_prob))]
  #  all_sv_timing$timing[all_sv_data_row] = as.character(sv_timing$timing[i])
  #}
  for (i in 1:nrow(all_sv_data)) {

	  if (grepl("_2", all_sv_data$id[i])) {
		  row_mapid = gsub("_2", "_1", all_sv_data$id[i])
	  } else {
		  row_mapid = all_sv_data$id[i]
	  }

	  svmap_index = match(row_mapid, svclone_dat$original_ID)
	  if (!is.na(svmap_index)) {

		  match_1 = sv_assignments$chr==svclone_dat$chr1[svmap_index] & sv_assignments$pos==svclone_dat$pos1[svmap_index] & sv_assignments$chr2==svclone_dat$chr2[svmap_index] & sv_assignments$pos2==svclone_dat$pos2[svmap_index]
		  match_2 = sv_assignments$chr==svclone_dat$chr2[svmap_index] & sv_assignments$pos==svclone_dat$pos2[svmap_index] & sv_assignments$chr2==svclone_dat$chr1[svmap_index] & sv_assignments$pos2==svclone_dat$pos1[svmap_index]

		  if (any(match_1) & !any(match_2)) {
			  selection = match_1
		  } else if (any(match_2) & !any(match_1)) {
			  selection = match_2
		  } else {
			  selection = NULL
		  }

		  if (!is.null(selection)) {
			  assign_index = which(selection)
			  if (length(assign_index)==1) {
			  	all_sv_data$cluster[i] = as.character(sv_assignments$cluster[assign_index])
			  	all_sv_data_probs[i, grepl("cluster", colnames(all_sv_data_probs))] = sv_assignments_prob[assign_index, grepl("cluster", colnames(sv_assignments_prob))]
			  	all_sv_timing$timing[i] = as.character(sv_timing$timing[assign_index])
  				all_sv_timing$prob_clonal_early[i] = sv_timing$prob_clonal_early[assign_index]
  				all_sv_timing$prob_clonal_late[i] = sv_timing$prob_clonal_late[assign_index]
  				all_sv_timing$prob_subclonal[i] = sv_timing$prob_subclonal[assign_index]
			  } else {
				# in this scenario multiple chromosome/position entries mapped to this entry in all_sv_data, we'll keep this ambiguous therefore as we cannot map the probabilities uniquely
				  print("Found multiple mapping entries")
			  	print(sv_assignments[assign_index,])
				  all_sv_data$cluster[i] = NA
			  	all_sv_data_probs[i, grepl("cluster", colnames(all_sv_data_probs))] = NA
  				all_sv_timing$timing[i] = NA
  				all_sv_timing$prob_clonal_early[i] = NA
  				all_sv_timing$prob_clonal_late[i] = NA
  				all_sv_timing$prob_subclonal[i] = NA
			  }
		  }
	  }
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

#' Estimate cluster sizes from cluster locations and provided mutation data
#' @param cluster_locations Locations of clusters in CCF
#' @param vcf A VCF file in ICGC PCAWG consensus format, read in
#' @param bb GenomicRanges object, output of loadBB, a loaded in segments file
#' @param purity The tumour purity estimate
#' @param sex Either male or female
#' @param is_wgd Specify TRUE is the tumour is estimated to have undergone a whole genome doubling
#' @param rho_snv Overdispersion parameter to be used when fitting beta-binomial
#' @param xmin Power adjustment parameter, higher values will make it more likely mutations are assigned to smaller clusters
#' @param deltaFreq Distance in CP to match mutation clusters
#' @param max_iters Maximum number of iterations of EM to allow (Default: 10)
#' @param min_snvs_assign_change Minimum number of SNVs to change between two iterations, fewer than this number stops EM (Default: 10)
#' @return A vector with the estimated cluster locations
estimate_cluster_size = function(cluster_locations, vcf, bb, purity, sex, is_wgd, rho_snv, xmin, deltaFreq, max_iters=10, min_snvs_assign_change=10) {

  if (length(vcf)==0) {
	  return(rep(0, length(cluster_locations)))
  }

  # estimate_cluster_size(clusters$ccf, vcf_indel, bb, purity, sex, is_wgd, rho_snv, xmin, deltaFreq)
  alt_count = getAltCount(vcf)
  wt_count = getTumorDepth(vcf)-alt_count
  
  n_muts = length(alt_count)
  clusters = data.frame(cluster=1:length(cluster_locations), location=cluster_locations, n_ssms=NA)
  clusters$proportion = cluster_locations * purity

  if (nrow(clusters)==1) {
	  return(n_muts)
  }

  probs = as.data.frame(matrix(NA, ncol=nrow(clusters), nrow=n_muts))
  colnames(probs) = paste0("prob_cluster_", 1:length(cluster_locations))
  
  # ad-hoc establishment of multiplicity values to calculate starting probabilities
  overlap = findOverlaps(vcf, bb)
  total_cn = bb$total_cn[subjectHits(overlap)]

  # take only those mutations that overlap with copy number (i.e. that have a ccf value)
  alt_count = alt_count[queryHits(overlap)]
  wt_count = wt_count[queryHits(overlap)]
  
  mcn = mutationBurdenToMutationCopyNumber(alt_count/(alt_count+wt_count), total_cn, purity)
  mult = round(mcn)
  mult[mult==0] = 1
  
  # assign to most likely cluster without cluster sizes
  counter = 1
  assignments = list()
  cluster_sizes = list()
  for (i in 1:n_muts) {
    res = getClustLL2(NA, mult[i], total_cn[i], alt_count[i], wt_count[i], clusters$location, purity)
    # res = dtrbetabinom(dpin$mut.count[i],dpin$mut.count[i]+dpin$WT.count[i], ifelse(clusters$location==1, 1-.Machine$double.eps, clusters$location), rho=0, xmin=pmin(dpin$mut.count[i],0))# + .Machine$double.eps), ncol=length(whichStates)
    res = res-max(res)
    res = exp(res)
    res = res / sum(res)
    probs[i,] = res
  }
  assignments[[counter]] = probs
  clusters$n_ssms = colSums(probs, na.rm=T)
  cluster_sizes[[counter]] = clusters
  counter = counter+1
  
  running = T
  while(running) {
    
    # e step - assign using current cluster sizes
    MCN <- computeMutCn(vcf, bb, cluster_sizes[[counter-1]], purity, gender=sex, isWgd=is_wgd, rho=rho_snv, n.boot=0, xmin=xmin, deltaFreq=deltaFreq)
    probs = get_probs(cluster_sizes[[counter-1]], MCN, vcf)[,3:(nrow(cluster_sizes[[counter-1]])+2)]
    assignments[[counter]] = probs
    
    # m step - update cluster sizes
    clusters$n_ssms = colSums(probs, na.rm=T)
    cluster_sizes[[counter]] = clusters
    
    # check stop condition
    if (counter==max_iters | mean(abs(cluster_sizes[[counter-1]]$n_ssms - cluster_sizes[[counter]]$n_ssms)) < min_snvs_assign_change) {
      if (counter==max_iters) { print(paste0("Reached max iters for sample ", samplename)) }
      running = F
    } else {
      assignments[[counter]] = probs
      cluster_sizes[[counter]] = clusters
      counter = counter+1
    }
  }
  return(cluster_sizes[[counter-1]]$n_ssms)
}

#' Obtain a mutation to cluster assignment probability table
get_probs = function(final_clusters, MCN, vcf_snv) {
  
  n_subclones = nrow(final_clusters)-1
  if (n_subclones==0) {
    r = t(t(sapply(MCN$D$pAllSubclones, function(x) 0)))
  } else if (n_subclones==1) {
    r = matrix(unlist(sapply(MCN$D$pAllSubclones, function(x) if(length(x)!=0) x else rep(NA, n_subclones))))
  } else {
    r = matrix(unlist(lapply(MCN$D$pAllSubclones, function(x) if (is.null(x)) { rep(NA, n_subclones) } else { x })), ncol=n_subclones, byrow=T)
  }
  snv_assignments_prob = data.frame(chr=as.character(seqnames(vcf_snv)), 
                                    pos=as.numeric(start(vcf_snv)), 
                                    clone=1-rowSums(r), 
                                    r, stringsAsFactors=F)
  
  if (n_subclones==0) {
    snv_assignments_prob = snv_assignments_prob[,1:3]
  }
  
  # set cluster number in the header
  colnames(snv_assignments_prob) = c("chr", "pos", paste0("cluster_", final_clusters$cluster))
  return(snv_assignments_prob)
}
