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


#' Create a custom temp SV vcf file that goes into computeMutCn. Using the dummy vcf means
#' the computeMutCn does not need to be modified for the different use case that are SVs.
create_dummy_sv_vcf = function(svclone_output) {
  
  #' Get preferred breakpoints -> thtat are mapped to cna <- add chr/pos of selected to parse_sv_data output
  #' 
  #' Format into vcf <- load a template
  #' 
  #' return(vcf)
  
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
  return(assignment_ll)
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
  
  cluster_counts = table(best_cluster)
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
assign_moritz = function(MCN, clusters, purity) {
  best_cluster = sapply(MCN$D$CNF, function(x) if (is.na(x)) NA else which.min(abs(x-clusters$proportion)))
  cluster_counts = table(best_cluster)
  clusters_new_2 = data.frame(clusters$cluster, sapply(clusters$cluster, function(x) cluster_counts[[as.character(x)]]), clusters$proportion, clusters$ccf)
  colnames(clusters_new_2) = colnames(clusters)
  
  mcn = mutationBurdenToMutationCopyNumber(burden=MCN$D$altCount / (MCN$D$altCount + MCN$D$wtCount), cellularity=purity, normalCopyNumber=rep(2, nrow(MCN$D)), totalCopyNumber=MCN$D$MajCN + MCN$D$MinCN)
  ccf = MCN$D$CNF / purity
  plot_data_2 = data.frame(mcn=mcn, ccf=ccf, cluster=factor(best_cluster, levels=rev(unique(sort(clusters$cluster)))))
  return(list(plot_data=plot_data_2, clusters_new=clusters_new_2))
}

########################################################################
# PCAWG11 Summary table
########################################################################

get_clusters_entry = function(clusters, assignments_table, indel_assignments=NULL, sv_assignments=NULL, min_clonal_ccf=0.9, max_clonal_ccf=1.1) {
  assignments = table(assignments_table$cluster)
  if (!is.null(indel_assignments)) { indel_assignments = table(indel_assignments$cluster) }
  if (!is.null(sv_assignments)) { sv_assignments = table(sv_assignments$cluster) }
  
  total_muts = sum(assignments)
  if (total_muts < 100) {
    threshold = total_muts*FRAC_SNVS_CLUSTER
  } else {
    threshold = 30
  }
  
  # kept_clusters = names(assignments)[assignments > (total_muts*FRAC_SNVS_CLUSTER) & assignments > MIN_NUM_SNVS_CLUSTER]
  kept_clusters = names(assignments)[assignments > threshold]
  superclones_to_merge = names(assignments)[clusters$location > max_clonal_ccf & clusters$no.of.mutations < total_muts*FRAC_SNVS_CLUSTER & rep(length(kept_clusters) > 1, length(kept_clusters))]
  
  # kept_clusters = names(assignments)
  
  # Count the number of subclones and SNVs assigned
  num_clonal = 0
  num_subclonal = 0
  num_subclones = 0
  num_superclones = 0
  num_superclonal = 0
  indel_clonal = 0
  indel_subclonal = 0
  indel_superclonal = 0
  sample_entry$sv_clonal = 0
  sample_entry$sv_subclonal = 0
  sample_entry$sv_superclonal = 0
  cluster_locations = c()
  cluster_sizes = c()
  for (cluster in rev(kept_clusters)) {
    if (clusters[clusters$cluster.no==cluster,]$location > min_clonal_ccf) {
      # Clonal
      num_clonal = num_clonal + assignments[cluster]
      if (!is.null(indel_assignments)) { indel_clonal = indel_clonal + indel_assignments[cluster] }
      if (!is.null(sv_assignments)) { sv_clonal = sv_clonal + sv_assignments[cluster] }
      
      if (clusters[clusters$cluster.no==cluster,]$location > max_clonal_ccf & ! cluster %in% superclones_to_merge) {
        # Superclonal
        num_superclones = num_superclones + 1
        num_superclonal = num_superclonal + assignments[cluster]
        if (!is.null(indel_assignments)) { indel_superclonal = indel_superclonal + indel_assignments[cluster] }
      }
    } else {
      # Subclonal
      num_subclonal = num_subclonal + assignments[cluster]
      num_subclones = num_subclones + 1
      if (!is.null(indel_assignments)) { indel_subclonal = indel_subclonal + indel_assignments[cluster] }
      if (!is.null(sv_assignments)) { sv_subclonal = sv_subclonal + sv_assignments[cluster] }
    }
    
    #' TODO: indel and sv here required ?
    if (! cluster %in% superclones_to_merge) {
      cluster_locations = c(cluster_locations, clusters[clusters$cluster.no==cluster,]$location)
      
      temp_cluster_size = assignments[cluster]
      if (!is.null(indel_assignments)) { temp_cluster_size = temp_cluster_size + indel_assignments[cluster] }
      if (!is.null(sv_assignments)) { temp_cluster_size = temp_cluster_size + sv_assignments[cluster] }
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
get_summary_table_entry = function(samplename, summary_table, cluster_info, snv_assignment_table, indel_assignment_table=NULL, sv_assignment_table=NULL) {
  sample_entry = summary_table[summary_table$samplename==samplename, ,drop=F]
  # temp because using old summary table
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
  colnames(res) = c("num_subclones", "num_clonal", "num_subclonal", "num_superclones", "num_superclonal", "indel_clonal", "indel_subclonal", "indel_superclonal")
  res$frac_clonal = round(res$num_clonal / (res$num_subclonal+res$num_clonal), 3)
  sample_entry[, colnames(sample_entry) %in% colnames(res)] = res
  
  # Append cluster details
  res = lapply(clust_details, function(x) data.frame(paste(round(x[[1]], 4), collapse=";"), paste(x[[2]], collapse=";")))
  res = do.call(rbind, res)
  colnames(res) = c("cluster_locations", "cluster_sizes")
  sample_entry[, colnames(sample_entry) %in% colnames(res)] = res
  return(sample_entry)
}

