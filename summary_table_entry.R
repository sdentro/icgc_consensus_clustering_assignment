args = commandArgs(T)
samplename = args[1]

source("~/repo/moritz_mut_assignment/MutationTime.R")

load(paste0(samplename, "_assignment.RData"))


##############################################
#' Get summary table status
get_clusters_entry = function(clusters, assignments_table, indel_assignments=NULL, min_clonal_ccf=0.9, max_clonal_ccf=1.1) {
  assignments = table(assignments_table$cluster)
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
  indel_subclonal=0
  indel_superclonal = 0
  cluster_locations = c()
  cluster_sizes = c()
  for (cluster in rev(kept_clusters)) {
    if (clusters[clusters$cluster.no==cluster,]$location > min_clonal_ccf) {
           # Clonal
           num_clonal = num_clonal + assignments[cluster]
           if (!is.null(indel_assignments)) { indel_clonal = indel_clonal + indel_assignments[cluster] }

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
    }

    if (! cluster %in% superclones_to_merge) {
        cluster_locations = c(cluster_locations, clusters[clusters$cluster.no==cluster,]$location)

        temp_cluster_size = assignments[cluster]
        if (!is.null(indel_assignments)) { temp_cluster_size = temp_cluster_size + indel_assignments[cluster] }
        cluster_sizes = c(cluster_sizes, temp_cluster_size)
    }
  }
  return(list(clust_stats=list(num_subclones, num_clonal, num_subclonal, num_superclones, num_superclonal, indel_clonal, indel_subclonal, indel_superclonal), clust_details=list(cluster_locations, cluster_sizes)))
}
FRAC_SNVS_CLUSTER = 0.01

sample_entry = summary_table[summary_table$samplename==samplename, ,drop=F]
# temp because using old summary table
sample_entry$indel_clonal = 0
sample_entry$indel_subclonal = 0
sample_entry$indel_superclonal = 0


colnames(clusters_new_2) = c("cluster.no", "no.of.muts", "proportion", "location")
res = get_clusters_entry(clusters_new_2, plot_data_2)
clust_stats = lapply(list(res), function(x) x$clust_stats)
# Saving these for later in the script to be appended to the table
clust_details = lapply(list(res), function(x) x$clust_details)
res = as.data.frame(matrix(unlist(clust_stats), ncol=8, byrow=T))
colnames(res) = c("num_subclones", "num_clonal", "num_subclonal", "num_superclones", "num_superclonal", "indel_clonal", "indel_subclonal", "indel_superclonal")
res$frac_clonal = round(res$num_clonal / (res$num_subclonal+res$num_clonal), 3)

sample_entry[, colnames(sample_entry) %in% colnames(res)] = res

#### Cluster details
res = lapply(clust_details, function(x) data.frame(paste(round(x[[1]], 4), collapse=";"), paste(x[[2]], collapse=";")))
res = do.call(rbind, res)
colnames(res) = c("cluster_locations", "cluster_sizes")
sample_entry[, colnames(sample_entry) %in% colnames(res)] = res

write.table(sample_entry, file=paste0(samplename, "_summary_table_entry.txt"), quote=F, sep="\t", row.names=F)


