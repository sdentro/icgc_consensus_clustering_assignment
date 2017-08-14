args = commandArgs(T)
infile = args[1]
outdir = args[2]

load(infile)
samplename = paste(unlist(stringr::str_split(basename(infile), "_"))[1:3], collapse="_")
write.table(final_pcawg11_output$final_clusters, file=paste0(outdir, samplename, "_subclonal_structure.txt"), row.names=F, sep="\t", quote=F)
write.table(final_pcawg11_output$snv_assignments, file=paste0(outdir, samplename,"_mutation_assignments.txt"), row.names=F, sep="\t", quote=F)
