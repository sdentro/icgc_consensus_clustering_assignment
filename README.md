# ICGC PCAWG-11 consensus mutation assignment

### Overview

This repo contains code used to assign mutations (SNV, indel and SV) to consensus mutation clusters from subclonal reconstruction methods.

### Dependencies

 * VariantAnnotation
 * VGAM
 * ggplot2
 * gridExtra
 * grid

### How to run

### Output - run.R

This script runs the assignments

 * [samplename]_assignment.RData : All output stored in an RData archive (requires loading of dependencies)
 * [samplename]_pcawg11_output.RData : All PCAWG-11 output (no dependencies required)
 * [samplename]_summary_table_entry.txt : A PCAWG-11 summary table entry for this sample
 * [samplename]_final_assignment.png : A figure showing the data and assignment, with comparison to regular binomial assignment

### Output - run_final_output.R

This script produces the output that is shared PCAWG-wide

 
### Output - run_final_mutccf.R

This script produces per mutation raw CCF estimates

### Output - run_pcawg11.R

Produces PCAWG-11 formatted output for internal use

 * [samplename]_subclonal_structure.txt
 * [samplename]_mutation_assignments.txt
 * [samplename]_multiplicity.txt