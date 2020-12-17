# PhD Dissertation
Dissertation related scripts including GTs list generation, differential expression, conservation, and more

This repo houses scripts for each of the major sections of the dissertation project, executed by shell wrappers in the following order:
 
  - [gts_list](./gts_list): generation of the list of human glycosyltransferases (GTs)
 
  - [diff_exp](./diff_exp): summary analysis of DE for the GTs and enrichment of GTs and families in cancer

  - [conservation](./conservation): identification of GT orthologs in mouse and analysis of conserved expression profiles

  - [multi](./multi): integration and merging of multimodal data, including miRNA DE, scRNA cell type specificity, literature mining

  - [glycans](./other): identification of potentially impacted glycans and the top affected residues


  All python and shell code and analysis was initially developed and performed on the MGPC server hosted by GW SMHS.
  
  For more information, go to: https://smhs.gwu.edu/mgpc/  
  
  Additional R scripts (not currently shared here) were executed in a separate local environment, primarily for the purposes of generating images. 
  
  __A note on coding style__: code presented herein was initially developed for internal use, and as such has not been optimized. There exists ample opportunity for improvement with respect to defining reusable functions instead of copying/pasting redundant code blocks, and streamlining loops.
