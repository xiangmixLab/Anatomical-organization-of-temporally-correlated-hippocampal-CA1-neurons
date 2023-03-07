# Anatomical-organization-of-temporally-correlated-hippocampal-CA1-neurons
Code repository for the publication "Anatomical organization of temporally correlated neural calcium activity in the hippocampal CA1 region"
*MATLAB 2017A and beyond required to run the codes
- cluster_data_prepare: code to calculate the cluster data used in the study. 
- cluster_code: code for temporal and anatomical cluster calculation
  - anatomical_clustering_with_DBSCAN: anatomical cluster calculation, main function is DBSCAN_region_quantify_022422.m
  - Kmean_based_consensus_clustering: temporal cluster calculation, main function is cluster_determine_by_suoqin_NMF_firstPeakCoph_022422.m
  - cluster_calculation_wrapper: wrapper functions of the clustering codes used in cluster_data_prepare


Some of the codes are still under arrangement. Will indicate here when finished
