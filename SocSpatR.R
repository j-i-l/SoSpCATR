compare_soc_spat=function(node_locations, soc_membership){
	#node_locations = matrix of x,y coordinates for each node
	#soc_membership = membership vector from social partitioning - same order as node_locations
	#returns:  homogeneity, completeness, v-measure
	# see also:
	#     | V-Measure: A conditional entropy-based external cluster evaluation measure
      #	| Andrew Rosenberg and Julia Hirschberg, 2007
      #	| http://acl.ldc.upenn.edu/D/D07/D07-1043.pdf

	require(sabre)
	#Perform k-means clustering on the provided node locations.
	spat_membership=kmeans(node_locations,length(unique(soc_membership)),iter.max=300,algorithm="Lloyd")$cluster
	return(vmeasure(soc_membership,spat_membership))
}






