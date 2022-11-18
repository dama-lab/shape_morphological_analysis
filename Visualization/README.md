Tentative list of  Functions to generate Paper

- Ma, D., Cardoso, M. J., Zuluaga, M. A., Modat, M., Powell, N. M.,  Wiseman, F. K., Cleary, J. O., Sinclair, B., Harrison, I. F., Siow, B.,  Popuri, K., Lee, S., Matsubara, J. A., Sarunic, M. V, Beg, M. F.,  Tybulewicz, V. L. J., Fisher, E. M. C., Lythgoe, M. F., & Ourselin,  S. (2020). **Substantially thinner internal granular layer and  reduced molecular layer surface in the cerebellum of the Tc1 mouse model of Down Syndrome – a comprehensive morphometric analysis with active  staining contrast-enhanced MRI**. NeuroImage, 117271. https://doi.org/https://doi.org/10.1016/j.neuroimage.2020.117271
- Ma, D., Cardoso, M. J., Zuluaga, M. A., Modat, M., Powell, N.,  Wiseman, F., Tybulewicz, V., Fisher, E., Lythgoe, M. F., & Ourselin, S. (2015). **Grey Matter Sublayer Thickness Estimation in the Mouse Cerebellum**. In Medical Image Computing and Computer Assisted Intervention 2015 (pp. 644–651). https://doi.org/10.1007/978-3-319-24574-4_77

## Figure 5 - Surface morphometry
surface_morphometry();## Figure 6 - zscape:
tc1_zscape();



 Function list

 ===== Plot HelperFunction ======
 - pcshowRegistration <<==<< pcShowRegistrationDemo 
   Show registration performance (before v.s. after)
 - groupwise_purkinje_thickness

 	show thickness difference on purkinje layer
 - surface_morphometry
 - statColormap: Generate colormap for different statistical plot
 - individualSegmentationQC: Present QC of individual Segmentation
 - parcellatedSurThick: 

 - extractVol: get label volumes from segmentation matrix
 - parcellatedSurThick: get mean surface thickness for each parcellated label from vertex information

 ===== point Cloud Registration HelperFunctions ======
 - pcRegistration: single registration for a pair of point cloud
 - pcRegistrationGroup: 
 - pcRegistrationGlobalICP: Global ICP registration using the point cloud class from https://www.geo.tuwien.ac.at/downloads/pg/pctools/publish/globalICPHelp/globalICP.html  


 - pcDistance: find the pointwise distance (signed scalar) between two rigidly aligned point clouds
 - pcReorder: reorder point cloud based on the index of the other point cloud

 ======= vertex based analysis
 - vertexDistance: find the pointwise distance (signed scalar) between two rigidly aligned point clouds
 - vertexReorder: reorder point cloud based on the index of the other vertex points
 - vertexNorm: find the vertex norm of a surface
 - surfaceDisplace2Distance: convert surface displacement to distance

 ======= preprocessing
 - average_nonrigid_groupmean: Average individual measurements to groupwise average
 - groupAverageCorticalSublayerSegmentation: Segment cortical sublayers (x2) from AMBMC

 - cortical_sublayer_segmentation_individual 
 - cortical_sublayer_segmentation_groupMean

 ===== Other Processing HelperFunctions ======
 - extractStruct(vol, seg, labelNo): extract a single structure from the vol with segmentation