# Shape & Morphological Analysis and Rendering Toolkits (SMART)
Shape and Morphological Analysis (e.g. Cortical Shape, Surface Area, etc.)

Functionality includes:

- Laplacian-based structural thickness estimation (e.g. cerebral/cerebellar cortical thickness and hippocampal subfield thickness)
- Groupwise surface-based statistical-mapping (for shape metrics such as thickness, surface area, and structural volume)

It utilizes the "**MASMAT**" (Multi Atlas Segmentation and Morphometric Analysis Toolkit) to preprocess the groupwise  analysis pipeline: https://github.com/dama-lab/multi-atlas-segmentation

# QuickCheck

<img src="https://github.com/dama-lab/mouse-brain-atlas/blob/master/docs/segmentation_qc_all_28_original_aspect_ratio.png" alt="QuickCheck demo" style="zoom:50%;" />
<img src="https://github.com/dama-lab/mouse-brain-atlas/blob/master/docs/pipeline_results_step1.png" alt="Pipeline demo" style="zoom:50%;" />
<img src="https://ars.els-cdn.com/content/image/1-s2.0-S1053811920307576-gr5_lrg.jpg"  alt="Cortical Morpometry" style="zoom:50%;" />

# Folder structure:

## Processing

> Contain Image Data processing and Morphological Analysis functions

## Visualization

> Contain Shape and thickness visualization functions

## HelperFunctions

> Various help functions

### Laplacian Thickness

> To Calculate Laplacian-Based Cortical ThicknessEstimation

## Toolbox

> External Toolboxes

### iso2mesh

> To construct surface mesh from binary segmentation labels

# Reference:

  - Ma, D., Cardoso, M. J., Zuluaga, M. A., Modat, M., Powell, N. M., Wiseman, F. K., Cleary, J. O., Sinclair, B., Harrison, I. F., Siow, B., Popuri, K., Lee, S., Matsubara, J. A., Sarunic, M. V, Beg, M. F., Tybulewicz, V. L. J., Fisher, E. M. C., Lythgoe, M. F., & Ourselin, S. (2020). Substantially thinner internal granular layer and reduced molecular layer surface in the cerebellum of the Tc1 mouse model of Down Syndrome – a comprehensive morphometric analysis with active staining contrast-enhanced MRI. NeuroImage, 117271. https://doi.org/https://doi.org/10.1016/j.neuroimage.2020.117271
  - Ma, D., Cardoso, M. J., Zuluaga, M. A., Modat, M., Powell, N., Wiseman, F., Tybulewicz, V., Fisher, E., Lythgoe, M. F., & Ourselin, S. (2015). Grey Matter Sublayer Thickness Estimation in the Mouse Cerebellum. In Medical Image Computing and Computer Assisted Intervention 2015 (pp. 644–651). https://doi.org/10.1007/978-3-319-24574-4_77
