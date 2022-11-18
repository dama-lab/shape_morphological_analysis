# Laplacian-based Thickness Estimation Function List:

## - Main function: `thickness_map_3d`

> Calculate laplacian-based structural (e.g. cortical) thickness based on segmented labels files stored in nifti format

## - Internal function: `LaplacianThickness`

> Calculate the laplacian-based structural (e.g. cortical) thickness based on segmented labels

### - Internal Step 1: `LaplaceField`

> Calculate laplace gradient from the segmented labels

### - Internal step 2: `laplace_thickness`

> Calculate layer thickness map from the laplace gradient field



