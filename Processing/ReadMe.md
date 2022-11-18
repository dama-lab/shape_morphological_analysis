# Batch processing script

## `nonrigidToGroupAverageSpace.sh`

> bash script to propagate labels/thickness nifti files between groupwise average space to individual space

This script resamples the deformation field generated from the NiftyReg-based [GroupWise Registration (GWR)](https://github.com/dama-lab/multi-atlas-segmentation/tree/master/GWR) pipeline available here: https://github.com/dama-lab/multi-atlas-segmentation/tree/master/GWR

## `cortical_thickness_individual.m` / `cortical_sublayer_thickness.m`

> Matlab script to batch calculate cortical thickness