#/bin/bash

## label propagation
source ~/Code/multi-atlas-segmentation/MASHelperFunctions.sh > /dev/null 2>&1

# target_dir="/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/00_cerebellum_N4_with_mask"
# target_list="/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/TargetList/tc1_269456-ob_c"
# target_list="/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/TargetList/TargetList.txt"
# targetmask_dir="/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/01_mask_groupwise_f3d2/mask_manual_corrected"

target_dir="/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/00_groupwise_average"
target_list="/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/TargetList/average_nonrigid_it_10_nonan"
targetmask_dir="/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/00_groupwise_average/mask_manual"

atlas_dir="/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/Atlas/AMBMC-cere-larger-view"
result_dir="/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/04_cotical_layers"

## Prepare input dataset
cere_extract_dir="/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/00_cerebellar_extract"
targetmask_dil1_dir="/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/01_mask_groupwise_f3d2/mask_manual_corrected_dil1"

mkdir -p $cere_extract_dir
mkdir -p $targetmask_dil1_dir

for target_id in $(cat $target_list); do
	# extract cerebellar
	seg_maths $target_dir/$target_id -mul $targetmask_dir/$target_id $cere_extract_dir/$target_id.nii.gz
	# dilate input mask
	seg_maths $targetmask_dir/$target_id -dil 1 $targetmask_dil1_dir/$target_id.nii.gz
done

# label propagation
mas_parcellation_batch -T $cere_extract_dir -t $target_list -A $atlas_dir -r $result_dir -M $targetmask_dil1_dir -e local

