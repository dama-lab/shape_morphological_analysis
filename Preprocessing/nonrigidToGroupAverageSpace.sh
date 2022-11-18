#/bin/bash

## Idea brainstorming
## [Shape]: 
# 1. non-rigidly propagate groupwise average labels to native space
# [Shape]:
#   2. rigidly propagate back to groupwise space
# [Thickness]:
#   <Option 1 good idea, been used> 
#       2. calculate thickness (both gran and mol) in the native space
#       3. non-rigidly propagate thickness back to groupwise space
#   <Option 2 bad idea>
#       2. Convert sublayer label to surface
#       3. Calculate surface-based thickness (closest points or registration-based)
#       4. rigidly+KNN or Non-rigidly register surface back to groupwise space
##

rootDir='/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum'
groupAvgDir=$rootDir'/00_groupwise_average'
transformationDir=$groupAvgDir'/transformations/'
N4Dir=$rootDir/'00_cerebellum_N4_with_mask'
parcellationDir=$rootDir/'02_parcellation'
purkinjeDir=$rootDir/'07_Purkinje/round2_3_remove_sulci/4_remove_sulci'
tissueLabelDir=$rootDir/'03_tissue_label'
thicknessDir=$rootDir/'06_thickness/groupMeanBackpropagate'

resampledRootDir=$groupAvgDir/'resample_nrr_10'
# resampledN4Dir=$resampledRootDir'/N4'
# resampledParcDir=$resampledRootDir'/parcellation/'
# resampledPurkinjeDir=$resampledRootDir'/purkinje/'

maskDir=$groupAvgDir/mask_manual

# targetList=$rootDir'/TargetList/TargetList.txt'
# targetList=$rootDir'/TargetList/TargetList-22.txt'
targetList=$rootDir'/TargetList/TargetList_All.txt'

##
seg_files="sublayerParcellationQC" #"thicknessCortical thicknessGranular thicknessMolecular" #"sublayerParcellation" # "sublayerParcellationQC" #"parcellation" #"N4" # Wm1Gm2RemoveSulci" # "sublayerParcellation" # "displacement" # "thicknessMolecular" # "thicknessCortical" # "granularWM" # "GM_mask_remove_sulci" # "granualar" "N4 parcellation purkinje"
labelFusionFiles= # "GM_mask_remove_sulci" # "purkinje"

groupAvgName="average_nonrigid_it_10"
groupAvg=$groupAvgDir/$groupAvgName"_nonan.nii.gz"
groupAvgMask="$maskDir/${groupAvgName}_nonan.nii.gz"
groupAvgGmRemoveSulci="$$maskDir/GM_mask_remove_sulci"
groupAvgPurkinje=$groupAvgDir/'segmentations/gmLayerSegKmean_purkinje_lcc_manual.nii.gz'
groupAvgGranularWM=$groupAvgDir/'segmentations/gmLayerSegKmean_granularPlusWMManual.nii.gz'
groupAvgGranMoleParc=$groupAvgDir/'segmentations/gmLayerSegKmean_GranMoleParc.nii.gz'
groupAvgGranMoleParcQC=$groupAvgDir/'segmentations/gmLayerSegKmean_GranMoleParc_QC.nii.gz'
groupAvgWm1Gm2RemoveSulci=$groupAvgDir/'segmentations/average_nonrigid_it_10_WM1_GM2_remove_sulci_manual.nii.gz'
groupAvgParcellation=$resampledRootDir/'average_nonrigid_it_10_lobuleParcellationRecolor.nii.gz'

## Setup flag
labelPropagationFlag=0
labelBackPropagation=1
labelMergeFlag=0
labelFusionFlag=0
qc_flag=0

## source file for quickcheck generator
if [[ $qc_flag -eq 1 ]]; then
    source $HOME/Code/multi-atlas-segmentation/MASHelperFunctions.sh > /dev/null 2>&1
fi

##
for seg_file in $seg_files; do
    interType=0
    ## determine variables
    if [[ $seg_file == "N4" ]]; then
        # resample individual images
        floatDir=$N4Dir
        # for back label propagation
        groupSegFile=$groupAvg
        interType=3
    elif [[ $seg_file == "parcellation" ]]; then
        # label propagation (mapping and NN resampling)
        floatDir=$parcellationDir
        groupSegFile=$groupAvgParcellation
    elif [[ $seg_files == "purkinje" ]]; then
        # propagate purkinje layer
        floatDir=$purkinjeDir
        # for back label propagation
        groupSegFile=$groupAvgPurkinje
    #elif [[ $seg_files == "GM_mask_remove_sulci" ]]; then
        # for label propagation: put below
        # for back label propagation: not yet needed
    elif [[ $seg_file == "granularWM" ]]; then
        # for back label propagation
        groupSegFile=$groupAvgGranularWM
    elif [[ $seg_file == "thicknessCortical" ]]; then
        floatDir=$thicknessDir/'thickness_map'
    elif [[ $seg_file == "thicknessGranular" ]]; then
        # floatDir=$thicknessDir/'thickness_granular'
        # defined as WMD_ratio
        # floatDir=$thicknessDir/'WMD'
        floatDir=$thicknessDir/'thickness_WMD'
    elif [[ $seg_file == "thicknessMolecular" ]]; then
        # floatDir=$thicknessDir/'thickness_molecular'
        # defined as WMD_ratio
        # floatDir=$thicknessDir/'CSFD'
        floatDir=$thicknessDir/'thickness_CSFD'
    elif [[ $seg_file == "sublayerParcellation" ]]; then
        groupSegFile=$groupAvgGranMoleParc
    elif [[ $seg_file == "sublayerParcellationQC" ]]; then
        groupSegFile=$groupAvgGranMoleParcQC
    elif [[ $seg_file == "Wm1Gm2RemoveSulci" ]]; then
        groupSegFile=$groupAvgWm1Gm2RemoveSulci
    fi

    # Initialize variables
    atlas_no=0 # for Label Merge
    seg_file_1=""
    seg_file_n=""
    for target in $(cat $targetList); do
        ############ transform/resample parcellations to the groupsise space ###############
        if [[ $labelPropagationFlag -eq 1 ]]; then
            echo ">>>>>> propagating <$seg_file> for $target ..."
            nrrcpp="$transformationDir/nrr_10/nrr_cpp_${target}_it10.nii.gz"
            affmat="$transformationDir/aff_10/aff_mat_${target}_it10.txt"
            labelPropDir=$resampledRootDir/$seg_file
            mkdir -p $labelPropDir
            if [[ $seg_files == "displacement" ]]; then
                # convert control points to displacement field (still contain the initial affine)
                reg_transform -ref $groupAvg -disp $nrrcpp $labelPropDir/$target.nii.gz
                ## Attempt to remove affine transformation first (doesn't seem to work)
                # mkdir -p $labelPropDir/invAff
                # # Step 1: inverse the affine transformation
                # reg_transform -invAff $affmat $labelPropDir/invAff/$target.txt
                # # Step 2: combine the control points with the inverse affine
                # reg_transform -comp $nrrcpp $affmat $labelPropDir/$target.nii.gz -ref $groupAvg
            elif [[ $seg_files == "GM_mask_remove_sulci" ]]; then
                # propagate GM_remove_sulci
                floatDir=$tissueLabelDir/$seg_files/$target
                reg_resample -inter $interType  -voff\
                    -flo $floatDir/GM_with_sulci.nii.gz -ref $groupAvg -cpp $nrrcpp \
                    -res $labelPropDir/$target.nii.gz
            else
                reg_resample -inter $interType  -voff\
                    -flo $floatDir/$target.nii.gz -ref $groupAvg -cpp $nrrcpp \
                    -res $labelPropDir/$target.nii.gz
            fi
            # Quickcheck
            if [[ qc_flag -eq 1 ]]; then
                mkdir -p $labelPropDir/QuickCheck
                mas_quickcheck $groupAvg $labelPropDir/$target.nii.gz $labelPropDir/QuickCheck $target 
            fi
        fi

        ############ back-transform/resample parcellations to the groupsise space ###############
        if [[ $labelBackPropagation == 1 ]]; then
            echo "back-propagating <$seg_files> for $target ..."
            backpropDir=$resampledRootDir/backprop_$seg_file
            mkdir -p $backpropDir
            nrrcppbackward="$transformationDir/nrr_10/nrr_cpp_${target}_it10_backward.nii.gz"
            reg_resample -inter $interType  -voff\
                -flo $groupSegFile -ref $N4Dir/$target.nii.gz -cpp $nrrcppbackward \
                -res $backpropDir/$target.nii.gz
            if [[ qc_flag -eq 1 ]]; then
                mkdir -p $backpropDir/QuickCheck
                mas_quickcheck $N4Dir/$target.nii.gz $backpropDir/$target.nii.gz $backpropDir/QuickCheck $target
            fi
        fi

        ############# Preparing Label Merge ###############
        if [[ $labelMergeFlag -eq 1 ]]; then
            # prepare command line
            if [[ $atlas_no -eq 0 ]]; then
                seg_file_1="$resampledRootDir/$seg_file/$target.nii.gz"
            else
                seg_file_n="$seg_file_n $resampledRootDir/$seg_file/$target.nii.gz"
            fi
            let atlas_no+=1
        fi
    done
    let atlas_no-=1

    ############# Label Merge ###############
    # merge 4D label to prepare for label fusion:
    if [[ $labelMergeFlag -eq 1 ]]; then
        if [[ $seg_files == "displacement" ]]; then
            mergeDim=5
        else
            mergeDim=4
        fi
        echo " >>>>>> merging $seg_file into ${mergeDim}D volume ..."
        merged_4d_file=$resampledRootDir/$seg_file/4D.nii.gz
        seg_maths $seg_file_1 -merge $atlas_no $mergeDim $seg_file_n $merged_4d_file
    fi
done

############ Label Fusion #############
if [[ $labelFusionFlag -eq 1 ]]; then
    for labelFusionFile in $labelFusionFiles; do
        echo "     label fusion <$labelFusionFile> ....."
        groupAvgSegFile="$resampledRootDir/${groupAvgName}_${labelFusionFile}.nii.gz"
        seg_LabFusion -unc -v 1 -STEPS 5 8 $groupAvg $resampledRootDir/$seg_file/4D.nii.gz \
            -in "$resampledRootDir/$labelFusionFile/4D.nii.gz" \
            -out "$groupAvgSegFile" -mask "$groupAvgMask" #_mask_manual.nii.gz"
        # crop only regions within the mask
        seg_maths "$groupAvgSegFile" -mul "$groupAvgMask" "$groupAvgSegFile"
        # QuickCheck
        if [[ qc_flag -eq 1 ]]; then
            mkdir -p "$resampledRootDir"/QuickCheck
            mas_quickcheck "$groupAvg" "$groupAvgSegFile" "$resampledRootDir"/QuickCheck $groupAvgName
        fi
    done
fi



## Original stepwise setup
# ############ transform/resample parcellations to the groupsise space ###############
# if [[ $labelPropagationFlag -eq 1 ]]; then
#     for target in $(cat $targetList); do
#         echo "resample $target ..."
#         nrrcpp="$transformationDir/nrr_10/nrr_cpp_${target}_it10.nii.gz"
#         nrrcppbackward="$transformationDir/nrr_10/nrr_cpp_${target}_it10_backward.nii.gz"
        
#         for seg_file in $seg_files; do
#             echo "        raw $seg_file ..."
#             mkdir -p $resampledRootDir/$seg_file
#             if [[ $seg_file == "N4" ]]; then
#                 # resample individual images
#                 floatDir=$N4Dir
#             elif [[ $seg_file == "parcellation" ]]; then
#                 # label propagation (mapping and NN resampling)
#                 floatDir=$parcellationDir
#             elif [[ $seg_files == "purkinje" ]]; then
#                 # propagate purkinje layer
#                 floatDir=$purkinjeDir
#             fi
#             reg_resample \
#                 -flo $floatDir/$target.nii.gz -ref $groupAvg -cpp $nrrcpp \
#                 -inter 3 -res $resampledRootDir/$seg_file/$target.nii.gz -voff
#         done
#     done
# fi

# ############# Label Merge ###############
# if [[ $labelMergeFlag -eq 1 ]]; then
#     atlas_no=0
#     for seg_file in $seg_files; do
#         # prepare command line
#         if [[ $seg_file == "N4" ]]; then
#             resampledDir=$resampledN4Dir
#         elif [[ $seg_file == "parcellation" ]]; then
#             resampledDir=$resampledParcDir
#         elif [[ $seg_files == "purkinje" ]]; then
#             resampledDir=$resampledPurkinjeDir
#         fi


#         for target in $(cat $targetList); do
#             if [[ $atlas_no -eq 0 ]]; then
#                 N4_1="$resampledN4Dir/$target.nii.gz"
#                 parcellation_1="$resampledDir/$target.nii.gz"
#             else
#                 N4_n="$N4_n $resampledN4Dir/$target.nii.gz"
#                 parcellation_n="$parcellation_n $resampledDir/$target.nii.gz"
#             fi
#             let atlas_no+=1
#         done
#         let atlas_no-=1


#         # prepare 4D label for label fusion:
    
#         echo "merging $seg_file into 4D volume ..."
#         merged_4d_file=$resampledRootDir/$seg_file/4D.nii.gz
#         merge_cmd="seg_maths \$${seg_file}_1 -merge $atlas_no 4 \$${seg_file}_n $merged_4d_file"
#         eval $merge_cmd
#     done
# fi

# ############ Label Fusion #############
# if [[ $labelFusionFlag -eq 1 ]]; then
#     seg_LabFusion -unc -v 1 -in "$resampledRootDir/$labelFusionFile/4D.nii.gz" -out "$groupAvg$labelFusionFile" -mask "$groupAvgMask" #_mask_manual.nii.gz"
#     # crop only regions within the mask
#     seg_maths "$groupAvgParcellation" -mul "$groupAvgMask" "$groupAvgParcellation"
# fi


