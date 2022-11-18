        thicknessFile = fullfile(thicknessDir, target);
        
        %%
        reslice_nii([rawVolNii,'.nii.gz'],[rawVolNii,'_reslice.nii.gz']);
        nii = load_nii([rawVolNii,'_reslice.nii.gz']);
        NiftiInfo = niftiinfo([rawVolNii,'_reslice.nii.gz']);
        NiftiInfo.Datatype = 'double';
        nii.img = single(thickness);
        save_nii(nii,[thicknessFile,'.nii.gz']);
        %%
        niftiwrite(thickness, thicknessFile);
        NiftiInfo = niftiinfo(thicknessFile);
        NiftiInfo.PixelDimensions = volNiftiInfo.PixelDimensions;
        NiftiInfo.TransformName = volNiftiInfo.TransformName;
        NiftiInfo.Transform = volNiftiInfo.Transform;
        NiftiInfo.raw = NiftiInfo.raw;
        NiftiInfo.Datatype = double;
        NiftiInfo.BitsPerPixel = 64;
        NiftiInfo.Version = 'NIfTI2';
        %%
        niftiwrite(thickness, thicknessFile,NiftiInfo,'Compressed',true); 
        
        
        %%
        %load_nifti_hdr([thicknessFile,'.nii.gz']);
        volNiftiInfo.raw.sizeof_hdr = 559903;

        
        reslice_nii([segVolNii,'.nii.gz'],[segVolNii,'_reslice.nii.gz'])
        load_nii([segVolNii,'_reslice.nii.gz'])
        
        