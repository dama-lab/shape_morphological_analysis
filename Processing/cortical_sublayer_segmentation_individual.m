function cortical_sublayer_segmentation_individual()
    
    %% define target list
    data_dir = 'C:\Users\madam\OneDrive\Data\MRI\Downs-Tc1-Ben-Sinclair\cerebellum\shrink_fov_20150108\20150304\07_Purkinje\5_remove_patch';
    dirlist = dir(data_dir);
    targetlist = dirlist(3:end);
    
    parfor target_no = 1:length(targetlist)
        
        target_id = targetlist(target_no).name; % 'tc1_269455-ob_c.nii.gz';

        %% test alternative methods (Laplacian ~ thickness)
        addpath(genpath('C:\Users\madam\Dropbox\Documents\SFU\Tools\my_code\layer_thickness\cortical_thickness\'));
        addpath(genpath('C:\Users\madam\OneDrive\Study\PhD\Codes\MATLAB'));

        %% load Purkinji layer
        dir_purkinji = 'C:\Users\madam\OneDrive\Data\MRI\Downs-Tc1-Ben-Sinclair\cerebellum\shrink_fov_20150108\20150304\07_Purkinje\round2_3_remove_sulci\4_remove_sulci';
        purkinji_nii = load_untouch_nii(fullfile(dir_purkinji,target_id));
        purkinji = purkinji_nii.img;

        nii = purkinji_nii;

        %% %%%%%%%%%%%% WM to purkinji
        %% load WM
        dir_WM = 'C:\Users\madam\OneDrive\Data\MRI\Downs-Tc1-Ben-Sinclair\cerebellum\shrink_fov_20150108\20150304\03_tissue_label\WM_mask';
        WM_nii = load_untouch_nii(fullfile(dir_WM,target_id));
        WM = WM_nii.img;

        %% load GM (with_sulci)
        dir_GM_with_sulci = 'C:\Users\madam\OneDrive\Data\MRI\Downs-Tc1-Ben-Sinclair\cerebellum\shrink_fov_20150108\20150304\03_tissue_label\GM1_WM2_Boundary3_manual_group_f3d2';
        GM_with_sulci_nii = load_untouch_nii(fullfile(dir_GM_with_sulci,target_id));
        GM_with_sulci = GM_with_sulci_nii.img == 2;

        %% prepare Laplacien
        granular_layer = WM + 2*purkinji;
        %% Laplacian
        [p_gran,cortical] = laplace_equation_3d(granular_layer,2,0,1);
        % granular layer
        p_gra = imfill((p_gran>0.04)*1.* GM_with_sulci);


        %% calculate molacular layer from granullar layer
        se = strel('sphere',1);
        p_gra_dil = ~imdilate(p_gra,se);
        p_mol = (GM_with_sulci .* ~purkinji) .* ~imfill(imdilate(p_gra,se) .* ~purkinji);

    %     %% granular layer dilation - purkinji layer - WM - GM_with_sulci
    %     p_gra_new = p_gra;
    %     for i = 1:2
    %         % Ref 1: https://www.mathworks.com/help/images/ref/imdilate.html#d117e141904
    %         % Ref 2: https://www.mathworks.com/help/images/ref/strel.html#d117e267529
    %         se = strel('sphere',1);
    %         p_gra_dil = imdilate(p_gra_new,se);
    %         p_gra_new = p_gra_dil .* (GM_with_sulci) .* (~purkinji);
    %     end
    %     p_gra_new = imfill(p_gra_new);

        %% %%%%%%%%%%% pial to purkinji
    %     %% load GM
    %     dir_GM = 'C:\Users\madam\OneDrive\Data\MRI\Downs-Tc1-Ben-Sinclair\cerebellum\shrink_fov_20150108\20150304\03_tissue_label\GM1_WM2_removing_sulci';
    %     GM_nii = load_untouch_nii(fullfile(dir_GM,target_id));
    %     pial = GM_nii.img == 0;
    %     %% prepare Laplacian
    %     molecular_layer = pial + 2*purkinji;
    %     %% Laplacian
    %     [p_mole,cortical] = laplace_equation_3d(molecular_layer,2,0,1);
    %     % molacular layer
    %     p_mol = (p_mole>0.04)*1.* GM_with_sulci;

        %% parcellation
        dir_parcellation = 'C:\Users\madam\OneDrive\Data\MRI\Downs-Tc1-Ben-Sinclair\cerebellum\shrink_fov_20150108\20150304\02_parcellation';
        parcellate_nii = load_untouch_nii(fullfile(dir_parcellation,target_id));
        parcellate = double(parcellate_nii.img);
        parcellate_gra = parcellate .* (p_gra .* (parcellate>1));
        parcellate_mol = parcellate .* (p_mol .* (parcellate>1));

        %% save nii
        output_dir = 'C:\Users\madam\OneDrive\Data\MRI\Downs-Tc1-Ben-Sinclair\cerebellum\shrink_fov_20150108\20150304\09_sublayer_pacellation\';
        type_list = {'granular_layer','molecular_layer','granular_parcellation','molecular_parcellation'};
        for type_id = 1:length(type_list)
            subdir = fullfile(output_dir,type_list{type_id});
            if ~exist(subdir,'dir'); mkdir(subdir); end
        end
        % granular layer
        save_nii_with_hdr(nii,p_gra,fullfile(output_dir,type_list{1},target_id));
        save_nii_with_hdr(nii,p_mol,fullfile(output_dir,type_list{2},target_id));
        save_nii_with_hdr(nii,parcellate_gra,fullfile(output_dir,type_list{3},target_id));
        save_nii_with_hdr(nii,parcellate_mol,fullfile(output_dir,type_list{4},target_id));

        %% convex % not working yet
    %     %% use find to get the x,y,z, index
    %     [y,x,z] = ind2sub(size(nii.img),find(nii.img));
    %     
    %     %% use convhull to find the convexed mask
    %     k = convhull(x,y,z);
    %     
    %     %% convert convhull back to binary mask: poly2mask (ind2vec)
    %     [circXY(:,1),circXY(:,2)] = pol2cart(linspace(0,2*pi,50)', 1);
    %     sqXY = [-1 -1;1 -1;1 1;-1 1; -1 -1];
    %     C = {[sqXY*5 ones(5,1)] % Start with a small square
    %     [circXY*40 ones(50,1)*30] % Blend to a large circle
    %     [sqXY*20 ones(5,1)*65] % Blend to a large square
    %     [circXY*10 ones(50,1)*99]}; % Blend to a small circle
    %     BW = blendedPolymask(C,x,y,z);
    end

end

function save_nii_with_hdr(nii,p,output_path)
    % nii: input nii structure
    % p:   3d volume to be saved
    nii.hdr.dime.datatype=64;
    nii.hdr.dime.bitpix=64;
    nii.img=p;
    save_untouch_nii(nii,output_path);
end