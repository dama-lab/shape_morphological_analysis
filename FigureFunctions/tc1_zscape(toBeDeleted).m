function tc1_zscape()
    % add "fdr_bh.m" to system path
    % addpath(genpath('C:\Users\madam\Dropbox\Documents\SFU\Tools\my_code\matlab\'));
    data_dir = 'C:\Users\madam\OneDrive\Data\Brain_MRI\UCL\Tc1_Cerebellum\08_stats\stats_20191020\remove_missegmented_subjects';
    % statDir = '/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/08_stats/stats_20191020/';
    data_dir = statDir;

    normalization = 'Raw_All';
    
    csv_filelist{1} = 'thickness_cortical_on_purkinji.csv';
    csv_filelist{2} = 'thickness_granular_WMD_on_purkinje.csv';
    csv_filelist{3} = 'thickness_molecular_CSFD_on_purkinje.csv';
    csv_filelist{4} = 'volume_full_raw.csv';
    csv_filelist{5} = 'volume_granular_parcellation.csv';
    csv_filelist{6} = 'volume_molecular_parcellation.csv';
    csv_filelist{7} = 'surface_area_full_Purkinje.csv';
    csv_filelist{8} = 'surface_area_granuler_Purkinje.csv';
    csv_filelist{9} = 'surface_area_molecular_Purkinje.csv';
    %csv_filelist{10} = 'CS_surface_area';
    
    %% TIV (for accurately calculating)
    tiv_file = fullfile(data_dir,'TIV_BV.csv');
%     tiv_cell = readcell(tiv_file);
%     tiv = cellfun(@str2double, tiv_cell(2:end,3));
    tivArray = readmatrix(tiv_file);
    tiv = tivArray(1:end,3);
    
    %% Read Measure for sort from file (need to calculate beforehand)
    measure_for_sort_file = fullfile(data_dir,'avg_cortical_thickness_for_sort.csv');
    measure_for_sort = csvread(measure_for_sort_file);
        
    %% surface area flag
    surface_area_flag = 1;
    if surface_area_flag == 1
        %% calculate full cortical surface area
        % thickness
        thickness_path{1} = fullfile(data_dir, csv_filelist{1}); % full thickness
        thickness_path{2} = fullfile(data_dir, csv_filelist{2}); % /-WMD 
        thickness_path{3} = fullfile(data_dir, csv_filelist{3}); % /-CSFD
        % volume
        volume_path{1} = fullfile(data_dir, csv_filelist{4}); % full volume
        volume_path{2} = fullfile(data_dir, csv_filelist{5}); % granuler volume
        volume_path{3} = fullfile(data_dir, csv_filelist{6}); % molecular volume
        % convert to the surface area
        surface_area_path{1} = fullfile(data_dir,'surface_area_full_Purkinje.csv');
        surface_area_path{2} = fullfile(data_dir,'surface_area_granuler_Purkinje.csv');
        surface_area_path{3} = fullfile(data_dir,'surface_area_molecular_Purkinje.csv');
        
        for pair_id = 2:3
            [surface_area,surface_area_cell] = calculate_surface_area(volume_path{pair_id},thickness_path{pair_id},surface_area_path{pair_id});
        end
        return
    end
    
    %%
    for id = 1:length(csv_filelist)
        csv_filename = csv_filelist{id};
        [~,filename,~] = fileparts(csv_filename);
        csv_path = fullfile(data_dir,csv_filename);
        csv_cell = readcell(csv_path);
        
        %% Remove outlier -TG22
        remove_outlier = 1;
        if remove_outlier == 1
            csv_cell = csv_cell([1:21,23:end],:);
            tiv = tiv([1:21,23:end],:);
            measure_for_sort = measure_for_sort([1:21,23:end],:);
        end

        %% prepare Zscape plot visualization
        zscape_rootdir = fullfile(data_dir,'zscape/-TG22',char(normalization));
        zscape_dir = fullfile(zscape_rootdir,filename);
        if ~exist(zscape_dir,'dir'); mkdir(zscape_dir); end

        %% prepare variables
        subject_list = csv_cell(3:end,2);
        feature_list = csv_cell(1,3:end);
        classification = csv_cell(3:end,1);
        diag_classes = unique(classification);

        feature = cell2mat(csv_cell(3:end,3:end));
        
        %% TIV through GLM
        % extract the control group only
        idx_ctl = (classification=="WT");
        tiv_ctl = tiv(idx_ctl==1);
        %
        feat_ctl= feature(idx_ctl==1,:);
        
        zscape_group = classification;
        % zscape_group = repmat("WT",1,21);
        
        %% GLM model fitting
        if normalization == "GLM"
            [~, feature] = GLM(feature, tiv, zscape_group,'WT');
        elseif normalization == "divide"
            feature = feature./tiv;
        end
%         % has bug here
%         mdl = fitlm(tiv_ctl,feat_ctl,'linear');
%         feat_pred = predict(mdl,tiv);
%         % calculate residual
%         feature = feature - feat_pred;
        
        %% calculate mean full cortical thickness on purkinji (only once)       
        calculate_sort_order = 0;
        if calculate_sort_order == 1
            mean_thickness = mean(feature,2);
            csvwrite(measure_for_sort_file,mean_thickness);
            % starting from 2019a, use below instead
            % writematrix(mean_thickness,fullfile(data_dir,'avg_cortical_thickness_for_sort.csv'));
        end

        
        %% unpaired t-test with fDR = 0.05
        [h, p, ci] = deal([],[],[]); 
        clear stats;
        for f = 1:length(feature_list)
            feature_wt = feature((classification == "WT"),f);
            feature_tg = feature((classification == "TG"),f);
            [h(f),p(f),ci(f,:),stats(f)] = ttest2(feature_wt,feature_tg);
        end
        fdr_q = 0.05;
        [h_adj, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,fdr_q,'pdep','yes');
        % make adj_p 3-decimal float
        trunc3 = @(x) sprintf('%0.3f',x);
        adj_p_round = arrayfun(trunc3,adj_p,'UniformOutput',false);
        % add star after significant
        add_signif_star = @(x) sprintf('%s*',x);
        add_spaces = @(x) sprintf('%s  ',x);
        adj_p_round(adj_p < 0.05) = cellfun(add_signif_star, adj_p_round(adj_p < 0.05),'UniformOutput',false);
        adj_p_round(adj_p >= 0.05) = cellfun(add_spaces, adj_p_round(adj_p >= 0.05),'UniformOutput',false);


        %% Prepare Zscape input
        % Construct DataInfo
        DataInfo.dataMatrix = feature;
        DataInfo.feature = feature_list;
        DataInfo.subjectName = subject_list;
        DataInfo.classification = classification;
        DataInfo.group_metrics = adj_p_round;
        DataInfo.measure_for_sort.data = measure_for_sort;
        DataInfo.measure_for_sort.name = 'full cortical thickness';
        
        % Construct DisplaySpec
        DisplaySpec.displayMode = zeros(length(subject_list),1);
        DisplaySpec.displayGroup = {'WT','TG'};
        DisplaySpec.refDict = {'WT','WT';'TG','WT'};
        % PlotSpec
        DisplaySpec.plotSpec.Title = [];
        DisplaySpec.plotSpec.zscape_modes = {'sorted','zscape'};
        DisplaySpec.plotSpec.threshold_values = 0;
        DisplaySpec.plotSpec.zscoreMinMax = [-4,4];
        DisplaySpec.plotSpec.maxWinSize=[10,10,800,800];
        DisplaySpec.plotSpec.cellWidth=30;
        DisplaySpec.plotSpec.sortGapScale=3;
        DisplaySpec.plotSpec.sortHightScale=0.5;

        %% plot zscape
        zscape_output = zscape_plot(DataInfo,DisplaySpec,zscape_dir);
        
        %% move file up
        movefile(fullfile(zscape_dir,'zscape_sorted_0.png'),fullfile(zscape_rootdir,[filename,'.png']));
        rmdir(zscape_dir, 's');
    end
end
    
function [feature,csv_cell] = extract_feature_from_csv(csv_path)
    csv_cell = readCsvToCell(csv_path);
    feature = cellfun(@str2double,csv_cell(3:end,3:end));
end
    
function [surface_area,csv_cell] = calculate_surface_area(full_volume_path,thickness_path,surface_area_path)
        % full volume
        [volume, csv_cell] = extract_feature_from_csv(full_volume_path);
        % full thickness
        thickness = extract_feature_from_csv(thickness_path);
        % calculate the surface area
        surface_area = volume ./ thickness;
        % convert matrix to cell array
        csv_cell(3:end,3:end) = num2cell(surface_area);
        % save_surface_area_file
        if exist('surface_area_path','var')
            writetable(cell2table(csv_cell),surface_area_path,'WriteVariableNames',false);
        end
        
end