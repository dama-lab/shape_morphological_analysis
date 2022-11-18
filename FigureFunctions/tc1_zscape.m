function tc1_zscape()
    %% calculate zscapes
    
    % (May need to remove Cop), volume analysis doesn't make sense
    
    plot_opt = load_opts();
    v2struct(plot_opt);

    %%
    % add "fdr_bh.m" to system path
    % addpath(genpath('C:\Users\madam\Dropbox\Documents\SFU\Tools\my_code\matlab\'));
    % data_dir = 'C:\Users\madam\OneDrive\Data\Brain_MRI\UCL\Tc1_Cerebellum\08_stats\stats_20191020\remove_missegmented_subjects';
    % statDi        r = '/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/08_stats/stats_20191020/remove_missegmented_subjects';
    % statDir = '/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/08_stats/stats_20191020/';
    %    for all 28 subjects
    statDir = '/home/dma73/Data/Brain_MRI/UCL/Tc1_Cerebellum/08_stats/stats_20200413/';
    
    
    data_dir = statDir;

    normalization = 'GLM'; % 'Raw_All'; % "divide"; % 
% 
%     csv_filelist{1} = 'area2thicknessRatio_moleAreaGranThick.csv';
%     csv_filelist{2} = 'area2thicknessRatio_granAreaMoleThick.csv';
%     csv_filelist{3} = 'area2thicknessRatio_full.csv';
%     csv_filelist{4} = 'area2thicknessRatio_granular.csv';
%     csv_filelist{5} = 'area2thicknessRatio_molecular.csv';

%     csv_filelist{1} = 'thickness2AreaRatio_granThickMoleArea.csv';
%     csv_filelist{2} = 'thickness2AreaRatio_moleThickGranArea.csv';    
%     csv_filelist{3} = 'thickness2areaRatio_full.csv';
%     csv_filelist{4} = 'thickness2areaRatio_granular.csv';
%     csv_filelist{5} = 'thickness2areaRatio_molecular.csv';
    
%     
%     
    csv_filelist{1} = 'volume_full_raw.csv';
    csv_filelist{2} = 'volume_granular_parcellation.csv';
    csv_filelist{3} = 'volume_molecular_parcellation.csv';
    
    csv_filelist{4} = 'thickness_cortical_on_purkinje.csv';
    csv_filelist{5} = 'thickness_granular_WMD_on_purkinje.csv';
    csv_filelist{6} = 'thickness_molecular_CSFD_on_purkinje.csv';

    csv_filelist{7} = 'surface_area_full_Purkinje.csv';
    csv_filelist{8} = 'surface_area_granular_Purkinje.csv';
    csv_filelist{9} = 'surface_area_molecular_Purkinje.csv';
    %csv_filelist{10} = 'CS_surface_area';
    
    %% TIV (for accurately calculating)
%     tiv_file = fullfile(data_dir,'TIV_BV.csv');
% %     tiv_cell = readcell(tiv_file);
% %     tiv = cellfun(@str2double, tiv_cell(2:end,3));
%     tivArray = readmatrix(tiv_file);
%     tiv = tivArray(1:end,3);
    
    %% normalize/regress out cerebellar volume
    tivTable = readTIVcsv(TIVCsv);
    tiv_raw  = tivTable.TIV;
    
    %% Read Measure for sort from file (need to calculate beforehand)
    measure_for_sort_file = fullfile(data_dir,'avg_cortical_thickness_for_sort.csv');
    measure_for_sort_raw = csvread(measure_for_sort_file);
        
    %% surface area flag
    surface_area_flag = 0;
    if surface_area_flag == 1
        %% calculate full cortical surface area
        % volume
        volume_path{1} = fullfile(data_dir, csv_filelist{1}); % full volume
        volume_path{2} = fullfile(data_dir, csv_filelist{2}); % granuler volume
        volume_path{3} = fullfile(data_dir, csv_filelist{3}); % molecular volume
        
        % thickness
        thickness_path{1} = fullfile(data_dir, csv_filelist{4}); % full thickness
        thickness_path{2} = fullfile(data_dir, csv_filelist{5}); % /-WMD 
        thickness_path{3} = fullfile(data_dir, csv_filelist{6}); % /-CSFD
        
        % convert to the surface area
        surface_area_path{1} = fullfile(data_dir,'surface_area_full_Purkinje.csv');
        surface_area_path{2} = fullfile(data_dir,'surface_area_granuler_Purkinje.csv');
        surface_area_path{3} = fullfile(data_dir,'surface_area_molecular_Purkinje.csv');
        
        %%
        for pair_id = 1:3
            [surface_area,surface_area_cell] = calculate_surface_area(volume_path{pair_id},thickness_path{pair_id},surface_area_path{pair_id});
        end
        return
    end
          
    %% statistical comparison + plot zscape
    t_array = [];
    p_array = [];
    adj_p_array = [];
    refMeanArray= [];
    refStdArray = [];
    nMatrix = []; % sample size calculation for power analysis
    
    for id = 1:length(csv_filelist)
        csv_filename = csv_filelist{id};
        [~,filename,~] = fileparts(csv_filename);
        csv_path = fullfile(data_dir,csv_filename);
        csv_cell_raw = readcell(csv_path);

        %% Donot need to remove outlier -TG22 anymore
        remove_outlier = 0;
        if remove_outlier == 1
            csv_cell = csv_cell_raw([1:21,23:end],:);
            tiv = tiv_raw([1:21,23:end],:);
            measure_for_sort_raw = measure_for_sort_raw([1:21,23:end],:);
        end

        %% select subset
        expId = 1;
        subsets = 14:14;
        subsetId = 0;
        for subset = subsets
            subsetId = subsetId + 1;
            wt_idx = 1:14; % 14
            tg_idx = 15:28; % 14
            % permute through all possible combinations of subsample pick
            wtSubidxs = nchoosek(wt_idx, subset);
            tgSubidxs = nchoosek(tg_idx, subset);
            expNo = size(wtSubidxs,1)*size(tgSubidxs,1);
            %%
            for wtSub = 1:size(wtSubidxs,1)
                wtSubidx = wtSubidxs(wtSub,:);
                for tgSub = 1:size(tgSubidxs,1)
                    tgSubidx = tgSubidxs(wtSub,:);
                
                    fprintf('morph %d/%d; subset: %d/%d; exp %d/%d; ', ...
                        id, length(csv_filelist), subsetId, length(subsets), expId, expNo);
                    expId = expId + 1;

                    %% Alternative repeated random pick implementation
                    % fix random state. Ref: https://www.mathworks.com/help/matlab/math/random-integers.html
                    % rng(0,'twister');
                    % subset=11;
                    % wtSubidx = datasample(wt_idx, subset, 'Replace', false);
                    % tgSubidx = datasample(tg_idx, subset, 'Replace', false);

                    csv_cell = csv_cell_raw([1,2,wtSubidx+2,tgSubidx+2],:);
                    tiv = tiv_raw([wtSubidx,tgSubidx],:);
                    measure_for_sort = measure_for_sort_raw([wtSubidx,tgSubidx]);
                    
                    %% prepare variables
                    subject_list = csv_cell(3:end,2);
                    % convert subject list to string if it's in numerical format
                    subject_list = cellfun(@num2str,subject_list,'UniformOutput',false);
                    feature_list = csv_cell(1,3:end);
                    classification = csv_cell(3:end,1);
                    diag_classes = unique(classification);

                    feature = cell2mat(csv_cell(3:end,3:end));

                    %% investigate the distribution of the reference group 
                    refMeanArray = [refMeanArray, mean(feature(1:subset,:), 1)'];
                    refStdArray = [refStdArray, std(feature(1:subset,:), [], 1)'];
                    
                    %% Normalize TIV through GLM
                    % extract the control group only
                    idx_ctl = (classification=="wt");
                    tiv_ctl = tiv(idx_ctl==1);
                    %
                    feat_ctl= feature(idx_ctl==1,:);

                    zscape_group = classification;
                    % zscape_group = repmat("WT",1,21);

                    %% GLM model fitting
                    if normalization == "GLM"
                        [~, feature] = GLM(feature, tiv, zscape_group,'wt');
                    elseif normalization == "divide"
                        feature = feature./tiv;
                        fprintf('\n');
                    else
                        fprintf('\n');
                    end
            %         % has bug here
            %         mdl = fitlm(tiv_ctl,feat_ctl,'linear');
            %         feat_pred = predict(mdl,tiv);
            %         % calculate residual
            %         feature = feature - feat_pred;

                    %% sample size power analysis calculation
                    % Ref: https://www.mathworks.com/help/stats/sampsizepwr.html#namevaluepairarguments
                    powerAnalysis_flag = 1;
                    if powerAnalysis_flag == 1
                        %% confirm reference mean shall == 0, std != 1
                        ctrMean = mean(feature(idx_ctl,:),1);
                        ctrStd = std(feature(idx_ctl,:),0,1);
                        
                        %% zscore
                        featureZscore = (feature-ctrMean)./ctrStd;
                        
                        %% power analysis
                        % reference: mean=0, std=1
                        ctrMeanZ = mean(featureZscore(idx_ctl,:),1);
                        stdMeanZ = std(featureZscore(idx_ctl,:),0,1);
                        % test
                        testMeanZ = mean(featureZscore(idx_ctl~=1,:),1);
                        sampleNoTg = sum(idx_ctl~=1);
                        pwrDesire = 0.9; %0.70:0.05:0.95;
                        nArray = nan(size(testMean));
                        
                        for structId = 1:length(testMean)
                            %% Compute the power of the test % = 1
                            pwr = sampsizepwr('t2', [ctrMeanZ(structId), stdMeanZ(structId)], testMeanZ(structId),[],sampleNoTg);
                            %% calculate sample size for effect size of 0.95
                            nArray(structId) = sampsizepwr('t2', [ctrMeanZ(structId), stdMeanZ(structId)], testMeanZ(structId), pwrDesire, []);
                            %% calculate sample size of both group for effect size of 0.95
                            % [n1, n2] = sampsizepwr('t2', [ctrMeanZ(structId), stdMeanZ(structId)], testMeanZ(structId), pwrDesire, [], 'ratio', 1)
                        end
                        nMatrix = [nMatrix;nArray];
                    end
            
                    %% calculate mean full cortical thickness on purkinji (only once)       
                    calculate_sort_order = 0;
                    if calculate_sort_order == 1
                        mean_thickness = mean(feature,2);
                        % csvwrite(measure_for_sort_file,mean_thickness);
                        % starting from 2019a, use below instead
                        % measure_for_sort_file = fullfile(data_dir,'avg_cortical_thickness_for_sort.csv')
                        writematrix(mean_thickness,measure_for_sort_file);
                    end


                    %% unpaired t-test 
                    [h, p, t, ci] = deal([],[],[],[]); 
                    clear stats;
                    for f = 1:length(feature_list)
                        feature_wt = feature((classification == "wt"),f);
                        feature_tg = feature((classification == "tg"),f);
                        [h(f),p(f),ci(f,:),stats(f)] = ttest2(feature_tg, feature_wt);
                        t(f) = stats(f).tstat;
                    end
                    p_array = [p_array; p];
                    t_array = [t_array; t];

                    %% Multiple Comparison Correction with fDR = 0.05
                    fdr_correct=0;
                    if fdr_correct==1
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

                        adj_p_array = [adj_p_array; adj_p];
                    end

                    %% zscape
                    zscape_flag = 0;
                    if zscape_flag == 1
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
                        DisplaySpec.displayGroup = {'wt','tg'};
                        DisplaySpec.refDict = {'wt','wt';'tg','wt'};
                        % PlotSpec
                        DisplaySpec.plotSpec.Title = [];
                        DisplaySpec.plotSpec.zscape_modes = {'sorted','zscape'};
                        DisplaySpec.plotSpec.threshold_values = 0;
                        DisplaySpec.plotSpec.zscoreMinMax = [-6,6];
                        DisplaySpec.plotSpec.maxWinSize=[10,10,800,800];
                        DisplaySpec.plotSpec.cellWidth=30;
                        DisplaySpec.plotSpec.sortGapScale=3;
                        DisplaySpec.plotSpec.sortHightScale=0.5;


                        %% prepare Zscape plot visualization
                        % zscape_rootdir = fullfile(data_dir,'zscape',char(normalization));
                        zscape_rootdir = fullfile(data_dir,'zscape',sprintf('%s_%d_wt%d_tg% d', normalization, subset, wtSub, tgSub));
                        zscape_dir = fullfile(zscape_rootdir,filename);
                        if ~exist(zscape_dir,'dir'); mkdir(zscape_dir); end
                        
                        %% plot zscape
                        zscape_output = zscape_plot(DataInfo,DisplaySpec,zscape_dir);

                        %% move file up
                        movefile(fullfile(zscape_dir,'zscape_sorted_0.png'),fullfile(zscape_rootdir,[filename,'.png']));
                        rmdir(zscape_dir, 's');
                    end
                end
            end
        end
    end
    
    %% plot power analysis sample
    figure;
    boxplot(nMatrix');
    
    %% 
    titles = {'volume (full)', 'volume (granular)', 'volume (molecular)'; ...
             'thickness (full)', 'thickness (granular)', 'thickness (molecular)'; ...
             {'surface','area (full)'}, {'surface','area (granular)'}, {'surface','area (molecular)'}; ...
             };
    figure;
    for morphId = 1:size(nMatrix,1)
        subplot(3,3,morphId);
        bar(nMatrix(morphId,:));
        hold on;
        xlims = xlim;
        plot([xlims(1),xlims(2)],[15, 15],'k--')
        ylim([0,50]);
        xticks([]);
        title(titles{morphId});
        if ismember(morphId, [1,4,7])
            ylabel('sample size');
        end
        if ismember(morphId, [7,8,9])
            xlabel('Structures');
        end
    end
    
    
    %% plot distribution of mean/std
    figure;plot(refMeanArray(1,:)','.');
    
    
    %% HDR correct p_array
    fdr_q = 0.05;
    [h_adj, crit_p, adj_ci_cvrg, adj_p_array]=fdr_bh(p_array,fdr_q,'pdep','yes');

    %% Reshape (9 morph_metrics * '196' experiment * 15 structures )
    p_vol = reshape(p_array,length(csv_filelist), expNo, size(adj_p_array,2));
    adj_p_vol = reshape(adj_p_array,length(csv_filelist), expNo, size(adj_p_array,2));
    t_vol = reshape(t_array,length(csv_filelist), expNo, size(t_array,2));
    
    %% plot distribution of permuted p-value
    
    figure; 
    set(gcf,'Position',[0,0,2560,1220]);
    for morph = 1:size(adj_p_vol,1)
        subplot(3,3,morph);
        plot(squeeze(adj_p_vol(morph,:,:))','.');
        % violin(adj_p_vol(:,:,strut));
        % boxplot(squeeze(adj_p_vol(:,:,struct))');
        
        set(gca, 'XLim', [0,size(adj_p_vol,3)+1]);
        xlims = xlim;
        set(gca, 'XTick', xlims(1):xlims(2));
        
        hold on;
        plot(xlim,[0.05,0.05],'k--');
        ylim([0,1]);
    end
    
    %% plot distribution of permuted t-value
    figure; 
    set(gcf,'Position',[0,0,2560,1220]);
    for morph = 1:size(t_vol,1)
        subplot(3,3, morph);
        plot(squeeze(t_vol(morph,:,:))','.');
        % violin(t_vol(:,:,strut));
        % boxplot(t_vol(:,:,struct)');
        
        set(gca, 'XLim', [0,size(t_vol,3)+1]);
        xlims = xlim;
        set(gca, 'XTick', xlims(1):xlims(2));
        
        hold on;
        plot(xlim,[0,0],'k--');
    end

end
    
% function [feature,csv_cell] = extract_feature_from_csv(csv_path)
%     csv_cell = readCsvToCell(csv_path);
%     %% for removing_missegemnted_subjects
%     % feature = cellfun(@str2double,csv_cell(3:end,3:end));
%     %% for all subjects
%     feature = cellfun(@str2double,csv_cell(2:end,4:18));
% end
%     
% function [surface_area,csv_cell] = calculate_surface_area(full_volume_path,thickness_path,surface_area_path)
%         % full volume
%         [volume, csv_cell] = extract_feature_from_csv(full_volume_path);
%         % full thickness
%         thickness = extract_feature_from_csv(thickness_path);
%         % calculate the surface area
%         surface_area = volume ./ thickness;
%         % convert matrix to cell array
%         csv_cell(3:end,3:end) = num2cell(surface_area);
%         % save_surface_area_file
%         if exist('surface_area_path','var')
%             writetable(cell2table(csv_cell),surface_area_path,'WriteVariableNames',false);
%         end
%         
% end