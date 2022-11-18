function fig = vol_thickness_correlation_all_transpose()
%     fig = vol_thickness_correlation_all(...
%     thicknessFullCsv, surfaceAreaFullCsv, ...
%     thicknessGranCsv, surfaceAreaGranCsv, ...
%     thicknessMoleCsv, surfaceAreaMoleCsv, ...
%     normTIVflag, TIVCsv, figTitle, figFname)

%% load variables
plot_opt = load_opts();
v2struct(plot_opt);

% figFname = 'All_cortical_surfaceArea_vs_thickness_Normalize_TIV_all_28_transpose';
%figFname = 'All_cortical_surfaceArea_vs_thickness_Normalize_TIV_transpose_9-15';
normTIVflag = 1;
normalizeType = {'raw','GLM','Divide'};
%% plot correlation between volume and thickness
% by default normalize to TIV
if ~exist('normTIVflag','var'); normTIVflag = 1; end

%%%%%%%%%%%%%%%%%%%%
% for old compatibility with cell name in old csv (WT/TG)
color.WT = 'b';
color.TG = 'r';
% Text location of statistical resutls
% R-square
r2_positions.WT = [0.1,0.9];
r2_positions.TG = [0.6,0.9];
% P-value of 
p_positions.WT = [0.1,0.1];
p_positions.TG = [0.6,0.1];

%%%%%%%%%%%%%%%%%%%%
% new name (wt/tg)
color.wt = 'b';
color.tg = 'r';
% Text location of statistical resutls
% R-square
r2_positions.wt = [0.1,0.9];
r2_positions.tg = [0.6,0.9];
% P-value of 
p_positions.wt = [0.1,0.1];
p_positions.tg = [0.6,0.1];

%% load structural surface area
surfarea = [];
thickness = [];

% load structural surface area
[surfarea.full,targetList,structList,structLabels,subjectGroups] = readMeasureCsv(surfaceAreaFullCsv);
[surfarea.mole,~,~,~,~] = readMeasureCsv(surfaceAreaMoleCsv);
[surfarea.gran,~,~,~,~] = readMeasureCsv(surfaceAreaGranCsv);

% load structural thickness
[thickness.full,~,~,~] = readMeasureCsv(thicknessFullCsv);
[thickness.mole,~,~,~] = readMeasureCsv(thicknessMoleCsv);
[thickness.gran,~,~,~] = readMeasureCsv(thicknessGranCsv);

%% normalize TIV if normTIVflag == 1 (GLM) or 2 (divide)
layers = {'full','mole','gran'};
for layerId = 1:length(layers)
    if normTIVflag ==1
        TIVtable = readTIVcsv(TIVCsv);
        % Choose to use cerebellar as normalization effect
        TIV = TIVtable.Cerebellum;
        layer = layers{layerId};
        fprintf('\nnormalizing TIV (GLM) for surface area of: %s ...', layer);
        surfarea.(layer) = GLM(surfarea.(layer), TIV, subjectGroups,'wt');
        fprintf('\nnormalizing TIV (GLM) for thickness of:    %s ...', layer);
        thickness.(layer) = GLM(thickness.(layer), TIV, subjectGroups,'wt');
    elseif normTIVflag == 2
        TIV = readTIVcsv(TIVCsv);
        layer = layers{layerId};
        fprintf('\nnormalizing TIV (divide) for: %s ...', layer);
        surfarea.(layer) = surfarea.(layer)./(TIV.^(1/2));
        thickness.(layer) = thickness.(layer)./(TIV.^(1/3));
    end
end


%% plot correlation for all structures, color-coding the group (red: WT; blue: TG)
groups = unique(subjectGroups);
groupNo = length(groups);
structNo = size(surfarea.full,2);

fig = figure('Visible','on');
corrTypes = 5;
rows = structNo; %7; %8; % structNo;
columns = corrTypes;
p_corrs = nan(rows,columns,groupNo);

colTitles = {{'Full Cortical',''},{'Molecular Layer',''},{'Granular Layer',''}, ...
         {' Molecular Area  ';'vs. Granular Thickness'},{' Molecular Thickness  ';'vs. Granular Area'}};


t = tiledlayout(rows,columns,'TileSpacing','none','Padding','compact');
axId = 0;

LabelList = 1:structNo; % 1:structNo %9:15 % 1:8 % structNo = 15
for s = LabelList
    for c=1:corrTypes
        axes(s,c) = nexttile;
        sAreaBoth = [];
        thickBoth = [];
        for g = 1:groupNo
            % determine points
            group = groups{g};
            subjId = strcmp(subjectGroups, group);

            %% determine which two are correlated
            if (c==1) % row 1: full cortical gray matter
                sArea = surfarea.full(subjId,s);
                thick = thickness.full(subjId,s);
            elseif (c==2) % row 2: molecular layer
                sArea = surfarea.mole(subjId,s);
                thick = thickness.mole(subjId,s);
            elseif (c==3) % row 3: granular layer
                sArea = surfarea.gran(subjId,s);
                thick = thickness.gran(subjId,s);
            elseif (c==4) % row 4: molecular surfacearea vs. granular thickness 
                sArea = surfarea.mole(subjId,s);
                thick = thickness.gran(subjId,s);
            elseif (c==5) % row 5: granular surfacearea vs. molecular thickness
                sArea = surfarea.gran(subjId,s);
                thick = thickness.mole(subjId,s);
            end
            
            %% prepare for the correlation analysis
            sAreaBoth = [sAreaBoth;sArea];
            thickBoth = [thickBoth;thick];
            %% scatter plot
            scatPlt.(group) = plot(thick, sArea, '.', 'MarkerSize', 10, 'Color', color.(group));
            hold on;
            if normTIVflag == 1
                xlim([-4.9, 4.9]);
                ylim([-4.9, 4.9]);
            end
            
            % remove xticklabels / y ticklabels for the facets
            if c > 0
                yticklabels(gca,[]);
            end
            if s < structNo
                xticklabels(gca,[]);
            end
            % Add title
            if axId < columns
                title(colTitles{c},'Interpreter','latex');
            end
            % Add ylabel
            if ~mod(axId, columns)
                ylabel({sprintf('\\textbf{%s}',structList{s});'Surface Area'},'Interpreter','latex');
            end
            % Add ylabel
            if (axId >= (rows-1)*columns)
                xlabel('Thickness','Interpreter','latex');
            end
        end
        
        %% Linear fit for subjects from both groups
        coefs = polyfit(thickBoth,sAreaBoth,1);
        slope = coefs(1);
        intercept = coefs(2);

        % plot linear fit
        % line_plot = lsline; % (urgly)
        if normTIVflag == 1
            x = [-4, 4];
        else
            x = [min(thick),max(thick)];
        end
        y = slope*x + intercept;
        plot(x,y,'--','LineWidth',1); %,'Color', color.(group));
        hold on;

        % compute significance (p) of correlation of coefficient
        [R,P,RL,RU] = corrcoef(thick,sArea);
        % the off-diagnal entries of P is the p-vallue for the variable pair 
        p_corrs(s,c,g) = P(2,1); % = P(1,2);

        % Calculate R^2
        sAreaLinpred = polyval(coefs, thick);
        [b,bint,r,rint,stats] = regress(thick,[ones(size(sArea)),sArea]);
        r2 = stats(1);

        % display R^2
        r2_position = r2_positions.(group);
        r2_text = text(r2_position(1), r2_position(2), ...
            sprintf('$R^2$=%0.2f',r2), 'Color', color.(group), ...
            'FontSize', 10,'Units','normalized','Visible','on','Interpreter','latex');
   
        %% Add legend to the second last plot
        if axId == (length(LabelList * corrTypes) - 1)
            % corr_legend = legend([scatPlt.WT,scatPlt.TG],{'WT','TG'},'Location','southoutside','Orientation','horizontal','Interpreter','latex');
            corr_legend = legend([scatPlt.wt,scatPlt.tg],{'wt','tg'},'Location','southoutside','Orientation','horizontal','Interpreter','latex');
        end
         
        axId = axId + 1;
        
    end
end


%% update display p-value for exch axe
% Multiple Comparison Correction with False Discovery Rate (FDR) =0.05
fdr_q = 0.1;
[adj_h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_corrs,fdr_q,'pdep','yes');

%%
for s = 1:8 % 9:15 % 1:8 % structNo = 15
    for c=1:corrTypes
        for g = 1:groupNo
            group = groups{g};
            p_corr = adj_p(s,c,g);
            % bold P if <=0.05 (indicating significance)
            if p_corr <=0.05
                p_text = sprintf('\\textbf{$p$=%0.3f}',p_corr);
            else % Don't display p if not significant
                p_text = ''; sprintf('$p$=%0.3f',p_corr);
            end
            p_position = p_positions.(group);
            p_text = text(axes(s,c), p_position(1), p_position(2), ...
                p_text, 'color', color.(group),'FontSize', 10, ...
                'Units','normalized','Visible','on','Interpreter','latex');
        end
    end
end

%% Add Overall figure title
figTitle = "Cortical surface area v.s. thickness"; % (cont'd)";
if exist('figTitle','var')
    sgtitle(sprintf('\\textbf{%s}',figTitle),'Interpreter','latex');
end

%% Adjust Figure size
% % for first 8 structures
fig.Position = [0, 0, 1000, 1300];
% % for first 7 structures
fig.Position = [0, 0, 1000, 1100];
% % for all 15 structure
fig.Position = [0, 0, 1000, 2500];
%% save figure
figFname = sprintf('surfaceAreaVsThickness_all_28_transpose_%s_%d-%d', ...
                    normalizeType{normTIVflag+1},LabelList(1),LabelList(end));
%%
figPath = fullfile(resultDir,figFname);
if exist('figFname','var')
    MatlabVersion = version;
    MatlabVersion =MatlabVersion(end-5:end-2);
    if str2double(MatlabVersion) >= 2020
        t.Units = 'inches';
        t.OuterPosition = [0 0 10 25];
        exportgraphics (fig, [figPath,'.png'],'Resolution', 200);
    else
        print(fig, figPath, '-dpng', '-r200');
    end
end


%% alternative save figure (will contrain gray area)
% figPath = fullfile(resultDir,figFname);
% if exist('figFname','var')
%     export_fig(figPath,'-m1');
% end


end