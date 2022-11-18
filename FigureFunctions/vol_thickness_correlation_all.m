function fig = vol_thickness_correlation_all(...
    thicknessFullCsv, surfaceAreaFullCsv, ...
    thicknessGranCsv, surfaceAreaGranCsv, ...
    thicknessMoleCsv, surfaceAreaMoleCsv, ...
    normTIVflag, TIVCsv, figTitle, figFname)

%% plot correlation between volume and thickness
% by default normalize to TIV
if ~exist('normTIVflag','var'); normTIVflag = 1; end
color.WT = 'b';
color.TG = 'r';


%% load structural surface area
surfarea = [];
thickness = [];

% load structural surface area
[surfarea.full,targetList,structList,structLabels,subjectGroups] = readMeasureCsv(surfaceAreaFullCsv);
[surfarea.gran,~,~,~,~] = readMeasureCsv(surfaceAreaGranCsv);
[surfarea.mole,~,~,~,~] = readMeasureCsv(surfaceAreaMoleCsv);

% load structural thickness
[thickness.full,~,~,~] = readMeasureCsv(thicknessFullCsv);
[thickness.gran,~,~,~] = readMeasureCsv(thicknessGranCsv);
[thickness.mole,~,~,~] = readMeasureCsv(thicknessMoleCsv);

% normalize TIV if normTIVflag ==1
layers = {'full','gran','mole'};
for layerId = 1:length(layers)
    if normTIVflag ==1
        TIV = readTIVcsv(TIVCsv);
        layer = layers{layerId};
        surfarea.(layer) = GLM(surfarea.(layer), TIV, subjectGroups,'WT');
        thickness.(layer) = GLM(thickness.(layer), TIV, subjectGroups,'WT');
    end
end


%% plot correlation for all structures, color-coding the group (red: WT; blue: TG)
groups = unique(subjectGroups);
groupNo = length(groups);
structNo = size(surfarea.full,2);

fig = figure('Visible','off');
corrTypes = 5;
rows = corrTypes;
columns = structNo;
ylabels = {{'Full Cortical','surface Area'},{'Molecular Layer','Surface Area'},{'Granular Layer','Surface Area'}, ...
    {'Molecular Layer','Surface Area'},{'Granular Layer','Surface Area'}};
xlabels = {'Thickness',{'Granular Layer','Thickness'},{'Molecular Layer','Thickness'}};
tiledlayout(rows,columns);
axId = 0;

% Text location of statistical resutls
% R-square
r2_positions.WT = [0.15,0.9];
r2_positions.TG = [0.55,0.9];
% P-value of 
p_positions.WT = [0.15,0.1];
p_positions.TG = [0.55,0.1];

for c=1:corrTypes
    for s = 1:structNo
        nexttile;
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
            elseif (c==4) % row 4: molecular surfacearea vs. granular surface 
                sArea = surfarea.mole(subjId,s);
                thick = thickness.gran(subjId,s);
            elseif (c==5) % row 5: granular surfacearea vs. molecular surface
                sArea = surfarea.gran(subjId,s);
                thick = thickness.mole(subjId,s);
            end
           
            
            %% scatter plot
            scatPlt.(group) = plot(thick, sArea, '.', 'MarkerSize', 10, 'Color', color.(group));
            hold on;
            if normTIVflag == 1
                xlim([-4.9, 4.9]);
                ylim([-4.9, 4.9]);
            end
            
            % Linear fit
            coefs = polyfit(thick,sArea,1);
            slope = coefs(1);
            intercept = coefs(2);

            % plot linear fit
            % line_plot = lsline; % (urgly)
            if normTIVflag == 0
                x = [min(thick),max(thick)];
            else
                x = [-4, 4];
            end
            y = slope*x + intercept;
            plot(x,y,'--','LineWidth',1);
            hold on;

            % compute significance (p) of correlation of coefficient
            [R,P,RL,RU] = corrcoef(thick,sArea);
            % the off-diagnal entries of P is the p-vallue for the variable pair 
            p_corr = P(2,1); % = P(1,2);
            
            % display p-value
            p_position = p_positions.(group);
            p_text = text(p_position(1), p_position(2), ...
                {sprintf('p=%0.2f',p_corr)}, 'color', color.(group), ...
                'FontSize', 10,'Units','normalized','Visible','on');
            
            % Calculate R^2
            sAreaLinpred = polyval(coefs, thick);
            [b,bint,r,rint,stats] = regress(thick,[ones(size(sArea)),sArea]);
            r2 = stats(1);

            % display R^2
            r2_position = r2_positions.(group);
            r2_text = text(r2_position(1), r2_position(2), ...
                {sprintf('R^2=%0.2f',r2)}, 'color', color.(group), ...
                'FontSize', 10,'Units','normalized','Visible','on');

            % Add title
            if axId < columns
                title(structList{s});
            end
            % Add ylabel
            if ~mod(axId, columns)
                ylabel(ylabels{c});
            end
            % Add ylabel
            if axId < (rows-2)*columns
                xlabel(xlabels{1});
            elseif (axId >= (rows-2)*columns) && axId < (rows-1)*columns
                xlabel(xlabels{2});
            elseif (axId >= (rows-1)*columns)
                xlabel(xlabels{3});
            end

        end
            
        axId = axId + 1;
    end
end

%% Add legend
corr_legend = legend([scatPlt.WT,scatPlt.TG],groups,'Location','southoutside');

%% Add Overall figure title
if exist('figTitle','var')
    sgtitle(figTitle);
end

%% Adjust Figure size
fig.Position = [0, 0, 2500, 1000];
    
%% save figure
if exist('figFname','var')
    MatlabVersion = version;
    MatlabVersion =MatlabVersion(end-5:end-2);
        if str2double(MatlabVersion) >= 2020
        exportgraphics (fig, [figFname,'.png'],'Resolution', 600);
    else
        print(fig, figFname, '-dpng', '-r600');
    end
end

%% Close figures
close all;

end