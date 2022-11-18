function fig = vol_thickness_correlation(thicknessCsv, surfareaCsv, TIVCsv, normTIVflag,figTitle, figFname)
%% plot correlation between volume and thickness
% by default normalize to TIV
if ~exist('normTIVflag','var'); normTIVflag = 1; end
color.WT = 'b';
color.TG = 'r';

%% load structural measurements (thickness/surface_area)
[surfarea,targetList,structList,structLabels,subjectGroups] = readMeasureCsv(surfareaCsv);
[thickness,~,~,~] = readMeasureCsv(thicknessCsv);

groups = unique(subjectGroups);
groupNo = length(groups);
structNo = size(surfarea,2);


%% plot correlation for all structures, color-coding the group (red: WT; blue: TG)

%% normalize TIV if normTIVflag ==1
if normTIVflag ==1
    TIV = readTIVcsv(TIVCsv);
    surfarea2Plot = GLM(surfarea, TIV, subjectGroups,'WT');
    thickness2Plot = GLM(thickness, TIV, subjectGroups,'WT');
else
    surfarea2Plot = surfarea;
    thickness2Plot = thickness;
end

fig = figure;
rows = 3;
columns = 5;
ylabels = 'Surface Area';
xlabels = 'Thickness';
tiledlayout(rows,columns);
axId = 0;
for s = 1:structNo
    nexttile;
    for g = 1:groupNo
        % determine points
        group = groups{g};
        subjId = strcmp(subjectGroups, group);
        sArea = surfarea2Plot(subjId,s);
        thick = thickness2Plot(subjId,s);
       
        % scatter plot
        scatPlt.(group) = plot(thick, sArea, '.', 'MarkerSize', 20, 'Color', color.(group));
        hold on;
        if normTIVflag == 1
            xlim([-4, 4]);
            ylim([-4, 4]);
        end
        title(structList{s});
        
        
        % Linear fit
        coefs = polyfit(thick,sArea,1);
        slope = coefs(1);
        intercept = coefs(2);
        
        % plot linear fit
        % line_plot = lsline; % (urgly)
        if normTIVflag == 0
            x = [min(thick),max(thick)];
        else
            x = [-3.5, 3.5];
        end
        y = slope*x + intercept;
        plot(x,y,'--','LineWidth',2);
        hold on;
        
        % Calculate R^2
        sAreaLinpred = polyval(coefs, thick);
        [b,bint,r,rint,stats] = regress(thick,[ones(size(sArea)),sArea]);
        r2 = stats(1);
        
        % display R^2
        r2_positions.WT = [0.55,0.1];
        r2_positions.TG = [0.15,0.1];
        
        r2_position = r2_positions.(group);
        r2_text = text(r2_position(1), r2_position(2), ...
            {sprintf('R^2=%0.2f',r2)}, 'color', color.(group), ...
            'FontSize', 10,'Units','normalized','Visible','on');
%         set(r2_text,)
        
              
        % Add ylabel
        if ~mod(axId, columns)
            ylabel(ylabels);
        end
        % Add ylabel
        if axId >= (rows-1)*columns
            xlabel(xlabels);
        end

    end
    axId = axId + 1;
end

%% Add legend
corr_legend = legend([scatPlt.WT,scatPlt.TG],groups);

%% Add Overall figure title
if exist('figTitle','var')
    sgtitle(figTitle);
end

%% Adjust Figure size
fig.Position = [0, 0, 1300, 700];
    
%% save figure
if exist('figFname','var')
    MatlabVersion = version;
    MatlabVersion =MatlabVersion(end-5:end-2);
    if MatlabVersion >= 2020
        exportgraphics (fig, figFname,'Resolution', 600);
    else
        print(fig, figFname, '-dpng', '-r600');
    end
end

%% Close figures
close all;

end

