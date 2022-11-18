function fig = plot_tiv_bv_cerebellar_volume()
    %% 
    
    %% load options
    plot_opts = load_opts();
    v2struct(plot_opts);
    %%
    [volTable,~,group,~,~,~] = readTIVcsv(TIVCsv);
    %
    group = [repmat("WT",sum(strcmp(volTable.population,'wt')),1);
    repmat("Tc1",sum(strcmp(volTable.population,'tg')),1)];
    
    %% construct volume matrix
    volMatrix = [volTable.TIV,volTable.BV,volTable.Cerebellum];
    
    %% Figure plot
    close all;
    fig = figure;
    ha = tight_subplot(1,4,[0.1,0.05],[0.1,0.1],[0.05,0.05]);

    %% Resize figure
    fig.Position(3) = 700;
    fig.Position(4) = 180;
    r2_positions.WT = [0.6,0.05];
    r2_positions.Tc1 = [0.05,0.65];

    % boxplot
    titles = {'TIV','BV','Cerebellum'};
    for axId = 1:size(volMatrix,2)
        %% Box plot
        volArray = volMatrix(:,axId);
        boxplot(ha(axId),volArray, group,'BoxStyle','filled','Colors','br');
        title(ha(axId),titles{axId});
        %% Statistical T-test
        wtId = strcmp(group,"WildType");
        tgId = ~wtId;
        [h,p,ci,stats] = ttest2(volArray(wtId),volArray(tgId),'Tail','left');
        [h(axId),p(axId),ci(:,axId),stats(axId)] = ...
            ttest2(volArray(wtId),volArray(tgId),'Tail','left');
        %% adjust y lim
        if axId < 3
            ylim(ha(axId),[400,600]);
        end    
    end
    
    % Multiple comparison correction
    [adj_h, crit_p,adj_c_cvrp, adj_p] = fdr_bh(p,0.05);
    
    % scatter plot
    
    % align the verticle locations
    ha(end).Position(2) = ha(end-1).Position(2);
    ha(end).Position(4) = ha(end-1).Position(4);
    % for wild type
    groups = {'Tc1','WT'};
    markers = {'r.','b.'};
    colors = {'r','b'};
    for g = 1:length(groups)
        %% determine group
        gCurrent = groups{g};
        groupId = strcmp(group,gCurrent);
        
        %% get TIV/Cereebellum
        TIV = volMatrix(groupId, 1);
        Cerebellum = volMatrix(groupId, 3);

        scatPlt.(gCurrent) = plot(ha(end),Cerebellum,TIV,markers{g},'MarkerSize',10);
        hold on;

        xlim(ha(end),[45,70]);
        ylim(ha(end),[450,600]);
        xticks(ha(end),[]); % [50,65]);
        xticklabels(ha(end), []); % [50,65]);
        yticks(ha(end),[]); % [50,65]);
        yticklabels(ha(end), []); % [50,65]);
        xlabel(ha(end),'Cerebellum');
        ylabel(ha(end),'TIV');

        % Linear fit 
        coefs = polyfit(Cerebellum, TIV, 1);
        slope = coefs(1);
        intercept = coefs(2);
        x = [47,68];
        y = slope*x+intercept;
        plot(x,y,'--','LineWidth',1,'Color', colors{g});
        hold on;

        %% Calculate R^2
        % sAreaLinpred = polyval(coefs, thick);
        [b,bint,r,rint,stats] = regress(Cerebellum,[ones(size(TIV)),TIV]);
        r2 = stats(1);
        
        %% display R^2
        r2_position = r2_positions.(gCurrent);
        r2_text = text(r2_position(1), r2_position(2), ...
            sprintf('$R^2$=%0.2f',r2), 'Color', colors{g}, ...
            'FontSize', 10,'Units','normalized','Visible','on','Interpreter','latex');

    end

    
    %% Legend
    legend([scatPlt.Tc1, scatPlt.WT],groups,'Location','northwest','NumColumns',1);
    %% Title
    title(ha(end),'TIV vs. Cerebellum');
    %% save figure
    figFname = fullfile(resultDir,'/TIV_BV_Cerebellar_Volume.png');
    saveMethod = "native"; % "export_fig"; % 
    if strcmp(saveMethod, "export_fig")
        % preferred way
        export_fig(figFname, '-m2');
    elseif strcmp(saveMethod, "native")
        % Native way to save figure
        MatlabVersion = version;
        MatlabVersion =MatlabVersion(end-5:end-2);
        if str2double(MatlabVersion) >= 2020
            t.Units = 'inches';
            t.OuterPosition = [0 0 10 25];
            exportgraphics (fig, [figFname,'.png'],'Resolution', 200);
        else
            print(fig, figFname, '-dpng', '-r200');
        end
    end

    