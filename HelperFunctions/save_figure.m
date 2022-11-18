function save_figure(fig, figFname, saveMethod, resolution)
    %% save Figure
    % save methods:
    %   1) 'export_fig'
    %   2) 'native' (default)
    % resolution: default=300
    
    %% use native method by default
    if ~exist('saveMethod','var') || isempty(saveMethod)
        saveMethod = "native"; "export_fig";
    end
    if  ~exist('resolution','var') || isempty(resolution)
        resolution = 300;
    end
    
    %% 
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
            exportgraphics (fig, [figFname,'.png'],'Resolution', resolution);
        else
            print(fig, figFname, '-dpng', sprintf('-r%d',resolution));
        end
    end

