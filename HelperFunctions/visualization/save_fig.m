function save_fig(fig, figFname)
    %% Save figure
    
    saveMethod = "native";  %"export_fig";
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
            exportgraphics(fig, figFname, 'Resolution', 200);
        else
            print(fig, figFname, '-dpng', '-r200');
        end
    end