function [surfFig,TR, colormaps] = surfPlot(vertex,triedge, vertexData, facecolor, ...
    colormaps,clim, viewAngle, opts, ax)
    %% plot surface from vertex, triedge, and (optionally) vertex color intensity
    %  Input:
    %       vertex: nx3 vertex coordinate 
    %       tridge: nx3 or nx4 triangular edge
    %       vertexData: nx1 vertex intensity color
    
    % create triangular
    TR = triangulation(double(triedge(:,1:3)),double(vertex));
    
    if ~exist('facecolor','var') || isempty(facecolor)
        % grayColormap
        facecolor = [0.9 0.9 0.9]; 
    end
    if ~exist('axis','var')
        ax = gca;
    end
    %% plot after removing isolated nodes (beautiful!)
    axes(ax);
    if exist('vertexData','var') && ~isempty(vertexData)
        surfFig = trisurf(TR, vertexData,'EdgeColor','none');
        %% determing color shading
        if exist('opts','var') && isfield(opts,'shading')
            shading(ax,opts.shading);
        else % default: shading interp
            shading(ax,'interp');
        end
    else
        surfFig = trisurf(TR, 'EdgeColor','none','FaceColor',facecolor);
    end
    
    %% Setup view angle
    if exist('viewAngle','var')
        view(ax, viewAngle);
    end
    %% beautify figure
    % generate colormap if not defined
    if exist('vertexData','var') && (~exist('colormaps','var') || isempty(colormaps))
        colormaps = jet(length(unique(vertexData)));
    elseif exist('colormaps','var') && ~isempty(colormaps)
        colormap(ax, colormaps);
    end
    % beautify
    daspect(ax, [1 1 1]);
%     lh = camlight(ax,'headlight'); lh.Style = 'infinite';    
    ll = camlight(ax,'left'); ll.Style = 'infinite';
    lr = camlight(ax,'right'); lr.Style = 'infinite';
    axis(ax, 'vis3d', 'off'); 
    lighting(ax,'gouraud');
    material(surfFig,'shiny'); 
    % set color ImageBORG9400limits if specified
    if exist('clim','var'); caxis(ax,clim); end
end

