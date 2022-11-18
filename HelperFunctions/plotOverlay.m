function plotOverlay(Vol, ax)
    %% get the given or current axis
    if ~exist('ax','var') || isempty(ax)
        ax = gca;
    end
    
    %% show segmentation images
    ax = axes;
    imagesc(ax,Vol);
    colormap(ax,colormapReverse);
    linkprop([ha(haId),ax],'Position');

end