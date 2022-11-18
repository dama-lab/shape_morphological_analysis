function [cmaps, clim, tt,pval] = statColormap(statTest,stat)
    %% generate colormap for different statistics
    % Ref: SurfStatView
    
    switch statTest
        case "T" % T-statistics
            %% Colormap: [-1,1] = gray (0.9)
            cmaps = jet(1200);
            cmaps(501:700,:) = ones(200,3)*0.9;
            clim = [-6,6];
        case "P" % P-value (using Random Theory)
            if ~isfield(stat,'pval')
                error('no Pvalue presented ...')
            else
                pval = stat.pval;
            end
            %% [colormap] two-level (peek and cluster)
            cmaps=[zeros(1,3);
            zeros(127,1)   (0:126)'/127   ones(127,1); ones(1,3)*0.9;
            ones(127,1)    (0:126)'/127   zeros(127,1)];
        
            %% Significant peak value
            if ~isfield(pval,'thresh')
                pval.thresh=0.05;
            end
            %% clim (for caxis)
            clim = [0 255]*pval.thresh;
            signifpeak=pval.P<pval.thresh;
            %% convert p-value => Significant points (t-value)
            % t1: significant clusters
            if isfield(pval,'C')
                signifclus=pval.C<pval.thresh;
                t1=signifclus.*(1-signifpeak).*(127-pval.C/pval.thresh*126);
            else
                signifclus=0;
                t1=0;
            end
            % t2: significant peeks
            t2=signifpeak.*(255-pval.P/pval.thresh*126);
            % t3: non-significant regions 
            t3=(1-signifpeak).*(1-signifclus)*128;
            tt=(t1+t2+t3).*pval.mask*pval.thresh;
        case "Parcellation"
            % stat is a parcellation volume in this case
            %% tt
            tt = [];
            %% clim
            labels = unique(stat);
            clim = [min(labels),max(labels)];
            %% cmap
            cmaps = cmap3DSlicer;
            % cmaps = cmapITKSNAP;
    end
end

function cmaps = cmap3DSlicerReverse()
    cmaps = cmap3DSlicer();
    cmaps = [cmaps(:,3) cmaps(:,2) cmaps(:,1)];
end

function cmaps = cmap3DSlicerPermute()
    cmaps = cmap3DSlicer();
    cmaps = [cmaps(:,2) cmaps(:,3) cmaps(:,1)];
end

function cmaps = cmap3DSlicer()
    % Copied from Slicer3_2010_Brain_Labels / Slicer3_2010_Label_Colors.txt
    cmaps = [];
    cmaps(1,:) = [0,0,0];
    cmaps(2,:) = [231,197,108];
    cmaps(3,:) = [177,122,101];
    cmaps(4,:) = [205,92,92];
    cmaps(5,:) = [238,232,170];
    cmaps(6,:) = [0,191,255];
    cmaps(7,:) = [65,105,255];
    cmaps(8,:) = [250,128,114];
    cmaps(9,:) = [70,130,180];
    cmaps(10,:) = [211,211,211];
    cmaps(11,:) = [255,255,240];
    cmaps(12,:) = [144,238,144];
    cmaps(13,:) = [224,255,255];
    cmaps(14,:) = [216,191,216];
    cmaps(15,:) = [238,130,238];
    cmaps(16,:) = [255,215,0];
    cmaps(17,:) = [106,90,205];
    cmaps(18,:) = [221,160,221];
    cmaps(19,:) = [233,150,122];
    cmaps(20,:) = [165,42,42];
    cmaps(21,:) = [255,250,250];
    cmaps(22,:) = [147,112,219];
    cmaps = cmaps/255;
end

function cmaps = cmapITKSNAP()
    % (copied from ITKSNAP.label)
    cmaps = [];
    cmaps(1,:) = [0,0,0];
    cmaps(2,:) = [255,255,255];
    cmaps(3,:) = [255,17,0];
    cmaps(4,:) = [0,255,0];
    cmaps(5,:) = [0,0,255];
    cmaps(6,:) = [255,255,0];
    cmaps(7,:) = [0,255,255];
    cmaps(8,:) = [255,239,213];
    cmaps(9,:) = [0,0,205];
    cmaps(10,:) = [205,133,63];
    cmaps(11,:) = [210,180,40];
    cmaps(12,:) = [102,205,170];
    cmaps(13,:) = [0,0,128];
    cmaps(14,:) = [0,139,139];
    cmaps(15,:) = [46,139,87];
    cmaps(16,:) = [255,228,225];
    cmaps(17,:) = [106,90,205];
    cmaps(18,:) = [221,160,221];
    cmaps(19,:) = [233,150,122];
    cmaps(20,:) = [165,42,42];
    cmaps(21,:) = [255,250,250];
    cmaps(22,:) = [147,112,219];
    cmaps = cmaps/255;
end