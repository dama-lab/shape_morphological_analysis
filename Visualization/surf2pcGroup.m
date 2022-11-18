function [fig,axs] = surf2pcGroup(targetLists, groups, volDir, maskDir, colorLim, bg_color, iso2meshOpt, labelRemapTable, labelsToRemap, bgLabel)
    %%
    % volDir: contain volume to plot surface (with color intensity into)
    % maskDir (optional): contain layer .nii file to plot onto (e.g. purkinje layer) 
    % colorLim (optional): , e.g. [1.356,25.6447]
    % bgLabel: (optional): background label (default: 0) (recommend to preset bgLabel = NaN)
    
    % Use white ('w') background by default % Ref: https://www.mathworks.com/matlabcentral/answers/452951-how-make-the-background-of-pcshow-white-instead-of-black
    if ~exist('bg_color','var') || isempty(bg_color); bg_color = 'w'; end
    
    % by default, using iso2meshFlag to generate surface plot
    if ~exist('iso2meshOpt','var') || ~isstruct(iso2meshOpt) || ~isfield(iso2meshOpt,'iso2meshFlag')
        iso2meshOpt.iso2meshFlag = true;
    end

    % by default, set background label to zero
    if ~exist('bgLabel','var') || isempty(bgLabel); bgLabel = 0; end
    
    groupNum = size(targetLists,2);
    targetNum = size(targetLists,1);

    fig = figure;
    tiledlayout(groupNum,targetNum);
    colorLimRange = [];
    axs = []; % Prepare axes vector to rotatesimutaneously

    for groupId = 1:size(targetLists,2)
        for targetId = 1:size(targetLists,1)
            target = targetLists{targetId,groupId};
            % skip empty target
            if target == ""
                nexttile;
                sizeN = 1000;
                blackImg = zeros(sizeN,sizeN,3);
                imshow(blackImg);
                colorMin = min(colorLimRange(:,1));
                colorMax = max(colorLimRange(:,2));
                %% add colorbar
                if exist('colorLim','var') && ~isempty(colorLim)
                    caxis(axs(end),colorLim);
                else
                    caxis(axs(end),[colorMin,colorMax]);
                end
                colorbar;
                % %% create uniformly distributed random number betwen [colorMin,colorMax]
                % randImg = rand(sizeN)*(colorMax-colorMin)+colorMin;
                continue
            end
            fprintf('(%d/%d) %s - (%d/%d) %s\n',groupId,length(groups),groups{groupId}, ...
                    targetId, length(targetLists), target);

            % target = 'tc1_269455-ob_c.nii.gz';
            
            %% read volNii (contain color intensity into) 
            volNii = fullfile(volDir,target);
            vol = niftiread(volNii);
            % convert intensity type to single (has to be: uint8, single, or double)
            vol = single(vol); 

            %% check if maskDir not specified or same as volDir
            if ~exist('maskDir','var') || isempty(maskDir) || strcmp(maskDir, volDir)
                mask = (vol ~= bgLabel);
            else
                maskNii = fullfile(maskDir,target);
                mask = niftiread(maskNii);
            end
            % colormapSurf = uint8(niftiread(colormapNii));

            
            %% remap the label if "labelRemapTable", "labelsToRemap" are specified
            if exist('labelRemapTable','var') && exist('labelsToRemap','var')
                volNew = zeros(size(vol));
                % row number
                rowId = find(strcmpi(labelRemapTable.Label,target));
                labelRemapTable(rowId,:);
                for label = labelsToRemap(1):labelsToRemap(end)
                    labelStr = num2str(label);
                    newValue = labelRemapTable.(labelStr)(rowId,:);
                    volNew(vol==label) = newValue;
                end
                vol = volNew;
                % empty volNew
                volNew = [];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            convertOpt.convert_type = 'vol2pc';
            convertOpt.vol = vol;
            convertOpt.mask = mask;
            convertOpt.bgLabel = 0;
            convertOpt.iso2meshOpt = iso2meshOpt;
            convertOpt = mri_converter(convertOpt);
            convertOpt.iso2meshOpt = iso2meshOpt;
            pc = convertOpt.pc;
            %pc = pcdownsample(pc,'random',0.5);    
                

            %% Visualization
%             fig;
    %         subplot()
    %         tiledlayout(groupNum,targetNum);
            nexttile;
            pcshow(pc,'MarkerSize',50);
            axis off;
            set(gca,'color',bg_color); % no need if set axis off

            title(target,'Interpreter','none');
            % get current colormap limit
            display(caxis);
            colorLimRange = [colorLimRange; caxis];
            % set colormap limit to consistent rang
            % caxis(colorLim);
            colormap('jet');
            % colormap('hsv');
            axs = [axs, gca];
        end
    end

    
    %% Syncronize all subplot figures
    % Link xylim Ref: https://www.mathworks.com/help/matlab/ref/linkaxes.html
    % linkaxes(axs,'xy');
    % Link camera view. Ref: https://www.mathworks.com/help/matlab/ref/linkprop.html
    hlink = linkprop(axs,{'CameraPosition','CameraUpVector','CameraViewAngle',...
        'CameraViewAngleMode','CameraUpVectorMode','CLim','Colormap'}); % ,'View'
    
    %% set ColorLim
    colorMin = min(colorLimRange(:,1));
    colorMax = max(colorLimRange(:,2));
    display([colorMin, colorMax]);
    % re-color the surfaces with customized or auto-determined colorlim
    if exist('colorLim','var') && ~isempty(colorLim)
        caxis(axs(end),colorLim);
    else
        caxis(axs(end),[colorMin,colorMax]);
    end
    % colorbar(axs(end));
    
    % Set Camera view (rotation)
    axs(1).CameraViewAngle = 8; % scaling
    axs(1).View = [-90,50]; % Camera angle
    % Alternatively:
    % axs(1).CameraUpVector = [1,0,0.5];
    
    % Set figure color (will make title of target name invisible)
    set(gcf,'color',bg_color);
    
    %% To further update (add/remove) linked target/property:
    % Ref: https://www.mathworks.com/help/matlab/ref/linkprop.html?s_tid=srchtitle
    % removetarget
    
    % or removeprop
    
    
    %% To further manually adjust the camera angle
    % Adjust the rotate simutaniously
    % rotate3d on;
    % After the adjustment, turn the rotation off
    % rotate3d off;
    %     for i = 1:36
    %        camorbit(10,i,'data',[0 1 0])
    %        drawnow
    %     end

end