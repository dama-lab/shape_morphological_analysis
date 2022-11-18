function surfGroup = vol2surfGroup(targetLists, groups, maskDir, bgLabel, iso2meshOpt)
    %% convert volume (mask) to surface
    % maskDir (optional): contain layer .nii file to plot onto (e.g. purkinje layer) 
        
    % by default, using iso2meshFlag to generate surface plot
    if ~exist('iso2meshOpt','var') || ~isstruct(iso2meshOpt) || ~isfield(iso2meshOpt,'iso2meshFlag')
        iso2meshOpt.iso2meshFlag = true;
    end
    
    if ~isfield(iso2meshOpt,'downsampleRate'); iso2meshOpt.downsampleRate=1; end
    
    % by default, set background label to zero
    if ~exist('bgLabel','var') || isempty(bgLabel); bgLabel = 0; end
    
    groupNum = length(groups); % group number
    targetNum = size(targetLists,1); % target number in each group

    
    %% create a containers.map (Ref: https://www.mathworks.com/help/matlab/ref/containers.map.html)
    surfGroup = containers.Map;
    
    %%
    for groupId = 1:size(targetLists,2)
        for targetId = 1:size(targetLists,1)
            target = targetLists{targetId,groupId};
            %% skip empty target
            if target == ""
                continue
            end
            
            %%
            fprintf('(%d/%d) %s - (%d/%d) %s\n',groupId,length(groups),groups{groupId}, ...
                    targetId, length(targetLists), target);

            %% read mask nifti
            maskNii = fullfile(maskDir, target);
            mask = niftiread(maskNii);
            % convert to binary
            mask = (mask ~= bgLabel);
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%
%             convertOpt.convert_type = 'vol2surf';
%             convertOpt.mask = mask;
%             convertOpt.bgLabel = 0;
%             convertOpt.iso2meshOpt = iso2meshOpt;
%             convertOpt = mri_converter(convertOpt);
%             % store in a containers.map 
%             surfGroup(target) = convertOpt.surf;

            %% convert mask to surface points, and store in dict (containers.map)
            surfGroup(target) = vol2surface(mask,0,iso2meshOpt);
        end
    end
end