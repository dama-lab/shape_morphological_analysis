function [LCC, CC] = LargestConnectedComponent(vol)
    %% Largest connected compoinent
    CC = bwconncomp(vol);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    LCC = zeros(size(vol));
    LCC(CC.PixelIdxList{idx}) = 1;
    
    
    %% Alternative
    % Ref: https://www.mathworks.com/matlabcentral/answers/75784-how-to-isolate-and-display-the-largest-connected-component
    CC = bwlabeln(vol);
    



