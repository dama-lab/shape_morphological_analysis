function labelVolRemapped = labelRemap(labelVol)
    %% remap label numbers into sequential labels: [1,2,3,...]
    labels = unique(labelVol);
    
    labelVolRemapped = zeros(size(labelVol));
    
    for labelId = 1:length(labels)
        labelVolRemapped(labelVol==labels(labelId)) = labelId-1;
    end
        
