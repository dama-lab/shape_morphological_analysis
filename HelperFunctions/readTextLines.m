function lines = readTextLines(inputPath)
    % Ref: prof_DNN_script.m
    fileID = fopen(inputPath,'r');
    lines = textscan(fileID, '%s','Delimiter','\n');
    fclose(fileID);
    lines = lines{1};
end