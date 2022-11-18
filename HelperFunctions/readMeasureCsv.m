function [measure,targetList,structList,structLabels,groups] = readMeasureCsv(measureCsv)
%% load structural measurements (volume/thickness/surface_area)
measureCell = readcell(measureCsv);

% Hack: remove subject 22
remove_22 = 0;
if remove_22 == 1
    measureCell = measureCell([1:end-3,end-1:end],:);
end

structList = measureCell(1,3:end);
structLabels = measureCell(2,3:end); %cell2mat(measureCell(2,3:end));
targetList = measureCell(3:end,2);
groups = measureCell(3:end,1);
measure = cell2mat(measureCell(3:end,3:end));
end

