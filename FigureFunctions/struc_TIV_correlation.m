function fig = struct_TIV_correlation(measureCsv, TIVCsv)

%% load structural measurements (volume/thickness/surface_area)
[measure,targetList,structList,structLabels] = readMeasureCsv(measureCsv);

%% Load TIV volume
TIV = readTIVcsv(TIVCsv);






