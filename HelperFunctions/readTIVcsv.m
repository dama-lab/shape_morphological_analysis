function [TIVtable, TIVstruct, TIV, BV, group, subjectList] = readTIVcsv(TIVCsv)
%READTIVCSV
% [TIV, BV, group, subjectList]
%  = readTIVcsv(TIVCsv)
%   

TIVtable = readtable(TIVCsv);
TIV = TIVtable.TIV;
BV = TIVtable.BV;
group = TIVtable.population;
subjectList= TIVtable.file_name;


varLists = TIVtable.Properties.VariableNames;
for varId = 1:length(varLists)
    TIVstruct.(varLists{varId}) = TIVtable.(varLists{varId});
end

