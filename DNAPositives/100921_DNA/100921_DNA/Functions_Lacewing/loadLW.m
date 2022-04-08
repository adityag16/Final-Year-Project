function [time,out] = loadLW(file_data,file_data_all)
%LOADLW loads data from an experiment run on Lacewing

table1 = readtable(file_data, 'Delimiter', ',');
time = table1.TimeElapsed';

data_raw_all = importdata(file_data_all); % file_sata_all is the vschem_export file
out = data_raw_all';

end

