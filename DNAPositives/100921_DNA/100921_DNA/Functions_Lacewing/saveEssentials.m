function saveEssentials(folder_in,folder_out,filename)
% The function saves a mimised version of the data files to accomodate for
% our RAM requirements :-)

% Import data file
data = importdata([folder_in filesep filename '.mat']);

% Name of output file
filename_out = [folder_out filesep filename '_ltd'];

% Variables to save
time = data.time;
ampl_curve_lin = data.ampl_curve_lin;
pixel_on = data.pixel_on;

% Save
save(filename_out,'time','ampl_curve_lin','pixel_on');

end