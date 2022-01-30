close all; clear all; clc;
addpath('./Functions_Lacewing/');

%% Global options
target = 'George';
% Initialisation
n_exp=1;
options = {'onlyaverage','wells',1};
Ct = zeros(18,1);

%% Copy paste this section for each filename
% Options
options = {'onlyaverage','wells',1};
start_exp{n_exp} = [650 800]; % in s
end_exp{n_exp} = 3500; % in s
pH_init{n_exp} = 8; 
foldername{n_exp} = './Data';
filename{n_exp} = 'beta1.app.1e4';

% Execution
file_data = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_data_export.csv'];
file_data_all = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_vsChem_export.csv'];
[time{n_exp},out_smooth{n_exp},ampl_curve_V{n_exp},ampl_curve_pH{n_exp},...
    ampl_curve_lin{n_exp},ampl_curve_sigm{n_exp},ampl_curve_norm{n_exp},exp_comp{n_exp},...
    sigm_param{n_exp},pixel_on{n_exp},pixel_on_dna{n_exp},...
    xm{n_exp},xi{n_exp},xf{n_exp},Ct(n_exp),ampl{n_exp}]= ...
processDNADataLW(file_data,file_data_all,start_exp{n_exp},end_exp{n_exp},...
    pH_init{n_exp},options);
Ct(n_exp)
n_exp = n_exp+1;

%% Copy paste this section for each filename
% Options
options = {'onlyaverage','wells',1};
start_exp{n_exp} = [770 900]; % in s
end_exp{n_exp} = 2500; % in s
pH_init{n_exp} = 8; 
foldername{n_exp} = './Data';
filename{n_exp} = 'beta2.app.1e5';

% Execution
file_data = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_data_export.csv'];
file_data_all = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_vsChem_export.csv'];
[time{n_exp},out_smooth{n_exp},ampl_curve_V{n_exp},ampl_curve_pH{n_exp},...
    ampl_curve_lin{n_exp},ampl_curve_sigm{n_exp},ampl_curve_norm{n_exp},exp_comp{n_exp},...
    sigm_param{n_exp},pixel_on{n_exp},pixel_on_dna{n_exp},...
    xm{n_exp},xi{n_exp},xf{n_exp},Ct(n_exp),ampl{n_exp}]= ...
processDNADataLW(file_data,file_data_all,start_exp{n_exp},end_exp{n_exp},...
    pH_init{n_exp},options);
Ct(n_exp)
n_exp = n_exp+1;

%% Copy paste this section for each filename
% Options
options = {'onlyaverage','wells',1};
start_exp{n_exp} = [650 800]; % in s
end_exp{n_exp} = 3000; % in s
pH_init{n_exp} = 8; 
foldername{n_exp} = './Data';
filename{n_exp} = 'beta3.app.1e4.4realz2';

% Execution
file_data = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_data_export.csv'];
file_data_all = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_vsChem_export.csv'];
[time{n_exp},out_smooth{n_exp},ampl_curve_V{n_exp},ampl_curve_pH{n_exp},...
    ampl_curve_lin{n_exp},ampl_curve_sigm{n_exp},ampl_curve_norm{n_exp},exp_comp{n_exp},...
    sigm_param{n_exp},pixel_on{n_exp},pixel_on_dna{n_exp},...
    xm{n_exp},xi{n_exp},xf{n_exp},Ct(n_exp),ampl{n_exp}]= ...
processDNADataLW(file_data,file_data_all,start_exp{n_exp},end_exp{n_exp},...
    pH_init{n_exp},options);
Ct(n_exp)
n_exp = n_exp+1;

%% Copy paste this section for each filename
% Options
options = {'onlyaverage','wells',1};
start_exp{n_exp} = [650 800]; % in s
end_exp{n_exp} = 3500; % in s
pH_init{n_exp} = 8; 
foldername{n_exp} = './Data';
filename{n_exp} = 'beta5.app.1e5';

% Execution
file_data = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_data_export.csv'];
file_data_all = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_vsChem_export.csv'];
[time{n_exp},out_smooth{n_exp},ampl_curve_V{n_exp},ampl_curve_pH{n_exp},...
    ampl_curve_lin{n_exp},ampl_curve_sigm{n_exp},ampl_curve_norm{n_exp},exp_comp{n_exp},...
    sigm_param{n_exp},pixel_on{n_exp},pixel_on_dna{n_exp},...
    xm{n_exp},xi{n_exp},xf{n_exp},Ct(n_exp),ampl{n_exp}]= ...
processDNADataLW(file_data,file_data_all,start_exp{n_exp},end_exp{n_exp},...
    pH_init{n_exp},options);
Ct(n_exp)
n_exp = n_exp+1;

%% Copy paste this section for each filename
% Options
options = {'onlyaverage','wells',1};
start_exp{n_exp} = [650 800]; % in s
end_exp{n_exp} = 3000; % in s
pH_init{n_exp} = 8; 
foldername{n_exp} = './Data';
filename{n_exp} = 'gamma1.app.1e5';

% Execution
file_data = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_data_export.csv'];
file_data_all = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_vsChem_export.csv'];
[time{n_exp},out_smooth{n_exp},ampl_curve_V{n_exp},ampl_curve_pH{n_exp},...
    ampl_curve_lin{n_exp},ampl_curve_sigm{n_exp},ampl_curve_norm{n_exp},exp_comp{n_exp},...
    sigm_param{n_exp},pixel_on{n_exp},pixel_on_dna{n_exp},...
    xm{n_exp},xi{n_exp},xf{n_exp},Ct(n_exp),ampl{n_exp}]= ...
processDNADataLW(file_data,file_data_all,start_exp{n_exp},end_exp{n_exp},...
    pH_init{n_exp},options);
Ct(n_exp)
n_exp = n_exp+1;

%% Copy paste this section for each filename
% Options
options = {'onlyaverage','wells',1};
start_exp{n_exp} = [650 800]; % in s
end_exp{n_exp} = 4000; % in s
pH_init{n_exp} = 8; 
foldername{n_exp} = './Data';
filename{n_exp} = 'gamma2.app.1e4';

% Execution
file_data = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_data_export.csv'];
file_data_all = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_vsChem_export.csv'];
[time{n_exp},out_smooth{n_exp},ampl_curve_V{n_exp},ampl_curve_pH{n_exp},...
    ampl_curve_lin{n_exp},ampl_curve_sigm{n_exp},ampl_curve_norm{n_exp},exp_comp{n_exp},...
    sigm_param{n_exp},pixel_on{n_exp},pixel_on_dna{n_exp},...
    xm{n_exp},xi{n_exp},xf{n_exp},Ct(n_exp),ampl{n_exp}]= ...
processDNADataLW(file_data,file_data_all,start_exp{n_exp},end_exp{n_exp},...
    pH_init{n_exp},options);
Ct(n_exp)
n_exp = n_exp+1;

%% Copy paste this section for each filename
% Options
options = {'onlyaverage','wells',1};
start_exp{n_exp} = [650 800]; % in s
end_exp{n_exp} = 4000; % in s
pH_init{n_exp} = 8; 
foldername{n_exp} = './Data';
filename{n_exp} = 'gamma3.app.1e5';

% Execution
file_data = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_data_export.csv'];
file_data_all = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_vsChem_export.csv'];
[time{n_exp},out_smooth{n_exp},ampl_curve_V{n_exp},ampl_curve_pH{n_exp},...
    ampl_curve_lin{n_exp},ampl_curve_sigm{n_exp},ampl_curve_norm{n_exp},exp_comp{n_exp},...
    sigm_param{n_exp},pixel_on{n_exp},pixel_on_dna{n_exp},...
    xm{n_exp},xi{n_exp},xf{n_exp},Ct(n_exp),ampl{n_exp}]= ...
processDNADataLW(file_data,file_data_all,start_exp{n_exp},end_exp{n_exp},...
    pH_init{n_exp},options);
Ct(n_exp)
n_exp = n_exp+1;

%% Copy paste this section for each filename
% Options
options = {'onlyaverage','wells',1};
start_exp{n_exp} = [650 800]; % in s
end_exp{n_exp} = 4000; % in s
pH_init{n_exp} = 8; 
foldername{n_exp} = './Data';
filename{n_exp} = 'gamma5.app.1e4';

% Execution
file_data = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_data_export.csv'];
file_data_all = [foldername{n_exp} filesep filename{n_exp} filesep filename{n_exp} '_vsChem_export.csv'];
[time{n_exp},out_smooth{n_exp},ampl_curve_V{n_exp},ampl_curve_pH{n_exp},...
    ampl_curve_lin{n_exp},ampl_curve_sigm{n_exp},ampl_curve_norm{n_exp},exp_comp{n_exp},...
    sigm_param{n_exp},pixel_on{n_exp},pixel_on_dna{n_exp},...
    xm{n_exp},xi{n_exp},xf{n_exp},Ct(n_exp),ampl{n_exp}]= ...
processDNADataLW(file_data,file_data_all,start_exp{n_exp},end_exp{n_exp},...
    pH_init{n_exp},options);
Ct(n_exp)
n_exp = n_exp+1;

%%
save('../Outputs/Out_Dengue_LW.mat','filename',...
    'time','out_smooth','ampl_curve_V','ampl_curve_pH','ampl_curve_lin',...
    'ampl_curve_sigm','ampl_curve_norm',...
    'exp_comp','sigm_param','pixel_on','pixel_on_dna','xm','xi','xf','Ct',...
    'ampl');
disp('Save successful');