function [Vtemp,Vpeltier] = getPeltierLW(out,time,coord_peltier)
%getPeltierLW returns the voltage sent to the Peltier module for
%temperature regulation
% - out the array output.
% - time is the time vector for the experiment.

% Temperature voltage for regulation
% For whole array : 29:3:50, 32:3:56
% For temp array : 10:17, 11:end
Vtemp = mean(out(coord_peltier,:),1);

Vpeltier_reg = Vtemp(199);

% PID parameters
temp0 = 395;
tempPk = 20;
tempPi = 1;
tempPd = 10;

%Temperature error
temp_err = Vpeltier_reg - Vtemp;

% Regulation terms
tempP = tempPk*temp_err(2:end);
tempD = tempPd*(temp_err(2:end)-temp_err(1:end-1));
Vpeltier = [0,temp0 + tempP + tempD];

end

