% Author: Junran Chen
% Date: 2024-May-19
% Function: OCV GA's OBJ function

function [OBJ] = Junran_OCV_OBJ(Opt_Param)
%% Load data
load('input-data/Samsung30T_OCV.mat');
run param/params_Samsung30T.m
% p.OCP_Cathode   = Opt_Param(1:19);
volt_exp = meas.Voltage(1:1205,1);
time_exp = meas.Time(1:1205,1);
current_exp = meas.Current(1:1205,1);
I = -current_exp'/p.Area;
t = time_exp';
p.delta_t = t(2)-t(1);
V0 = 4.173;


%% OCV calculation -> Theta got from Weihan's paper
Theta_n = linspace(0.855, 0.008, 1205);
Theta_p = linspace(0.254, 0.915, 1205);
V_p = Junran_refPotentialCathode(p, Theta_p);
V_n = Junran_refPotentialAnode(p, Theta_n);
OCV_estim = V_p - V_n;
error = volt_exp' - OCV_estim;
OBJ = sqrt(mean(error.^2));  %RMSE
if ~isreal(OBJ)
    OBJ = 100000;
end

end