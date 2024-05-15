% Author: Junran Chen
% Date: 2024-May-15
% Function: Main function to use GA indentify EMs parameters

global Batt;
format long
Data = load("input-data/UDDS.mat");
Batt.RecordingTime          = time_exp';
Batt.I                      = -Data(1:13703,2);

