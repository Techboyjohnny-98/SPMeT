% Author: Junran Chen
% Date: 2024-May-15
% Function: Main function to use GA indentify OCV curves
clear all;
close all;


%% Parameter initialization
nvars = 19;

OCP_Cathode_1_init = 2.16216;   % OCP parameters for Cathode. [0.1, 2.5]
OCP_Cathode_1_lb = 0.1;
OCP_Cathode_1_ub = 2.5;
OCP_Cathode_2_init = 0.07645;   %[0, 2.5]
OCP_Cathode_2_lb = 0;
OCP_Cathode_2_ub = 2.5;
OCP_Cathode_3_init = 30.834;   %[0, 100]
OCP_Cathode_3_lb = 0;
OCP_Cathode_3_ub = 100;
OCP_Cathode_4_init = 54.4806;   %[0, 100]
OCP_Cathode_4_lb = 0;
OCP_Cathode_4_ub = 100;
OCP_Cathode_5_init = 2.1581;   %[0, 2.5]
OCP_Cathode_5_lb = 0;
OCP_Cathode_5_ub = 2.5;
OCP_Cathode_6_init = 52.294;   %[0, 100]
OCP_Cathode_6_lb = 0;
OCP_Cathode_6_ub = 100;
OCP_Cathode_7_init = 50.294;   %[0, 100]
OCP_Cathode_7_lb = 0;
OCP_Cathode_7_ub = 100;
OCP_Cathode_8_init = 0.14169;   %[0, 2.5]
OCP_Cathode_8_lb = 0;
OCP_Cathode_8_ub = 2.5;
OCP_Cathode_9_init = 11.0923;   %[0, 100]
OCP_Cathode_9_lb = 0;
OCP_Cathode_9_ub = 100;
OCP_Cathode_10_init = 19.8543;   %[0, 100]
OCP_Cathode_10_lb = 0;
OCP_Cathode_10_ub = 100;
OCP_Cathode_11_init = 0.2051;   %[0, 2.5]
OCP_Cathode_11_lb = 0;
OCP_Cathode_11_ub = 2.5;
OCP_Cathode_12_init = 1.4684;   %[0, 100]
OCP_Cathode_12_lb = 0;
OCP_Cathode_12_ub = 100;
OCP_Cathode_13_init = 5.4888;   %[0, 100]
OCP_Cathode_13_lb = 0;
OCP_Cathode_13_ub = 100;
OCP_Cathode_14_init = 0.2531;   %[0, 1]
OCP_Cathode_14_lb = 0;
OCP_Cathode_14_ub = 1;
OCP_Cathode_15_init = 0.56478;   %[0, 1]
OCP_Cathode_15_lb = 0;
OCP_Cathode_15_ub = 1;
OCP_Cathode_16_init = 0.1316;   %[0, 1]
OCP_Cathode_16_lb = 0;
OCP_Cathode_16_ub = 1;
OCP_Cathode_17_init = 0.02167;   %[0, 1]
OCP_Cathode_17_lb = 0;
OCP_Cathode_17_ub = 1;
OCP_Cathode_18_init = 0.525;   %[0, 1]
OCP_Cathode_18_lb = 0;
OCP_Cathode_18_ub = 1;
OCP_Cathode_19_init = 0.006;   %[0, 1]
OCP_Cathode_19_lb = 0;
OCP_Cathode_19_ub = 1;




InitialPopulation_Data = [
    OCP_Cathode_1_init,...
    OCP_Cathode_2_init,...
    OCP_Cathode_3_init,...
    OCP_Cathode_4_init,...
    OCP_Cathode_5_init,...
    OCP_Cathode_6_init,...
    OCP_Cathode_7_init,...
    OCP_Cathode_8_init,...
    OCP_Cathode_9_init,...
    OCP_Cathode_10_init,...
    OCP_Cathode_11_init,...
    OCP_Cathode_12_init,...
    OCP_Cathode_13_init,...
    OCP_Cathode_14_init,...
    OCP_Cathode_15_init,...
    OCP_Cathode_16_init,...
    OCP_Cathode_17_init,...
    OCP_Cathode_18_init,...
    OCP_Cathode_19_init,...
    ];

%% Parameter boundary
lb = [
    OCP_Cathode_1_lb,...
    OCP_Cathode_2_lb,...
    OCP_Cathode_3_lb,...
    OCP_Cathode_4_lb,...
    OCP_Cathode_5_lb,...
    OCP_Cathode_6_lb,...
    OCP_Cathode_7_lb,...
    OCP_Cathode_8_lb,...
    OCP_Cathode_9_lb,...
    OCP_Cathode_10_lb,...
    OCP_Cathode_11_lb,...
    OCP_Cathode_12_lb,...
    OCP_Cathode_13_lb,...
    OCP_Cathode_14_lb,...
    OCP_Cathode_15_lb,...
    OCP_Cathode_16_lb,...
    OCP_Cathode_17_lb,...
    OCP_Cathode_18_lb,...
    OCP_Cathode_19_lb
    ];

ub = [
    OCP_Cathode_1_ub,...
    OCP_Cathode_2_ub,...
    OCP_Cathode_3_ub,...
    OCP_Cathode_4_ub,...
    OCP_Cathode_5_ub,...
    OCP_Cathode_6_ub,...
    OCP_Cathode_7_ub,...
    OCP_Cathode_8_ub,...
    OCP_Cathode_9_ub,...
    OCP_Cathode_10_ub,...
    OCP_Cathode_11_ub,...
    OCP_Cathode_12_ub,...
    OCP_Cathode_13_ub,...
    OCP_Cathode_14_ub,...
    OCP_Cathode_15_ub,...
    OCP_Cathode_16_ub,...
    OCP_Cathode_17_ub,...
    OCP_Cathode_18_ub,...
    OCP_Cathode_19_ub
    ];
PopInitRange_Data      = [lb; ub];

%% GA population parameters
PopulationSize_Data    = 50;% fill here;
Generations_Data       = 100;% fill here;·

%% Run GA
% [x,fval,exitflag,output,population,score] = Junran_runGA(nvars,...
%                                                    lb(1:17),...
%                                                    ub(1:17),...
%                                                    PopInitRange_Data(:,1:17),...
%                                                    PopulationSize_Data,...
%                                                    Generations_Data,...
%                                                    InitialPopulation_Data(1:17));
% %% Run GA
[x,fval,exitflag,output,population,score] = Junran_runGA(nvars,...
                                                   lb,...
                                                   ub,...
                                                   PopInitRange_Data,...
                                                   PopulationSize_Data,...
                                                   Generations_Data,...
                                                   InitialPopulation_Data);


