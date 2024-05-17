% Author: Junran Chen
% Date: 2024-May-15
% Function: Main function to use GA indentify EMs parameters
clear all;
close all;
% global Batt;
% format long
% Data = load("input-data/UDDS.mat");
% Batt.RecordingTime          = Data.time_exp;
% Batt.I                      = -Data.current_exp;

%% Parameter initialization
nvars = 63;
% L_p_init = 54.6e-6;         % Thickness of negative electrode [m]
% L_n_init = 60.6e-6;         % Thickness of separator [m]
% c_s_p_max_init = 50778;     % Max concentration in cathode, [mol/m^3]
% c_s_n_max_init = 32095;     % Max concentration in anode, [mol/m^3]
% epsilon_s_p_init = 0.4322;  % Volume fraction in solid for pos. electrode
% epsilon_s_n_init = 0.4815;  % Volume fraction in solid for neg. electrode
% Area_init = 0.3828;         % Electrode current collector area [m^2]
% % k_p0_init = 6.82e-11;       % Reaction rate in pos. electrode
% % k_n0_init = 9.19e-11;       % Reaction rate in neg. electrode
% R_s_p_init = 7.08e-6;       % Radius of solid particles in positive electrode [m]
% R_s_n_init = 8.71e-6;       % Radius of solid particles in negative electrode [m]
% D_s_p0_init = 5.04e-14;     % Diffusion coeff for solid in pos. electrode, [m^2/s]
% D_s_n0_init = 2.98e-14;     % Diffusion coeff for solid in neg. electrode, [m^2/s]
% R_f_n_init = 0.0081;        % Resistivity of SEI layer, [Ohms*m^2]
% epsilon_e_p_init = 0.3867;  % Volume fraction in electrolyte for pos. electrode
% epsilon_e_n_init = 0.4373;  % Volume fraction in electrolyte for neg. electrode
% c_e_init = 1133;            % Fixed electrolyte concentration for SPM, [mol/m^3]
% L_s_init = 19e-6;           % Thickness of separator [m]
% t_plus_init = 0.2744;       % Transference number

L_p_init = 100e-6;         % Thickness of negative electrode [m]
L_n_init = 100e-6;         % Thickness of separator [m]
c_s_p_max_init = 5.1219e+04;     % Max concentration in cathode, [mol/m^3]
c_s_n_max_init = 2.4984e+04;     % Max concentration in anode, [mol/m^3]
epsilon_s_p_init = 0.5;  % Volume fraction in solid for pos. electrode
epsilon_s_n_init = 0.6;  % Volume fraction in solid for neg. electrode
Area_init = 1;         % Electrode current collector area [m^2]
% k_p0_init = 6.82e-11;       % Reaction rate in pos. electrode
% k_n0_init = 9.19e-11;       % Reaction rate in neg. electrode
R_s_p_init = 10e-6;       % Radius of solid particles in positive electrode [m]
R_s_n_init = 10e-6;       % Radius of solid particles in negative electrode [m]
D_s_p0_init = 1e-13;     % Diffusion coeff for solid in pos. electrode, [m^2/s]
D_s_n0_init = 3.9e-14;     % Diffusion coeff for solid in neg. electrode, [m^2/s]
R_f_n_init = 1e-3;        % Resistivity of SEI layer, [Ohms*m^2]
epsilon_e_p_init = 0.3;  % Volume fraction in electrolyte for pos. electrode
epsilon_e_n_init = 0.3;  % Volume fraction in electrolyte for neg. electrode
c_e_init = 1e3;            % Fixed electrolyte concentration for SPM, [mol/m^3]
L_s_init = 25e-6;           % Thickness of separator [m]
t_plus_init = 0.4;       % Transference number

OCP_Anode_1_init = 0.194;   % OCP parameters for Anode. [0.1, 0.2]
OCP_Anode_1_lb = 0.1;
OCP_Anode_1_ub = 0.2;
OCP_Anode_2_init = 1.5;     % [0.5, 1.6]
OCP_Anode_2_lb = 0.5;
OCP_Anode_2_ub = 1.6;
OCP_Anode_3_init = 120;     % [30, 150]
OCP_Anode_3_lb = 30;
OCP_Anode_3_ub = 150;
OCP_Anode_4_init = 0.0351;  % [0, 0.2]
OCP_Anode_4_lb = 0;
OCP_Anode_4_ub = 0.2;
OCP_Anode_5_init = 0.286;   % [0.1, 1]
OCP_Anode_5_lb = 0.1;
OCP_Anode_5_ub = 1;
OCP_Anode_6_init = 0.083;   % [0, 0.2]
OCP_Anode_6_lb = 0;
OCP_Anode_6_ub = 0.2;
OCP_Anode_7_init = 0.0045;  % [0, 0.2]
OCP_Anode_7_lb = 0;
OCP_Anode_7_ub = 0.2;
OCP_Anode_8_init = 0.849;   % [0.1, 1]
OCP_Anode_8_lb = 0.1;
OCP_Anode_8_ub = 1;
OCP_Anode_9_init = 0.119;   % [0, 0.2]
OCP_Anode_9_lb = 0;
OCP_Anode_9_ub = 0.2;
OCP_Anode_10_init = 0.035;  % [0, 0.2]
OCP_Anode_10_lb = 0;
OCP_Anode_10_ub = 0.2;
OCP_Anode_11_init = 0.9233;   % [0.1, 1]
OCP_Anode_11_lb = 0.1;
OCP_Anode_11_ub = 1;
OCP_Anode_12_init = 0.05;   % [0, 0.2]
OCP_Anode_12_lb = 0;
OCP_Anode_12_ub = 0.2;
OCP_Anode_13_init = 0.0147;  % [0, 0.2]
OCP_Anode_13_lb = 0;
OCP_Anode_13_ub = 0.2;
OCP_Anode_14_init = 0.5;   % [0.1, 1]
OCP_Anode_14_lb = 0.1;
OCP_Anode_14_ub = 1;
OCP_Anode_15_init = 0.034;   % [0, 0.2]
OCP_Anode_15_lb = 0;
OCP_Anode_15_ub = 0.2;
OCP_Anode_16_init = 0.102;  % [0, 0.2]
OCP_Anode_16_lb = 0;
OCP_Anode_16_ub = 0.2;
OCP_Anode_17_init = 0.194;   % [0.1, 1]
OCP_Anode_17_lb = 0.1;
OCP_Anode_17_ub = 1;
OCP_Anode_18_init = 0.142;   % [0, 0.2]
OCP_Anode_18_lb = 0;
OCP_Anode_18_ub = 0.2;
OCP_Anode_19_init = 0.022;  % [0, 0.2]
OCP_Anode_19_lb = 0;
OCP_Anode_19_ub = 0.2;
OCP_Anode_20_init = 0.9;   % [0.1, 1]
OCP_Anode_20_lb = 0.1;
OCP_Anode_20_ub = 1;
OCP_Anode_21_init = 0.0164;   % [0, 0.2]
OCP_Anode_21_lb = 0;
OCP_Anode_21_ub = 0.2;
OCP_Anode_22_init = 0.011;  % [0, 0.2]
OCP_Anode_22_lb = 0;
OCP_Anode_22_ub = 0.2;
OCP_Anode_23_init = 0.124;   % [0.1, 1]
OCP_Anode_23_lb = 0.1;
OCP_Anode_23_ub = 1;
OCP_Anode_24_init = 0.0226;   % [0, 0.2]
OCP_Anode_24_lb = 0;
OCP_Anode_24_ub = 0.2;
OCP_Anode_25_init = 0.0155;  % [0, 0.2]
OCP_Anode_25_lb = 0;
OCP_Anode_25_ub = 0.2;
OCP_Anode_26_init = 0.105;   % [0.1, 1]
OCP_Anode_26_lb = 0.1;
OCP_Anode_26_ub = 1;
OCP_Anode_27_init = 0.029;   % [0, 0.2]
OCP_Anode_27_lb = 0;
OCP_Anode_27_ub = 0.2;

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
    L_p_init,...
    L_n_init,...
    c_s_p_max_init,...
    c_s_n_max_init,...
    epsilon_s_p_init,...
    epsilon_s_n_init,...
    Area_init,...
    R_s_p_init,...
    R_s_n_init,...
    D_s_p0_init,...
    D_s_n0_init,...
    R_f_n_init,...
    epsilon_e_p_init,...
    epsilon_e_n_init,...
    c_e_init,...
    L_s_init,...
    t_plus_init,...
    OCP_Anode_1_init,...
    OCP_Anode_2_init,...
    OCP_Anode_3_init,...
    OCP_Anode_4_init,...
    OCP_Anode_5_init,...
    OCP_Anode_6_init,...
    OCP_Anode_7_init,...
    OCP_Anode_8_init,...
    OCP_Anode_9_init,...
    OCP_Anode_10_init,...
    OCP_Anode_11_init,...
    OCP_Anode_12_init,...
    OCP_Anode_13_init,...
    OCP_Anode_14_init,...
    OCP_Anode_15_init,...
    OCP_Anode_16_init,...
    OCP_Anode_17_init,...
    OCP_Anode_18_init,...
    OCP_Anode_19_init,...
    OCP_Anode_20_init,...
    OCP_Anode_21_init,...
    OCP_Anode_22_init,...
    OCP_Anode_23_init,...
    OCP_Anode_24_init,...
    OCP_Anode_25_init,...
    OCP_Anode_26_init,...
    OCP_Anode_27_init,...
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
    % low boundary
L_p_lb = 35e-6;         % Thickness of negative electrode [m]
L_n_lb = 35e-6;         % Thickness of separator [m]
c_s_p_max_lb = 48000;     % Max concentration in cathode, [mol/m^3]
c_s_n_max_lb = 20000;     % Max concentration in anode, [mol/m^3]
epsilon_s_p_lb = 0.35;  % Volume fraction in solid for pos. electrode
epsilon_s_n_lb = 0.4;  % Volume fraction in solid for neg. electrode
Area_lb = 0.3;         % Electrode current collector area [m^2]
% k_p0_lb = 1e-11;       % Reaction rate in pos. electrode
% k_n0_lb = 1e-11;       % Reaction rate in neg. electrode
R_s_p_lb = 1e-6;       % Radius of solid particles in positive electrode [m]
R_s_n_lb = 1e-6;       % Radius of solid particles in negative electrode [m]
D_s_p0_lb = 1e-14;     % Diffusion coeff for solid in pos. electrode, [m^2/s]
D_s_n0_lb = 1e-14;     % Diffusion coeff for solid in neg. electrode, [m^2/s]
R_f_n_lb = 0.001;        % Resistivity of SEI layer, [Ohms*m^2]
epsilon_e_p_lb = 0.27;  % Volume fraction in electrolyte for pos. electrode
epsilon_e_n_lb = 0.26;  % Volume fraction in electrolyte for neg. electrode
c_e_lb = 1;            % Fixed electrolyte concentration for SPM, [mol/m^3]
L_s_lb = 10e-6;           % Thickness of separator [m]
t_plus_lb = 0.25;       % Transference number
    % up boundary
L_p_ub = 100e-6;         % Thickness of negative electrode [m]
L_n_ub = 100e-6;         % Thickness of separator [m]
c_s_p_max_ub = 52000;     % Max concentration in cathode, [mol/m^3]
c_s_n_max_ub = 36000;     % Max concentration in anode, [mol/m^3]
epsilon_s_p_ub = 0.7;  % Volume fraction in solid for pos. electrode
epsilon_s_n_ub = 0.7;  % Volume fraction in solid for neg. electrode
Area_ub = 1.5;         % Electrode current collector area [m^2]
% k_p0_ub = 10e-11;       % Reaction rate in pos. electrode
% k_n0_ub = 9.19e-11;       % Reaction rate in neg. electrode
R_s_p_ub = 11e-6;       % Radius of solid particles in positive electrode [m]
R_s_n_ub = 11e-6;       % Radius of solid particles in negative electrode [m]
D_s_p0_ub = 10e-14;     % Diffusion coeff for solid in pos. electrode, [m^2/s]
D_s_n0_ub = 10e-14;     % Diffusion coeff for solid in neg. electrode, [m^2/s]
R_f_n_ub = 0.01;        % Resistivity of SEI layer, [Ohms*m^2]
epsilon_e_p_ub = 0.45;  % Volume fraction in electrolyte for pos. electrode
epsilon_e_n_ub = 0.5;  % Volume fraction in electrolyte for neg. electrode
c_e_ub = 1200;            % Fixed electrolyte concentration for SPM, [mol/m^3]
L_s_ub = 30e-6;           % Thickness of separator [m]
t_plus_ub = 0.43;       % Transference number

lb = [
    L_p_lb,...
    L_n_lb,...
    c_s_p_max_lb,...
    c_s_n_max_lb,...
    epsilon_s_p_lb,...
    epsilon_s_n_lb,...
    Area_lb,...
    R_s_p_lb,...
    R_s_n_lb,...
    D_s_p0_lb,...
    D_s_n0_lb,...
    R_f_n_lb,...
    epsilon_e_p_lb,...
    epsilon_e_n_lb,...
    c_e_lb,...
    L_s_lb,...
    t_plus_lb,...
    OCP_Anode_1_lb,...
    OCP_Anode_2_lb,...
    OCP_Anode_3_lb,...
    OCP_Anode_4_lb,...
    OCP_Anode_5_lb,...
    OCP_Anode_6_lb,...
    OCP_Anode_7_lb,...
    OCP_Anode_8_lb,...
    OCP_Anode_9_lb,...
    OCP_Anode_10_lb,...
    OCP_Anode_11_lb,...
    OCP_Anode_12_lb,...
    OCP_Anode_13_lb,...
    OCP_Anode_14_lb,...
    OCP_Anode_15_lb,...
    OCP_Anode_16_lb,...
    OCP_Anode_17_lb,...
    OCP_Anode_18_lb,...
    OCP_Anode_19_lb,...
    OCP_Anode_20_lb,...
    OCP_Anode_21_lb,...
    OCP_Anode_22_lb,...
    OCP_Anode_23_lb,...
    OCP_Anode_24_lb,...
    OCP_Anode_25_lb,...
    OCP_Anode_26_lb,...
    OCP_Anode_27_lb,...
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
    L_p_ub,...
    L_n_ub,...
    c_s_p_max_ub,...
    c_s_n_max_ub,...
    epsilon_s_p_ub,...
    epsilon_s_n_ub,...
    Area_ub,...
    R_s_p_ub,...
    R_s_n_ub,...
    D_s_p0_ub,...
    D_s_n0_ub,...
    R_f_n_ub,...
    epsilon_e_p_ub,...
    epsilon_e_n_ub,...
    c_e_ub,...
    L_s_ub,...
    t_plus_ub,...
    OCP_Anode_1_ub,...
    OCP_Anode_2_ub,...
    OCP_Anode_3_ub,...
    OCP_Anode_4_ub,...
    OCP_Anode_5_ub,...
    OCP_Anode_6_ub,...
    OCP_Anode_7_ub,...
    OCP_Anode_8_ub,...
    OCP_Anode_9_ub,...
    OCP_Anode_10_ub,...
    OCP_Anode_11_ub,...
    OCP_Anode_12_ub,...
    OCP_Anode_13_ub,...
    OCP_Anode_14_ub,...
    OCP_Anode_15_ub,...
    OCP_Anode_16_ub,...
    OCP_Anode_17_ub,...
    OCP_Anode_18_ub,...
    OCP_Anode_19_ub,...
    OCP_Anode_20_ub,...
    OCP_Anode_21_ub,...
    OCP_Anode_22_ub,...
    OCP_Anode_23_ub,...
    OCP_Anode_24_ub,...
    OCP_Anode_25_ub,...
    OCP_Anode_26_ub,...
    OCP_Anode_27_ub,...
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
PopulationSize_Data    = 5;% fill here;
Generations_Data       = 10;% fill here;·

%% Run GA
[x,fval,exitflag,output,population,score] = Junran_runGA(nvars,...
                                                   lb,...
                                                   ub,...
                                                   PopInitRange_Data,...
                                                   PopulationSize_Data,...
                                                   Generations_Data,...
                                                   InitialPopulation_Data);



