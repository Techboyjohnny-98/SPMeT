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
nvars = 17;
L_p_init = 7.716e-6;         % Thickness of negative electrode [m]
L_n_init = 8.564e-6;         % Thickness of separator [m]
c_s_p_max_init = 50778;     % Max concentration in cathode, [mol/m^3]
c_s_n_max_init = 32095;     % Max concentration in anode, [mol/m^3]
epsilon_s_p_init = 0.4322;  % Volume fraction in solid for pos. electrode
epsilon_s_n_init = 0.4815;  % Volume fraction in solid for neg. electrode
Area_init = 1;         % Electrode current collector area [m^2]
% k_p0_init = 6.82e-11;       % Reaction rate in pos. electrode
% k_n0_init = 9.19e-11;       % Reaction rate in neg. electrode
R_s_p_init = 7.08e-6;       % Radius of solid particles in positive electrode [m]
R_s_n_init = 8.71e-6;       % Radius of solid particles in negative electrode [m]
D_s_p0_init = 5.04e-14;     % Diffusion coeff for solid in pos. electrode, [m^2/s]
D_s_n0_init = 2.98e-14;     % Diffusion coeff for solid in neg. electrode, [m^2/s]
R_f_n_init = 0.0081;        % Resistivity of SEI layer, [Ohms*m^2]
epsilon_e_p_init = 0.3867;  % Volume fraction in electrolyte for pos. electrode
epsilon_e_n_init = 0.4373;  % Volume fraction in electrolyte for neg. electrode
c_e_init = 1133;            % Fixed electrolyte concentration for SPM, [mol/m^3]
L_s_init = 19e-6;           % Thickness of separator [m]
t_plus_init = 0.2744;       % Transference number
n_Li_s_init = 2.5;               % Total moles of lithium in solid phase [mol]
% 
    % Fast dfn
% L_p_init = 100e-6;         % Thickness of negative electrode [m]
% L_n_init = 100e-6;         % Thickness of separator [m]
% c_s_p_max_init = 5.1219e+04;     % Max concentration in cathode, [mol/m^3]
% c_s_n_max_init = 2.4984e+04;     % Max concentration in anode, [mol/m^3]
% epsilon_s_p_init = 0.5;  % Volume fraction in solid for pos. electrode
% epsilon_s_n_init = 0.6;  % Volume fraction in solid for neg. electrode
% Area_init = 1;         % Electrode current collector area [m^2]
% % k_p0_init = 6.82e-11;       % Reaction rate in pos. electrode
% % k_n0_init = 9.19e-11;       % Reaction rate in neg. electrode
% R_s_p_init = 10e-6;       % Radius of solid particles in positive electrode [m]
% R_s_n_init = 10e-6;       % Radius of solid particles in negative electrode [m]
% D_s_p0_init = 1e-13;     % Diffusion coeff for solid in pos. electrode, [m^2/s]
% D_s_n0_init = 3.9e-14;     % Diffusion coeff for solid in neg. electrode, [m^2/s]
% R_f_n_init = 1e-3;        % Resistivity of SEI layer, [Ohms*m^2]
% epsilon_e_p_init = 0.3;  % Volume fraction in electrolyte for pos. electrode
% epsilon_e_n_init = 0.3;  % Volume fraction in electrolyte for neg. electrode
% c_e_init = 1e3;            % Fixed electrolyte concentration for SPM, [mol/m^3]
% L_s_init = 25e-6;           % Thickness of separator [m]
% t_plus_init = 0.4;       % Transference number
    % A123 from J.C's 2012 paper
% L_p_init = 6.5e-5;         % Thickness of negative electrode [m]
% L_n_init = 2.88e-5;         % Thickness of separator [m]
% c_s_p_max_init = 1.035e+04;     % Max concentration in cathode, [mol/m^3]
% c_s_n_max_init = 2.948e+04;     % Max concentration in anode, [mol/m^3]
% epsilon_s_p_init = 0.4794;  % Volume fraction in solid for pos. electrode
% epsilon_s_n_init = 0.3812;  % Volume fraction in solid for neg. electrode
% Area_init = 0.3108;         % Electrode current collector area [m^2]
% % k_p0_init = 6.82e-11;       % Reaction rate in pos. electrode
% % k_n0_init = 9.19e-11;       % Reaction rate in neg. electrode
% R_s_p_init = 1.637e-7;       % Radius of solid particles in positive electrode [m]
% R_s_n_init = 3.6e-6;       % Radius of solid particles in negative electrode [m]
% D_s_p0_init = 1e-13;     % Diffusion coeff for solid in pos. electrode, [m^2/s]
% D_s_n0_init = 3.9e-14;     % Diffusion coeff for solid in neg. electrode, [m^2/s]
% R_f_n_init = 1e-3;        % Resistivity of SEI layer, [Ohms*m^2]
% epsilon_e_p_init = 0.52;  % Volume fraction in electrolyte for pos. electrode
% epsilon_e_n_init = 0.62;  % Volume fraction in electrolyte for neg. electrode
% c_e_init = 1e3;            % Fixed electrolyte concentration for SPM, [mol/m^3]
% L_s_init = 1.697e-6;           % Thickness of separator [m]
% t_plus_init = 0.2495;       % Transference number





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
%     n_Li_s
    ];
% load ./GA_results/May20_50_100.mat
% InitialPopulation_Data = x;
%% Parameter boundary
    % low boundary
L_p_lb = 1e-6;         % Thickness of negative electrode [m]
L_n_lb = 1e-6;         % Thickness of separator [m]
c_s_p_max_lb = 48000;     % Max concentration in cathode, [mol/m^3]
c_s_n_max_lb = 20000;     % Max concentration in anode, [mol/m^3]
epsilon_s_p_lb = 0.35;  % Volume fraction in solid for pos. electrode
epsilon_s_n_lb = 0.4;  % Volume fraction in solid for neg. electrode
Area_lb = 1;         % Electrode current collector area [m^2]
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
L_s_lb = 0.1e-6;           % Thickness of separator [m]
t_plus_lb = 0.25;       % Transference number
n_Li_s_lb = 0.1;           % % Total moles of lithium in solid phase [mol]
    % up boundary
L_p_ub = 100e-6;         % Thickness of negative electrode [m]
L_n_ub = 100e-6;         % Thickness of separator [m]
c_s_p_max_ub = 52000;     % Max concentration in cathode, [mol/m^3]
c_s_n_max_ub = 36000;     % Max concentration in anode, [mol/m^3]
epsilon_s_p_ub = 0.7;  % Volume fraction in solid for pos. electrode
epsilon_s_n_ub = 0.7;  % Volume fraction in solid for neg. electrode
Area_ub = 1;         % Electrode current collector area [m^2]
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
n_Li_s_ub = 3;          % Total moles of lithium in solid phase [mol]

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



