%% Params for Electrochemical Model
%   Created May 24, 2012 by Scott Moura
%
%   Most parameters are from DUALFOIL model
%   D_e is from Capiglia et al, for c_e = 1000 mol/m^3
%   Equilibrium potentials are from Bosch Klein TCST 2011
%   Electrode area chosen to correspond to 2.3 Ah cell
%   Thermal parameters are adopted from L. Zhang (HIT) 2014 10.1016/j.jpowsour.2014.07.110
%   NOTE: Paramters represent an illustrative cell


%% Geometric Params
% Thickness of each layer
p.L_n = 100e-6;     % Thickness of negative electrode [m]
p.L_s = 25e-6;     % Thickness of separator [m]
p.L_p = 100e-6;     % Thickness of positive electrode [m]

L_ccn = 25e-6;    % Thickness of negative current collector [m]
L_ccp = 25e-6;    % Thickness of negative current collector [m]

% Particle Radii
p.R_s_n = 10e-6;   % Radius of solid particles in negative electrode [m]
p.R_s_p = 10e-6;   % Radius of solid particles in positive electrode [m]

% Volume fractions
p.epsilon_s_n = 0.6;      % Volume fraction in solid for neg. electrode
p.epsilon_s_p = 0.5;      % Volume fraction in solid for pos. electrode

p.epsilon_e_n = 0.3;   % Volume fraction in electrolyte for neg. electrode
p.epsilon_e_s = 0.304;   % Volume fraction in electrolyte for separator
p.epsilon_e_p = 0.3;   % Volume fraction in electrolyte for pos. electrode

% make element to caclulate phi_{s} by Saehong Park 
p.epsilon_f_n = 1-p.epsilon_s_n-p.epsilon_e_n;  % Volume fraction of filler in neg. electrode
p.epsilon_f_p = 1-p.epsilon_s_p-p.epsilon_e_p;  % Volume fraction of filler in pos. electrode

epsilon_f_n = p.epsilon_f_n;  % Volume fraction of filler in neg. electrode
epsilon_f_p = p.epsilon_f_p;  % Volume fraction of filler in pos. electrode


% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]

% Mass densities
rho_sn = 1800;    % Solid phase in negative electrode [kg/m^3]
rho_sp = 5010;    % Solid phase in positive electrode [kg/m^3]
rho_e =  1324;    % Electrolyte [kg/m^3]
rho_f = 1800;     % Filler [kg/m^3]
rho_ccn = 8954;   % Current collector in negative electrode
rho_ccp = 2707;   % Current collector in positive electrode

% Compute cell mass [kg/m^2]
m_n = p.L_n * (rho_e*p.epsilon_e_n + rho_sn*p.epsilon_s_n + rho_f*epsilon_f_n);
m_s = p.L_s * (rho_e*p.epsilon_e_n);
m_p = p.L_p * (rho_e*p.epsilon_e_p + rho_sp*p.epsilon_s_p + rho_f*epsilon_f_p);
m_cc = rho_ccn*L_ccn + rho_ccp*L_ccp;

% Lumped density [kg/m^2]
p.rho_avg = m_n + m_s + m_p + m_cc;

%% Transport Params
% Diffusion coefficient in solid
p.D_s_n0 = 3.9e-14;  % Diffusion coeff for solid in neg. electrode, [m^2/s]
p.D_s_p0 = 1e-13;  % Diffusion coeff for solid in pos. electrode, [m^2/s]

p.brug = 1.5;       % Bruggeman porosity

% Conductivity of solid
p.sig_n = 100;    % Conductivity of solid in neg. electrode, [1/Ohms*m]
p.sig_p = 10;    % Conductivity of solid in pos. electrode, [1/Ohms*m]

% Miscellaneous
p.t_plus = 0.4;       % Transference number
p.Faraday = 96485.33289;    % Faraday's constant, [Coulumbs/mol]
p.Area = 1;           % Electrode current collector area [m^2]

%% Kinetic Params
p.R = 8.314472;       % Gas constant, [J/mol-K]

p.alph = 0.5;         % Charge transfer coefficients

p.R_f_n = 1e-3;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_p = 0;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_c = 0;         % Contact Resistance/Current Collector Resistance, [Ohms-m^2]

% Nominal Reaction rates
p.k_n0 = 1e-5;  % Reaction rate in neg. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
p.k_p0 = 3e-7; % Reaction rate in pos. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]

%% Thermodynamic Params

% Ambient Temperature
p.T_amb = 298.15; % [K]

% Activation Energies
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
% All units are [J/mol]
p.E.kn = 37.48e3;
p.E.kp = 39.57e3;
p.E.Dsn = 42.77e3;
p.E.Dsp = 18.55e3;
p.E.De = 37.04e3;
p.E.kappa_e = 34.70e3;

% Reference temperature
p.T_ref = 298.15; %[K]

% Heat transfer parameters
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
p.C1 = 62.7;    % [J/K]
p.C2 = 4.5;     % [J/K]
p.h12 = 10; %1.9386; % [W/K]
p.h2a = 21.45;  % [W/K]

%% Aging submodel parameters

%   SEI Layer Growth model
%   Adopted from Ramadass et al (2004) [Univ of South Carolina]
%   "Development of First Principles Capacity Fade Model for Li-Ion Cells"
%   DOI: 10.1149/1.1634273
%   NOTE: These parameters have NOT been experimentally validated by eCAL

p.kappa_P = 1;      % [S/m] conductivity of side rxn product
p.M_P = 7.3e1;      % [kg/mol] molecular weight of side rxn product
p.rho_P = 2.1e3;    % [kg/m^3] mass density of side rxn product
p.i0s = 0; %1.5e-6;     % [A/m^2] exchange current density of side rxn
p.Us = 0.4;         % [V] reference potential of side rxn

%% Concentrations
% Maxima based on DUALFOIL 
% line 588 in DUALFOIL Fortran code

p.c_s_n_max = 3.6e3 * 372 * 1800 / p.Faraday;   % Max concentration in anode, [mol/m^3]
%p.c_s_n_max = 3.6e3 * 372 * 2260 / p.Faraday;   % Max concentration in anode, [mol/m^3]

%p.c_s_p_max = 3.6e3 * 247 * 5010 / p.Faraday;    % Max concentration in cathode, [mol/m^3]
p.c_s_p_max = 3.6e3 * 274 * 5010 / p.Faraday;    % Max concentration in cathode, [mol/m^3]

p.n_Li_s = 2.5; %2.781;        % Total moles of lithium in solid phase [mol]
p.c_e = 1e3;              % Fixed electrolyte concentration for SPM, [mol/m^3]

%% Cutoff voltages
p.volt_max = 3.6; %4.1113; %4.7;
p.volt_min = 2; %2.6;
%% OCP 
p.OCP_Anode     = [0.194,...   %
                    1.5,...     %
                    120,...
                    0.0351,...
                    0.286,...
                    0.083,... 
                    0.0045,...
                    0.849,...
                    0.119,...
                    0.035,...
                    0.9233,...
                    0.05,...
                    0.0147,...
                    0.5,...
                    0.034,...
                    0.102,...
                    0.194,...
                    0.142,...
                    0.022,...
                    0.9,...
                    0.0164,...
                    0.011,...
                    0.124,...
                    0.0226,...
                    0.0155,...
                    0.105,...
                    0.029];


p.OCP_Cathode   = [2.16216,...
                    0.07645,...
                    30.834,...
                    54.4806,...
                    2.1581,...
                    52.294,...
                    50.294,...
                    0.14169,...
                    11.0923,...
                    19.8543,...
                    0.2051,...
                    1.4684,...
                    5.4888,...
                    0.2531,...
                    0.56478,...
                    0.1316,...
                    0.02167,...
                    0.525,... 
                    0.006];



