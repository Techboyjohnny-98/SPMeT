%% Reference Potential for Neg Electrode: Unref(theta_n)
%   Created July 12, 2011 by Scott Moura

function [Uref,varargout] = Junran_refPotentialAnode(p,theta)

if(~isreal(theta))
    beep;
    error('dfn:err','Complex theta_n');
%     pause;
end
% if theta<0
%     theta = 0;
% end
%% From Weihan's paper

% Uref = OCP_Anode_1...
%      + OCP_Anode_2*exp(+OCP_Anode_3*theta) ...
%      + OCP_Anode_4*tanh((theta+OCP_Anode_5)/OCP_Anode_6) ... 
%      + OCP_Anode_7*tanh((theta+OCP_Anode_8)/OCP_Anode_9) ...
%      + OCP_Anode_10*tanh((theta+OCP_Anode_11)/OCP_Anode_12) ...
%      + OCP_Anode_13*tanh((theta+OCP_Anode_14)/OCP_Anode_15) ...
%      + OCP_Anode_16*tanh((theta+OCP_Anode_17)/OCP_Anode_18) ...
%      + OCP_Anode_19*tanh((theta+OCP_Anode_20)/OCP_Anode_21) ...
%      + OCP_Anode_22*tanh((theta+OCP_Anode_23)/OCP_Anode_24) ...
%      + OCP_Anode_25*tanh((theta+OCP_Anode_26)/OCP_Anode_27);

%%


% % Polynomail Fit
% Uref = ppvalFast(p.Uppn,theta);
OCP_Anode_1 = 0.1379;
OCP_Anode_2 = 0.7526;
OCP_Anode_3 = -35.61;
OCP_Anode_4 = -0.0153;
OCP_Anode_5 = -0.6142;
OCP_Anode_6 = 0.0156;
OCP_Anode_7 = -0.1312;
OCP_Anode_8 = -0.3173;
OCP_Anode_9 = 0.0721;
OCP_Anode_10= -0.1212;
OCP_Anode_11= -0.2120;
OCP_Anode_12= 0.0940;
OCP_Anode_13= -0.1291;
OCP_Anode_14= -0.4524;
OCP_Anode_15= 0.1584;
OCP_Anode_16= -0.1099;
OCP_Anode_17= -0.3976;
OCP_Anode_18= 0.1596;
OCP_Anode_19= -0.1083;
OCP_Anode_20= -0.4246;
OCP_Anode_21= 0.1539;
OCP_Anode_22= -0.1543;
OCP_Anode_23= -0.4003;
OCP_Anode_24= 0.0985;
OCP_Anode_25= 0.7192;
OCP_Anode_26= -0.3684;
OCP_Anode_27= 0.1573;
% DUALFOIL: MCMB 2528 graphite (Bellcore) 0.01 < x < 0.9
Uref = OCP_Anode_1...
     + OCP_Anode_2*exp(+OCP_Anode_3*theta) ...
     + OCP_Anode_4*tanh((theta+OCP_Anode_5)/OCP_Anode_6) ... 
     + OCP_Anode_7*tanh((theta+OCP_Anode_8)/OCP_Anode_9) ...
     + OCP_Anode_10*tanh((theta+OCP_Anode_11)/OCP_Anode_12) ...
     + OCP_Anode_13*tanh((theta+OCP_Anode_14)/OCP_Anode_15) ...
     + OCP_Anode_16*tanh((theta+OCP_Anode_17)/OCP_Anode_18) ...
     + OCP_Anode_19*tanh((theta+OCP_Anode_20)/OCP_Anode_21) ...
     + OCP_Anode_22*tanh((theta+OCP_Anode_23)/OCP_Anode_24) ...
     + OCP_Anode_25*tanh((theta+OCP_Anode_26)/OCP_Anode_27);

% Gradient of OCP wrt theta
if(nargout >= 2)

%     % Polynomial Fit
%     dUref = ppvalFast(p.dUppn,theta);
%     varargout{1} = dUref / p.c_s_n_max;

dUref = +OCP_Anode_2*(OCP_Anode_3/p.c_s_n_max)*exp(+OCP_Anode_3*theta)  ...
 +(OCP_Anode_4/(OCP_Anode_6*p.c_s_n_max))*((cosh((theta+OCP_Anode_5)/OCP_Anode_6)).^(-2)) ...
 +(OCP_Anode_7/(p.c_s_n_max*OCP_Anode_9))*((cosh((theta+OCP_Anode_8)/OCP_Anode_9)).^(-2)) ...
 +(OCP_Anode_10/(p.c_s_n_max*OCP_Anode_12))*((cosh((theta+OCP_Anode_11)/OCP_Anode_12)).^(-2)) ...
 +(OCP_Anode_13/(p.c_s_n_max*OCP_Anode_15))*((cosh((theta+OCP_Anode_14)/OCP_Anode_15)).^(-2)) ...
 +(OCP_Anode_16/(p.c_s_n_max*OCP_Anode_18))*((cosh((theta+OCP_Anode_17)/OCP_Anode_18)).^(-2)) ...
 +(OCP_Anode_19/(p.c_s_n_max*OCP_Anode_21))*((cosh((theta+OCP_Anode_20)/OCP_Anode_21)).^(-2)) ...
 +(OCP_Anode_22/(p.c_s_n_max*OCP_Anode_24))*((cosh((theta+OCP_Anode_23)/OCP_Anode_24)).^(-2)) ...
 +(OCP_Anode_25/(p.c_s_n_max*OCP_Anode_27))*((cosh((theta+OCP_Anode_26)/OCP_Anode_27)).^(-2));
varargout{1} = dUref;

end

% Gradient of OCP wrt temperature
if(nargout >= 3)
    
    dUdT = 0;
    varargout{2} = dUdT;
    
end
% 
