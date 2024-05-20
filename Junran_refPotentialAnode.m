%% Reference Potential for Neg Electrode: Unref(theta_n)
%   Created July 12, 2011 by Scott Moura

function [Uref,varargout] = Junran_refPotentialAnode(p,theta)

if(~isreal(theta))
    beep;
    error('dfn:err','Complex theta_n');
%     pause;
end
% %% From Weihan's paper
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

%%

% 
% % % Polynomail Fit
% % Uref = ppvalFast(p.Uppn,theta);
% OCP_Anode_1 = p.OCP_Anode(1);
% OCP_Anode_2 = p.OCP_Anode(2);
% OCP_Anode_3 = p.OCP_Anode(3);
% OCP_Anode_4 = p.OCP_Anode(4);
% OCP_Anode_5 = p.OCP_Anode(5);
% OCP_Anode_6 = p.OCP_Anode(6);
% OCP_Anode_7 = p.OCP_Anode(7);
% OCP_Anode_8 = p.OCP_Anode(8);
% OCP_Anode_9 = p.OCP_Anode(9);
% OCP_Anode_10 = p.OCP_Anode(10);
% OCP_Anode_11 = p.OCP_Anode(11);
% OCP_Anode_12 = p.OCP_Anode(12);
% OCP_Anode_13 = p.OCP_Anode(13);
% OCP_Anode_14 = p.OCP_Anode(14);
% OCP_Anode_15 = p.OCP_Anode(15);
% OCP_Anode_16 = p.OCP_Anode(16);
% OCP_Anode_17 = p.OCP_Anode(17);
% OCP_Anode_18 = p.OCP_Anode(18);
% OCP_Anode_19 = p.OCP_Anode(19);
% OCP_Anode_20 = p.OCP_Anode(20);
% OCP_Anode_21 = p.OCP_Anode(21);
% OCP_Anode_22 = p.OCP_Anode(22);
% OCP_Anode_23 = p.OCP_Anode(23);
% OCP_Anode_24 = p.OCP_Anode(24);
% OCP_Anode_25 = p.OCP_Anode(25);
% OCP_Anode_26 = p.OCP_Anode(26);
% OCP_Anode_27 = p.OCP_Anode(27);
% % DUALFOIL: MCMB 2528 graphite (Bellcore) 0.01 < x < 0.9
% Uref = OCP_Anode_1...
%      + OCP_Anode_2*exp(-OCP_Anode_3*theta) ...
%      + OCP_Anode_4*tanh((theta-OCP_Anode_5)/OCP_Anode_6) ... 
%      - OCP_Anode_7*tanh((theta-OCP_Anode_8)/OCP_Anode_9) ...
%      - OCP_Anode_10*tanh((theta-OCP_Anode_11)/OCP_Anode_12) ...
%      - OCP_Anode_13*tanh((theta-OCP_Anode_14)/OCP_Anode_15) ...
%      - OCP_Anode_16*tanh((theta-OCP_Anode_17)/OCP_Anode_18) ...
%      - OCP_Anode_19*tanh((theta-OCP_Anode_20)/OCP_Anode_21) ...
%      - OCP_Anode_22*tanh((theta-OCP_Anode_23)/OCP_Anode_24) ...
%      + OCP_Anode_25*tanh((theta-OCP_Anode_26)/OCP_Anode_27);
% 
% % Gradient of OCP wrt theta
% if(nargout >= 2)
% 
% %     % Polynomial Fit
% %     dUref = ppvalFast(p.dUppn,theta);
% %     varargout{1} = dUref / p.c_s_n_max;
% 
% dUref = -OCP_Anode_2*(OCP_Anode_3/p.c_s_n_max)*exp(-OCP_Anode_3*theta)  ...
%  +(OCP_Anode_4/(OCP_Anode_6*p.c_s_n_max))*((cosh((theta-OCP_Anode_5)/OCP_Anode_6)).^(-2)) ...
%  -(OCP_Anode_7/(p.c_s_n_max*OCP_Anode_9))*((cosh((theta-OCP_Anode_8)/OCP_Anode_9)).^(-2)) ...
%  -(OCP_Anode_10/(p.c_s_n_max*OCP_Anode_12))*((cosh((theta-OCP_Anode_11)/OCP_Anode_12)).^(-2)) ...
%  -(OCP_Anode_13/(p.c_s_n_max*OCP_Anode_15))*((cosh((theta-OCP_Anode_14)/OCP_Anode_15)).^(-2)) ...
%  -(OCP_Anode_16/(p.c_s_n_max*OCP_Anode_18))*((cosh((theta-OCP_Anode_17)/OCP_Anode_18)).^(-2)) ...
%  -(OCP_Anode_19/(p.c_s_n_max*OCP_Anode_21))*((cosh((theta-OCP_Anode_20)/OCP_Anode_21)).^(-2)) ...
%  -(OCP_Anode_22/(p.c_s_n_max*OCP_Anode_24))*((cosh((theta-OCP_Anode_23)/OCP_Anode_24)).^(-2)) ...
%  +(OCP_Anode_25/(p.c_s_n_max*OCP_Anode_27))*((cosh((theta-OCP_Anode_26)/OCP_Anode_27)).^(-2));
% varargout{1} = dUref;
% 
% end
% 
% % Gradient of OCP wrt temperature
% if(nargout >= 3)
%     
%     dUdT = 0;
%     varargout{2} = dUdT;
%     
end
% 
