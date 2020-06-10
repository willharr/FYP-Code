% Thermal Conduction Model with variable boundary conditions
% Will Harradence
% Imperial Aeronautics 2019/20
% FYP 

%% Setup and Options
clear
clc
close all

load m_dot_resultsR1E1.mat
load m_dot_resultsR1E2.mat
load m_dot_resultsR1E3.mat
load m_dot_resultsR1E4.mat
load m_dot_resultsR1E5.mat

T_amb = 293;
P_amb = 101325;

%tank conditions
T_T = T_amb; %tank temp
P_T = coolPropNos('P','T',T_amb,'Q',0); %tank pressure
h_T = coolPropNos('H','P',P_T,'Q',0); %tank enthalpy
beta_T = 0; %quality of 0, all liquid in tank outflow

%Coolant Loop inlet conditions
lineLoss = 10e5; %10 bar pressure loss in feed system, ESTIMATE
h_1 = h_T; %ideal line losses are adiabatic
P_1 = P_T - lineLoss; %pressure at entrance to coolant loop
T_1 = coolPropNos('T','P',P_1,'H',h_T); %K
beta_1 = coolPropNos('Q','P',P_1,'H',h_T);

rhoL = coolPropNos('D','T',T_1,'Q',0);
rhoG = coolPropNos('D','T',T_1,'Q',1);
rho = coolPropNos('D','T',T_1,'Q',beta_1);
xi = beta_1;
x0 = 0.7;

alph0 = x0*(rho/rhoG);
alphi = xi*(rho/rhoL);

radius = [0.02 0.276 0.4085 0.8775 1.34];

%for j = 1:length(m_dot_resultsR1E1(:,1))

    %for i = 1:length(m_dot_resultsR1E1(1,:))

        n_channels = 1000;
        A = ((2*pi*radius(5))/(n_channels*0.5))^2;

        G = m_dot_resultsR1E5(6,1)/(n_channels*A);

        part1 = ((1-x0)^2)/(rhoL*(1-alph0)) + (x0^2)/(rhoG*alph0);
        part2 = ((1-xi)^2)/(rhoL*(1-alphi)) + (xi^2)/(rhoG*alphi);

        deltaP_R1E5 = (G^2)*(part1 - part2);

        %clear A part1 part2 G
    %end
%end

function out = coolPropNos(out_type,in1_type,in1,in2_type,in2)

    out = py.CoolProp.CoolProp.PropsSI(out_type,in1_type,in1,in2_type,in2,'N2O');
    
end
