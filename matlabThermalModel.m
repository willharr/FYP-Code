%Thermal Conduction Model with MATLAB
% Will Harradence
% Imperial Aeronautics 2019/20

close all
clear
clc

%% Material Properties
% 
% %nozzle gas properties
% gamma = 1.13; %CEA
% mach = 1;
% 
% molar_mass = 0.032;%kg/mol, from CEA mass fractions of products
% molar_mass_i = convmass(molar_mass,'kg','lbm');
% 
% Dia = 0.05;
% Dia_i = convlength(Dia,'m','in');
% 
% r_throat = 0.15;
% 
% Pr = (4*gamma)/(9*gamma-5); %[1]
% r = Pr^0.33; %[1]
% 
% T_ns = 3300; %K
% T_ns_i = convtemp(T_ns,'K','F');
% 
% mu_i = (46.6*10e-10)*molar_mass_i^(0.5)*convtemp(T_ns,'K','F')^0.6; %[1]
% P_throat = 16.7e5; %Pa
% P_throat_i = convpres(P_throat,'Pa','psi');
% 
% c_star = 1406; %m/s
% c_star_i = convvel(c_star,'m/s','ft/s');
% 
% Cp_g = 1727; %J/kg*K
% Cp_g_i = 4.125e-1; %BTU/lb*R

%wall material properties - COPPER
C_w = 376; %J/kg*K 
k = 413; %conduction constant

%Nitrous properties
C_nos = 72; %J/mol*K
molar_nos = 0.044; %kg/mol
C_nos = C_nos/0.044; %J/kg*K

%nitrous temp
T_c = 300;

%gravity
g = 9.81;
g_i = convvel(9.81,'m/s','ft/s');

%Wall conditions
T_w_g = convtemp(800,'K','R');


%% Modelling
model = createpde('thermal','transient');

geom = decsg([3,4,0,0,0.1,0.1,0,0.05,0.05,0]');

geometryFromEdges(model,geom);

thermalProperties(model,'ThermalConductivity',k,'MassDensity',8960,'SpecificHeat',C_w);

pdegplot(model,'EdgeLabels','on')

coolVal = @coolantFlux;
nozzleVal = @nozzleFlux;
thermalBC(model,'Edge',2,'HeatFlux',coolVal);
thermalBC(model,'Edge',3,'HeatFlux',0);
thermalBC(model,'Edge',1,'HeatFlux',0);
thermalBC(model,'Edge',4,'HeatFlux',nozzleVal);

thermalIC(model,500);

tlist = 0:0.001:5;

generateMesh(model);

results = solve(model,tlist);

T = results.Temperature;

[qx,qy] = evaluateHeatFlux(results);

figure(2)
pdeplot(model,'XYData',T(:,end),'Contour','on',...
                     'FlowData',[qx(:,end),qy(:,end)],'ColorMap','hot')
                 
%% Nozzle Side Heat Transfer
function q_dot = nozzleFlux(~,state)
    
    %gravity
    g = 9.81;
    g_i = convvel(9.81,'m/s','ft/s');
    
    %nozzle gas properties
    gamma = 1.13; %CEA
    mach = 1;

    molar_mass = 0.032;%kg/mol, from CEA mass fractions of products
    molar_mass_i = convmass(molar_mass,'kg','lbm');

    Dia = 0.05;
    Dia_i = convlength(Dia,'m','in');
    r_throat = 1; %inches

    Pr = (4*gamma)/(9*gamma-5); %[1]
    r = Pr^0.33; %[1]

    T_ns = 3300; %K
    T_ns_i = convtemp(T_ns,'K','R');

    mu_i = (46.6*10e-10)*molar_mass_i^(0.5)*T_ns_i^0.6; %[1]
    P_throat = 16.7e5; %Pa
    P_throat_i = convpres(P_throat,'Pa','psi');

    c_star = 1406; %m/s
    c_star_i = convvel(c_star,'m/s','ft/s');

    Cp_g = 1727; %J/kg*K
    Cp_g_i = 4.125e-1; %BTU/lb*R

    T_aw_g = T_ns*((1+r*((gamma-1)/2)*(mach^2))/(1+((gamma-1)/2)*(mach^2)));

    %the following 2 lines are in imperial units (due to empirical
    %constants)
    T_w_g = convtemp(state.u,'K','R');
    sigma_i = (((0.5*(T_w_g/T_ns_i)*(1+((gamma-1)/2)*(mach^2))+0.5)^(0.68))*(1+((gamma-1)/2)*mach^2)^0.12)^(-1); %Bartz
    h_g_i = ((0.026/(Dia_i^0.2))*(((mu_i^0.2)*Cp_g_i)/(Pr^0.6))*((g_i*P_throat_i)/c_star_i)^0.8)*((Dia_i/r_throat)^0.1)*(4.5^0.9)*sigma_i; %Bartz - CHECK constant

    h_g = h_g_i/2.941e6; %convert heat transfer coeff to SI from BTU/in^2*s*F CHECK
    q_dot = h_g*(T_aw_g-state.u);
    
end
%% Coolant side heat transfer

function q_dot = coolantFlux(~,state)
T_aw_c = 300; %K, low speed flow
T_c = 300;
T_c_i = convtemp(300,'K','R');

C_nos = 72; %J/mol*K
molar_nos = 0.044; %kg/mol
C_nos = C_nos/0.044; %J/kg*K
Cp_nos_i = C_nos/4192.11; %BTU/lbm*R
k_nos = 0.0915;

mu_nos_i = (46.6e-10)*(molar_nos^0.5)*(T_c^0.6);%lbm/in*s [2]
mu_nos = mu_nos_i/17.83;

Pr = (mu_nos*C_nos)/k_nos;

d = 0.2; %in, coolant channel diameter

mdot = 0.1; %kg/s
mdot_i = convmass(mdot,'kg','lbm');
G_i = mdot_i/(pi*d^2*4); %lbm/in^2*s

h_c_i = ((0.029*Cp_nos_i*(mu_nos_i^0.2))/(Pr^(2/3)))*((G_i^0.8)/(d^0.2))*(T_c_i/convtemp(state.u,'K','R'))^0.55;

h_c = h_c_i/(2.941e6); 
q_dot = h_c*(T_aw_c-state.u);
end

%{
References

[1]Bartz
[2]Huzel and Huang
%}