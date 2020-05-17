% Thermal Conduction Model with variable boundary conditions
% Will Harradence
% Imperial Aeronautics 2019/20
% FYP 

%% Setup and Options
close all
clear
clc

load n2o %loads nitrous virial equations data file - OBSOLETE WITH COOLPROP

animate = 0; % option to plot animated change in temperature

%% Material Properties

%wall material properties - COPPER
C_w = 376; %J/kg*K 
k = 413; %conduction constant

%Nitrous properties
C_nos = heatCapLiquid(n2o); %J/kg*K
molar_nos = 0.044; %kg/mol

%nitrous saturation temp
T_c = 300;
n2o.T = T_c;

%gravity
g = 9.81;
g_i = convvel(9.81,'m/s','ft/s');

%% Environmental Conditions

T_amb = 293;
P_amb = 101535;

%% Nozzle Flow Properties
%gets variation in flow properties through a supersonic con-di nozzle with
%axial distance, given a set of input parameters and assuming a conical
%expansion section. Combustion properties are given by NASA tool CEAOnline

nozzle.Pcc = 35e5; %combustion chamber pressure
nozzle.gamma = 1.1726; %specific heat ratio
nozzle.Patm = 101325; %ambient (exit) pressure
nozzle.Rcc = 0.020; %combustion chamber diameter
nozzle.m_dot = 0.150; %kg/s, total mass flow rate
nozzle.m_molar = 30.53e-3; %kg/mol of reaction products
nozzle.Tcc = 3252; %combustion chamber stagnation temperature
nozzle.c_star = 1400; %characteristic velocity
nozzle.g = 9.81; %gravity

nozzle = nozzleGeometry(nozzle);
nozzle.M_array = unique(nozzle.M_array); %variation of Mach with x
nozzle.x_array = unique(nozzle.x_array); %array of x positions
nozzle.A_array = unique(nozzle.A_array); %variation of area with x
nozzle.y_array = unique(nozzle.y_array,'stable'); %variation of radius with x
save nozzle 

%% Geometry Definition
% Currently only supports rectangular 2D wall geometry, representing a full
% wall cross section at an arbitrary radial position

geom_w = 0.002; %through wall thickness, [m]
geom_x = [nozzle.x_array(end) nozzle.x_array(end) nozzle.x_array(1) nozzle.x_array(1)]; %list of x_coordinates of verticies
geom_y = [geom_w/2 -geom_w/2 -geom_w/2 geom_w/2]; %list of y_coordinates of verticies
geom_circ = [1 0 nozzle.Rt+nozzle.R nozzle.R];
geom_mat = [3 4 geom_x geom_y]'; %see below
%{
formats geometry matrix for the decsg function, specifying:
3 - rectangular
4 - number of edges
geom_x - list of x-coordinates of verticies
geom_y - list of y-coordinates of verticies
%}

%Coolant Channel Geometry
channel_h = 0.005; %height (radial)
channel_d = 0.005; %depth (circumferential)
channel_A = channel_h*channel_d; %area

%% Modelling

model = createpde('thermal','transient'); %generate pde object

geom = decsg(geom_mat); %generate geometry
geometryFromEdges(model,geom); %add geometry to model

thermalProperties(model,'ThermalConductivity',k,'MassDensity',8960,'SpecificHeat',C_w); %apply thermal properties to model, in this case a copper wall

pdegplot(model,'EdgeLabels','on') %plot geometry with labelled edges

%set BC function handles - probably not needed, could just put handles
%straight into BC function calls
coolVal = @coolantFlux; 
nozzleVal = @nozzleFlux;

% Set boundary conditions

%coolant side
%thermalBC(model,'Edge',4,'HeatFlux',coolVal); %supercritical coolant
thermalBC(model,'Edge',4,'Temperature',T_c); %nucleate boiling coolant

%side walls - insulated
thermalBC(model,'Edge',3,'HeatFlux',0);
thermalBC(model,'Edge',1,'HeatFlux',0);

%hot gas side
thermalBC(model,'Edge',2,'HeatFlux',nozzleVal);

%Set initial condition
thermalIC(model,500); %500K initial guess 

%Set timesteps (for plotting/data collection only, not related to
%simulation)
t_fin = 2;
t_step = 0.001;
n_steps = t_fin/t_step;
tlist = 0:t_step:t_fin; %tells the solver at what times you want to collect data

generateMesh(model); %generate mesh

results = solve(model,tlist); %solve

T_model = results.Temperature; %extract temperature data

[n_mesh, ~] = size(T_model);

%% Coolant Evaluation
% NB using beta for volumetric flow quality to avoid confusion with heat,
% Q, however CoolProp uses Q for quality.

[~, q_dot_c] = evaluateHeatFlux(results,nozzle.x_array,ones(1,length(nozzle.x_array))*(geom_w/2),tlist(end));
Q_c = trapz(nozzle.x_array,(q_dot_c'.*nozzle.y_array*pi*2));

%tank conditions
T_T = T_amb; %tank temp
P_T = py.CoolProp.CoolProp.PropsSI('P','T',T_amb,'Q',0,'N2O'); %tank pressure
h_T = py.CoolProp.CoolProp.PropsSI('H','P',P_T,'Q',0,'N2O'); %tank enthalpy
beta_T = 0;

%inlet conditions
lineLoss = 10e5; %10 bar pressure loss in feed system

h_1 = h_T; %ideal line losses are adiabatic
P_1 = P_T - lineLoss; %pressure at entrance to coolant loop
T_1 = py.CoolProp.CoolProp.PropsSI('T','P',P_1,'H',h_T,'N2O');
beta_1 = py.CoolProp.CoolProp.PropsSI('Q','P',P_1,'H',h_T,'N2O');

%cooling segment
h_vap_1 = coolPropNos('H','P',P_1,'Q',1)-coolPropNos('H','P',P_1,'Q',0);
quality_limit = 0.5;
cool_cap = (quality_limit*h_vap_1)+coolPropNos('H','P',P_1,'Q',0) - h_1;
m_dot_cool = Q_c/cool_cap;







%% Plotting 
%------------------------------------------------
figure
pdeplot(model,'XYData',T_model(:,end),'Contour','on','ColorMap','hot')
daspect([5 1 1])
hold on
plot(nozzle.x_array,(nozzle.y_array/10)-3e-3,'k',[nozzle.x_array(1) nozzle.x_array(end)],[-3e-3 -3e-3],'-.b')
xlabel('Axial Position from Throat [m]')
ylabel('Radial Position from Centerline [m]')
axis('tight')
set(gca,'FontSize',14)
grid on
box on
              
figure %convergence check
hold on 
plot(T_model(floor(n_mesh/2),:),'b');
plot(T_model(floor(n_mesh/4),:),'r');
plot(T_model(floor((3*n_mesh)/4),:),'k');
hold off
grid on
box on
xlabel('Time [ms]')
ylabel('Temperature')
set(gca,'FontSize',14)

[qx,qy] = evaluateHeatFlux(results);
figure
pdeplot(model,'XYData',results.Temperature(:,end), ...
                     'Contour','on',...
                     'FlowData',[qx(:,end),qy(:,end)], ...
                     'ColorMap','hot')

if animate == 1
for i = 1:t_fin/t_step
    figure(3)
    pdeplot(model,'XYData',T_model(:,i),'Contour','on',...
                         'FlowData',[qx_history(:,end),qy_history(:,end)],'ColorMap','hot')
    pause(0.05)
    i
end   
end

%% Nozzle Side Heat Transfer
function q_dot = nozzleFlux(location,state)
    %% Description
    %{
        This function generates a heat flux for a given position in the
        nozzle and wall temperature at that position, using the Bartz
        method for rocket nozzle heat transfer.
 
        Properties with the appendix "_i" are in imperial units, which are
        given in the comments.
    %}
    
    %% Get nozzle properties from saved file
    % The matlab pde solver will not accept boundary condition functions
    % with more input arguments than [location, state], so in order to pass
    % information from the main script to this function, it is saved in the
    % main script and loaded here. 
    
    load nozzle.mat

    g_i = convvel(nozzle.g,'m/s','ft/s'); %gravity
    
    %% Nozzle and hot-side gas properties
    gamma = nozzle.gamma; %CEA
    mach = interp1(nozzle.x_array,nozzle.M_array,location.x); %gets local mach number from nozzle calculation and mesh position
    
    m_molar_i = convmass(nozzle.m_molar,'kg','lbm'); %lbm/mol, from CEA mass fractions of products

    Dt = nozzle.Rt*2; %diameter of throat
    Dt_i = convlength(Dt,'m','in');
    R = convlength(nozzle.R,'m','in'); %inches, axial radius of throat section
    
    A_local = interp1(nozzle.x_array,nozzle.A_array,location.x);
    A_ratio = nozzle.At/A_local;
    
    Pr = (4*gamma)/(9*gamma-5); %[1] Prandtl number
    recovery = Pr^0.33; %[1], recovery factor
    
    Tcc = nozzle.Tcc;
    Tcc_i = convtemp(nozzle.Tcc,'K','R');% Rankine - stagnation temperature, assumed to be equal to combustion chamber temp/flame temp

    mu_i = (46.6*10e-10)*m_molar_i^(0.5)*Tcc_i^0.6; %[1] 

    Pcc_i = convpres(nozzle.Pcc,'Pa','psi');

    c_star_i = convvel(nozzle.c_star,'m/s','ft/s'); % ft/s, characteristic velocity

    Cp_g = 1727; %J/kg*K
    Cp_g_i = 4.125e-1; %BTU/lb*R

    T_aw_n = Tcc*((1+recovery*((gamma-1)/2)*(mach^2))/(1+((gamma-1)/2)*(mach^2))); % K - adiabatic wall temp
    
    %% Get nozzle-side wall temp from PDE
    T_w_n = convtemp(state.u,'K','R'); %Rankine - nozzle-side wall temp
    
    %% Calculate heat transfer coefficient
    sigma_i = (((0.5*(T_w_n/Tcc_i)*(1+((gamma-1)/2)*(mach^2))+0.5)^(0.68))*(1+((gamma-1)/2)*mach^2)^0.12)^(-1); %Bartz
    h_g_i = ((0.026/(Dt_i^0.2))*(((mu_i^0.2)*Cp_g_i)/(Pr^0.6))*((g_i*Pcc_i)/c_star_i)^0.8)*((Dt_i/R)^0.1)*(A_ratio^0.9)*sigma_i; %Bartz - CHECK constant
    
    %% Convert to metric
    h_g = h_g_i*2.941e6; %convert heat transfer coeff to SI from BTU/in^2*s*F CHECK
    
    %% Calculate heat flux
    q_dot = h_g*(T_aw_n-state.u);
    
end

%% Coolant side heat transfer

function q_dot = coolantFlux(~,state)
%% Description
%{
    This coolant side heat flux script uses and equation from Huzel and
    Huang to determine the heat flux to a coolant fluid when the coolant is
    supercritical. 
%}

T_aw_c = 300; %K, adiabatic wall temp for relatively low speed flow
T_c = 300; % coolant temperature
T_c_i = convtemp(300,'K','R');

C_nos = 72; %J/mol*K - Specific heat capacity of coolant, nitrous oxide in this case
molar_nos = 0.044; %kg/mol - molar mass of nitrous oxide
C_nos = C_nos/0.044; %J/kg*K - convert to per-mass value
Cp_nos_i = C_nos/4192.11; %BTU/lbm*R
k_nos = 0.0915; % Thermal conductivity of nitrous

mu_nos_i = (46.6e-10)*(molar_nos^0.5)*(T_c^0.6); %lbm/in*s, dynamic viscosity of nitrous [2]
mu_nos = mu_nos_i/17.83;

Pr = (mu_nos*C_nos)/k_nos;

d = 0.2; %in, coolant channel diameter

mdot = 0.1; %kg/s
mdot_i = convmass(mdot,'kg','lbm');
G_i = mdot_i/(pi*d^2*4); %lbm/in^2*s

h_c_i = ((0.029*Cp_nos_i*(mu_nos_i^0.2))/(Pr^(2/3)))*((G_i^0.8)/(d^0.2))*(T_c_i/convtemp(state.u,'K','R'))^0.55;

h_c = h_c_i*(2.941e6); 

q_dot = h_c*(T_aw_c-state.u);
end

function out = coolPropNos(out_type,in1_type,in1,in2_type,in2)

    out = py.CoolProp.CoolProp.PropsSI(out_type,in1_type,in1,in2_type,in2,'N2O');
    
end

%{
References

[1]Bartz (1959)
[2]Huzel and Huang (1992)

%}