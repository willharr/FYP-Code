% Thermal Conduction Model with variable boundary conditions
% Will Harradence
% Imperial Aeronautics 2019/20
% FYP 

%% Setup and Options

load n2o %loads nitrous virial equations data file - OBSOLETE WITH COOLPROP
load looper 
load looper2
load materials
load materials_list
load wall_thickness
%{
    format of Materials table:

    Material (str)
    Melting Temp 1
    Max Service Temp 2
    Conductivity 3 
    Heat Capacity 4 
    Density 5
%}

animate = 0; % option to plot animated change in temperature
accurate_geom = 1;
plotting_general = 0;

%% Environmental Conditions

T_amb = 293;
P_amb = 101325;
%gravity
g = 9.81;
g_i = convvel(9.81,'m/s','ft/s');

%% Material Properties

%wall material properties - COPPER
C_w = materials(4,looper); %J/kg*K %heat capacity
k = materials(3,looper); %conduction constant
rho_w = materials(5,looper);

%Nitrous properties
C_nos = coolPropNos('C','T',T_amb,'Q',0); %J/kg*K
molar_nos = 0.044; %kg/mol


%% Line losses

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

%nitrous saturation temp
T_c = coolPropNos('Tcrit','T',T_1,'P',P_1);

%% Nozzle Flow Properties
%gets variation in flow properties through a supersonic con-di nozzle with
%axial distance, given a set of input parameters and assuming a conical
%expansion section. Combustion properties are given by NASA tool CEAOnline

nozzle.Pcc = 200e5; %combustion chamber pressure
nozzle.gamma = 1.1726; %specific heat ratio
nozzle.Patm = P_amb; %ambient (exit) pressure
nozzle.Rcc = 1.34; %combustion chamber diameter
nozzle.m_dot = 760.22; %kg/s, total mass flow rate
nozzle.m_molar = 30.53e-3; %kg/mol of reaction products
nozzle.Tcc = 3252; %combustion chamber stagnation temperature
nozzle.c_star = 1321; %characteristic velocity
nozzle.g = 9.81; %gravity

nozzle.t_wall = wall_thickness(looper2); %wall thickness

nozzle = nozzleGeometry(nozzle); %Call nozzle geometry function to get a conical nozzle shape

nozzle.M_array = unique(nozzle.M_array); %variation of Mach with x
nozzle.x_array = unique(nozzle.x_array); %array of x positions
nozzle.A_array = unique(nozzle.A_array); %variation of area with x
nozzle.y_array = unique(nozzle.y_array,'stable'); %variation of radius with x
save nozzle 

%% Geometry Definition
% Currently only supports rectangular 2D wall geometry, representing a full
% wall cross section at an arbitrary radial position

%Rectangular Geometry
geom_w = nozzle.t_wall; %through wall thickness, [m]
geom_x = [nozzle.x_array(end) nozzle.x_array(end) nozzle.x_array(1) nozzle.x_array(1)]; %list of x_coordinates of verticies
geom_y = [geom_w/2 -geom_w/2 -geom_w/2 geom_w/2]; %list of y_coordinates of verticies

%Accurate Nozzle Geometry - OPTIONAL

%Collation of  geometry matrix
geom_mat = [3 4 geom_x geom_y]'; %see below
%{
formats geometry matrix for the decsg function, specifying:
3 - rectangular
4 - number of edges
geom_x - list of x-coordinates of verticies
geom_y - list of y-coordinates of verticies
%}

%Coolant Channel Geometry - CURRENTLY NOT USED
channel_h = 0.005; %height (radial)
channel_d = 0.005; %depth (circumferential)
channel_A = channel_h*channel_d; %area

%% Modelling

model = createpde('thermal','transient'); %generate pde object

if accurate_geom == 0 %Option to form 2D model as rectantgle
    geom = decsg(geom_mat); %generate geometry
elseif accurate_geom == 1 %Option to form 2D model as actual nozzle shape
    geom = decsg(nozzle.geom, nozzle.setFunc, nozzle.nameList);
end

geometryFromEdges(model,geom); %add geometry to model

thermalProperties(model,'ThermalConductivity',k,'MassDensity',rho_w,'SpecificHeat',C_w); %apply thermal properties to model, in this case a copper wall

%if edges_plot_presolve == 1
    % Edges plot
%     figure
%     pdegplot(model,'EdgeLabels','on') %plot geometry with labelled edge
%     axis('image')
%     hold on
%     plot(nozzle.outer_x_array,nozzle.outer_y_array,'xr')
%     plot(nozzle.x_array,nozzle.y_array,'xg')
%     hold off
%end


%set BC function handles - probably not needed, could just put handles
%straight into BC function calls
coolVal = @coolantFlux; 
nozzleVal = @nozzleFlux;

% Set boundary conditions

if accurate_geom == 0
    
    %coolant side
    %thermalBC(model,'Edge',4,'HeatFlux',coolVal); %supercritical coolant
    thermalBC(model,'Edge',4,'Temperature',T_c); %nucleate boiling coolant

    %side walls - insulated
    thermalBC(model,'Edge',3,'HeatFlux',0);
    thermalBC(model,'Edge',1,'HeatFlux',0);

    %hot gas side
    thermalBC(model,'Edge',2,'HeatFlux',nozzleVal);

elseif accurate_geom == 1
    
    %coolant side
    %thermalBC(model,'Edge',4,'HeatFlux',coolVal); %supercritical coolant
    thermalBC(model,'Edge',7,'Temperature',T_c); %nucleate boiling coolant
    thermalBC(model,'Edge',8,'Temperature',T_c);
    thermalBC(model,'Edge',14,'Temperature',T_c);
    thermalBC(model,'Edge',9,'Temperature',T_c);
    thermalBC(model,'Edge',10,'Temperature',T_c);

    %side walls - insulated
    thermalBC(model,'Edge',2,'HeatFlux',0);
    thermalBC(model,'Edge',3,'HeatFlux',0);

    %hot gas side
    thermalBC(model,'Edge',1,'HeatFlux',nozzleVal);
    thermalBC(model,'Edge',12,'HeatFlux',nozzleVal);
    thermalBC(model,'Edge',4,'HeatFlux',nozzleVal);
    
    %internal - no bc: edges 5, 6, 11, 13 

end

%Set initial condition
thermalIC(model,500); %500K initial guess 

%Set timesteps (for plotting/data collection only, not related to
%simulation)
t_fin = 5;
t_step = 0.1;
n_steps = t_fin/t_step;
tlist = 0:t_step:t_fin; %tells the solver at what times you want to collect data

generateMesh(model); %generate mesh

results = solve(model,tlist); %solve

T_model = results.Temperature; %extract temperature data

[n_mesh, ~] = size(T_model);

%% Coolant Evaluation
% NB using beta for mass vapour quality to avoid confusion with heat (Q) or
% position (x) however CoolProp uses Q for quality. Beta is often used for
% VOLUMETRIC flow quality instead.

%Get heat flux results along the cooled surface in x and y
if accurate_geom == 0
    [q_dotx, q_doty] = evaluateHeatFlux(results,nozzle.x_array,ones(1,length(nozzle.x_array))*(geom_w/2),tlist(end)); %Gets 2D heat flux on coolant surface at steady state
elseif accurate_geom == 1
    [q_dotx, q_doty] = evaluateHeatFlux(results,nozzle.outer_x_array,nozzle.outer_y_array,tlist(end));
end

%pre-allocate net heat flux array
q_dot_c = zeros(1,length(q_dotx));

for j = 1:length(q_dotx)
    q_dot_c(j) = norm([q_dotx(j),q_doty(j)]); %get magnitude of heat flux
    
    if isnan(q_dot_c(j)) %checks for NAN in the array and sets as average of either side, will break if two nan's in a row
        if j == length(q_dotx) 
            q_dot_c(j) = q_dot_c(j-1);
        elseif j == 1
            q_dot_c(j) = q_dot_c(j+1);
        else
            q_dot_c(j) = mean([q_dot_c(j-1) q_dot_c(j+1)]);
        end
    end
    
end

%multiplies 2D heat flux by local nozzle circumference (assuming
%axisymmetric conditions) and integrates to get total heat rate across
%nozzle wall
if accurate_geom == 0
    Q_c = trapz(nozzle.x_array,((q_dot_c.*nozzle.y_array).*pi.*2)); %generates 3D (total surface) heat flux accounting for changing diameter of nozzle.
elseif accurate_geom == 1
    Q_c = trapz(nozzle.outer_x_array,((q_dot_c.*nozzle.outer_y_array).*pi.*2)); %generates 3D (total surface) heat flux accounting for changing diameter of nozzle.
end

%cooling segment
h_vap_1 = coolPropNos('H','P',P_1,'Q',1)-coolPropNos('H','P',P_1,'Q',0); %unit?
quality_limit = 0.7;
cool_cap = (quality_limit*h_vap_1)+coolPropNos('H','P',P_1,'Q',0) - h_1; %
m_dot_cool = Q_c/cool_cap;

P_2 = P_1; %isobaric heat addition - IDEALISATION
h_2 = coolPropNos('H','P',P_2,'Q',quality_limit);

%% Plotting 
%------------------------------------------------
if plotting_general == 1
    
    % Edges plot
    figure
    pdegplot(model,'EdgeLabels','on') %plot geometry with labelled edge
    axis('image')
    hold on
    plot(nozzle.outer_x_array,nozzle.outer_y_array,'xr')
    plot(nozzle.x_array,nozzle.y_array,'xg')
    hold off
    
    %Results Heat Map
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

    %results plot 
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

    %P-H Diagram----------------------------------------

    n = 1000;
    T = linspace(270,309.559,n-1);
    T = [T 309.56];
    width = 18;
    height = 15;

    P_arr = linspace(88000,7.245e+06,n);
    for i = 1:n

        H_L(i) = py.CoolProp.CoolProp.PropsSI('H','P',P_arr(i),'Q',0,'N2O');
        H_G(i) = py.CoolProp.CoolProp.PropsSI('H','P',P_arr(i),'Q',1,'N2O');

    end

    m = 100;
    H_arr = linspace(0,5e5,m);
    P_arr2 = linspace(88000,7.244e+06,m);

    T_HP = zeros(m,m);

    for i =1:m
       for j = 1:m
        T_HP(i,j) = py.CoolProp.CoolProp.PropsSI('T','H',H_arr(j),'P',P_arr2(i),'N2O');
       end
    end

    figure % P-H diagram
    plot(H_L/1e3,P_arr/1e5,'-b','LineWidth',3)
    hold on
    plot(H_G/1e3,P_arr/1e5,'-r','LineWidth',3)
    plot([h_T h_1]/1e3,[P_T P_1]/1e5,'--k')
    plot([h_1 h_2]/1e3,[P_1 P_2]/1e5,'--k')
    contour(H_arr/1e3,P_arr2/1e5,T_HP,[200 220 240 260 280 300 309.57 320 340],'-k','ShowText','on','LabelSpacing',800)
    hold off
    grid on
    box on
    axis([0,500,0,75])
    xlabel('Enthalpy [kJ/kgK]')
    ylabel('Pressure [bar]')
    legend('Saturated Liquid','Saturated Vapour','Isotherms','Cooling Loop','orientation','horizontal','Location','northoutside')
    set(gca,'FontSize',16)
    set(gcf,'units','centimeters','position',[5,5,width,height])

end

%% Data Collection
load m_dot_results
load q_dot_results 
load max_temp_results

max_temp = max(results.Temperature(:,end));

m_dot_results(looper,looper2) = m_dot_cool;
max_temp_results(looper,looper2) = max_temp;
q_dot_results(looper,looper2) = max(q_dot_c);

save('m_dot_results.mat','m_dot_results')
save('q_dot_results.mat','q_dot_results')
save('max_temp_results.mat','max_temp_results')



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