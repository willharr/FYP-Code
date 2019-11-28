%% FYP Regenerative Tank Sim
% Will Harradence
% Imperial Aeronautics 2019/20

%% Initialisation

clear
clc
close all
load n2o.mat

% Tank dimensions - cylindrical + hemispherical endcaps
tank.radius = 0.075;
tank.length = 0.75;
tank.vol = tank.length*pi*tank.radius^2 + (4/6)*pi*tank.radius^3;
ullage = 0.15;

% Engine setup
Pcc = 10e5; %fixed for now

Ainj = 0.00001; %arbitrary
Ninj = 1; %arbitrary
Kinj = 2; %from ref [1]
D_loss = Kinj/((Ninj*Ainj)^2);

% Initial Conditions
temp_init = 300;
n2o.T = temp_init;

tank.m_l = (1-ullage)*tank.vol*densityLiquid(n2o);
tank.m_v = ullage*tank.vol*densityVap(n2o);
tank.m = tank.m_l+tank.m_v;

% Loop Conditions
timestep = 0.01;
i = 1;
dm_v = 0.001;

%% Calculation
while(tank.m_l > 0.01) 
    mdot
    
    dQ = dm_v*heatVaporisation(n2o);

    dT = -dQ/(tank.m_l*heatCapLiquid(n2o));

    n2o.T = n2o.T + dT;

    n2o.rho_l = densityLiquid(n2o);
    n2o.rho_v = densityVap(n2o);
    n2o.Pvap = vapPressure(n2o);
    
    if n2o.Pvap < Pcc
        disp('Tank pressure less than combustion chamber pressure!')
        break
    end
    
    m_dot = sqrt((2*n2o.rho_l*(n2o.Pvap - Pcc))/D_loss);

    tank.m_l_old = tank.m_l - m_dot*timestep;

    tank.m = tank.m - m_dot*timestep;

    tank.m_l = (tank.vol-(tank.m/n2o.rho_v))/((1/n2o.rho_l)-(1/n2o.rho_v));
    
    if tank.m_l_old < tank.m_l
        disp('Burnout!')
        break
    end

    dm_v = tank.m_l_old - tank.m_l;
    tank.m_v = tank.m_v + dm_v;
    
    store.tank_T(i) = n2o.T;
    store.tank_m(i) = tank.m;
    store.tank_ml(i) = tank.m_l;
    store.tank_mv(i) = tank.m_v;
    store.time(i) = i*timestep;
    
    i = i+1;

end

figure(1)
plot(store.time,store.tank_T,'r')
ylabel('temp')
grid on

figure(2)
plot(store.time,store.tank_m,'b',store.time,store.tank_ml,'r',store.time,store.tank_mv,'k')
ylabel('mass')
grid on




%% References

% [1] Aspire Space Technical Paper "Modelling the nitrous run tank emptying"
