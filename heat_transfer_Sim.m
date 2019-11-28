%% FYP Heat Transfer Sim
% Will Harradence
% Imperial Aeronautics 2019/20

%% Notes
% *_i denotes a value in Imperial units, generally BTU/ft/lbm unless
% otherwise specified. All other variables use SI units (J/m/kg etc.)

clear 
clc
close all

%% Grid Setup

L = 1; %m
W = 1; %m

N = 20;
M = 22;

dx = L/(N);
dy = W/(M);

x1 = (0:(N-1))*dx; %horizontal vector of x positions from 0-L
y1 = (0:(M-1))'*dy; %vetical vector of y positions from 0-W

%plotting arrays
x = repmat(x1,M,1); %list of x-pos for each y-pos
y = repmat(y1,1,N); %list of y-pos for each x-pos

%temp array
temp_init = 500;

u = temp_init*ones(M,N);
du = u.*0;


%% Material Properties

%nozzle gas properties
gamma = 1.13; %CEA
mach = 1;

molar_mass = 0.032;%kg/mol, from CEA mass fractions of products
molar_mass_i = convmass(molar_mass,'kg','lbm');

Dia = 0.05;
Dia_i = convlength(Dia,'m','in');

r_throat = 1;

Pr = (4*gamma)/(9*gamma-5); %[1]
r = Pr^0.33; %[1]

T_ns = 3300; %K
T_ns_i = convtemp(T_ns,'K','F');

mu_i = (46.6*10e-10)*molar_mass_i^(0.5)*convtemp(T_ns,'K','F')^0.6; %[1]
P_throat = 16.7e5; %Pa
P_throat_i = convpres(P_throat,'Pa','psi');

c_star = 1406; %m/s
c_star_i = convvel(c_star,'m/s','ft/s');

Cp_g = 1727; %J/kg*K
Cp_g_i = 4.125e-1; %BTU/lb*R

%wall material properties - COPPER
C_w = 376; %J/kg*K 
k = 413; %conduction constant

%Nitrous properties
C_nos = 72; %J/mol*K
molar_nos = 0.044; %kg/mol
C_nos = C_nos/0.044; %J/kg*K

%% Timestep
tfin = 1;

dtmax = 0.25*dx*dy/k;
disp(dtmax)
dt = input('Input timestep: ');
%dt = dtmax/2;

const = k*dt/dx;

%% Intial Conditions

u(1,:) = temp_init;
u(M,:) = temp_init;

%% Calculation Loop
t = 0;
while t<tfin
    
    t = t+dt
    
    % Conduction through wall
    for i = 2:M-1
        for j = 2:N-1
            
            du(i,j) = const*(u(i+1,j)+u(i-1,j)-2*u(i,j));
            
        end
    end
    
    du(1,:) = du(2,:);
    du(M,:) = du(M-1,:);
    
    % Nozzle Side Heat Transfer
    T_aw_g = T_ns*((1+r*((gamma-1)/2)*(mach^2))/(1+((gamma-1)/2)*(mach^2)));
    
    %the following 2 lines are in imperial units (due to empirical
    %constants)
    sigma_i = (((0.5*(u(1,N/2)/T_ns_i)*(1+((gamma-1)/2)*(mach^2))+0.5)^(0.68))*(1+((gamma-1)/2)*mach^2)^0.12)^(-1);
    h_g_i = ((0.026/Dia_i^(0.2))*(((mu_i^0.2)*Cp_g_i)/Pr^0.6)*((P_throat_i/c_star_i)^0.8)*(Dia_i/r_throat)^0.1)*(4.5^0.9)*sigma_i;
    
    h_g = h_g_i/2.044e4; %convert heat transfer coeff to SI
    q_g = h_g*(T_aw_g-u(1,1));
    Q_g = q_g*C_w*dx;
    
    du(1,:) = du(1,:) + Q_g ; %set nozzle side temp change w.r.t. heat transfer
    
    % Coolant side heat transfer
    T_aw_c = 300;
    h_c = (2.15e-8)*5;
    
    q_c = h_c*(T_aw_c-u(M,N/2));
    Q_c = q_c*C_nos*dx;
    
    du(M,:) = du(M,:) + Q_c;

    
    %Update temp matrix
    u = u + du;
    
    %Enforce boundary conditions
    %u(1,:) = 2500;
    %u(M,:) = temp_init;
    
    u(:,1) = u(:,2);
    u(:,N) = u(:,N-1);
    
    %plot
    figure(1)
    surf(x,y,u)
%     figure(2)
%     subplot(3,1,1)
%     bar([0,1],[du(1,1),du(M,1)])
%     xlabel('du')
%     subplot(3,1,2)
%     bar([0,1],[q_g,q_c])
%     xlabel('q')
%     subplot(3,1,3)
%     bar([0,1],[u(1,1),u(M,1)])
%     xlabel('u')
    drawnow

end

