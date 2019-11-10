%% FYP Heat Transfer Sim
% Will Harradence
% Imperial Aeronautics 2019/20

clear 
clc
close all

%% Grid Setup

L = 0.1;
W = 0.05;

N = 100;
M = 100;

dx = L/(N);
dy = W/(M);

x1 = (0:(N-1))*dx; %horizontal vector of x positions from 0-L
y1 = (0:(M-1))'*dy; %vetical vector of y positions from 0-W

%plotting arrays
x = repmat(x1,M,1); %list of x-pos for each y-pos
y = repmat(y1,1,N); %list of y-pos for each x-pos

%temp array
u = 300*ones(M,N);
du = u;


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

%% Timestep
tfin = 1;

dtmax = 0.25*dx*dy/k; 
dt = input('Input timestep: ');
%dt = dtmax/2;

C = k*dt/(dx*dy);

%% Intial Conditions

u(1,:) = 400;
u(M,:) = 400;

%% Calculation Loop
t = 0;

while t<tfin
    
    t = t+dt
    
    
    T_aw = T_ns_i*((1+r*((gamma-1)/2)*(M^2))/(1+((gamma-1)/2)*(M^2)));
    sigma = (((0.5*(u(1,N/2)/T_ns_i)*(1+((gamma-1)/2)*(mach^2))+0.5)^(0.68))*(1+((gamma-1)/2)*mach^2)^0.12)^(-1);
    h_i = ((0.026/Dia_i^(0.2))*(((mu_i^0.2)*Cp_g_i)/Pr^0.6)*((P_throat_i/c_star_i)^0.8)*(Dia_i/r_throat)^0.1)*(4.5^0.9)*sigma;
    
    h = h_i/2.044e4;
    q = h*(T_aw-u(1,N/2));
    u(1,:) = q*C_w*W^2;
    disp(u(1,1))
    
    for i = 2:M-1
        for j = 2:N-1
            
            du(i,j) = C*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j));
            
        end
    end
    
    u = u + du;

    %u(1,:) = 2500;
    u(M,:) = 200;
    
    u(:,1) = u(:,2);
    u(:,N) = u(:,N-1);
    
    surf(x,y,u)
    drawnow

end

