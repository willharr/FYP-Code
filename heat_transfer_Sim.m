%% FYP Heat Transfer Sim
% Will Harradence
% Imperial Aeronautics 2019/20

clear 
clc
close all

%% Grid Setup

L = 0.1;
W = 0.05;

N = 20;
M = 20;

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
molar_mass = 0.032; %kg/mol, from CEA mass fractions of products
Dia = 0.05;
r_throat = 1;

Pr = (4*gamma)/(9*gamma-5); %[1]
r = Pr^0.33; %[1]
T_ns = 3300; %K
mu = (46.6*10e-10)*molar_mass^(0.5)*convtemp(T_ns,'K','R')^0.6; %[1]
P_throat = 16.7e5; %bar
c_star = 1406; %m/s
Cp_g = 1727 %J/kg*K

%wall material properties - COPPER
C_w = 376; %J/kg*K 
k = 5; %conduction constant

%% Timestep
tfin = 1;

dtmax = 0.25*dx*dy/k; 
%dt = input('Input timestep: ');
dt = dtmax/2;

C = k*dt/(dx*dy);

%% Intial Conditions

u(1,:) = 200;
u(M,:) = 200;

%% Calculation Loop
t = 0;
hold on
while t<tfin
    
    t = t+dt
    
    
    T_aw = T_ns*((1+r*((gamma-1)/2)*(M^2))/(1+((gamma-1)/2)*(M^2)));
    sigma = (((0.5*(u(1,N/2)/T_ns)*(1+((gamma-1)/2)*mach^2)+0.5)^(0.68))*(1+((gamma-1)/2)*mach^2)^0.12)^(-1);
    h = ((0.026/Dia^(0.2))*(((mu^0.2)*Cp_g)/Pr^0.6)*((P_throat/c_star)^0.8)*(Dia/r_throat)^0.1)*(4.5^0.9)*sigma;
    
    q = h*(T_aw-u(1,N/2));
    u(1,:) = q*C_w;
    
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

