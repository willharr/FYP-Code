%% Nozzle Configuration and Properties
% Will Harradence
% Imperial Aeronautics 2019/20

%% Setup
clear
clc
close all

%% Inputs

type = 'conical';

Pcc = 35e5; %combustion chamber pressure
gamma = 1.1726; %specific heat ratio
Patm = 101325; %ambient (exit) pressure
Rcc = 0.020; %combustion chamber diameter
m_dot = 0.150; %kg/s, total mass flow rate
m_molar = 30.53e-3; %kg/mol of reaction products
R_gas = 8.314/m_molar; %specific gas constant for exhaust gases
Tcc = 3252;
c_star = 1400;
g = 9.81;

%% General Properties

exp_ratio = ( ((2/(gamma+1))^(1/(gamma-1))) * ((Pcc/Patm)^(1/gamma)) )*(sqrt( ((gamma+1)/(gamma-1)) * (1-((Patm/Pcc)^((gamma-1)/gamma))) )^-1); %expansion ratio, Ae/At
Pthroat = Pcc*(2/(gamma+1))^(gamma/(gamma-1)); %pressure at throat

Acc = pi*Rcc^2;
At = m_dot*(sqrt(Tcc)/Pcc)*sqrt(R_gas/gamma)*((gamma+1)/2)^((gamma+1)/(2*(gamma-1)));
Ae = exp_ratio*At;

Rt = sqrt(At/pi);
Re = sqrt(Ae/pi);
Rcc = sqrt(Acc/pi);
R = Rt;


%% Conical Nozzle - Radiused
% Generates cross sectional plot of nozzle contour and array of nozzle
% areas with respect to x
if isequal(type,'conical')
    
alpha = 15; %deg, nozzle diverging half-angle
beta = 30; %deg, nozzle converging half-angle

x_radIn = -R*sind(beta);
x_radOut = R*sind(alpha);
x_inlet = x_radIn - (Rcc-R*(1-cosd(beta))-Rt)/tand(beta);
x_outlet = x_radOut + (Re-R*(1-cosd(alpha))-Rt)/tand(alpha);

n = 100;
x_arr1 = linspace(x_inlet,x_radIn,n);
x_arr2 = linspace(x_radIn,x_radOut,n);
x_arr3 = linspace(x_radOut,x_outlet,n);
x_array = [x_arr1 x_arr2 x_arr3];

y_con = zeros(1,3*n);
A_con = zeros(1,3*n);

for i = 1:n
    
    y_con(i) = tand(180-beta)*(x_arr1(i)-x_inlet) + Rcc;
    y_con(i+n) = -sqrt((R^2) - (x_arr2(i)^2)) + Rt+R;
    y_con(i+2*n) = tand(alpha)*(x_arr3(i)-x_radOut) + Rt+R*(1-cosd(alpha));

end

for i = 1:3*n

    
    A_con(i) = pi.*y_con(i).^2;
    
    if x_array(i) <= 0
        [M_con(i), T_con(i), P_con(i), rho_con(i) ,~] = flowisentropic(gamma,A_con(i)/At,'sub');
    else
        [M_con(i), T_con(i), P_con(i), rho_con(i) ,~] = flowisentropic(gamma,A_con(i)/At,'sup');
    end
    
    V_con(i) = M_con(i)*sqrt(gamma*R_gas*T_con(i));
    
end

end
%% Aerospike nozzle
%Primary resource: (1) Rei Masuda. Optimal design of annular aerospike engine nozzle. Ann Arbor: California State University, Long Beach; 2002

if isequal(type,'spike')
    
tol = 0.0001;

A_throat_check = (c_star*m_dot)/(Pcc*g);

if abs(A_throat_check - At) > tol
    
    disp('Error: Throat area mismatch. Consistency of analysis cannot be guaranteed')
end

%A_throat = A_throat_check;
A_throat = At;

M_exit = sqrt((((Pthroat/Patm)^((gamma-1)/gamma))-1)*(2/(gamma-1))); %exit mach number of the flow

M_spike = linspace(1,M_exit,n);

theta_throat = sqrt((gamma+1)/(gamma-1))*atan(sqrt(((gamma-1)/(gamma+1))*(M_exit^2 - 1))) - atan(sqrt(M_exit^2 - 1)); %effectively nu(M_exit), the PM Angle for M_exit

%A_exit = A_throat * (M_exit*((2/(gamma+1))*(1+((gamma-1)/2)*(M_exit^2)))^((-gamma+1)/(2*(gamma-1))))^(-1); %effective exit area of spike
A_exit = Ae;

r_exit = sqrt(A_exit/pi);

epsilon = A_exit/A_throat;

for i = 1 : n
    
    nu_spike(i) = sqrt((gamma+1)/(gamma-1))*atan(sqrt(((gamma-1)/(gamma+1))*(M_spike(i)^2 - 1))) - atan(sqrt(M_spike(i)^2 - 1)); %general PM angle equation

    mu_spike(i) = asin(1/M_spike(i)); %mach angle equation
    
    alpha_spike(i) = theta_throat - nu_spike(i) + mu_spike(i); %arbitraty relation for clarity
    
    beta_spike(i) = (1/M_spike(i))*((2/gamma+1)*(1+(gamma-1/2)*(M_spike(i)^2)))^((gamma-1)/2*(gamma-1)); %equal to A/At, inverse isentropic area ratio
    
    r_spike(i) = r_exit*sqrt(1 - ((beta_spike(i)/epsilon)*cos(theta_throat-nu_spike(i))));

    l(i) = (r_exit-r_spike(i))/sin(alpha_spike(i));
end

%l_tot = (r_exit-r_spike(n)/sin(alpha_spike(n)))*cos(alpha_spike(n));

x = linspace(0,l(n),n);

figure
plot(x,r_spike)
grid on
axis('image')

end
%% Heat Transfer Coefficient
%gravity
g = 9.81;
g_i = convvel(9.81,'m/s','ft/s');

%nozzle gas properties

molar_mass = 0.032;%kg/mol, from CEA mass fractions of products
molar_mass_i = convmass(molar_mass,'kg','lbm');

Dia_i = convlength(2*Rt,'m','in');
r_throat = convlength(R,'m','in'); %inches

Pr = (4*gamma)/(9*gamma-5); %[1]
r = Pr^0.33; %[1]

T_ns = Tcc; %K
T_ns_i = convtemp(T_ns,'K','R');

mu_i = (46.6*10e-10)*molar_mass_i^(0.5)*T_ns_i^0.6; %[1]

P_i = convpres(Pcc,'Pa','psi');

c_star = 1406; %m/s
c_star_i = convvel(c_star,'m/s','ft/s');

Cp_g = 1727; %J/kg*K
Cp_g_i = 4.125e-1; %BTU/lb*R

T_w_g = 500;

for i = 1:3*n
    
    T_aw_g(i) = T_ns*((1+r*((gamma-1)/2)*(M_con(i)^2))/(1+((gamma-1)/2)*(M_con(i)^2)));
    sigma_i(i) = (((0.5*(T_w_g/T_ns_i)*(1+((gamma-1)/2)*(M_con(i)^2))+0.5)^(0.68))*(1+((gamma-1)/2)*M_con(i)^2)^0.12)^(-1); %Bartz
    h_g_i(i) = ((0.026/(Dia_i^0.2))*(((mu_i^0.2)*Cp_g_i)/(Pr^0.6))*((g_i*P_i)/c_star_i)^0.8)*((Dia_i/r_throat)^0.1)*((At/A_con(i))^0.9)*sigma_i(i); %Bartz - CHECK constant
    h_g(i) = h_g_i(i)*2.941e6;
    
    h_g_simple(i) = (rho_con(i)*V_con(i))^0.8; % equation 4.11 in Huzel and Huang, very simple correlation
    
end

%% Plotting

figure % Conical nozzle contour
plot(x_array,y_con,'k',[x_array(1) x_array(end)],[0 0],'-.b')
grid on
axis('image')

figure % nozzle flow properties with x
plot(x_array,M_con,'k',x_array,T_con,'r',x_array,P_con,'b');
legend('M','T','P')
grid on

figure % nozzle heat transfer coeff with x
plot(x_array,h_g)
grid on

figure
plot(x_array,h_g_simple)
grid on


