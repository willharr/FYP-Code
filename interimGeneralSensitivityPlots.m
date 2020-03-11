%% FYP Engine Parameter Sensitivity Analysis
clear 
clc
close all

%% Defaults
% Based on A-1 sample engine from Huzel and Huang
R_def = 0.815;
gamma_def = 1.22;
T_def = 3588;
Pcc_def = 6.89e6; %bar
Pa_def = 101325;
exp_def = 14;
g = 9.81;
At = 0.314193;
n = 20;
Pe = ones(1,n)*101325;

%gamma loop: 1 - 2
R = R_def*ones(1,n);
T = T_def*ones(1,n);
pcc = Pcc_def*ones(1,n);
pa = Pa_def*ones(1,n);
exp = exp_def*ones(1,n);

gamma = linspace(gamma_def*0.5,gamma_def*1.5,n);
for i = 1:n  
    
    gamma_struct.c_star(i) = sqrt(g*gamma(i)*R(i)*T(i))/(gamma(i)*sqrt((2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1))));

    gamma_struct.Cf(i) = sqrt(((2*gamma(i)^2)/(gamma(i)-1))*((2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1)))*(1-(Pe(i)/pcc(i))^((gamma(i)-1)/gamma(i))))+exp(i)*((Pe(i)-pa(i))/pcc(i));

    gamma_struct.F(i) = pcc(i)*At*gamma_struct.Cf(i);

    gamma_struct.Isp(i) = (gamma_struct.c_star(i)*gamma_struct.Cf(i))/g;

end

% R loop: 100 - 400

%R = R_def*ones(1,100);
T = T_def*ones(1,n);
pcc = Pcc_def*ones(1,n);
pa = Pa_def*ones(1,n);
exp = exp_def*ones(1,n);
gamma = gamma_def*ones(1,n);

R = linspace(R_def*0.5,R_def*1.5,n);

for i = 1:n  
    
    R_struct.c_star(i) = sqrt(g*gamma(i)*R(i)*T(i)) / (gamma(i)*sqrt((2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1)) ));

    R_struct.Cf(i) = sqrt(((2*gamma(i)^2)/(gamma(i)-1))*((2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1)))*(1-(pa(i)/pcc(i))^((gamma(i)-1)/gamma(i))))+exp(i)*((Pe(i)-pa(i))/pcc(i));

    R_struct.F(i) = pcc(i)*At*R_struct.Cf(i);

    R_struct.Isp(i) = (R_struct.c_star(i)*R_struct.Cf(i))/g;

end

% T loop 2500 - 6500 k
R = R_def*ones(1,n);
%T = T_def*ones(1,100);
pcc = Pcc_def*ones(1,n);
pa = Pa_def*ones(1,n);
exp = exp_def*ones(1,n);
gamma = gamma_def*ones(1,n);

T = linspace(T_def*0.5,T_def*1.5,n);

for i = 1:n  
    
    T_struct.c_star(i) = sqrt(g*gamma(i)*R(i)*T(i))/(gamma(i)*sqrt((2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1))));

    T_struct.Cf(i) = sqrt(((2*gamma(i)^2)/(gamma(i)-1))*((2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1)))*(1-(Pe(i)/pcc(i))^((gamma(i)-1)/gamma(i))))+exp(i)*((Pe(i)-pa(i))/pcc(i));

    T_struct.F(i) = pcc(i)*At*T_struct.Cf(i);

    T_struct.Isp(i) = (T_struct.c_star(i)*T_struct.Cf(i))/g;

end

% Pcc loop 1 - 300 bar
R = R_def*ones(1,n);
T = T_def*ones(1,n);
%pcc = Pcc_def*ones(1,100);
pa = Pa_def*ones(1,n);
exp = exp_def*ones(1,n);
gamma = gamma_def*ones(1,n);

pcc = linspace(Pcc_def*0.5,Pcc_def*1.5,n);

for i = 1:n  
    
    Pcc_struct.c_star(i) = sqrt(g*gamma(i)*R(i)*T(i))/(gamma(i)*sqrt((2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1))));

    Pcc_struct.Cf(i) = sqrt(((2*gamma(i)^2)/(gamma(i)-1))*((2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1)))*(1-(Pe(i)/pcc(i))^((gamma(i)-1)/gamma(i))))+exp(i)*((Pe(i)-pa(i))/pcc(i));

    Pcc_struct.F(i) = pcc(i)*At*Pcc_struct.Cf(i);

    Pcc_struct.Isp(i) = (Pcc_struct.c_star(i)*Pcc_struct.Cf(i))/g;

end

% Pa loop 0 - 1 bar
R = R_def*ones(1,n);
T = T_def*ones(1,n);
pcc = Pcc_def*ones(1,n);
%Pa = Pa_def*ones(1,100);
exp = exp_def*ones(1,n);
gamma = gamma_def*ones(1,n);

pa = linspace(0,1e5,n);

for i = 1:n
    
    Pa_struct.c_star(i) = sqrt(g*gamma(i)*R(i)*T(i))/(gamma(i)*sqrt((2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1))));

    Pa_struct.Cf(i) = sqrt(((2*gamma(i)^2)/(gamma(i)-1))*((2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1)))*(1-(Pe(i)/pcc(i))^((gamma(i)-1)/gamma(i))))+exp(i)*((Pe(i)-pa(i))/pcc(i));

    Pa_struct.F(i) = pcc(i)*At*Pa_struct.Cf(i);

    Pa_struct.Isp(i) = (Pa_struct.c_star(i)*Pa_struct.Cf(i))/g;

end

%exp loop 3 - 100;
R = R_def*ones(1,n);
T = T_def*ones(1,n);
pcc = Pcc_def*ones(1,n);
pa = Pa_def*ones(1,n);
%exp = exp_def*ones(1,100);
gamma = gamma_def*ones(1,n);

exp = linspace(exp_def*0.5,exp_def*1.5,n);

for i = 1:n 
    
    exp_struct.c_star(i) = sqrt(g*gamma(i)*R(i)*T(i))/(gamma(i)*sqrt((2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1))));

    exp_struct.Cf(i) = sqrt(((2*gamma(i)^2)/(gamma(i)-1))*((2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1)))*(1-(Pe(i)/pcc(i))^((gamma(i)-1)/gamma(i))))+exp(i)*((Pe(i)-pa(i))/pcc(i));

    exp_struct.F(i) = pcc(i)*At*exp_struct.Cf(i);

    exp_struct.Isp(i) = (exp_struct.c_star(i)*exp_struct.Cf(i))/g;

end

x_arr = linspace(50,150,n);
x_arr2 = linspace(50,100,n);

figure(1)
hold on
plot(x_arr,gamma_struct.F)
plot(x_arr,R_struct.F,'dk')
plot(x_arr,T_struct.F,'sk')
plot(x_arr,Pcc_struct.F)
plot(x_arr2,Pa_struct.F)
plot(x_arr,exp_struct.F)
legend('\gamma','R','T','P_{cc}','P_{a}','r_{exp}')
hold off 
grid on
xlabel('% Change')
ylabel('Thust [N]')


figure(2)
hold on
plot(x_arr,gamma_struct.Isp)
plot(x_arr,R_struct.Isp)
plot(x_arr,T_struct.Isp,'sk')
plot(x_arr,Pcc_struct.Isp)
plot(x_arr2,Pa_struct.Isp)
plot(x_arr,exp_struct.Isp)
legend('\gamma','R','T','P_{cc}','P_{a}','r_{exp}')
hold off 
xlabel('% Change')
ylabel('I_{sp} [m/s]')
grid on

figure(3)
hold on
plot(x_arr,gamma_struct.Cf)
plot(x_arr,R_struct.Cf,'dk')
plot(x_arr,T_struct.Cf,'sk')
plot(x_arr,Pcc_struct.Cf)
plot(x_arr2,Pa_struct.Cf)
plot(x_arr,exp_struct.Cf)
legend('\gamma','R','T','P_{cc}','P_{a}','r_{exp}')
hold off 
xlabel('% Change')
ylabel('Coefficient of Thrust')
grid on

figure(4)
hold on
plot(x_arr,gamma_struct.c_star)
plot(x_arr,R_struct.c_star,'dk')
plot(x_arr,T_struct.c_star,'sk')
plot(x_arr,Pcc_struct.c_star,'o')
plot(x_arr2,Pa_struct.c_star,'x')
plot(x_arr,exp_struct.c_star)
legend('\gamma','R','T','P_{cc}','P_{a}','r_{exp}')
hold off 
xlabel('% Change')
ylabel('C star')
grid on



