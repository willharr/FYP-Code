clear 
clc
close all

load n2o.mat

n = 1000;

T = linspace(270,309.559,n-1);

T = [T 309.56];


for i = 1:n
   
    n2o.T = T(i);
    heatVap(i) = heatVaporisation(n2o);
    Pvap(i) = vapPressure(n2o);
    CLiquid(i) = heatCapLiquid(n2o);
    CVap(i) = heatCapVap(n2o);
    rhoL(i) = densityLiquid(n2o);
    conduct(i) = conductLiquid(n2o);
    
end

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

T_arr = linspace(183,309,m);

for i = 1:m
    
    P_phase(i) = py.CoolProp.CoolProp.PropsSI('P','T',T_arr(i),'Q',0,'N2O'); 
    
end

width = 18;
height = 15;

figure
plot(T_arr,P_phase/1e5,'-b','LineWidth',3)
hold on
plot([182,400],[n2o.P_crit/1e5,n2o.P_crit/1e5],'--k','LineWidth',2)
plot([n2o.T_crit,n2o.T_crit],[0,100],'--k','LineWidth',2)
grid on
box on
axis('tight')
xlabel('Temperature [K]')
ylabel('Presure [bar]')
legend('Saturation Line','Critical Lines')
set(gca,'FontSize',16)
set(gcf,'units','centimeters','position',[5,5,width,height])


figure
plot(H_L/1e3,P_arr/1e5,'-b','LineWidth',3)
hold on
plot(H_G/1e3,P_arr/1e5,'-r','LineWidth',3)
contour(H_arr/1e3,P_arr2/1e5,T_HP,[200 220 240 260 280 300 309.57 320 340],'-k','ShowText','on','LabelSpacing',800)
hold off
grid on
box on
axis([0,500,0,75])
xlabel('Enthalpy [kJ/kgK]')
ylabel('Pressure [bar]')
legend('Saturated Liquid','Saturated Vapour','Isotherms','orientation','horizontal','Location','northoutside')
set(gca,'FontSize',16)
set(gcf,'units','centimeters','position',[5,5,width,height])

figure
plot(T,heatVap/1e3,'LineWidth',3)
grid on
xlabel('Temperature [K]')
ylabel('Enthalpy of Vaporisation [kJ/kg]')
set(gca,'FontSize',16)
set(gcf,'units','centimeters','position',[5,5,width,height])

figure
plot(T,Pvap/1e5,'LineWidth',3)
grid on
xlabel('Temperature [K]')
ylabel('Vapour Pressure [bar]')
set(gca,'FontSize',16)
set(gcf,'units','centimeters','position',[5,5,width,height])

figure
plot(T,CLiquid,'LineWidth',3)
hold on
plot(T,CVap,'LineWidth',3)
hold off
grid on
xlabel('Temperature [K]')
ylabel('Specific Heat Capacity [J/kg*K]')
set(gca,'FontSize',16)
set(gcf,'units','centimeters','position',[5,5,width,height])

figure
plot(T,rhoL,'LineWidth',3)
grid on
xlabel('Temperature [K]')
ylabel('Density [kg/m^3]')
set(gca,'FontSize',16)
set(gcf,'units','centimeters','position',[5,5,width,height])

figure
plot(T,conduct,'LineWidth',3)
grid on
xlabel('Temperature [K]')
ylabel('Conductivity [W/mK]')
set(gca,'FontSize',16)
set(gcf,'units','centimeters','position',[5,5,width,height])