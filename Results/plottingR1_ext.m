% Plotting Script
% Will Harradence
% Imperial Aeronautics 2019/20
% FYP 

clear
clc
close all

load extrapolatedR1.mat
load m_dot_resultsR1E1.mat
load q_dot_resultsR1E1.mat
load max_temp_resultsR1E1.mat
load materials.mat
load colours.mat

fontsize = 10;
fonttype = 'Times New Roman';
width = 9;
height = 9;

t_wall_1 = [
0.002
0.003
0.004];

%% Material Comparison

%Combined temperature and mass flow rate plot E2
figure
yyaxis left
plot(t_ext_2,max_temp_resultsR1E2ext(2,:),'Color',colours(1,:))
hold on
plot(t_ext_2,max_temp_resultsR1E2ext(3,:),'-','Color',colours(2,:))
plot(t_ext_2,max_temp_resultsR1E2ext(4,:),'-','Color',colours(3,:))
plot(t_ext_2,max_temp_resultsR1E2ext(5,:),'-','Color',colours(4,:))
plot(t_ext_2,max_temp_resultsR1E2ext(6,:),'-','Color',colours(5,:))
plot(t_ext_2,max_temp_resultsR1E2ext(7,:),'-','Color',colours(6,:))
ylabel('Hot Side Wall Temperature [K]')
yyaxis right
plot(t_ext_2,m_dot_resultsR1E2ext(2,:),'--','Color',colours(1,:))
plot(t_ext_2,m_dot_resultsR1E2ext(3,:),'--','Color',colours(2,:))
plot(t_ext_2,m_dot_resultsR1E2ext(4,:),'--','Color',colours(3,:))
plot(t_ext_2,m_dot_resultsR1E2ext(5,:),'--','Color',colours(4,:))
plot(t_ext_2,m_dot_resultsR1E2ext(6,:),'--','Color',colours(5,:))
plot(t_ext_2,m_dot_resultsR1E2ext(7,:),'--','Color',colours(6,:))
plot([0.005,0.035],[8,8],'--k')
hold off
grid on
box on
legend('Inconel 625','\beta Ti','Al 6063','Cu-Zr Alloy','Be I250','Al 7068','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
ylabel('Coolant Mass Flow Rate [kg/s]')
set(gcf,'units','centimeters','position',[5,5,width,height])
saveas(gcf,'R1_twoaxis','epsc')


%==========================================================================

% Just mass flow rate E2
figure
plot(t_ext_2(6:end),m_dot_resultsR1E2ext(2,6:end),'-','Color',colours(1,:))
hold on
plot(t_ext_2(6:end),m_dot_resultsR1E2ext(3,6:end),'-','Color',colours(2,:))
plot(t_ext_2(6:end),m_dot_resultsR1E2ext(4,6:end),'-','Color',colours(3,:))
plot(t_ext_2(6:end),m_dot_resultsR1E2ext(5,6:end),'-','Color',colours(4,:))
plot(t_ext_2(6:end),m_dot_resultsR1E2ext(6,6:end),'-','Color',colours(5,:))
plot(t_ext_2(6:end),m_dot_resultsR1E2ext(7,6:end),'-','Color',colours(6,:))
plot([0,0.035],[8,8],'--k')
plot([0,0.035],[0.8,0.8],'-.k')
%===============
plot(t_ext_2(1:6),m_dot_resultsR1E2ext(2,1:6),'--','Color',colours(1,:))
plot(t_ext_2(1:6),m_dot_resultsR1E2ext(3,1:6),'--','Color',colours(2,:))
plot(t_ext_2(1:6),m_dot_resultsR1E2ext(4,1:6),'--','Color',colours(3,:))
plot(t_ext_2(1:6),m_dot_resultsR1E2ext(5,1:6),'--','Color',colours(4,:))
plot(t_ext_2(1:6),m_dot_resultsR1E2ext(6,1:6),'--','Color',colours(5,:))
plot(t_ext_2(1:6),m_dot_resultsR1E2ext(7,1:6),'--','Color',colours(6,:))
hold off
grid on
box on
legend('Inconel 625','\beta Ti','Al 6063','Cu-Zr Alloy','Be I250','Al 7068','Propulsion','10% Propulsion','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
ylabel('Coolant Mass Flow Rate [kg/s]')
axis([0 0.035 0 100])
set(gcf,'units','centimeters','position',[5,5,width,height])
saveas(gcf,'R1E2mdotvstwall','epsc')

% Just mass flow rate E3
figure
plot(t_ext_3(6:end),m_dot_resultsR1E3ext(2,6:end),'-','Color',colours(1,:))
hold on
plot(t_ext_3(6:end),m_dot_resultsR1E3ext(3,6:end),'-','Color',colours(2,:))
plot(t_ext_3(6:end),m_dot_resultsR1E3ext(4,6:end),'-','Color',colours(3,:))
plot(t_ext_3(6:end),m_dot_resultsR1E3ext(5,6:end),'-','Color',colours(4,:))
plot(t_ext_3(6:end),m_dot_resultsR1E3ext(6,6:end),'-','Color',colours(5,:))
plot(t_ext_3(6:end),m_dot_resultsR1E3ext(7,6:end),'-','Color',colours(6,:))
plot([0,0.06],[21.57,21.57],'--k')
plot([0,0.06],[2.157,2.157],'-.k')
%==============
plot(t_ext_3(1:6),m_dot_resultsR1E3ext(2,1:6),'--','Color',colours(1,:))
plot(t_ext_3(1:6),m_dot_resultsR1E3ext(3,1:6),'--','Color',colours(2,:))
plot(t_ext_3(1:6),m_dot_resultsR1E3ext(4,1:6),'--','Color',colours(3,:))
plot(t_ext_3(1:6),m_dot_resultsR1E3ext(5,1:6),'--','Color',colours(4,:))
plot(t_ext_3(1:6),m_dot_resultsR1E3ext(6,1:6),'--','Color',colours(5,:))
plot(t_ext_3(1:6),m_dot_resultsR1E3ext(7,1:6),'--','Color',colours(6,:))
hold off
grid on
box on
legend('Inconel 625','\beta Ti','Al 6063','Cu-Zr Alloy','Be I250','Al 7068','Propulsion','10% Propulsion','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
ylabel('Coolant Mass Flow Rate [kg/s]')
axis([0 0.052 0 200])
set(gcf,'units','centimeters','position',[5,5,width,height])
saveas(gcf,'R1E3mdotvstwall','epsc')

%===========================================================================

%just max temp E2
figure
plot(t_ext_2(6:end),max_temp_resultsR1E2ext(2,6:end),'-','Color',colours(1,:))
hold on
plot(t_ext_2(6:end),max_temp_resultsR1E2ext(3,6:end),'-','Color',colours(2,:))
plot(t_ext_2(6:end),max_temp_resultsR1E2ext(4,6:end),'-','Color',colours(3,:))
plot(t_ext_2(6:end),max_temp_resultsR1E2ext(5,6:end),'-','Color',colours(4,:))
plot(t_ext_2(6:end),max_temp_resultsR1E2ext(6,6:end),'-','Color',colours(5,:))
plot(t_ext_2(6:end),max_temp_resultsR1E2ext(7,6:end),'-','Color',colours(6,:))
%===========================
plot(t_ext_2(1:6),max_temp_resultsR1E2ext(2,1:6),'--','Color',colours(1,:))
plot(t_ext_2(1:6),max_temp_resultsR1E2ext(3,1:6),'--','Color',colours(2,:))
plot(t_ext_2(1:6),max_temp_resultsR1E2ext(4,1:6),'--','Color',colours(3,:))
plot(t_ext_2(1:6),max_temp_resultsR1E2ext(5,1:6),'--','Color',colours(4,:))
plot(t_ext_2(1:6),max_temp_resultsR1E2ext(6,1:6),'--','Color',colours(5,:))
plot(t_ext_2(1:6),max_temp_resultsR1E2ext(7,1:6),'--','Color',colours(6,:))
%===========================
plot(0.01,materials(2,2),'x','Color',colours(2,:));
plot(0.01,materials(2,3),'x','Color',colours(3,:));
plot(0.01,materials(2,4),'x','Color',colours(4,:));
plot(0.01,materials(2,5),'x','Color',colours(5,:));
plot(0.01,materials(2,6),'x','Color',colours(6,:));
plot(0.01,materials(1,2),'d','Color',colours(2,:));
plot(0.01,materials(1,3),'d','Color',colours(3,:));
plot(0.01,materials(1,4),'d','Color',colours(4,:));
plot(0.01,materials(1,5),'d','Color',colours(5,:));
plot(0.01,materials(1,6),'d','Color',colours(6,:));
ylabel('Hot Side Wall Temperature [K]')
hold off
grid on
box on
legend('Inconel 625','\beta Ti','Al 6063','Cu-Zr Alloy','Be I250','Al 7068','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
set(gcf,'units','centimeters','position',[5,5,width,height])
saveas(gcf,'R1E2tempvstwallmaterials','epsc')


%% Engine Size

% Beryllium data
figure
plot(t_wall_1,max_temp_resultsR1E1(6,:),'Color',colours(1,:))
hold on
plot(t_ext_2(6:end),max_temp_resultsR1E2ext(6,6:end),'-','Color',colours(2,:))
plot(t_ext_3(6:end),max_temp_resultsR1E3ext(6,6:end),'-','Color',colours(3,:))
plot(t_ext_4(12:end),max_temp_resultsR1E4ext(6,12:end),'-','Color',colours(4,:))
plot(t_ext_5(24:end),max_temp_resultsR1E5ext(6,24:end),'-','Color',colours(5,:))
%===================
plot(t_ext_2(1:6),max_temp_resultsR1E2ext(6,1:6),'--','Color',colours(2,:))
plot(t_ext_3(1:6),max_temp_resultsR1E3ext(6,1:6),'--','Color',colours(3,:))
plot(t_ext_4(1:12),max_temp_resultsR1E4ext(6,1:12),'--','Color',colours(4,:))
plot(t_ext_5(1:24),max_temp_resultsR1E5ext(6,1:24),'--','Color',colours(5,:))
ylabel('Hot Side Wall Temperature [K]')
hold off
grid on
box on
legend('Engine 1','Engine 2','Engine 3','Engine 4','Engine 5','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
axis([0 0.08 0 3000])
set(gcf,'units','centimeters','position',[5,5,width,height])
saveas(gcf,'R1Betempvstwallengines','epsc')

figure
semilogy(t_wall_1,m_dot_resultsR1E1(6,:),'Color',colours(1,:))
hold on
semilogy(t_ext_2(6:end),m_dot_resultsR1E2ext(6,6:end),'-','Color',colours(2,:))
semilogy(t_ext_3(6:end),m_dot_resultsR1E3ext(6,6:end),'-','Color',colours(3,:))
semilogy(t_ext_4(12:end),m_dot_resultsR1E4ext(6,12:end),'-','Color',colours(4,:))
semilogy(t_ext_5(24:end),m_dot_resultsR1E5ext(6,24:end),'-','Color',colours(5,:))
%===========================
semilogy(t_ext_2(1:6),m_dot_resultsR1E2ext(6,1:6),'--','Color',colours(2,:))
semilogy(t_ext_3(1:6),m_dot_resultsR1E3ext(6,1:6),'--','Color',colours(3,:))
semilogy(t_ext_4(1:12),m_dot_resultsR1E4ext(6,1:12),'--','Color',colours(4,:))
semilogy(t_ext_5(1:24),m_dot_resultsR1E5ext(6,1:24),'--','Color',colours(5,:))
%=============================
semilogy([0 0.15],[0.15 0.15],'-.','Color',colours(1,:))
semilogy([0 0.15],[8 8],'-.','Color',colours(2,:))
semilogy([0 0.15],[21.57 21.57],'-.','Color',colours(3,:))
semilogy([0 0.15],[84 84],'-.','Color',colours(4,:))
semilogy([0 0.15],[760 760],'-.','Color',colours(5,:))
ylabel('Coolant Mass Flow Rate [kg/s]')
hold off
grid on
box on
legend('Engine 1','Engine 2','Engine 3','Engine 4','Engine 5','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
axis([0 0.07 0.1 10000])
set(gcf,'units','centimeters','position',[5,5,width,height])
saveas(gcf,'R1Bemdotvstwallengines','epsc')

figure
plot(t_wall_1,q_dot_resultsR1E1(6,:),'Color',colours(1,:))
hold on
plot(t_ext_2(6:end),q_dot_resultsR1E2ext(6,6:end),'-','Color',colours(2,:))
plot(t_ext_3(6:end),q_dot_resultsR1E3ext(6,6:end),'-','Color',colours(3,:))
plot(t_ext_4(12:end),q_dot_resultsR1E4ext(6,12:end),'-','Color',colours(4,:))
plot(t_ext_5(24:end),q_dot_resultsR1E5ext(6,24:end),'-','Color',colours(5,:))
%======================
plot(t_ext_2(1:6),q_dot_resultsR1E2ext(6,1:6),'--','Color',colours(2,:))
plot(t_ext_3(1:6),q_dot_resultsR1E3ext(6,1:6),'--','Color',colours(3,:))
plot(t_ext_4(1:12),q_dot_resultsR1E4ext(6,1:12),'--','Color',colours(4,:))
plot(t_ext_5(1:24),q_dot_resultsR1E5ext(6,1:24),'--','Color',colours(5,:))
%======================
ylabel('Heat Flux Through Wall [J/s]')
hold off
grid on
box on
legend('Engine 1','Engine 2','Engine 3','Engine 4','Engine 5','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
set(gcf,'units','centimeters','position',[5,5,width,height])
%saveas(gcf,'.eps','epsc');

