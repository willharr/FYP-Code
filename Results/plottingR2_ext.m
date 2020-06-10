% Plotting Script
% Will Harradence
% Imperial Aeronautics 2019/20
% FYP 

clear
clc
close all

load extrapolatedR2.mat
load m_dot_resultsR2E1.mat
load q_dot_resultsR2E1.mat
load max_temp_resultsR2E1.mat
load materials.mat
load colours.mat

fontsize = 10;
fonttype = 'Times New Roman';
% position = [left bottom width height]
width = 15;
height = 15;

t_ext_1 = [
0.002
0.003
0.004];

%% Material Comparison

%Combined temperature and mass flow rate plot E2
figure
yyaxis left
plot(t_ext_2,max_temp_resultsR2E2ext(:,1),'Color',colours(1,:))
hold on
plot(t_ext_2,max_temp_resultsR2E2ext(:,2),'-','Color',colours(2,:))
plot(t_ext_2,max_temp_resultsR2E2ext(:,3),'-','Color',colours(3,:))
plot(t_ext_2,max_temp_resultsR2E2ext(:,4),'-','Color',colours(4,:))
plot(t_ext_2,max_temp_resultsR2E2ext(:,5),'-','Color',colours(5,:))
plot(t_ext_2,max_temp_resultsR2E2ext(:,6),'-','Color',colours(6,:))
ylabel('Hot Side Wall Temperature [K]')
yyaxis right
plot(t_ext_2,m_dot_resultsR2E2ext(:,1),'--','Color',colours(1,:))
plot(t_ext_2,m_dot_resultsR2E2ext(:,2),'--','Color',colours(2,:))
plot(t_ext_2,m_dot_resultsR2E2ext(:,3),'--','Color',colours(3,:))
plot(t_ext_2,m_dot_resultsR2E2ext(:,4),'--','Color',colours(4,:))
plot(t_ext_2,m_dot_resultsR2E2ext(:,5),'--','Color',colours(5,:))
plot(t_ext_2,m_dot_resultsR2E2ext(:,6),'--','Color',colours(6,:))
plot([0.005,0.035],[8,8],'--k')
hold off
grid on
box on
legend('Inconel 625','\beta Ti','Al 6063','Cu-Zr Alloy','Be I250','Al 7068','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
ylabel('Coolant Mass Flow Rate [kg/s]')
set(gcf,'units','centimeters','position',[5,5,width,height])
saveas(gcf,'R2_twoaxis','png')

%=============================================================================

% Just mass flow rate E2
figure
plot(t_ext_2(6:end),m_dot_resultsR2E2ext(6:end,1),'-','Color',colours(1,:))
hold on
plot(t_ext_2(6:end),m_dot_resultsR2E2ext(6:end,2),'-','Color',colours(2,:))
plot(t_ext_2(6:end),m_dot_resultsR2E2ext(6:end,3),'-','Color',colours(3,:))
plot(t_ext_2(6:end),m_dot_resultsR2E2ext(6:end,4),'-','Color',colours(4,:))
plot(t_ext_2(6:end),m_dot_resultsR2E2ext(6:end,5),'-','Color',colours(5,:))
plot(t_ext_2(6:end),m_dot_resultsR2E2ext(6:end,6),'-','Color',colours(6,:))
plot([0,0.035],[8,8],'--k')
plot([0,0.035],[0.8,0.8],'-.k')
%============================
plot(t_ext_2(1:6),m_dot_resultsR2E2ext(1:6,1),'--','Color',colours(1,:))
plot(t_ext_2(1:6),m_dot_resultsR2E2ext(1:6,2),'--','Color',colours(2,:))
plot(t_ext_2(1:6),m_dot_resultsR2E2ext(1:6,3),'--','Color',colours(3,:))
plot(t_ext_2(1:6),m_dot_resultsR2E2ext(1:6,4),'--','Color',colours(4,:))
plot(t_ext_2(1:6),m_dot_resultsR2E2ext(1:6,5),'--','Color',colours(5,:))
plot(t_ext_2(1:6),m_dot_resultsR2E2ext(1:6,6),'--','Color',colours(6,:))
hold off
grid on
box on
legend('Inconel 625','\beta Ti','Al 6063','Cu-Zr Alloy','Be I250','Al 7068','Propulsion','10% Propulsion','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
ylabel('Coolant Mass Flow Rate [kg/s]')
axis([0 0.035 0 8.5])
set(gcf,'units','centimeters','position',[5,5,width,height])
saveas(gcf,'R2E2mdotvstwall','png')

% Just mass flow rate E3
figure
plot(t_ext_3(6:end),m_dot_resultsR2E3ext(6:end,1),'-','Color',colours(1,:))
hold on
plot(t_ext_3(6:end),m_dot_resultsR2E3ext(6:end,2),'-','Color',colours(2,:))
plot(t_ext_3(6:end),m_dot_resultsR2E3ext(6:end,3),'-','Color',colours(3,:))
plot(t_ext_3(6:end),m_dot_resultsR2E3ext(6:end,4),'-','Color',colours(4,:))
plot(t_ext_3(6:end),m_dot_resultsR2E3ext(6:end,5),'-','Color',colours(5,:))
plot(t_ext_3(6:end),m_dot_resultsR2E3ext(6:end,6),'-','Color',colours(6,:))
plot([0,0.06],[21.57,21.57],'--k')
plot([0,0.06],[2.157,2.157],'-.k')
%============================
plot(t_ext_3(1:6),m_dot_resultsR2E3ext(1:6,1),'--','Color',colours(1,:))
plot(t_ext_3(1:6),m_dot_resultsR2E3ext(1:6,2),'--','Color',colours(2,:))
plot(t_ext_3(1:6),m_dot_resultsR2E3ext(1:6,3),'--','Color',colours(3,:))
plot(t_ext_3(1:6),m_dot_resultsR2E3ext(1:6,4),'--','Color',colours(4,:))
plot(t_ext_3(1:6),m_dot_resultsR2E3ext(1:6,5),'--','Color',colours(5,:))
plot(t_ext_3(1:6),m_dot_resultsR2E3ext(1:6,6),'--','Color',colours(6,:))
hold off
grid on
box on
legend('Inconel 625','\beta Ti','Al 6063','Cu-Zr Alloy','Be I250','Al 7068','Propulsion','10% Propulsion','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
ylabel('Coolant Mass Flow Rate [kg/s]')
axis([0 0.052 0 21.9])
set(gcf,'units','centimeters','position',[5,5,width,height])
saveas(gcf,'R2E3mdotvstwall','png')

%===========================================================================

%just max temp E2
figure
plot(t_ext_2(6:end),max_temp_resultsR2E2ext(6:end,1),'-','Color',colours(1,:))
hold on
plot(t_ext_2(6:end),max_temp_resultsR2E2ext(6:end,2),'-','Color',colours(2,:))
plot(t_ext_2(6:end),max_temp_resultsR2E2ext(6:end,3),'-','Color',colours(3,:))
plot(t_ext_2(6:end),max_temp_resultsR2E2ext(6:end,4),'-','Color',colours(4,:))
plot(t_ext_2(6:end),max_temp_resultsR2E2ext(6:end,5),'-','Color',colours(5,:))
plot(t_ext_2(6:end),max_temp_resultsR2E2ext(6:end,6),'-','Color',colours(6,:))
%=================================
plot(0.05,0,'xk')
plot(0.05,0,'dk')
%=======
plot(t_ext_2(1:6),max_temp_resultsR2E2ext(1:6,1),'--','Color',colours(1,:))
plot(t_ext_2(1:6),max_temp_resultsR2E2ext(1:6,2),'--','Color',colours(2,:))
plot(t_ext_2(1:6),max_temp_resultsR2E2ext(1:6,3),'--','Color',colours(3,:))
plot(t_ext_2(1:6),max_temp_resultsR2E2ext(1:6,4),'--','Color',colours(4,:))
plot(t_ext_2(1:6),max_temp_resultsR2E2ext(1:6,5),'--','Color',colours(5,:))
plot(t_ext_2(1:6),max_temp_resultsR2E2ext(1:6,6),'--','Color',colours(6,:))
%=================================
plot(0.032,materials(2,1),'x','Color',colours(1,:));
plot(0.032,materials(1,1),'d','Color',colours(1,:));
plot(0.034,materials(2,2),'x','Color',colours(2,:));
plot(0.034,materials(1,2),'d','Color',colours(2,:));
plot(0.003,materials(2,3),'x','Color',colours(3,:));
plot(0.003,materials(1,3),'d','Color',colours(3,:));
plot(0.008,materials(2,4),'x','Color',colours(4,:));
plot(0.008,materials(1,4),'d','Color',colours(4,:));
plot(0.006,materials(2,5),'x','Color',colours(5,:));
plot(0.006,materials(1,5),'d','Color',colours(5,:));
plot(0.002,materials(2,6),'x','Color',colours(6,:));
plot(0.002,materials(1,6),'d','Color',colours(6,:));
ylabel('Hot Side Wall Temperature [K]')
hold off
grid on
box on
legend('Inconel 625','\beta Ti','Al 6063','Cu-Zr Alloy','Be I250','Al 7068','Max Service Temp','Melting Temp','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
axis([0 0.035 0 3000])
set(gcf,'units','centimeters','position',[5,5,width,height])
saveas(gcf,'R2E2tempvstwallmaterials','png')


%% Engine Size

% Beryllium data
figure
plot(t_ext_1,max_temp_resultsR2E1(:,5),'Color',colours(1,:))
hold on
plot(t_ext_2(6:end),max_temp_resultsR2E2ext(6:end,5),'-','Color',colours(2,:))
plot(t_ext_3(6:end),max_temp_resultsR2E3ext(6:end,5),'-','Color',colours(3,:))
plot(t_ext_4(12:end),max_temp_resultsR2E4ext(12:end,5),'-','Color',colours(4,:))
plot(t_ext_5(24:end),max_temp_resultsR2E5ext(24:end,5),'-','Color',colours(5,:))
%=========================
plot([0 0.07],[materials(2,5) materials(2,5)],'-.k');
plot([0 0.07],[materials(1,5) materials(1,5)],':k');
%==========================
plot(t_ext_2(1:6),max_temp_resultsR2E2ext(1:6,5),'--','Color',colours(2,:))
plot(t_ext_3(1:6),max_temp_resultsR2E3ext(1:6,5),'--','Color',colours(3,:))
plot(t_ext_4(1:12),max_temp_resultsR2E4ext(1:12,5),'--','Color',colours(4,:))
plot(t_ext_5(1:24),max_temp_resultsR2E5ext(1:24,5),'--','Color',colours(5,:))
ylabel('Hot Side Wall Temperature [K]')
hold off
grid on
box on
legend('Engine 1','Engine 2','Engine 3','Engine 4','Engine 5','Max Service Temp','Melting Temp','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
axis([0 0.07 0 3000])
set(gcf,'units','centimeters','position',[5,5,width,height])
saveas(gcf,'R2Betempvstwallengines','png')

figure
semilogy(t_ext_1,m_dot_resultsR2E1(:,5),'Color',colours(1,:))
hold on
semilogy(t_ext_2(6:end),m_dot_resultsR2E2ext(6:end,5),'-','Color',colours(2,:))
semilogy(t_ext_3(6:end),m_dot_resultsR2E3ext(6:end,5),'-','Color',colours(3,:))
semilogy(t_ext_4(12:end),m_dot_resultsR2E4ext(12:end,5),'-','Color',colours(4,:))
semilogy(t_ext_5(24:end),m_dot_resultsR2E5ext(24:end,5),'-','Color',colours(5,:))
%===============================
semilogy(t_ext_2(1:6),m_dot_resultsR2E2ext(1:6,5),'--','Color',colours(2,:))
semilogy(t_ext_3(1:6),m_dot_resultsR2E3ext(1:6,5),'--','Color',colours(3,:))
semilogy(t_ext_4(1:12),m_dot_resultsR2E4ext(1:12,5),'--','Color',colours(4,:))
semilogy(t_ext_5(1:24),m_dot_resultsR2E5ext(1:24,5),'--','Color',colours(5,:))
%==============================
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
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','off','FontName', fonttype)
xlabel('Wall Thickness [m]')
axis([0 0.07 0.01 1000])
set(gcf,'units','centimeters','position',[5,5,width,height])
saveas(gcf,'R2Bemdotvstwallengines','png')

figure
plot(t_ext_1,q_dot_resultsR2E1(:,5),'Color',colours(1,:))
hold on
plot(t_ext_2(6:end),q_dot_resultsR2E2ext(6:end,5),'-','Color',colours(2,:))
plot(t_ext_3(6:end),q_dot_resultsR2E3ext(6:end,5),'-','Color',colours(3,:))
plot(t_ext_4(12:end),q_dot_resultsR2E4ext(12:end,5),'-','Color',colours(4,:))
plot(t_ext_5(24:end),q_dot_resultsR2E5ext(24:end,5),'-','Color',colours(5,:))
%=================================
plot(t_ext_2(1:6),q_dot_resultsR2E2ext(1:6,5),'--','Color',colours(2,:))
plot(t_ext_3(1:6),q_dot_resultsR2E3ext(1:6,5),'--','Color',colours(3,:))
plot(t_ext_4(1:12),q_dot_resultsR2E4ext(1:12,5),'--','Color',colours(4,:))
plot(t_ext_5(1:24),q_dot_resultsR2E5ext(1:24,5),'--','Color',colours(5,:))
%=================================
ylabel('Heat Flux Through Wall [J/s]')
hold off
grid on
box on
legend('Engine 1','Engine 2','Engine 3','Engine 4','Engine 5','Location','southoutside','Numcolumns',2)
set(gca,'FontSize',fontsize,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on','FontName', fonttype)
xlabel('Wall Thickness [m]')
set(gcf,'units','centimeters','position',[5,5,width,height])
%saveas(gcf,'R2 vs','png')

