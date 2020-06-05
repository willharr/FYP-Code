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

t_wall_1 = [
0.002
0.003
0.004];

%% Material Comparison


%Combined temperature and mass flow rate plot E2
figure
yyaxis left
plot(t_wall_2,max_temp_resultsR1E2(2,:),'Color',colours(1,:))
hold on
plot(t_wall_2,max_temp_resultsR1E2(3,:),'-','Color',colours(2,:))
plot(t_wall_2,max_temp_resultsR1E2(4,:),'-','Color',colours(3,:))
plot(t_wall_2,max_temp_resultsR1E2(5,:),'-','Color',colours(4,:))
plot(t_wall_2,max_temp_resultsR1E2(6,:),'-','Color',colours(5,:))
plot(t_wall_2,max_temp_resultsR1E2(7,:),'-','Color',colours(6,:))
ylabel('Hot Side Wall Temperature [K]')
yyaxis right
plot(t_wall_2,m_dot_resultsR1E2(2,:),'--','Color',colours(1,:))
plot(t_wall_2,m_dot_resultsR1E2(3,:),'--','Color',colours(2,:))
plot(t_wall_2,m_dot_resultsR1E2(4,:),'--','Color',colours(3,:))
plot(t_wall_2,m_dot_resultsR1E2(5,:),'--','Color',colours(4,:))
plot(t_wall_2,m_dot_resultsR1E2(6,:),'--','Color',colours(5,:))
plot(t_wall_2,m_dot_resultsR1E2(7,:),'--','Color',colours(6,:))
plot([0.005,0.035],[8,8],'--k')
hold off
grid on
box on
legend('Inconel 625','\beta Ti','Al 6063','Cu-Zr Alloy','Be I250','Al 7068','Location','best')
set(gca,'FontSize',14,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on')
xlabel('Wall Thickness [m]')
ylabel('Coolant Mass Flow Rate [kg/s]')
%saveas(gcf,'.eps','epsc');

% Just mass flow rate E2
figure
plot(t_wall_2,m_dot_resultsR1E2(2,:),'--','Color',colours(1,:))
hold on
plot(t_wall_2,m_dot_resultsR1E2(3,:),'--','Color',colours(2,:))
plot(t_wall_2,m_dot_resultsR1E2(4,:),'--','Color',colours(3,:))
plot(t_wall_2,m_dot_resultsR1E2(5,:),'--','Color',colours(4,:))
plot(t_wall_2,m_dot_resultsR1E2(6,:),'--','Color',colours(5,:))
plot(t_wall_2,m_dot_resultsR1E2(7,:),'--','Color',colours(6,:))
plot([0.005,0.035],[8,8],'--k')
hold off
grid on
box on
legend('Inconel 625','\beta Ti','Al 6063','Cu-Zr Alloy','Be I250','Al 7068','Location','best')
set(gca,'FontSize',14,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on')
xlabel('Wall Thickness [m]')
ylabel('Coolant Mass Flow Rate [kg/s]')
%saveas(gcf,'.eps','epsc');

%just max temp E2
figure
plot(t_wall_2,max_temp_resultsR1E2(2,:),'Color',colours(1,:))
hold on
plot(t_wall_2,max_temp_resultsR1E2(3,:),'-','Color',colours(2,:))
plot(t_wall_2,max_temp_resultsR1E2(4,:),'-','Color',colours(3,:))
plot(t_wall_2,max_temp_resultsR1E2(5,:),'-','Color',colours(4,:))
plot(t_wall_2,max_temp_resultsR1E2(6,:),'-','Color',colours(5,:))
plot(t_wall_2,max_temp_resultsR1E2(7,:),'-','Color',colours(6,:))
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
legend('Inconel 625','\beta Ti','Al 6063','Cu-Zr Alloy','Be I250','Al 7068','Location','best')
set(gca,'FontSize',14,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on')
xlabel('Wall Thickness [m]')
%saveas(gcf,'.eps','epsc');


%% Engine Size

% Beryllium data
figure
plot(t_wall_1,max_temp_resultsR1E1(6,:),'Color',colours(1,:))
hold on
plot(t_wall_2,max_temp_resultsR1E2(6,:),'-','Color',colours(2,:))
plot(t_wall_3,max_temp_resultsR1E3(6,:),'-','Color',colours(3,:))
plot(t_wall_4,max_temp_resultsR1E4(6,:),'-','Color',colours(4,:))
plot(t_wall_5,max_temp_resultsR1E5(6,:),'-','Color',colours(5,:))
ylabel('Hot Side Wall Temperature [K]')
hold off
grid on
box on
legend('Engine 1','Engine 2','Engine 3','Engine 4','Engine 5','Location','best')
set(gca,'FontSize',14,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on')
xlabel('Wall Thickness [m]')
%saveas(gcf,'.eps','epsc');

figure
plot(t_wall_1,m_dot_resultsR1E1(6,:),'Color',colours(1,:))
hold on
plot(t_wall_2,m_dot_resultsR1E2(6,:),'-','Color',colours(2,:))
plot(t_wall_3,m_dot_resultsR1E3(6,:),'-','Color',colours(3,:))
plot(t_wall_4,m_dot_resultsR1E4(6,:),'-','Color',colours(4,:))
plot(t_wall_5,m_dot_resultsR1E5(6,:),'-','Color',colours(5,:))
plot([0 0.15],[0.15 0.15],'--','Color',colours(1,:))
plot([0 0.15],[8 8],'--','Color',colours(2,:))
plot([0 0.15],[21.57 21.57],'--','Color',colours(3,:))
plot([0 0.15],[84 84],'--','Color',colours(4,:))
plot([0 0.15],[760 760],'--','Color',colours(5,:))
ylabel('Coolant Mass Flow Rate [kg/s]')
hold off
grid on
box on
legend('Engine 1','Engine 2','Engine 3','Engine 4','Engine 5','Location','best')
set(gca,'FontSize',14,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on')
xlabel('Wall Thickness [m]')
%saveas(gcf,'.eps','epsc');

figure
plot(t_wall_1,q_dot_resultsR1E1(6,:),'Color',colours(1,:))
hold on
plot(t_wall_2,q_dot_resultsR1E2(6,:),'-','Color',colours(2,:))
plot(t_wall_3,q_dot_resultsR1E3(6,:),'-','Color',colours(3,:))
plot(t_wall_4,q_dot_resultsR1E4(6,:),'-','Color',colours(4,:))
plot(t_wall_5,q_dot_resultsR1E5(6,:),'-','Color',colours(5,:))
ylabel('Heat Flux Through Wall [J/s]')
hold off
grid on
box on
legend('Engine 1','Engine 2','Engine 3','Engine 4','Engine 5','Location','best')
set(gca,'FontSize',14,'YMinorTick','on', 'YMinorGrid','on','XMinorTick','on', 'XMinorGrid','on')
xlabel('Wall Thickness [m]')
%saveas(gcf,'.eps','epsc');

