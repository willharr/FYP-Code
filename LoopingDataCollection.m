% Looping Data Collection Script
% Will Harradence
% Imperial Aeronautics 2019/20
% FYP 

close all
clear
clc

items = 3;
item_list = ["max_temp_","m_dot_","q_dot_"];

delete max_temp_results.mat
delete q_dot_results.mat
delete m_dot_results.mat

%Reset loop counter
load looper
looper = 1;
save('looper.mat','looper')

%Reset second loop counter
load looper2
looper2 = 1;
save('looper2.mat','looper2')


load materials_list %get list of material names
wall_thickness = [
0.05
0.055
0.06
0.065
0.07
0.075
0.08
0.085
0.09
0.095
0.1
0.105
0.11
0.115
0.12
0.125
0.13
];

save('wall_thickness.mat','wall_thickness')

loops = length(materials_list); %set number of material loops

loops2 = length(wall_thickness);

%generate results arrays
m_dot_results = zeros(loops,loops2); 
save('m_dot_results.mat','m_dot_results')

q_dot_results = zeros(loops,loops2);
save('q_dot_results.mat','q_dot_results')

max_temp_results = zeros(loops,loops2);
save('max_temp_results.mat','max_temp_results')

%% Calculation Loops

for counter = 1:loops %materials loop
    
    for counter2 = 1:loops2
    tic 
    save('loopVariables.mat','counter','counter2','loops2','loops','materials_list')
    save('looper.mat','looper')
    save('looper2.mat','looper2') 
    
    
    matlabThermalModelLoop
    clear
    
    load loopVariables
    load m_dot_results
    load q_dot_results
    load max_temp_results
    load wall_thickness
    load looper
    load looper2
    
    looper2 = looper2+1;
    fprintf('%s at t_wall = %.3f : \n Mdot = %.4f \n Qdot = %.4f \n Tmax = %.2f \n',materials_list(counter),wall_thickness(counter2),m_dot_results(counter,counter2),q_dot_results(counter,counter2),max_temp_results(counter,counter2))
    toc
    end
    looper2 = 1;
    looper = looper+1;
end

m_dot_results = m_dot_results';
q_dot_results = q_dot_results';
max_temp_results = max_temp_results';

save('\Results\Run2\q_dot_resultsR2E1.mat','q_dot_results')
save('\Results\Run2\m_dot_resultsR2E1.mat','m_dot_results')
save('\Results\Run2\max_temp_resultsR2E1.mat','max_temp_results')
