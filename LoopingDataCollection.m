% Looping Data Collection Script
% Will Harradence
% Imperial Aeronautics 2019/20
% FYP 

close all
clear
clc

%Reset loop counter
load looper
looper = 1;
save('looper.mat','looper')
load materials_list

loops = length(materials_list);
m_dot_results = zeros(1,loops);
save('m_dot_results.mat','m_dot_results')

q_dot_results = zeros(1,loops);
save('q_dot_results.mat','q_dot_results')

max_temp_results = zeros(1,loops);
save('max_temp_results.mat','max_temp_results')

for counter = 1:loops
    
    save('loopVariables.mat','counter','loops','materials_list')
    save('looper.mat','looper')
    
    matlabThermalModelLoop
    clear
    
    load loopVariables
    load m_dot_results
    load q_dot_results
    load max_temp_results
    load looper
    
    looper = looper+1;
    fprintf('%s : \n Mdot = %.4f \n Qdot = %.4f \n Tmax = %.2f \n',materials_list(counter),m_dot_results(counter),q_dot_results(counter),max_temp_results(counter))
    
end
