% Data Extrapolation for Thin Walls Script
% Will Harradence
% Imperial Aeronautics 2019/20
% FYP 

clear
clc

load m_dot_resultsR1E2.mat
load m_dot_resultsR1E3.mat
load m_dot_resultsR1E4.mat
load m_dot_resultsR1E5.mat

load max_temp_resultsR1E2.mat
load max_temp_resultsR1E3.mat
load max_temp_resultsR1E4.mat
load max_temp_resultsR1E5.mat

load q_dot_resultsR1E2.mat
load q_dot_resultsR1E3.mat
load q_dot_resultsR1E4.mat
load q_dot_resultsR1E5.mat


t_wall_2 = 0.006:0.002:0.032;
t_wall_3 = 0.006:0.002:0.052;
t_wall_4 = 0.012:0.002:0.044;
t_wall_5 = 0.05:0.005:0.13;

for i = 1:7
    
    T_ext2 = interp1(t_wall_2,max_temp_resultsR1E2(i,:),(0.001:0.001:0.005),'linear','extrap');
    T_ext3 = interp1(t_wall_3,max_temp_resultsR1E3(i,:),(0.001:0.001:0.005),'linear','extrap');
    T_ext4 = interp1(t_wall_4,max_temp_resultsR1E4(i,:),(0.001:0.001:0.011),'linear','extrap');
    T_ext5 = interp1(t_wall_5,max_temp_resultsR1E5(i,:),(0.001:0.002:0.045),'linear','extrap');
    
    max_temp_resultsR1E2ext(i,:) = [T_ext2 max_temp_resultsR1E2(i,:)]; 
    max_temp_resultsR1E3ext(i,:) = [T_ext3 max_temp_resultsR1E3(i,:)]; 
    max_temp_resultsR1E4ext(i,:) = [T_ext4 max_temp_resultsR1E4(i,:)]; 
    max_temp_resultsR1E5ext(i,:) = [T_ext5 max_temp_resultsR1E5(i,:)]; 
    
    M_ext2 = interp1(t_wall_2,m_dot_resultsR1E2(i,:),(0.001:0.001:0.005),'linear','extrap');
    M_ext3 = interp1(t_wall_3,m_dot_resultsR1E3(i,:),(0.001:0.001:0.005),'linear','extrap');
    M_ext4 = interp1(t_wall_4,m_dot_resultsR1E4(i,:),(0.001:0.001:0.011),'linear','extrap');
    M_ext5 = interp1(t_wall_5,m_dot_resultsR1E5(i,:),(0.001:0.002:0.045),'linear','extrap');
    
    m_dot_resultsR1E2ext(i,:) = [M_ext2 m_dot_resultsR1E2(i,:)]; 
    m_dot_resultsR1E3ext(i,:) = [M_ext3 m_dot_resultsR1E3(i,:)]; 
    m_dot_resultsR1E4ext(i,:) = [M_ext4 m_dot_resultsR1E4(i,:)]; 
    m_dot_resultsR1E5ext(i,:) = [M_ext5 m_dot_resultsR1E5(i,:)]; 
    
    Q_ext2 = interp1(t_wall_2,q_dot_resultsR1E2(i,:),(0.001:0.001:0.005),'linear','extrap');
    Q_ext3 = interp1(t_wall_3,q_dot_resultsR1E3(i,:),(0.001:0.001:0.005),'linear','extrap');
    Q_ext4 = interp1(t_wall_4,q_dot_resultsR1E4(i,:),(0.001:0.001:0.011),'linear','extrap');
    Q_ext5 = interp1(t_wall_5,q_dot_resultsR1E5(i,:),(0.001:0.002:0.045),'linear','extrap');
    
    q_dot_resultsR1E2ext(i,:) = [Q_ext2 q_dot_resultsR1E2(i,:)]; 
    q_dot_resultsR1E3ext(i,:) = [Q_ext3 q_dot_resultsR1E3(i,:)]; 
    q_dot_resultsR1E4ext(i,:) = [Q_ext4 q_dot_resultsR1E4(i,:)]; 
    q_dot_resultsR1E5ext(i,:) = [Q_ext5 q_dot_resultsR1E5(i,:)];
end

t_ext_2 = [(0.001:0.001:0.005) t_wall_2];
t_ext_3 = [(0.001:0.001:0.005) t_wall_3];
t_ext_4 = [(0.001:0.001:0.011) t_wall_4];
t_ext_5 = [(0.001:0.002:0.045) t_wall_5];

clear m_dot_resultsR1E2.mat
clear m_dot_resultsR1E3.mat
clear m_dot_resultsR1E4.mat
clear m_dot_resultsR1E5.mat

clear max_temp_resultsR1E2.mat
clear max_temp_resultsR1E3.mat
clear max_temp_resultsR1E4.mat
clear max_temp_resultsR1E5.mat

clear q_dot_resultsR1E2.mat
clear q_dot_resultsR1E3.mat
clear q_dot_resultsR1E4.mat
clear q_dot_resultsR1E5.mat
clear i 

save extrapolatedR1
