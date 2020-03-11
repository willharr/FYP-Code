function [c] = conductLiquid(fluid)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if fluid.T >= fluid.T_crit
    
   disp('ERROR: T out of range') 
   return
    
end

Tr = fluid.T/fluid.T_crit;
Tr1 = 1-Tr;

b1 = fluid.b_conductLiquid(1);
b2 = fluid.b_conductLiquid(2);
b3 = fluid.b_conductLiquid(3);
b4 = fluid.b_conductLiquid(4);

c = b1*(1 + b2*(Tr1^(1/3))...
+ b3*(Tr1^(2/3)) + b4*(Tr1^(1)));

c = c/1000; %convert to W/m*K
end

