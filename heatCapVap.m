function [c] = heatCapVap(fluid)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Tr = fluid.T/fluid.T_crit;
Tr1 = 1-Tr;

b1 = fluid.b_heatCapVap(1);
b2 = fluid.b_heatCapVap(2);
b3 = fluid.b_heatCapVap(3);
b4 = fluid.b_heatCapVap(4);
b5 = fluid.b_heatCapVap(5);

c = b1*(1 + b2*(Tr1^(-2/3))...
+ b3*(Tr1^(-1/3)) + b4*(Tr1^(1/3)) + b5*(Tr1^(2/3)));

c = c*1000; %convert to J/kg*K
end

