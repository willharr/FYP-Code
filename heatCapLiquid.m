function [c] = heatCapLiquid(fluid)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Tr = fluid.T/fluid.T_crit;
Tr1 = 1-Tr;

b1 = fluid.b_heatCapLiquid(1);
b2 = fluid.b_heatCapLiquid(2);
b3 = fluid.b_heatCapLiquid(3);
b4 = fluid.b_heatCapLiquid(4);
b5 = fluid.b_heatCapLiquid(5);

c = b1*(1 + b2*(Tr1^(-1))...
+ b3*(Tr1^(1)) + b4*(Tr1^(2)) + b5*(Tr1^(3)));

c = c*1000; %convert to J/kg*K
end

