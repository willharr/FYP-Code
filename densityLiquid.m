function [rho] = densityLiquid(fluid)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Tr = fluid.T/fluid.T_crit;
Tr1 = 1-Tr;

b1 = fluid.b_densityLiquid(1);
b2 = fluid.b_densityLiquid(2);
b3 = fluid.b_densityLiquid(3);
b4 = fluid.b_densityLiquid(4);

rho = fluid.rho_crit*exp(b1*(Tr1^(1/3)) + b2*(Tr1^(2/3))...
+ b3*(Tr1^(1)) + b4*(Tr1^(4/3)));

end

