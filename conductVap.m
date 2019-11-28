function [c] = conductVap(fluid)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Tr = fluid.T/fluid.T_crit;
Tr1 = 1-Tr;

b1 = fluid.b_conductVap(1);
b2 = fluid.b_conductVap(2);
b3 = fluid.b_conductVap(3);
b4 = fluid.b_conductVap(4);

c = exp(b1 + b2*(Tr1^(-2/3))...
+ b3*(Tr1^(-1/3)) + b4*(Tr1^(1/3)) + b5*(Tr1^(2/3)));

c = c/1000; % convert to W/m*K
end

