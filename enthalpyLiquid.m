function [h] = enthalpyLiquid(fluid)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Tr = fluid.T/fluid.T_crit;
Tr1 = 1-Tr;

b1 = fluid.b_enthalpyLiquid(1);
b2 = fluid.b_enthalpyLiquid(2);
b3 = fluid.b_enthalpyLiquid(3);
b4 = fluid.b_enthalpyLiquid(4);
b5 = fluid.b_enthalpyLiquid(5);

h = b1*(Tr1^(0)) + b2*(Tr1^(1/3))...
+ b3*(Tr1^(2/3)) + b4*(Tr1^(1)) + b5*(Tr1^(4/3));

h = 1000*h;
end

