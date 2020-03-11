function [p] = vapPressure(fluid)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if fluid.T >= fluid.T_crit
    
   disp('ERROR: T out of range') 
   return
    
end

Tr = fluid.T/fluid.T_crit;
Tr1 = 1-Tr;

b1 = fluid.b_vapPressure(1);
b2 = fluid.b_vapPressure(2);
b3 = fluid.b_vapPressure(3);
b4 = fluid.b_vapPressure(4);


p = fluid.P_crit*exp((Tr^-1)*(b1*Tr1 + b2*(Tr1^(3/2))...
+ b3*(Tr1^(5/2)) + b4*(Tr1^(5))));

end

