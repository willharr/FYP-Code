function [h] = heatVaporisation(fluid)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

h = enthalpyVap(fluid)-enthalpyLiquid(fluid);

h = 1000*h; %convert to J/kg
end

