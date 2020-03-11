function [nozzle] = nozzleGeometry(nozzle)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %% Inputs

    type = 'conical';

    Pcc = nozzle.Pcc; %combustion chamber pressure
    gamma = nozzle.gamma; %specific heat ratio
    Patm = nozzle.Patm; %ambient (exit) pressure
    Rcc = nozzle.Rcc; %combustion chamber diameter
    m_dot = nozzle.m_dot; %kg/s, total mass flow rate
    m_molar = nozzle.m_molar; %kg/mol of reaction products
    R_gas = 8.314/m_molar; %specific gas constant for exhaust gases
    Tcc = nozzle.Tcc;
    c_star = nozzle.c_star;
    g = nozzle.g;

    %% General Properties

    exp_ratio = ( ((2/(gamma+1))^(1/(gamma-1))) * ((Pcc/Patm)^(1/gamma)) )*(sqrt( ((gamma+1)/(gamma-1)) * (1-((Patm/Pcc)^((gamma-1)/gamma))) )^-1); %expansion ratio, Ae/At
    Pthroat = Pcc*(2/(gamma+1))^(gamma/(gamma-1)); %pressure at throat

    Acc = pi*Rcc^2;
    At = m_dot*(sqrt(Tcc)/Pcc)*sqrt(R_gas/gamma)*((gamma+1)/2)^((gamma+1)/(2*(gamma-1)));
    Ae = exp_ratio*At;

    Rt = sqrt(At/pi);
    Re = sqrt(Ae/pi);
    Rcc = sqrt(Acc/pi);
    R = Rt;


    %% Conical Nozzle - Radiused
    % Generates cross sectional plot of nozzle contour and array of nozzle
    % areas with respect to x
    if isequal(type,'conical')

        alpha = 15; %deg, nozzle diverging half-angle
        beta = 30; %deg, nozzle converging half-angle

        x_radIn = -R*sind(beta);
        x_radOut = R*sind(alpha);
        x_inlet = x_radIn - (Rcc-R*(1-cosd(beta))-Rt)/tand(beta);
        x_outlet = x_radOut + (Re-R*(1-cosd(alpha))-Rt)/tand(alpha);

        n = 20;
        x_arr1 = linspace(x_inlet,x_radIn,n);
        x_arr2 = linspace(x_radIn,x_radOut,n);
        x_arr3 = linspace(x_radOut,x_outlet,n);
        x_array = [x_arr1 x_arr2 x_arr3];

        y_con = zeros(1,3*n);
        A_con = zeros(1,3*n);

        for i = 1:n

            y_con(i) = tand(180-beta)*(x_arr1(i)-x_inlet) + Rcc;
            y_con(i+n) = -sqrt((R^2) - (x_arr2(i)^2)) + Rt+R;
            y_con(i+2*n) = tand(alpha)*(x_arr3(i)-x_radOut) + Rt+R*(1-cosd(alpha));

        end

        for i = 1:3*n


            A_con(i) = pi.*y_con(i).^2;

            if x_array(i) <= 0
                [M_con(i), T_con(i), P_con(i), rho_con(i) ,~] = flowisentropic(gamma,A_con(i)/At,'sub');
            else
                [M_con(i), T_con(i), P_con(i), rho_con(i) ,~] = flowisentropic(gamma,A_con(i)/At,'sup');
            end

            V_con(i) = M_con(i)*sqrt(gamma*R_gas*T_con(i));

        end

    end
    
    nozzle.A_array = A_con;
    nozzle.M_array = M_con;
    nozzle.Rt = Rt;
    nozzle.R = R;
    nozzle.x_array = x_array;
    nozzle.At = At;
    nozzle.y_array = y_con;
end

