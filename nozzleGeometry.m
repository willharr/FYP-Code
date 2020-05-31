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
    tw = nozzle.t_wall;

    %% General Properties

    exp_ratio = ( ((2/(gamma+1))^(1/(gamma-1))) * ((Pcc/Patm)^(1/gamma)) )*(sqrt( ((gamma+1)/(gamma-1)) * (1-((Patm/Pcc)^((gamma-1)/gamma))) )^-1); %expansion ratio, Ae/At
    Pthroat = Pcc*(2/(gamma+1))^(gamma/(gamma-1)); %pressure at throat

    Acc = pi*Rcc^2;
    At = m_dot*(sqrt(Tcc)/Pcc)*sqrt(R_gas/gamma)*((gamma+1)/2)^((gamma+1)/(2*(gamma-1)));
    Ae = exp_ratio*At;

    Rt = sqrt(At/pi); %throat radius
    Re = sqrt(Ae/pi);
    Rcc = sqrt(Acc/pi);
    R = Rt; %throat curvature radius


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
        x_array = [x_arr1 x_arr2(2:end) x_arr3(2:end)];

        y_con = zeros(1,(3*n)-2);
        A_con = zeros(1,(3*n)-2);

        for i = 1:n

            y_con(i) = tand(180-beta)*(x_array(i)-x_inlet) + Rcc;
            
        end
        
        for i = 1:n-1
            y_con(i+n) = -sqrt((R^2) - (x_array(i+n)^2)) + Rt+R;
            y_con(i-1+(2*n)) = tand(alpha)*(x_array(i-1+(2*n))-x_radOut) + Rt+R*(1-cosd(alpha));
        end

        for i = 1:(3*n)-2
            
            A_con(i) = pi.*y_con(i).^2;

            if x_array(i) <= 0
                [M_con(i), T_con(i), P_con(i), rho_con(i) ,~] = flowisentropic(gamma,A_con(i)/At,'sub');
            else
                [M_con(i), T_con(i), P_con(i), rho_con(i) ,~] = flowisentropic(gamma,A_con(i)/At,'sup');
            end

            V_con(i) = M_con(i)*sqrt(gamma*R_gas*T_con(i));

        end

    end
    
    % y_con breakdown: 1->n is converging, n+1->2*n-1 is radius, 2*n->3*n-2
    y_inlet = y_con(1);
    y_radIn = y_con(n);
    y_radOut = y_con((n*2)-1);
    y_outlet = y_con(end);
    
    %% Nozzle Geometry For Thermal Model
    
    % full nozzle
        tcb = tw*cosd(90-beta);
        tca = tw*cosd(90-alpha);
        tsb = tw*sind(90-beta);
        tsa = tw*sind(90-alpha);

        rect_con = [3 4 x_inlet+tcb x_radIn+tcb x_radIn x_inlet y_inlet+tsb y_radIn+tsb y_radIn y_inlet];
        rect_div = [3 4  x_radOut-tca x_outlet-tca x_outlet x_radOut y_radOut+tsa y_outlet+tsa y_outlet y_radOut];
        circ_inner = [1 0 Rt+R R 0 0 0 0 0 0];
        circ_outer = [1 0 Rt+R R-tw 0 0 0 0 0 0];
        rect_rad = [2 4 0 x_radOut 0 x_radIn Rt+R y_radOut 0 y_radIn];

        gd = [rect_con' rect_div' circ_inner' circ_outer' rect_rad'];

        ns = char('rcon','rdiv','ci','co','rrad');
        ns = ns';

        sf = '(rcon+rdiv) + (rrad*(ci-co))';
    
    %% Outer Surface Points array
    tw = tw-(0.1*tw);
    tcb = tw*cosd(90-beta);
    tca = tw*cosd(90-alpha);
    tsb = tw*sind(90-beta);
    tsa = tw*sind(90-alpha);
    
    out_xarr1 = x_arr1+tcb;
    out_xarr2 = linspace(x_radIn+tcb,x_radOut-tca,n);
    out_xarr2 = out_xarr2(2:end);
    out_xarr3 = x_arr3(2:end)-tca;
    
    out_yarr1 = y_con(1:n)+tsb;
    out_yarr3 = y_con((2*n):end)+tsa;
    for i = 1:n-1
        out_yarr2(i) = -sqrt(((R-tw)^2) - (out_xarr2(i)^2)) + Rt+R;
    end
    
    out_xarray = [out_xarr1 out_xarr2 out_xarr3];
    out_yarray = [out_yarr1 out_yarr2 out_yarr3];
    
    %% Outputs
    nozzle.geomBounds = [x_inlet, x_radIn, x_radOut, x_outlet; y_inlet, y_radIn, y_radOut, y_outlet];
    nozzle.outerGeomBounds = [x_inlet+tcb, x_radIn+tcb, x_radOut-tca, x_outlet-tca, y_inlet+tsb, y_radIn+tsb, y_radOut+tsa, y_outlet+tsa];
    nozzle.angles = [alpha beta];
    nozzle.geom = gd;
    nozzle.nameList = ns;
    nozzle.setFunc = sf;
    nozzle.S = sum(pi.*(y_con.^2));
    nozzle.A_array = A_con;
    nozzle.M_array = M_con;
    nozzle.Rt = Rt;
    nozzle.R = R;
    nozzle.x_array = x_array;
    nozzle.outer_x_array = out_xarray;
    nozzle.outer_y_array = out_yarray;
    nozzle.At = At;
    nozzle.y_array = y_con;
end

