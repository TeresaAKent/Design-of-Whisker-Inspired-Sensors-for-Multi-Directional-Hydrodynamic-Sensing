% Convert the Drag Force to rotations in x and y
function [pi_w]=get_phi_from_drag2(F_d, Theta,M_arm,lb)
    %Define spring variables
    p = 0.00041;
    t = 0;
    l = 0.001;
    w = 0.00025;
    h = 0.0001;
    n = 6;
    E = 200.0*10^9; %% Young's modulus
    v = 0.35; 
    Gu = 77.2*10^9;
    Ia = w*h*h*h/12;
    Ib = w*w*w*h/12;
    %Ib=Ia;
    Ja = w*h*h*h*(1/3-0.21*h/w*(1-h*h*h*h/w/w/w/w/12));
    Jb = w*w*w*h*(1/3-0.21*w/h*(1-w*w*w*w/h/h/h/h/12));
    %Jb=Ja;

    %%%Define Acquiring Values%%%
    theta = zeros(1);
    phi = zeros(1);
    delta = zeros(1);
    M = zeros(1);
    T = zeros(1);
    Fz = zeros(1);

    %%%Equations for C Matrix%%%
    C = zeros(3);
    C(1,1) = 8*p/(E*Ia)+6*l/(Gu*Jb);
    C(1,2) = 0;
    C(2,1) = 0;
    C(1,3) = -32*p^2/(E*Ia)-48*l*p/(Gu*Jb);
    C(3,1) = C(1,3);
    C(2,2) = 8*p/Gu/Ja+13*l/E/Ib;
    C(2,3) = 0;
    C(3,2) = 0;
    C(3,3) = 512*p^3/3/E/Ia+6*l^2/Gu/Ja+4*l^3/E/Ib+230*l/Gu/Jb;


    M_0=F_d*cos(Theta)*M_arm;
    T_0=F_d*sin(Theta)*M_arm;
    denom=1/C(1,1)+1/C(2,2)+(lb^2-lb*C(3,1)/C(1,1))/C(3,3);
    %denom=1/C(1,1)+1/C(2,2);
    pi_w=M_0/2*1/denom;
%     options = optimoptions(@fsolve,'Algorithm','trust-region','Display','iter','UseParallel',true,'OptimalityTolerance',1.0000e-3);
%     fun=@root2d;
%     x0=[-m2;-fz;-pi_w;0]
%     x = fsolve(@(x)nonlinear_function(x,M_0,C(1,1),C(1,3),C(3,1),C(3,3),C(2,2),lb),x0,options);
%     pi_w=x(3)
    
    
end

