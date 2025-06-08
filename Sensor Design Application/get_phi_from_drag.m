% Convert the Drag Force to rotations in x and y
function [pi_w]=get_phi_from_drag(F_d, Theta,M_arm,lb)
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
    Ib=Ia;
    Ja = w*h*h*h*(1/3-0.21*h/w*(1-h*h*h*h/w/w/w/w/12));
    Jb = w*w*w*h*(1/3-0.21*w/h*(1-w*w*w*w/h/h/h/h/12));
    Jb=Ja;

    %%%Define Acquiring Values%%%
    theta = zeros(1);
    phi = zeros(1);
    delta = zeros(1);
    M = zeros(1);
    T = zeros(1);
    Fz = zeros(1);

    %%%Equations for C Matrix%%%
    C = zeros(3);
    C(1,1) = 1/E/Ia*p*(2+n)+1/Gu/Jb*(n*(n-1)*t+2*n*l);
    C(1,2) = 0;
    C(2,1) = 0;
    C(1,3) = -1/E/Ia*p*p*(2+2*n+n*n/2)-1/Gu/Jb*p*(1/6*n*(4*n+7)*(n-1)*t+n*(n+2)*l);
    C(3,1) = C(1,3);
    C(2,2) = 1/Gu/Ja*p*(2+n)+1/E/Ib*((2*n+1)*l+n*(n-1)*t);
    C(2,3) = -1/Gu/Ja*p*t*n/2+1/E/Ib*(-n*(n-1)*t*t/2-n*l*t);
    C(3,2) = C(2,3);
    C(3,3) = 1/E/Ia*p*p*p*(1/3*(2+n)^3)+1/Gu/Ja*p*(1/6*n*(2*n-1)*(n-1)*t*t+n*(n-1)*l*t+n*l*l)+1/E/Ib*(1/6*n*n*(n-1)^2*t*t*t+2/3*n*l*l*l+n*(n-1)*l*l*t+(2/3*n*n*n-n*n+1/3*n)*l*t*t)+1/Gu/Jb*p*p*(1/2*n*(n-1)*(n*n+3*n+3)*t+(2/3*n*n*n+2*n*n+7/3*n)*l);
    M_0=F_d*cos(Theta)*M_arm;
    T_0=F_d*sin(Theta)*M_arm;
    [pi_w,m2,fz]=linearSolve(M_0,C(1,1),C(1,3),C(3,1),C(3,3),lb);
%     options = optimoptions(@fsolve,'Algorithm','trust-region','Display','iter','UseParallel',true,'OptimalityTolerance',1.0000e-3);
%     fun=@root2d;
%     x0=[-m2;-fz;-pi_w;0]
%     x = fsolve(@(x)nonlinear_function(x,M_0,C(1,1),C(1,3),C(3,1),C(3,3),C(2,2),lb),x0,options);
%     pi_w=x(3)
    
    
end

