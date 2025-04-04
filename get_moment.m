%%%Bottom plate displacement to moment by serpentine spring system
function [Mx, My] = get_moment(theta_x, theta_y, dbx, dby, lb)
    
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
    % Ib=Ia;
    Ja = w*h*h*h*(1/3-0.21*h/w*(1-h*h*h*h/w/w/w/w/12));
    Jb = w*w*w*h*(1/3-0.21*w/h*(1-w*w*w*w/h/h/h/h/12));
    % Jb=Ja;

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

    %%Calculate Mx, My
    C1 = inv(C);
    dis = [theta_x; theta_y; dbx; dby];
    C2 = [C1(1,2) C1(1,1) C1(1,3) 0; C1(2,2) C1(2,1) C1(2,3) 0; C1(3,2) C1(3,1) C1(3,3) 0; C1(1,1) C1(1,2) 0 C1(1,3); C1(2,1) C1(2,2) 0 C1(2,3); C1(3,1) C1(3,2) 0 C1(3,3)];
    LB = [0 2 0 2 0 2*lb;2 0 2*lb 0 2 0];
    M = LB*C2*dis;
    Mx = M(1,1);
    My = M(2,1);
end