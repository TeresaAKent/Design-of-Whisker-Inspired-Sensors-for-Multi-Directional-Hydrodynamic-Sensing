%%%Bottom plate displacement to moment by serpentine spring system
function [Mx, My] = get_momentw(theta_x, theta_y, dbx, dby, lb, w)
    
    %Define spring variables
    p = 0.00041; % pitch
    t = 0; % Spring thickness
    l = 0.001; % Spring length
    h = 0.0001; % Spring height
    n = 6; % number of turns
    E = 200.0*10^9; %% Young's modulus
    v = 0.35; 
    Gu = 77.2*10^9;
    Ia = w*h*h*h/12;
    Ib = w*w*w*h/12;
    Ib=Ia;
    Ja = w*h*h*h*(1/3-0.21*h/w*(1-h*h*h*h/w/w/w/w/12));
    Jb = w*w*w*h*(1/3-0.21*w/h*(1-w*w*w*w/h/h/h/h/12));
    Jb=Ja;


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
    
    %%Calculate Mx, My
    C1 = inv(C);
    dis = [theta_x; theta_y; dbx; dby];
    C2 = [C1(1,2) C1(1,1) C1(1,3) 0; C1(2,2) C1(2,1) C1(2,3) 0; C1(3,2) C1(3,1) C1(3,3) 0; C1(1,1) C1(1,2) 0 C1(1,3); C1(2,1) C1(2,2) 0 C1(2,3); C1(3,1) C1(3,2) 0 C1(3,3)];
    LB = [0 2 0 2 0 2*lb;2 0 2*lb 0 2 0];
    M = LB*C2*dis;
    Mx = M(1,1);
    My = M(2,1);
end