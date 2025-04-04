%%%%Actual Calculation
%clear all;
%%%Whisker displacement in spherical coordinate to bottom plate anular and
%%%positional displacement

%%Parameter Definition
lw = 0.01; %Whisker Length
lb = 0.0005; %Bottom plate edge length/2
lr = 0.001; %Root length

Rm = 0.0005; %Magnet Radius
Lm = 0.001; %Magnet z-length/2

ds = 0.00245; %distance from sensor to center of the magnet

%%Calculation
theta_w = 0;
pi_index_Num = 1;
theta_index_Num = 1;

pi_angle = zeros(1);
theta_angle = zeros(1);

theta_x = zeros(1);
theta_y = zeros(1);
dbx = zeros(1);
dby = zeros(1);
Mx = zeros(1);
My = zeros(1);
PSx = zeros(1);
PSy = zeros(1);
PSz = zeros(1);
BMx = zeros(1);
BMy = zeros(1);
BMz = zeros(1);
Bx = zeros(1);
By = zeros(1);
Bz = zeros(1);

%%Bfield Calculation

for theta_w = 0 : 0.5*pi  : 1*pi
    pi_index_Num = 1;
    theta_angle(theta_index_Num, 1) = theta_w;
    for pi_w = 0 : 5*pi/180 : 10*pi/180
        pi_angle(pi_index_Num, 1) = pi_w;
        [theta_x(pi_index_Num, theta_index_Num), theta_y(pi_index_Num, theta_index_Num)] = get_cart_angle(pi_w, theta_w, lw);
        [dbx(pi_index_Num, theta_index_Num), dby(pi_index_Num, theta_index_Num)] = get_plate_disp(pi_w, theta_w, lb);
        [Mx(pi_index_Num, theta_index_Num), My(pi_index_Num, theta_index_Num)] = get_moment(theta_x(pi_index_Num, theta_index_Num), theta_y(pi_index_Num, theta_index_Num), dbx(pi_index_Num, theta_index_Num), dby(pi_index_Num, theta_index_Num), lb);

        [PSx(pi_index_Num, theta_index_Num), PSy(pi_index_Num, theta_index_Num), PSz(pi_index_Num, theta_index_Num)] = senP_from_mag(pi_w, theta_w, lr, Lm, ds);
        [B1 B2 B3] = CylMag(PSx(pi_index_Num, theta_index_Num),PSy(pi_index_Num, theta_index_Num),PSz(pi_index_Num, theta_index_Num),Rm,Lm);
        
        BMx(pi_index_Num, theta_index_Num) = B1;
        BMy(pi_index_Num, theta_index_Num) = B2;
        BMz(pi_index_Num, theta_index_Num) = B3;
        
        [Bx(pi_index_Num, theta_index_Num), By(pi_index_Num, theta_index_Num), Bz(pi_index_Num, theta_index_Num)] = magB_to_senB(BMx(pi_index_Num, theta_index_Num), BMy(pi_index_Num, theta_index_Num), BMz(pi_index_Num, theta_index_Num), pi_w, theta_w);
        pi_index_Num = pi_index_Num+1;
    end
    theta_index_Num = theta_index_Num+1;
end


%%%%Functions Definitions%%%%

% Get the Drag Force
function [M_d]=get_drag_moment(rho,v,C_d,A,M_arm)
    F_d = 1/2*rho*v^2*C_d*A;
    M_d = F_d*M_arm
end

% Convert the Drag Force to rotations in x and y
function [pi_w]=get_phi_from_drag(M_d, Theta)
    %Define spring variables
    p = 0.0003;
    t = 0;
    l = 0.0005;
    w = 0.000140;
    h = 0.0005;
    n = 12;
    E = 3.2*10^9; %% Young's modulus
    v = 0.35; 
    Gu = E/2/(1+v);
    Ia = w*h*h*h/12;
    Ib = w*w*w*h/12;
    Ja = w*h*h*h*(1/3-0.21*h/w*(1-h*h*h*h/w/w/w/w/12));
    Jb = w*w*w*h*(1/3-0.21*w/h*(1-w*w*w*w/h/h/h/h/12));

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
    
    fun=@root2d;
    x0=[0,0,0];
    x = fsolve(fun,x0);
    pi_w=x(3);
    
end

function F=root2d(x,C11,C13,C31,C33,lb,M_d)
    F(1)=2*x(1)+2*lb*x(2)-M_d;
    F(2)=C11*x(1)+C13*x(2)-x(3);
    F(3)=C31*x(1)+C33*x(2)-lb*sin(x(3));
end


function [theta_x, theta_y] = get_cart_angle(pi_w, theta_w, lw)
    theta_x = asin(lw*sin(pi_w)*sin(theta_w)/sqrt((lw*sin(pi_w)*sin(theta_w))^2+(lw*cos(pi_w))^2));
    theta_y = asin(lw*sin(pi_w)*cos(theta_w)/sqrt((lw*sin(pi_w)*cos(theta_w))^2+(lw*cos(pi_w))^2));
end

function [db1, db2] = get_plate_disp(pi_w, theta_w, lb)
    W = [0 0 -cos(theta_w); 0 0 -sin(theta_w); cos(theta_w) sin(theta_w) 0];
    R = eye(3) + sin(pi_w)*W + 2*sin(pi_w/2)^2*W*W;
    x1 = [lb, 0, 0];
    y1 = [0, lb, 0];
    D1 = R*x1';
    D2 = R*y1';
    db1 = D1(3,1);
    db2 = D2(3,1);
end

%%%Bottom plate displacement to moment by serpentine spring system

function [Mx, My] = get_moment(theta_x, theta_y, dbx, dby, lb)
    
    %Define spring variables
    p = 0.0003;
    t = 0;
    l = 0.0005;
    w = 0.000140;
    h = 0.0005;
    n = 12;
    E = 3.2*10^9; %% Young's modulus
    v = 0.35; 
    Gu = E/2/(1+v);
    Ia = w*h*h*h/12;
    Ib = w*w*w*h/12;
    Ja = w*h*h*h*(1/3-0.21*h/w*(1-h*h*h*h/w/w/w/w/12));
    Jb = w*w*w*h*(1/3-0.21*w/h*(1-w*w*w*w/h/h/h/h/12));

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
    
    %%Calculate Mx, My
    C1 = inv(C);
    dis = [theta_x; theta_y; dbx; dby];
    C2 = [C1(1,2) C1(1,1) C1(1,3) 0; C1(2,2) C1(2,1) C1(2,3) 0; C1(3,2) C1(3,1) C1(3,3) 0; C1(1,1) C1(1,2) 0 C1(1,3); C1(2,1) C1(2,2) 0 C1(2,3); C1(3,1) C1(3,2) 0 C1(3,3)];
    LB = [0 2 0 2 0 2*lb;2 0 2*lb 0 2 0];
    M = LB*C2*dis;
    Mx = M(1,1);
    My = M(2,1);
end

%%%Whisker root displacement to Magnetic flux density

%%Angle to position
function [PSx, PSy, PSz] = senP_from_mag(pi_w, theta_w, lr, lm, ds)
    d1 = lr + lm;
    d2 = lr + 2*lm + ds;
    PSx = d2*sin(pi_w)*cos(theta_w);
    PSy = d2*sin(pi_w)*sin(theta_w);
    PSz = d1 - d2*cos(pi_w);
end

%%%Magentic fulx density from sensor position
%Functions for calculations
function y = P1(K, epsilon, k)
    y = K-2/(1-k^2)*(K-epsilon);
end

function y = P2(K, Pi, gamma)
    y = -gamma/(1-gamma^2)*(Pi-K) - 1/(1-gamma^2)*(gamma^2*Pi-K);
end

%Bx, By, Bz Calculation
function [BMx, BMy, BMz] = CylMag(x,y,z,R,L)
    m0 = 4*pi*10^(-7);
    Br = 1.4;
    %Ms = 800000;
    Ms = 1073400;
    
    lo = sqrt(x^2+y^2);
    
    if lo == 0
    B_lo = 0;
    B_z = m0*Ms/2*((z+L)/sqrt((z+L)^2+R^2) - (z-L)/sqrt((z-L)^2+R^2));
    
    BMx = B_lo;
    BMy = B_lo;
    BMz = B_z;
    else

    %lo is not 0 case

    ksi = [z+L, z-L];
    alpha = [1/sqrt(ksi(1,1)^2+(lo+R)^2) , 1/sqrt(ksi(1,2)^2+(lo+R)^2)];
    beta = [alpha(1,1)*ksi(1,1) , alpha(1,2)*ksi(1,2)];
    gamma = (lo-R)/(lo+R);
    k = [sqrt((ksi(1,1)^2+(lo-R)^2)/(ksi(1,1)^2+(lo+R)^2)) , sqrt((ksi(1,2)^2+(lo-R)^2)/(ksi(1,2)^2+(lo+R)^2))];

    B_lo = 0;
    B_z = 0;

    K = ellipticK(1-k(1,1)^2);
    epsilon = ellipticE(1-k(1,1)^2);
    Pi = ellipticPi(1-gamma^2, 1-k(1,1)^2);

    B_lo = B_lo + m0*Ms*R/pi*alpha(1,1)*P1(K,epsilon,k(1,1));
    B_z = B_z + m0*Ms*R/pi/(lo+R)*beta(1,1)*P2(K,Pi,gamma);

    K = ellipticK(1-k(1,2)^2);
    epsilon = ellipticE(1-k(1,2)^2);
    Pi = ellipticPi(1-gamma^2, 1-k(1,2)^2);

    B_lo = B_lo - m0*Ms*R/pi*alpha(1,2)*P1(K,epsilon,k(1,2));
    B_z = B_z - m0*Ms*R/pi/(lo+R)*beta(1,2)*P2(K,Pi,gamma);
    
    BMx = B_lo * x/sqrt(x^2+y^2);
    BMy = B_lo * y/sqrt(x^2+y^2);
    BMz = B_z;
    
    end
end


%%Convert B vector from magnet coordinate to sensor coordinate
function [BSx, BSy, BSz] = magB_to_senB(BMx, BMy, BMz, pi_w, theta_w)
    Bsm = zeros(3);
    W = [0 0 -cos(theta_w); 0 0 -sin(theta_w); cos(theta_w) sin(theta_w) 0];
    R = eye(3) + sin(pi_w)*W + 2*sin(pi_w/2)^2*W*W;
    Bsm(:,1) = R*[1;0;0];
    Bsm(:,2) = R*[0;1;0];
    Bsm(:,3) = R*[0;0;1];
    b = [BMx, BMy, BMz]';
    a = Bsm*b;
    BSx = a(1,1);
    BSy = a(2,1);
    BSz = a(3,1);
end


