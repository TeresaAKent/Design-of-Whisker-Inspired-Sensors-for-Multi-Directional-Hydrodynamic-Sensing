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