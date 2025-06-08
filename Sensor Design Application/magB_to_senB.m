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
