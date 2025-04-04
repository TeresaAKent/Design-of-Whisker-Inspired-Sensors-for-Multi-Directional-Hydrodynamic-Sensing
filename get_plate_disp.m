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
