function [F_d]=get_drag_moment(rho,v,C_d,A,M_arm)
    F_d = 1/2*rho*v^2*C_d*A*sign(v);
end

