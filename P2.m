%%%Magentic fulx density from sensor position
%Functions for calculations
function y = P2(K, Pi, gamma)
    y = -gamma/(1-gamma^2)*(Pi-K) - 1/(1-gamma^2)*(gamma^2*Pi-K);
end

