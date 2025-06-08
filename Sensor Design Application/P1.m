%%%Magentic fulx density from sensor position
%Functions for calculations
function y = P1(K, epsilon, k)
    y = K-2/(1-k^2)*(K-epsilon);
end

