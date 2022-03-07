function res = DdexpSO3(x,z)
% Directional derivative of the dexp at x in the direction of z 
%
% INPUTS:    x       : 3-vector, element of so(3)
%            z       : 3-vector, element of so(3)
%
% OUTPUTS:   res     : resulting derivative 
%% Script
hatx = hat(x);
hatz = hat(z);
phi = norm(x);

beta = sin(phi/2)^2/((phi/2)^2);
alpha = sin(phi)/phi;

res = 0.5*beta*hatz...
    + 1/phi^2*(1-alpha)*(hatx*hatz+hatz*hatx)...
    + 1/phi^2*(alpha-beta)*(x'*z)*hatx...
    + 1/phi^2*(beta/2-3/phi^2*(1-alpha))*(x'*z)*hatx*hatx;
end