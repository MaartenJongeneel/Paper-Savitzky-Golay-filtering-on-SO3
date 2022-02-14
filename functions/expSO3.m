function res = expSO3(a)
% Computes the exponential mapping on SO(3)
%
% INPUTS:    a       : 3 vector, isomorphism to element of so(3)
%
% OUTPUTS:   res     : element of SO(3)
%
%% Exponential on SO(3), Rodrigues formula
phi = norm(a);
if phi ~=0
    res = eye(3) + sin(phi)/phi*hat(a) + (1-cos(phi))/phi^2*hat(a)*hat(a);
else
    res = eye(3);
end
end