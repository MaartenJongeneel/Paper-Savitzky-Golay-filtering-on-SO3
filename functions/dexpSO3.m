function res = dexpSO3(a)
% Computes the right trivialized tangent d exp on SO(3)
%
% INPUTS:    a       : 3 vector, isomorphism to element of so(3)
%
% OUTPUTS:   res     : diff exponential of a 
%
%% Closed form expression of right trivialized tangent d exp on SO(3)
phi = norm(a);
ahat = hat(a);
beta = sin(phi/2)^2/((phi/2)^2);
alpha = sin(phi)/phi;
res = eye(3)+0.5*beta*ahat+1/(phi^2)*(1-alpha)*ahat*ahat;
end