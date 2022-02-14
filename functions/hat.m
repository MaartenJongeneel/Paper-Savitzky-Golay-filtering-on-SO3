function res = hat(vec)
% Take the 3- or 6-vector representing an isomorphism of so(3) or se(3) and
% writes this as element of so(3) or se(3). 
%
% INPUTS:    vec     : 3- or 6-vector. Isomorphism of so(3) or se(3)
%
% OUTPUTS:   res     : element of so(3) or se(3)
%
%% Hat operator
if length(vec) == 3
    res = [0, -vec(3), vec(2); vec(3), 0, -vec(1); -vec(2), vec(1), 0];
elseif length(vec) == 6
    skew = [0, -vec(6), vec(5); vec(6), 0, -vec(4); -vec(5), vec(4), 0];
    v = [vec(1);vec(2);vec(3)];
    res = [skew, v; zeros(1,4)];
end

