function res = vee(mat)
% Takes an element of so(3) or se(3) and returns its isomorphism in R^n.
%
% INPUTS:    mat     : element of so(3) or se(3)
%
% OUTPUTS:   res     : 3- or 6-vector. Isomorphism of so(3) or se(3)
%
%% Vee operator

xi1 = (mat(3,2)-mat(2,3))/2;
xi2 = (mat(1,3)-mat(3,1))/2;
xi3 = (mat(2,1)-mat(1,2))/2;

if length(mat) == 3
   res = [xi1; xi2; xi3];
elseif length(mat) == 4
   res = [mat(1:3,4);xi1;xi2;xi3]; 
end


end
