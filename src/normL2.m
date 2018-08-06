function res = normL2(u,D)
% NORML2 calculates the L2-norm of u.
%
% Input:
%   u    solution vector (Np-by-1)
%   D    matrix containing metric information (Np-by-Np)
%
% Output:
%   res  computed L2 norm
%
% See also norm4box
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

res = sqrt(u'*D*u);

end