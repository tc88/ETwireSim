function [i,j,k] = canonical2idx(msh,n)
% CANONICAL2IDX translates the given canonical index n into the indices
% (i,j,k) for the cartesian coordinate system.
%
% Input:
%   msh     struct as defined by src/msh.txt
%           required fields: My,Mz
%   n       can be a scalar or a vector of canonical indices
%
% Output:
%   i       returns the x-index for the given canonical index/indices
%   j       returns the y-index for the given canonical index/indices
%   k       returns the z-index for the given canonical index/indices
%
% See also coords2idx, idx2coords, idx2canonical, canonical4box
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

k = 1+floor((n-1)/msh.Mz);
n2D = 1+mod(n-1,msh.Mz);
j = 1+floor((n2D-1)/msh.My);
i = 1+mod(n2D-1,msh.My);

end