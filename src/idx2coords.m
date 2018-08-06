function [x,y,z] = idx2coords(msh,n)
% IDX2COORDS returns the coordinates of the given canonical index.
%
% Input:
%   msh     struct as defined by src/msh.txt
%           required fields: x,y,z,My,Mz
%   n       canonical index (scalar or vector)
%
% Output:
%   [x,y,z] coordinates of given input indices
%           (scalar, row vector or column vector, according to input)
%
% See also coords2idx, idx2canonical, canonical2idx, canonical4box
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if isrow(n), nIsRow = true; else, nIsRow = false; end

% check input format
if nIsRow
    xmesh = msh.x; ymesh = msh.y; zmesh = msh.z;
else
    xmesh = msh.x'; ymesh = msh.y'; zmesh = msh.z';
end

% assign output
[i,j,k] = canonical2idx(msh,n);
x = xmesh(i);
y = ymesh(j);
z = zmesh(k);

end