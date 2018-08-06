function dists = dist2line(msh,x,y)
% DIST2LINE computes the orthogonal distance of each point in a given mesh
% to a line in z-direction defined by x and y.
%
% Input:
%   msh     struct as defined by src/msh.txt
%           required fields: np,x,y,z,Mx,My,Mz
%   x       x-coordinate of reference line (optional, default: 0)
%   y       y-coordinate of reference line (optional, default: 0)
%
% Output:
%   dists   vector of orthogonal distances to specified line (np-by-1)
%
% See also distFIT
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 3, y = 0; end
if nargin < 2, x = 0; end

% get i- and j-index of reference line
[~,iLine] = ismember(x,msh.x);
[~,jLine] = ismember(y,msh.y);
if isempty(iLine), error('1D subdomain is not part of xmesh'); end
if isempty(jLine), error('1D subdomain is not part of ymesh'); end

% compute distance to reference line
dists = zeros(msh.np,1);
for ipn = 1:msh.np
    [~,~,kk] = canonical2idx(msh,ipn);
    ipnCenter = idx2canonical(msh,iLine,jLine,kk);
    dists(ipn) = distFIT(msh,ipn,ipnCenter);
end

end