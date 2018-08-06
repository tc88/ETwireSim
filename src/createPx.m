function Px = createPx(np)
% CREATEPX creates the grid differential operator Px
%
% Input:
%   np  number of primary points in the mesh
%
% See also computePwire
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

rows = [1:np,1:(np-1)];
columns = [1:np,2:np];
values = [-ones(1,np),ones(1,np-1)];
Px = sparse(rows,columns,values,np,np);

end