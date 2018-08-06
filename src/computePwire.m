function Pwire = computePwire(wire)
% CREATEPWIRE creates the differential operator Pwire.
%
% Input:
%   wire    struct as defined by src/wire.txt
%           required fields: N,N1D
%
% Output:
%   Pwire   cell array of 1D wire derivative matrix
%           (cell array of length wire.N; each matrix of size N1D-by-N1D)
%
% See also computeDSdWire, computeKwire, computeMwire, computeQwire,
% computeR13, computeR31, createPx
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

N1D = wire.N1D;
Pwire = cell(wire.N,1);
for i = 1:wire.N
    rows = [1:N1D,1:(N1D-1)];
    columns = [1:N1D,2:N1D];
    values = [-ones(1,N1D),ones(1,N1D-1)];
    Pwire{i} = sparse(rows,columns,values,N1D,N1D);
end

end