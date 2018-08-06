function Mwire = computeMwire(wire,alphaWire)
% COMPUTEMWIRE generates the wire mass matrix.
%
% Input:
%   wire        struct as defined by src/wire.txt
%               required fields: N,N1D,L
%   alphaWire   material property of wire (scalar)
%
% Output:
%   Mwire   cell array of wire mass matrices
%           (cell array of length wire.N; each matrix of size N1D-by-N1D)
%
% See also computeDSdWire, computeKwire, computePwire, computeQwire,
% computeR13, computeR31
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

Mwire = cell(wire.N,1);
for i=1:wire.N
    Mwire{i} = spdiags([alphaWire./wire.L(i,:)';0],0,wire.N1D,wire.N1D);
end

end