function DSdWire = computeDSdWire(wire)
% COMPUTEDSDWIRE computes the dual edge lengths on the wire.
%
% Input:
%   wire      struct as defined by src/wire.txt
%             required fields: N,N1D,L
%
% Output:
%   DSdWire   cell array of dual 1D length matrices
%             (cell array of length wire.N; each matrix of size N1D-by-N1D)
%
% See also computeKwire, computeMwire, computePwire, computeQwire,
% computeR13, computeR31
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% calculate 1D length matrix for each wire
DSdWire = cell(wire.N,1);
for i = 1:wire.N
    sdWire = [wire.L(i,1)/2, (wire.L(i,1:(end-1))+wire.L(i,2:(end)))/2, wire.L(i,end)/2]';
    DSdWire{i} = spdiags(sdWire,0,wire.N1D,wire.N1D);
end

end