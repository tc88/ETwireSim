function Kwire = computeKwire(wire,Mwire)
% COMPUTEKWIRE generates the wire matrix. This method includes wire.N wires
% that each consist of N1D 1D-points whereas each 1D-point is connected to
% the 3D grid using 4 points in the neighborhood. Note that N1D can be
% different for different wires.
%
% Input:
%   wire            struct as defined by src/wire.txt
%                   required fields: N,Ps,R13,R31
%                   optional fields: select (default: ones(wire.N,1))
%   Mwire           conductance of the considered wires (wire.N-by-N1D)
%
% Output:
%   Kwire   wire stiffness matrix (np-by-np)
%
% See also computeDSdWire, computeMwire, computePwire, computeQwire,
% computeR13, computeR31
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if ~isfield(wire,'select'), wire.select = ones(wire.N,1); end

selectedWires = find(wire.select);
if isnumeric(wire.R31), np = size(wire.R31,1); else, np = size(wire.R31{selectedWires(1)},1); end
Kwire = sparse(np,np);
for i = 1:wire.N
    if ~wire.select(i), continue; end
    % extract R13 and R31 from cell if not directly stored as matrix
    if isnumeric(wire.R13), R13 = wire.R13; else, R13 = wire.R13{i}; end
    if isnumeric(wire.R31), R31 = wire.R31; else, R31 = wire.R31{i}; end
    MwireCurrent = Mwire{i};
    currentStamp = R31*wire.Ps{i}'*MwireCurrent*wire.Ps{i}*R13;
    Kwire = Kwire + currentStamp;
end

end