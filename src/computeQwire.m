function Qwire = computeQwire(phi3D,wire)
% COMPUTEQWIRE computes dissipated heat of wires.
%
% Input:
%   phi3D   3D solution of electric potential (np-by-1)
%   wire    struct as defined by src/wire.txt
%           required fields: N,N1D,select,Msigma,Ps,R13,R31
%
% Output:
%   Qwire   dissipated heat of wires (np-by-1)
%
% See also computeDSdWire, computeKwire, computeMwire, computePwire,
% computeR13, computeR31
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

np = size(phi3D,1);
Qwire = zeros(np,1);
QonWire = zeros(wire.N,wire.N1D);
for i = 1:wire.N
    if ~wire.select(i), continue; end
    Msigma1D = wire.Msigma{i};
    phi1D = wire.R13{i}*phi3D;
    Ps = wire.Ps{i};
    I = -Msigma1D*Ps*phi1D;
    QwireOnEdges = -I.*(Ps*phi1D);
    for j = 1:wire.N1D
        if j == 1
            QonWire(i,j) = 0.5*QwireOnEdges(j);
        elseif j == wire.N1D
            QonWire(i,j) = 0.5*QwireOnEdges(j-1);
        else
            QonWire(i,j) = 0.5*QwireOnEdges(j)+0.5*QwireOnEdges(j-1);
        end
    end
    QwireThis = wire.R31{i}*QonWire(i,:)';
    Qwire = Qwire + QwireThis;
end

end