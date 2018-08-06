function R31 = computeR31(np,wire,verbose)
% COMPUTER31 computes R31 which is the transpose of R13 for a symmetric
% coupling and just a mapping of the 1D points into the 3D space for a
% non-symmetric coupling.
%
% Input:
%   np          number of primary points in mesh
%   wire        struct as defined by src/wire.txt
%               required fields: N,N1D,idx,select
%   verbose     triggers console outputs and plots
%               (optional, default: [1 0])
%
% Output:
%   R31    projects from 1D solution to 3D solution
%          (cell array of length wire.N; each matrix of size np-by-N1D)
%
% See also computeDSdWire, computeKwire, computeMwire, computePwire,
% computeQwire, computeR13
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 3, verbose = [1,0]; end

if verbose(1), fprintf('calculation of R31 ...\n'); end

% assemble R31
R31 = cell(wire.N,1);
for i = 1:wire.N
    if ~wire.select(i), continue; end
    R31this = sparse(np,wire.N1D);
    for k = 1:wire.N1D
        R31this(wire.idx(i,k),k) = 1;
    end
    R31{i} = R31this;
end

if verbose(1), fprintf('calculation of R31 done!\n'); end

end