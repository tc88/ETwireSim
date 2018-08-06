function R13 = computeR13(msh,wire,verbose)
% COMPUTER13 calculates R13 matrix that is a representation of the integral
% phi1D = 1/|L|*\int_L phi ds, where L is the path given by the
% circumference of the square around the singular wire coupling points
% given by {wire.x,wire.y,wire.z}.
%
% Input:
%   msh        struct as defined by src/msh.txt
%              required fields: np,x,y,z,My,Mz
%   wire       struct as defined by src/wire.txt
%              required fields: N,N1D,x,y,z,select,ipnCpl
%              optional fields: cplCoeff
%                               (default: log(wire.r./log(wire.rCpl.val))
%   verbose    triggers console outputs and plots
%              (optional, default: [1 0])
%
% Output:
%   R13        projects from 3D solution to 1D solution at coupling points
%              (cell array of length wire.N; each matrix of size N1D-by-np)
%
% See also computeDSdWire, computeKwire, computeMwire, computePwire,
% computeQwire, computeR31
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 3, verbose = [1,0]; end
if ~isfield(wire,'cplCoeff'), wire.cplCoeff = log(wire.r)./(log(wire.rCpl.val)); end

% broadcast wire.cplCoeff for all 1D points
if size(wire.cplCoeff,2) == 1, wire.cplCoeff = repmat(wire.cplCoeff,1,wire.N1D); end

% find number of coupling points from wire.ipnCpl
Ncpl = size(wire.ipnCpl,3);
assert(Ncpl == 4, 'in the current version, Ncpl must be 4!');

if verbose(1), fprintf('calculation of R13 ...\n'); end

R13 = cell(wire.N,1);
for i = 1:wire.N
    if ~wire.select(i), continue; end
    R13this = sparse(wire.N1D,msh.np);
    for j = 1:wire.N1D
        ipnCpl = squeeze(wire.ipnCpl(i,j,:));

        if wire.rCpl.val(i,j) == 0
            % direct coupling of 1D with 3D point
            weights = 1;
        else
            % calculate 1D values along the path
            if ~isrow(ipnCpl), ipnCpl = ipnCpl'; end
            delta_s1D = distFIT(msh,[ipnCpl ipnCpl(1)]);
            % construct 1D dual mesh on rectangle
            if ~iscolumn(delta_s1D), delta_s1D = delta_s1D'; end
            delta_sd1D =  [(delta_s1D(1)      +delta_s1D(end))  /2; ...
                           (delta_s1D(1:end-1)+delta_s1D(2:end))/2];

            % calculate path length
            pathLength = sum(delta_s1D);

            % calculate weight for current 1D point
            weights = wire.cplCoeff(i,j)*delta_sd1D / pathLength;
        end
        % assemble R13
        R13this(j,ipnCpl) = weights;
    end
    R13{i} = R13this;
end

if verbose(1), fprintf('calculation of R13 done!\n'); end

end