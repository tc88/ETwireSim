function [res1,res2,res3] = coords2idx(msh,arg1,arg2,arg3)
% COORDS2IDX translates given coordinates into the canonical indexing
% scheme.
%
%   n = COORDS2IDX(msh,x,y,z) returns the canonical index of the given
%   coordinates. x, y and z must be scalars.
%   [I,J,K] = COORDS2IDX(msh,x,y,z) returns the indices of the given
%   coordinates. x, y and z must be scalars.
%   [n] = COORDS2IDX(msh,x,y,z,true) returns all permutations of the given
%   coordinates. x, y and z must be vectors (possibly of different length).
%   [...] = COORDS2IDX(msh,coords) returns the indices of the given
%   coordinates. coords must be a matrix of the form [xmesh ymesh zmesh]
%   where {x,y,z}mesh are column vectors.
%
%   See also idx2coords, idx2canonical, canonical2idx, canonical4box
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% extract from mesh and make mesh vectors row vectors
xmesh = msh.x; ymesh = msh.y; zmesh = msh.z;
xmesh = xmesh(:)'; ymesh = ymesh(:)'; zmesh = zmesh(:)';

if nargin == 4 || (nargin == 2 && size(arg1,2) == 3)
    % the indices of [x y z] coordinates are searched for
    allPermutations = false;
    if nargin == 2
        arg2 = arg1(:,2);
        arg3 = arg1(:,3);
        arg1 = arg1(:,1);
    end
    % if arguments are vectors, return indices of all permutations
    if nargin == 4 && (~isscalar(arg1) || ~isscalar(arg2) || ~isscalar(arg3))
        allPermutations = true;
    end
    % check input format
    columnInput = false;
    if iscolumn(arg1) && iscolumn(arg2) && iscolumn(arg3)
        columnInput = true;
        arg1 = arg1'; arg2 = arg2'; arg3 = arg3';
    elseif ~(isrow(arg1) && isrow(arg2) && isrow(arg3))
        error('Input vectors should either be all row or all column vectors!')
    end

    % [I,J,K] or the canonical index for the specified points is searched for
    [~,I] = ismembertol(arg1,xmesh);
    [~,J] = ismembertol(arg2,ymesh);
    [~,K] = ismembertol(arg3,zmesh);

    if ~all(I) || ~all(J) || ~all(K)
        error('The given coordinates are not part of the used discretization given by the mesh!');
    end

    % assign output
    if nargout == 0 || nargout == 1
        res1 = idx2canonical(msh,I,J,K,allPermutations);
    elseif nargout == 3
        res1 = I;
        res2 = J;
        res3 = K;
    else
        error('Wrong number of output arguments!')
    end

    % check output format
    if columnInput
        res1 = res1';
        if nargout == 3
            res2 = res2';
            res3 = res3';
        end
    end
else
    error('Wrong number or format of input arguments');
end

end