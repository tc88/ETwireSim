function canonicalIdx = idx2canonical(msh,i,j,k,allPermutations)
% IDX2CANONICAL translates given indices (i,j,k) into the canonical
% indexing scheme. One always has to define a single point, a line or a
% cube to be translated.
%
% Input:
%   msh                 struct as defined by src/msh.txt
%                       required fields: Mx,My,Mz
%   i                   can be a scalar or a vector for the x-index
%   j                   can be a scalar or a vector for the y-index
%   k                   can be a scalar or a vector for the z-index
%   allPermutations     set to true if all permutations of the input
%                       vectors shall be added to the output indices
%                       (optional, default: false)
%
% Output:
%   idx     a vector of all indices specified by the given indices
%           i,j and k in the canonical indexing scheme
%
% See also coords2idx, idx2coords, canonical2idx, canonical4box
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 5, allPermutations = false; end

if ~isrow(i), i = i'; end
if ~isrow(j), j = j'; end
if ~isrow(k), k = k'; end
if isempty(i), error('parameter i does not contain values!'); end
if isempty(j), error('parameter j does not contain values!'); end
if isempty(k), error('parameter k does not contain values!'); end

if (numel(i) == numel(j) && numel(j) == numel(k)) ...
        || (numel(i) == 1 && numel(j) == 1) ...
        || (numel(j) == 1 && numel(k) == 1) ...
        || (numel(i) == 1 && numel(k) == 1)
    canonicalIdx = 1+(i-1)*msh.Mx+(j-1)*msh.My+(k-1)*msh.Mz;
    return
else
    if ~allPermutations
        warning('idx2canonical should not be used for this function. Use canonical4box instead!');
    else
        i = unique(i); j = unique(j); k = unique(k);
        [ii,jj,kk] = ndgrid(i,j,k);
        canonicalIdx = 1+(ii(:)-1*msh.Mx+(jj(:)-1)*msh.My+(kk(:)-1)*msh.Mz);
        canonicalIdx = sort(canonicalIdx)';
    end
end

end