function dist = distFIT(msh,varargin)
% DISTFIT calculates the Euclidian distance between mesh points.
%
%   distFIT(msh,ipn) calculates the Euclidian distance between primary
%   mesh nodes given by the vector ipn. The distance is calculated between
%   all N+1 neighboring elements of ipn and returned as an output vector of
%   length N.
%
%   distFIT(msh,ipnFrom,ipnTo) calculates the Euclidian distance between
%   primary mesh nodes ipnFrom and primary mesh nodes ipnTo. ipnFrom must
%   be a vector of length N and ipnTo can be a scalar or also a vector of
%   length N. In any case, the return vector is of length N.
%
% See also dist2line
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% make sure input vectors are row vectors
for i = 1:length(varargin)
    if ~isrow(varargin{i})
        varargin{i} = varargin{i}';
    end
end

if length(varargin) == 1
    % compute all distances on a 1D line
    ipn = varargin{1};
    [x,y,z] = idx2coords(msh,ipn);
    xdiff = abs(diff(x));
    ydiff = abs(diff(y));
    zdiff = abs(diff(z));
elseif length(varargin) == 2
    % compute distance between ipnFrom and ipnTo
    ipnFrom = varargin{1};
    ipnTo = varargin{2};
    if numel(ipnFrom) > numel(ipnTo) && isscalar(ipnTo)
        ipnTo = repmat(ipnTo,1,length(ipnFrom));
    elseif numel(ipnFrom) ~= numel(ipnTo)
        error('wrong format of input data');
    end
    [x1,y1,z1] = idx2coords(msh,ipnFrom);
    [x2,y2,z2] = idx2coords(msh,ipnTo);
    xdiff = abs(x1-x2);
    ydiff = abs(y1-y2);
    zdiff = abs(z1-z2);
else
    error('wrong number of input arguments, see documentation for help');
end
dist = sqrt(xdiff.^2 + ydiff.^2 + zdiff.^2);

end