function idxBox = canonical4box(msh,box,type)
% CANONICAL4BOX returns all canonical indices that are specified by box
% indexing scheme. One always has to define a single point, a line or a
% cube to be translated.
%
% Depending on type, the box should be given as a vector of indices or
% coordinates as follows:
% type == 'ijk'   :    box = [imin imax jmin jmax kmin kmax] (default)
% type == 'xyz'   :    box = [xmin xmax ymin ymax zmin zmax]
%
% Input:
%   msh     struct as defined by src/msh.txt
%           required fields: x,y,z,Mx,My,Mz
%   box     box for which indices shall be returned
%           format: box = [imin imax jmin jmax kmin kmax]
%   type    chooses whether box is given using indices or coordinates. Must
%           be eiter 'ijk' or 'xyz'
%           (optional, default: 'ijk')
%
% Output:
%   idxBox  canonical indices of all entities inside the given box
%           (including its boundary)
%
% See also coords2idx, idx2coords, idx2canonical, canonical2idx
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 3, type = 'ijk'; end

% get ijk index from box input
switch type
    case 'ijk'
        imin = box(1); imax = box(2);
        jmin = box(3); jmax = box(4);
        kmin = box(5); kmax = box(6);
    case 'xyz'
        [~,imin] = ismembertol(box(1),msh.x);
        [~,imax] = ismembertol(box(2),msh.x);
        [~,jmin] = ismembertol(box(3),msh.y);
        [~,jmax] = ismembertol(box(4),msh.y);
        [~,kmin] = ismembertol(box(5),msh.z);
        [~,kmax] = ismembertol(box(6),msh.z);
    otherwise
        error('type must be either ''ijk'' or ''xyz''.');
end

% collect all canonical indices inside the specified box
i = unique(imin:imax);
j = unique(jmin:jmax);
k = unique(kmin:kmax);
idxBox = NaN*ones(1,length(i)*length(j)*length(k));
count = 1;
for kz = k
    for jy = j
        for ix = i
            idxBox(count) = 1+(ix-1)*msh.Mx+(jy-1)*msh.My+(kz-1)*msh.Mz;
            count = count + 1;
        end
    end
end
idxBox = sort(idxBox);

end