function res = norm4box(msh,u,box)
% NORM4BOX approximates the L2-norm on a subregion Omega of the entire
% domain D.
%
% Input:
%   msh     struct as defined by src/msh.txt
%           required fields: x,y,z,nx,ny,nz,Mx,My,Mz,DVd
%   u       solution vector (Np-by-1)
%   box     box on which the L2-norm shall be computed
%           format: box = [xmin xmax ymin ymax zmin zmax]
%
% Output:
%   res     computed L2 norm
%
% See also normL2
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% extract from input
[~,imin] = ismembertol(box(1),msh.x);
[~,imax] = ismembertol(box(2),msh.x);
[~,jmin] = ismembertol(box(3),msh.y);
[~,jmax] = ismembertol(box(4),msh.y);
[~,kmin] = ismembertol(box(5),msh.z);
[~,kmax] = ismembertol(box(6),msh.z);
DVd = msh.DVd;

% check whether coordinates are contained in mesh
if imin == 0 || imax == 0 || jmin == 0 || jmax == 0 || kmin == 0 || kmax == 0
    error('all coordinates given by box must be contained in the given mesh');
end

% divide dual volumes by 2/4/8 at boundary facet/edge/corner of box.
% no need to do this at boundary, since definition of dual mesh takes care
ipnBoxBndXmin = canonical4box(msh,[imin,imin,jmin,jmax,kmin,kmax]);
ipnBoxBndXmax = canonical4box(msh,[imax,imax,jmin,jmax,kmin,kmax]);
ipnBoxBndYmin = canonical4box(msh,[imin,imax,jmin,jmin,kmin,kmax]);
ipnBoxBndYmax = canonical4box(msh,[imin,imax,jmax,jmax,kmin,kmax]);
ipnBoxBndZmin = canonical4box(msh,[imin,imax,jmin,jmax,kmin,kmin]);
ipnBoxBndZmax = canonical4box(msh,[imin,imax,jmin,jmax,kmax,kmax]);
if imin ~= 1,      DVd(ipnBoxBndXmin,ipnBoxBndXmin) = 0.5*DVd(ipnBoxBndXmin,ipnBoxBndXmin); end
if imax ~= msh.nx, DVd(ipnBoxBndXmax,ipnBoxBndXmax) = 0.5*DVd(ipnBoxBndXmax,ipnBoxBndXmax); end
if jmin ~= 1,      DVd(ipnBoxBndYmin,ipnBoxBndYmin) = 0.5*DVd(ipnBoxBndYmin,ipnBoxBndYmin); end
if jmax ~= msh.ny, DVd(ipnBoxBndYmax,ipnBoxBndYmax) = 0.5*DVd(ipnBoxBndYmax,ipnBoxBndYmax); end
if kmin ~= 1,      DVd(ipnBoxBndZmin,ipnBoxBndZmin) = 0.5*DVd(ipnBoxBndZmin,ipnBoxBndZmin); end
if kmax ~= msh.nz, DVd(ipnBoxBndZmax,ipnBoxBndZmax) = 0.5*DVd(ipnBoxBndZmax,ipnBoxBndZmax); end

% calculate norm
ipnBox2eval = canonical4box(msh,[imin,imax,jmin,jmax,kmin,kmax]);
DVd2eval = DVd(ipnBox2eval,ipnBox2eval);
res = normL2(u(ipnBox2eval),DVd2eval);

end