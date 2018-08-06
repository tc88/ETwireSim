function [phiSol,Tsol] = solveCoupledET(msh,materials,idx,pots,T,time,wire,verbose)
% SOLVECOUPLEDET solves the linear electrothermal problem.
%
% Input:
%   msh         struct as defined by src/msh.txt
%               required fields: np,Sd
%   materials   struct as defined by src/materials.txt
%               required fields: Meps,Msigma,Mrhoc,Mlambda
%               optional fields: H11,H13
%   idx         struct as defined by src/idx.txt
%               required fields: elect.dir,elect.dof,therm.dir,therm.dof
%   pots        electric potential. Entries for degrees of freedom
%               (dofs) need to be NaN while all other entries are
%               interpreted as fixed potentials (Dirichlet conditions)
%               (np-by-1 for constant in time or np-by-nt)
%   T           struct for temperature data
%       .dir    vector of Dirichlet temperatures. DoFs need to be NaN while
%               all other entries are interpreted as fixed temperatures
%               (optional, np-by-1 for constant in time or np-by-nt)
%       .init   initial temperature vector (np-by-1)
%       .inf    ambient temperature (scalar)
%               (only required for convective boundary conditions)
%               (optional, default: 0)
%   time        vector of time steps
%   wire        wire struct as defined by src/wire.txt
%               required fields: N,N1D,L,Ps,R13,R31,sigma,lambda
%               optional fields: select (default: ones(wire.N,1))
%               (optional, default: [])
%   verbose     triggers console output and plots
%               (optional, default: [1 0])
%
% Output:
%   phiSol     electric potential solution (np-by-1)
%   Tsol       temperature solution (np-by-1)
%
% See also computeQwire
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 8,  verbose = [1 0];    end
if nargin < 7,  wire = [];          end

tstart = tic;
fprintf('running transient electrothermal simulation using FIT ...\n');

% initialize convective matrices
if ~isfield(materials,'H11'), materials.H11 = sparse(msh.np,msh.np); end
if ~isfield(materials,'H13'), materials.H13 = sparse(msh.np,1);      end

% if no ambient temperature is given, use zero
if ~isfield(T,'inf'), T.inf = 0; end

% if no Dirichlet temperatures are given, use dummy vector
if ~isfield(T,'dir'), T.dir = NaN*ones(msh.np,1); end

% extract variables from inputs
np = msh.np;
Sd = msh.Sd;
nt = length(time);
Meps = materials.Meps;
Mrhoc = materials.Mrhoc;
Mlambda = materials.Mlambda;
Msigma = materials.Msigma;

% extend Dirichlet vectors if a constant potential over time is given
if size(pots,2)  == 1, pots  = repmat(pots,1,nt);  end
if size(T.dir,2) == 1, T.dir = repmat(T.dir,1,nt); end

% initializations for time loop
phi.sol = zeros(np,nt);
phi.sol(idx.elect.dir,:) = pots(idx.elect.dir,:);
phi.dir = pots(idx.elect.dir,:);
Qwire = zeros(np,nt);
T.sol = zeros(np,nt);
T.sol(idx.therm.dir,:) = T.dir(idx.therm.dir,:);
T.sol(:,1) = T.init;

% thermal mass and stiffness matrix
Mth = Mrhoc;
Kth = Sd*Mlambda*Sd';

% convection matrix
H = materials.H11;
H13 = materials.H13;

% time loop
for idxTime = 2:nt
    % compute time step dt
    dt = time(idxTime) - time(idxTime-1);

    % print current simulation time
    if verbose(1), fprintf('Current simulation time: %d s\n',time(idxTime)); end

    % compute electrical mass and stiffness matrix
    Mel = Sd*Meps*Sd';
    Kel = Sd*Msigma*Sd';

    % if wires are present, add wire contribution to electric stiffness matrix
    if ~isempty(wire)
        wire.Ps = computePwire(wire);
        wire.Msigma = computeMwire(wire,wire.sigma);
        wire.DSd = computeDSdWire(wire);
        wire.Ksigma = computeKwire(wire,wire.Msigma);
        Kel = Kel + wire.Ksigma;
    end

    % split electrical system
    K11el = Kel(idx.elect.dof,idx.elect.dof);
    K12el = Kel(idx.elect.dof,idx.elect.dir);
    M11el = Mel(idx.elect.dof,idx.elect.dof);
    M12el = Mel(idx.elect.dof,idx.elect.dir);

    % set up electric system matrix and right hand side
    Ael = M11el/dt + K11el;
    rhsEl = M11el/dt*phi.sol(idx.elect.dof,idxTime-1) - (M12el/dt+K12el)*phi.dir(:,idxTime) + M12el/dt*phi.dir(:,idxTime-1);

    % solve electric system
    if verbose(1), fprintf('solving electrical system ...\n'); end
    phi.sol(idx.elect.dof,idxTime) = Ael\rhsEl;

    % if wires are present, add wire contribution to thermal stiffness matrix and compute Joule losses on the wires
    if ~isempty(wire)
        wire.Mlambda = computeMwire(wire,wire.lambda);
        wire.Klambda = computeKwire(wire,wire.Mlambda);
        Kth = Kth + wire.Klambda;
        Qwire(:,idxTime) = computeQwire(phi.sol(:,idxTime),wire);
    end

    % set up thermal system matrix and right hand side
    Ath = 1/dt*Mth+Kth+H;
    rhsTh = 1/dt*Mth*T.sol(:,idxTime-1) - Ath(:,idx.therm.dir)*T.dir(idx.therm.dir,idxTime) - H13*T.inf + Qwire(:,idxTime);

    % solve thermal system for current step of fix point iteration
    if verbose(1), fprintf('solving thermal system ...\n'); end
    T.sol(idx.therm.dof,idxTime) = Ath(idx.therm.dof,idx.therm.dof)\rhsTh(idx.therm.dof);
end

% return values
phiSol = phi.sol;
Tsol = T.sol;

fprintf('finished transient electrothermal simulation using FIT after %d seconds.\n',toc(tstart));

end