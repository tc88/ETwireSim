function resistor2D(refinement,verbose)
% RESISTOR2D shows the convergence of the FIT discretization error for a 2D
% rectangular resistor coupled to a lumped resistor at a 0D point.
%
% Note that the 2D-0D coupling is implemented on a 3D geometry, since FIT
% always requires three-dimensional geometries. However, the geometry and
% results are invariant in z-direction.
%
% The brick-shaped resistor is series-connected to an external resistor and
% a voltage source. The applied current is given by "I0prime". The outer
% boundary of an equivalent circular capacitor (of equal cross-sectional
% area) serves as reference electrode. At the centre point, a cylindrical
% wire with radius wire.r (much smaller than the domain size) is positioned.
% 
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

tstart = tic;
fprintf('running 2D resistor test case ...\n');

% name of the model that is also used for output files
modelname = 'resistor2D';

% load model and settings
data = load([modelname,lower(refinement),'.mat']);
Omega = data.Omega;                                                        % [m]    : definition of a subregion Omega in which 3D norm is evaluated
msh4refine = data.msh4refine;                                              % []     : mesh for different levels of refinement
materials4refine = data.materials4refine;                                  % []     : materials for different levels of refinement
wire4refine = data.wire4refine;                                            % []     : wire for different levels of refinement
Nz = data.Nz;                                                              % []     : number of points in z-direction (is required since FIT must be used in 3D)
r0 = data.r0;                                                              % [m]    : radius of an equivalent cylindrical resistor (with same cross-sectional area)
sigma = data.sigma;                                                        % [S/m]  : conductivity of resistor material
V0 = data.V0;                                                              % [V]    : externally applied voltage
R0prime = data.R0prime;                                                    % [Ohm*m]: resistance of the external resistor per unit length in z-direction
method = data.method;                                                      % []     : description of the different used methods
selectMethods = data.selectMethods;                                        % []     : numbers of selected methods as a vector
Nsamples2D = size(msh4refine,1);
maxMethods = max(selectMethods);

% getting wire radius
rWire = wire4refine{1,selectMethods(1)}.r;                                 % [m]   : wire radius
% check that wire radius is the same for all considered cases
for i = 1:Nsamples2D
    for m = selectMethods
        assert(wire4refine{i,m}.r == rWire);
    end
end

% calculate analytical references
u2Dana = @(I0prime,r) -I0prime/(2*pi*sigma)*log(r/r0);                     % [V]    : function for calculating the electric scalar potential
RintPrime = log(r0/rWire)/(2*pi*sigma);                                    % [Ohm*m]: resistance of the equivalent circular resistor per unit length in z-direction
I0prime = V0/(R0prime+RintPrime);                                          % [A/m]  : line current per unit length in z-direction
u0Dana = repmat(u2Dana(I0prime,rWire),Nz,1);                               % [V]    : analytical solution for the wire potential

% looping over different mesh sizes and methods
u0D = zeros(Nsamples2D,maxMethods,Nz);
epsL22D = zeros(Nsamples2D,maxMethods);
epsL20D = zeros(Nsamples2D,maxMethods);
h = zeros(Nsamples2D,maxMethods);
for i = 1:Nsamples2D
    for m = selectMethods
        msh = msh4refine{i,m};
        materials = materials4refine{i,m};
        wire = wire4refine{i,m};

        % extract from mesh
        np = msh.np;
        nz = msh.nz;
        DSvec = full(diag(msh.DS));
        DSvecXY = DSvec(1:2*np);
        DSvecXY(DSvecXY==0) = [];

        % compute required wire fields
        wire.cplCoeff = log(wire.r/r0)./log(wire.rCpl.val/r0);
        wire.select = 1;
        wire.R13 = computeR13(msh,wire,verbose);
        wire.L = diff(msh.z);
        if nz == 2
            wire.Gel = spdiags([1/R0prime/2;1/R0prime/2],0,nz,nz);
        else
            error('current implementation only works for nz == 2');
        end
        wire.R31 = computeR31(msh.np,wire,verbose);
        if verbose(1), fprintf('nr. 3D degrees of freedom: %d, hmax/hmin=%f\n',msh.np,max(DSvecXY)/min(DSvecXY)); end

        % index handling
        idx.elect.dir = [msh.ipnXmin,msh.ipnXmax,msh.ipnYmin,msh.ipnYmax];
        idx.elect.dir = sort(idx.elect.dir);
        idx.elect.dof = setdiff(1:np,idx.elect.dir);

        % assign boundary conditions
        pots = NaN*ones(np,1);
        for j = 1:length(idx.elect.dir)
            [~,~,kk] = canonical2idx(msh,idx.elect.dir(j));
            pots(idx.elect.dir(j)) = u2Dana(I0prime,distFIT(msh,idx.elect.dir(j),wire.idx(kk)'));
        end

        % compute distances to 1D subdomain
        dists = dist2line(msh,uniquetol(wire.x),uniquetol(wire.y));

        % solve system
        Knum = msh.Sd*materials.Msigma*msh.Sd'+wire.R31{1}*wire.Gel*wire.R13{1}; % Note: Ps reduces to a scalar = 1 for 2D-0D coupling
        qfit = wire.R31{1}*diag(wire.Gel)*V0-Knum(:,idx.elect.dir)*pots(idx.elect.dir);
        uAnaNum = u2Dana(I0prime,dists);
        Kbcs = Knum(idx.elect.dof,idx.elect.dof);
        qbcs = qfit(idx.elect.dof,1);
        ubcs = Kbcs\qbcs;

        % post-processing
        u2D = zeros(np,1);
        u2D(idx.elect.dof,1) = ubcs;
        u2D(idx.elect.dir,1) = pots(idx.elect.dir,1);
        u0D(i,m,:) = wire.R13{1}*u2D;

        % compute average and max of hXY
        if min(DSvecXY) < wire.r/2
            error('more than one element inside wire.r detected');
        end
        h(i,m) = mean(DSvecXY);

        % compute 2D error
        epsL22D(i,m) = norm4box(msh,u2D-uAnaNum,Omega)/norm4box(msh,uAnaNum,Omega);

        % compute 0D error
        DSd = computeDSdWire(wire);
        wire.DSd = DSd{1};
        epsL20D(i,m) = normL2(squeeze(u0D(i,m,:))-u0Dana,wire.DSd)/normL2(u0Dana,wire.DSd);
    end
end

figure(9); clf;

% plot 0D rel. error vs. 1/h
subplot(2,1,1);
loglog(1./h(:,selectMethods),epsL20D(:,selectMethods),'-x');
xlabel('1/h in 1/m');
ylabel('rel. error \epsilon_{L2}^{1D}');
legend({method(selectMethods).legend});

% plot 2D rel. error vs. 1/h
subplot(2,1,2);
loglog(1./h(:,selectMethods),epsL22D(:,selectMethods),'-x');
xlabel('1/h in 1/m');
ylabel('rel. error \epsilon_{L2}^{3D}');
legend({method(selectMethods).legend});
print([modelname,'results',refinement,'.pdf'],'-dpdf');

fprintf('finished 2D resistor test case after %d seconds.\n',toc(tstart));