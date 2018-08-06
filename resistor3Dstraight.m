function resistor3Dstraight(refinement,verbose)
% RESISTOR3DSTRAIGHT runs 1D-3D coupled thin wire simulations and compares
% result to analytical solution of a line current. Boundary conditions are
% set according to the analytical solution.
% A variable number of 3D points along the wire are used as 1D-3D coupling
% point. Thus, a study on the influence of the number of 1D points is
% possible.
% 
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

tstart = tic;
fprintf('running straight 3D resistor test case ...\n');

% name of the model that is also used for output files
modelname = 'resistor3Dstraight';

% load model and settings
data = load([modelname,refinement,'.mat']);
d = data.d;                                                                % [m]   : size of unit cube
r0 = data.r0;                                                              % [m]    : radius of an equivalent cylindrical resistor (with same cross-sectional area)
I0prime = data.I0prime;                                                    % [A/m] : line current source for manufactured solution
sigma = data.sigma;                                                        % [S/m] : conductivity of domain
method = data.method;                                                      % []    : description of the different used methods
selectMethods = data.selectMethods;                                        % []    : numbers of selected methods as a vector
NxyVecGlobal = data.NxyVecGlobal;                                          % []    : no. of points in x- and y-direction used for global grading
NxyVecLocal = data.NxyVecLocal;                                            % []    : no. of points in x- and y-direction used for local grading
Omega = data.Omega;                                                        % [m]   : definition of a subregion Omega in which 3D norm is evaluated
msh4refine = data.msh4refine;                                              % []    : mesh for different levels of refinement
materials4refine = data.materials4refine;                                  % []    : materials for different levels of refinement
wire4refine = data.wire4refine;                                            % []    : wire for different levels of refinement
Nsamples3D = size(msh4refine,1);
Nsamples1D = size(msh4refine,2);
maxMethods = max(selectMethods);

% getting wire radius
rWire = wire4refine{1,1,1}.r;                                              % [m]   : wire radius
% check that wire radius is the same for all considered cases
for i = 1:Nsamples3D
    for j = 1:Nsamples1D
        for m = selectMethods
            assert(wire4refine{i,j,m}.r == rWire);
        end
    end
end

% calculate analytical references
uAna3D = @(r,z)-I0prime*z/(2*pi*d*sigma).*log(r/r0);                       % [V]   : analytical solution in 3D
uAna1D = @(z)uAna3D(rWire,z);                                              % [V]   : analytical solution in 1D
u1DnormAnaL2exact = 1/sqrt(3)*uAna1D(1);                                   % [V]   : analytical 1D norm in L2
u1DnormAnaH1exact = uAna1D(1);                                             % [V]   : analytical 1D norm in H1
u3DnormAnaL2exact = I0prime/(d*sigma)*0.03709449913929426;                 % [V]   : integral solved by Mathematica, see resistor3DstraightAnaSol.nb

% initializations
h3D = zeros(Nsamples3D,Nsamples1D,maxMethods);
h1D = zeros(1,Nsamples1D);
u1DnormAnaL2grid = zeros(Nsamples3D,Nsamples1D,maxMethods);
u1DnormAnaH1grid = zeros(Nsamples3D,Nsamples1D,maxMethods);
u3DnormAnaL2grid = zeros(Nsamples3D,Nsamples1D,maxMethods);
u1DnormL2 = zeros(Nsamples3D,Nsamples1D,maxMethods);
u1DnormH1 = zeros(Nsamples3D,Nsamples1D,maxMethods);
u3DnormL2 = zeros(Nsamples3D,Nsamples1D,maxMethods);
u1DdiffNormL2 = zeros(Nsamples3D,Nsamples1D,maxMethods);
u1DdiffNormH1 = zeros(Nsamples3D,Nsamples1D,maxMethods);
u3DdiffNormL2 = zeros(Nsamples3D,Nsamples1D,maxMethods);
u1DnormRefL2 = zeros(1,1,maxMethods);
u1DnormRefH1 = zeros(1,1,maxMethods);
u3DnormRefL2 = zeros(1,1,maxMethods);

% start loop
tic;
for m = selectMethods
    for i = 1:Nsamples3D
		% choose Nxy depending on grading
		if strcmp(method(m).grading,'global')
			if i > length(NxyVecGlobal), continue; end
		elseif strcmp(method(m).grading,'local')
			if i > length(NxyVecLocal), continue; end
		else
			error('grading must be ''local'' or ''global''');
		end
        for j = 1:Nsamples1D
            msh = msh4refine{i,j,m};
            materials = materials4refine{i,j,m};
            wire = wire4refine{i,j,m};
            wire.idx = coords2idx(msh,[wire.x',wire.y',wire.z'])';
            for k = 1:4
                wire.ipnCpl(:,:,k) = coords2idx(msh,[wire.xCpl(:,:,k)',wire.yCpl(:,:,k)',wire.zCpl(:,:,k)'])';
            end

            if verbose(1), fprintf('solving problem using m=%d, Nx=Ny=%d, Nz=%d and N1D=%d\n',m,msh.ny,msh.nz,wire.N1D); end

            % extract from msh struct
            np = msh.np;

            % get edge lengths
            DSvec = full(diag(msh.DS));
            DSvec(DSvec==0) = [];

            % fill wire fields
            wire.R13 = computeR13(msh,wire,verbose);
            wire.L = diff(wire.z);
            wire.A = pi*wire.r^2;
            wire.sigma = wire.A*1e15*sigma;
            wire.Gel = 1./wire.L*wire.sigma;
            wire.Msigma1D = spdiags([wire.Gel';0],0,wire.N1D,wire.N1D);
            wire.R31 = computeR31(msh.np,wire,verbose);
            DSd = computeDSdWire(wire);
            wire.DSd = DSd{1};
            assert(all(wire.L==wire.L(1)));
            if verbose(1), fprintf('nr. 3D degrees of freedom: %d, hmax/hmin=%f\n',msh.np,max(DSvec)/min(DSvec)); end

            % assign boundary conditions
            pots = NaN*ones(np,1);
            dists = dist2line(msh,uniquetol(wire.x),uniquetol(wire.y));
            [~,~,zBnd] = idx2coords(msh,msh.ipnBnd);
            pots(msh.ipnBnd) = uAna3D(dists(msh.ipnBnd),zBnd');
            pots(wire.idx(1)) = uAna1D(min(msh.z));
            pots(wire.idx(end)) = uAna1D(max(msh.z));

            % index handling
            idx.elect.dir = find(~isnan(pots));
            idx.elect.dof = setdiff(1:np,idx.elect.dir);

            % solve system
            u3D = NaN*ones(np,1);
            u3D(idx.elect.dir) = pots(idx.elect.dir);
            Ps = createPx(wire.N1D);
            Knum = msh.Sd*materials.Msigma*msh.Sd'+wire.R31{1}*Ps'*wire.Msigma1D*Ps*wire.R13{1};
            Kbcs = Knum(idx.elect.dof,idx.elect.dof);
            qbcs = -Knum(idx.elect.dof,idx.elect.dir)*pots(idx.elect.dir);
            u3D(idx.elect.dof) = Kbcs\qbcs;
            u1D = wire.R13{1}*u3D;

            % calculate h3D and h1D
            h3D(i,j,m) = mean(DSvec);
            h1D(1,j) = wire.L(1);
            assert(all(wire.L == wire.L(1)),'only equidistant 1D discretizations are supported');

            % calculate L2 norm of absolute error of u1D
            u1Dana = uAna1D(wire.z');
            u1Ddiff = u1D-u1Dana;
            u1DdiffNormL2(i,j,m) = normL2(u1Ddiff,wire.DSd);
            u1DnormL2(i,j,m) = normL2(u1D,wire.DSd);
            u1DnormAnaL2grid(i,j,m) = normL2(u1Dana,wire.DSd);

            % calculate L2 norm of absolute error of derivative of u1D
            wire.DSinv = wire.Msigma1D/wire.sigma;
            u1DdiffNormH1(i,j,m) = normL2(Ps*u1Ddiff,wire.DSinv);
            u1DnormH1(i,j,m) = normL2(Ps*u1D,wire.DSinv);
            u1DnormAnaH1grid(i,j,m) = normL2(Ps*u1Dana,wire.DSinv);

            % calculate L2 norm of absolute error of u3D
            [~,~,zVec] = idx2coords(msh,(1:np)');
            u3Dana = uAna3D(dists,zVec);
            u3Ddiff = u3D-u3Dana;
            u3DdiffNormL2(i,j,m) = norm4box(msh,u3Ddiff,Omega);
            u3DnormL2(i,j,m) = norm4box(msh,u3D,Omega);
            u3DnormAnaL2grid(i,j,m) = norm4box(msh,u3Dana,Omega);
            % save fine grid solutions as reference
            if strcmp(method(m).grading,'global') && i==length(NxyVecGlobal) && j==1 ...
		    || strcmp(method(m).grading,'local' ) && i==length(NxyVecLocal ) && j==1
                u1DnormRefL2(1,1,m) = normL2(u1D,wire.DSd);
                u1DnormRefH1(1,1,m) = normL2(Ps*u1D,wire.DSinv);
                u3DnormRefL2(1,1,m) = norm4box(msh,u3D,Omega);
            end

            % print 1D solution
            if verbose(1), fprintf('phi1D = [%s]\n',sprintf('%d ',u1D')); end

            % export solution to .vtr file for Paraview
            if ismember(m,[2,6]) && i==Nsamples3D && j==1
                uParaview = u3D;
                uParaview(wire.idx) = u1D;
                fit_write_vtk(msh.x,msh.y,msh.z,sprintf('%sResultsMethod%d.vtr',modelname,m),{'phi',uParaview});
            end
        end
    end
end

% calculate errors
epsL21D = u1DdiffNormL2./u1DnormAnaL2grid;
epsH11D = u1DdiffNormH1./u1DnormAnaH1grid;
epsL23D = u3DdiffNormL2./u3DnormAnaL2grid;
deltaL21D = abs(u1DnormL2-u1DnormAnaL2exact)./u1DnormAnaL2exact;
deltaH11D = abs(u1DnormH1-u1DnormAnaH1exact)./u1DnormAnaH1exact;
deltaL23D = abs(u3DnormL2-u3DnormAnaL2exact)./u3DnormAnaL2exact;
DeltaL21D = abs(u1DnormL2-u1DnormRefL2)./u1DnormRefL2;
DeltaH11D = abs(u1DnormH1-u1DnormRefH1)./u1DnormRefH1;
DeltaL23D = abs(u3DnormL2-u3DnormRefL2)./u3DnormRefL2;

% select refinement to plot
i1D2plot = Nsamples3D;
j3D2plot = 1;

figure(1213); clf; hold off;

% plot 3D rel. error vs. 1/h3D for different meshes
subplot(3,1,1);
x2plot = 1./squeeze(h3D(:,j3D2plot,[1,2,3,4]));
y2plot = squeeze(epsL23D(:,j3D2plot,[1,2,3,4]));
loglog(x2plot,y2plot,'-x');
xlabel('1/h in 1/m');
ylabel('rel. error \epsilon_{L2}^{3D}');
legend(method([1,2,3,4]).legend,'Location','bestoutside');

% plot 3D rel. error vs. 1/h3D for different rCpl
subplot(3,1,2);
x2plot = 1./squeeze(h3D(:,j3D2plot,[2,5,6]));
y2plot = squeeze(epsL23D(:,j3D2plot,[2,5,6]));
loglog(x2plot,y2plot,'-x');
xlabel('1/h in 1/m');
ylabel('rel. error \epsilon_{L2}^{3D}');
legend(method([2,5,6]).legend,'Location','bestoutside');

% plot 1D rel. error vs. 1/h1D
subplot(3,1,3);
x2plot = 1./h1D;
y2plot = squeeze(deltaL21D(i1D2plot,:,selectMethods));
loglog(x2plot,y2plot,'-x');
xlabel('1/h in m');
ylabel('rel. error \delta_{L2}^{1D}');
legend(method(selectMethods).legend,'Location','bestoutside');

print([modelname,'Results',refinement,'.pdf'],'-dpdf');

fprintf('finished straight 3D resistor test case after %d seconds.\n',toc(tstart));