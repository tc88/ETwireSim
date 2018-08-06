function resistor3Dbent(refinement,verbose)
% RESISTOR3DBENT runs 1D-3D coupled simulations for a bent thin wire and
% compares result to the finest solution. Boundary conditions are set as
% homogeneous Dirichlet except for the wire end points.
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
fprintf('running bent 3D resistor test case ...\n');

% name of the model that is also used for output files
modelname = 'resistor3Dbent';

% load model and settings
data = load([modelname,refinement,'.mat']);
d = data.d;                                                                % [m]   : side length of domain (unit cube)
dPEC = data.dPEC;                                                          % [m]   : side length of PEC cubes at wire ends
Omega = data.Omega;                                                        % [m]   : definition of a subregion Omega in which 3D norm is evaluated
msh4geo = data.msh4geo;                                                    % []    : mesh for different levels of refinement without mesh lines at coupling points
msh4refine = data.msh4refine;                                              % []    : mesh for different levels of refinement with mesh lines at coupling points
materials4refine = data.materials4refine;                                  % []    : materials for different levels of refinement
wire4method = data.wire4refine;                                            % []    : wire for different levels of refinement
clear data;

Nsamples3D = size(msh4refine,1);
Nsamples1D = size(msh4refine,2);

% initializations
h3D = zeros(Nsamples3D,Nsamples1D);
h1D = zeros(Nsamples3D,1);
u1DnormL2 = zeros(Nsamples3D,Nsamples1D);
u1DnormH1 = zeros(Nsamples3D,Nsamples1D);
u3DnormL2 = zeros(Nsamples3D,Nsamples1D);
u1DnormRefL2 = zeros(1,1);
u1DnormRefH1 = zeros(1,1);
u3DnormRefL2 = zeros(1,1);

% start solving
for i = 1:Nsamples3D
    for j = 1:Nsamples1D
        % load msh and wire
        mshInit = msh4geo{i,j};
        wire = wire4method{i,j};
        materials = materials4refine{i,j};
        msh = msh4refine{i,j};

        if verbose(1)
            fprintf('solving problem using i=%d, Nz=%d and N1D=%d\n',i,wire.N1D+2,wire.N1D);
            fprintf('%f seconds have passed\n',toc);
            fprintf('nx4geo=%d   , ny4geo=%d   , nz4geo=%d\n'   ,mshInit.nx,mshInit.ny,mshInit.nz);
            fprintf('nx4method=%d, ny4method=%d, nz4method=%d\n',msh.nx    ,msh.ny    ,msh.nz);
        end

        % compute R13
        wire.cplCoeff = log(wire.r)./(log(wire.rCpl.val));
        wire.R13 = computeR13(msh,wire,verbose);

        np = msh.np;
        DSvec = full(diag(msh.DS));
        DSvec(DSvec==0) = [];
        if verbose(1), fprintf('nr. 3D degrees of freedom: %d, hmax/hmin=%f\n',msh.np,max(DSvec)/min(DSvec)); end

        % fill wire fields
        wire.idx = coords2idx(msh,[wire.x',wire.y',wire.z'])';
        wire.sigma = wire.A*1e15*unique(materials.sigma);
        wire.Gel = 1./wire.L*wire.sigma;
        wire.Msigma1D = spdiags([wire.Gel';0],0,wire.N1D,wire.N1D);
        wire.R31 = computeR31(msh.np,wire,verbose);
        DSd = computeDSdWire(wire);
        wire.DSd = DSd{1};

        % ensure that all inner wire nodes are far enough from PEC box
        assert(distFIT(msh,wire.idx(1:2))>dPEC);
        assert(distFIT(msh,wire.idx(end-1:end))>dPEC);

        % assign boundary conditions
        pots = NaN*ones(np,1);
        ipnPECstart = canonical4box(msh,[wire.x(1)-dPEC/2 wire.x(1)+dPEC/2 ...
                                         wire.y(1)-dPEC/2 wire.y(1)+dPEC/2 ...
                                         wire.z(1)-dPEC/2 wire.z(1)+dPEC/2] ...
                                         ,'xyz');
        ipnPECend = canonical4box(msh,[wire.x(end)-dPEC/2 wire.x(end)+dPEC/2 ...
                                       wire.y(end)-dPEC/2 wire.y(end)+dPEC/2 ...
                                       wire.z(end)-dPEC/2 wire.z(end)+dPEC/2] ...
                                       ,'xyz');
        pots(ipnPECstart) = 0;
        pots(ipnPECend) = 1;

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
        h3D(i,j) = mean(DSvec);
        h1D(i) = (wire.sParam(2)-wire.sParam(1))/d;

        % calculate L2 norm of 1D solution
        u1DnormL2(i,j) = normL2(u1D,wire.DSd);

        % calculate L2 norm of derivative of 1D solution
        wire.DSinv = wire.Msigma1D/wire.sigma;
        u1DnormH1(i,j) = normL2(Ps*u1D,wire.DSinv);

        % calculate L2 norm of absolute error of u3D
        u3DnormL2(i,j) = norm4box(msh,u3D,Omega);

        % print 1D solution
        if verbose(1), fprintf('phi1D = [%s]\n',sprintf('%d ',u1D')); end

        % save fine grid solutions as reference
        if i==Nsamples3D && j==1
            u1DnormRefL2(1,1) = normL2(u1D,wire.DSd);
            u1DnormRefH1(1,1) = normL2(Ps*u1D,wire.DSinv);
            u3DnormRefL2(1,1) = norm4box(msh,u3D,Omega);
            if verbose(1)
                fprintf('u1DnormRefL2 = %f\n',u1DnormRefL2);
                fprintf('u3DnormRefL2 = %f\n',u3DnormRefL2);
            end
        end

        if i==Nsamples3D
            % obtain u for Paraview visualization
            ipn4uParaview = coords2idx(msh,mshInit.x',mshInit.y',mshInit.z');
            ipnWire4mshInit = coords2idx(mshInit,wire.x,wire.y,wire.z);
            uParaview = u3D(ipn4uParaview);
            uParaview(ipnWire4mshInit) = u1D;

            % export data for Paraview
            filenameParaview = sprintf('%sResults.vtr',modelname);
            fit_write_vtk(mshInit.x,mshInit.y,mshInit.z,filenameParaview,{'phi',uParaview});

            % visualize 1D solution with respect to t
            figure(14); clf;
            plot(wire.sParam,u1D,'-x');
            xlabel('Wire parametrization s');
            ylabel('$$\overline{\varphi}_{h}$$ in V','Interpreter','Latex');
            print([modelname,'Results',refinement,'Phi.pdf'],'-dpdf');
        end
    end
end

% calculate errors
DeltaL21D = abs(u1DnormL2-u1DnormRefL2)./u1DnormRefL2;
DeltaH11D = abs(u1DnormH1-u1DnormRefH1)./u1DnormRefH1;
DeltaL23D = abs(u3DnormL2-u3DnormRefL2)./u3DnormRefL2;

figure(15); clf; hold off;

% plot 1D rel. error vs. 1/h
subplot(2,1,1);
x2plot = 1./squeeze(h3D(:,1));
y2plot = squeeze(DeltaL21D(:,1));
loglog(x2plot,y2plot,'-x');
xlabel('1/h in 1/m');
ylabel('rel. error \Delta_{L2}^{1D}');
legend('rCpl=1e-2/\kappa');

% plot 3D rel. error vs. 1/h
subplot(2,1,2);
x2plot = 1./squeeze(h3D(:,1));
y2plot = squeeze(DeltaL23D(:,1));
loglog(x2plot,y2plot,'-x');
xlabel('1/h in 1/m');
ylabel('rel. error \Delta_{L2}^{3D}');
legend('rCpl=1e-2/\kappa');
print([modelname,'Results',refinement,'Error.pdf'],'-dpdf');

fprintf('finished bent 3D resistor test case after %d seconds.\n',toc(tstart));