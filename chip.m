function chip(refinement,verbose)
% CHIP runs the transient electrothermal simulation of a microelectronic
% chip package.
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

tstart = tic;
fprintf('running chip test case ...\n');

% name of the model that is also used for output files
modelname = 'chip';

% load model and settings
data = load([modelname,refinement,'.mat']);
msh = data.msh;                                                            % []  : struct as defined by src/msh.txt
materials = data.materials;                                                % []  : struct as defined by src/materials.txt
T = data.T;                                                                % []  : temperature data
idx = data.idx;                                                            % []  : struct as defined by src/idx.txt
wire = data.wire;                                                          % []  : struct as defined by src/wire.idx
phiDir = data.phiDir;                                                      % [V] : Dirichlet potential data (np-by-1)
time = data.time;                                                          % [s] : time data (1-by-nt)
nt = length(time);

% obtain 1D-3D coupling matrices
wire.cplCoeff = log(wire.r)./(log(wire.rCpl.val));
wire.R13 = computeR13(msh,wire,verbose);
wire.R31 = computeR31(msh.np,wire,verbose);

% obtain edge lengths
dsVec = full(diag(msh.DS));
dsVec(dsVec==0) = [];

% print problem information
if verbose(1), fprintf('nr. 3D degrees of freedom: %d, hmax/hmin=%f\n',numel(idx.elect.dof),max(dsVec)/min(dsVec)); end

% solve problem
[phi3D,T3D] = solveCoupledET(msh,materials,idx,phiDir,T,time,wire,verbose);
phi1D = zeros(wire.N,wire.N1D,nt);
T1D = zeros(wire.N,wire.N1D,nt);
for i = 1:wire.N
    if ~wire.select(i), continue; end
    for k = 1:nt
        phi1D(i,:,k) = (wire.R13{i}*phi3D(:,k))';
        T1D(i,:,k) = (wire.R13{i}*T3D(:,k))';
    end
end

% visualize result
figure(1617); clf;
subplot(2,1,1);
plot(wire.sParam,1e3*phi1D(find(wire.select),:,end),'x-');
xlabel('Wire parametrization s');
ylabel('1D solution $$\overline{\varphi}_{h}(t_{0})$$ in mV','Interpreter','Latex');
subplot(2,1,2);
plot(wire.sParam,T1D(find(wire.select),:,end),'x-');
xlabel('Wire parametrization s');
ylabel('1D solution $$\overline{\mathbf{T}}(t_{0})$$ in K','Interpreter','Latex');
print([modelname,'Results',refinement,'.pdf'],'-dpdf');

% export results to Paraview
fit_write_vtk(msh.x,msh.y,msh.z,sprintf('%sResults.vtr',modelname),{'potential',phi3D(:,end);'Temperature',T3D(:,end)});

fprintf('finished chip test case after %d seconds.\n',toc(tstart));