%
% Global run using transport matrices
%
% Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
% and be put into the location 'NUMmodel/TMs'
% Simulations currently works with:
%  - MITgcm_2.8deg (low resolution; runs on a laptop)
%  - MITgcm_ECCO (higher resolution; requires more memory).
%
% Input:
%  p: parameter structure from parametersGlobal
%  sim: (optional) simulation to use for initial conditions
%  bCalcAnnualAverages: increases the simulation time by a factor 2-3
%
% Output:
%  sim: structure with simulation results
%
function sim = simulateGlobal(p, sim, bCalcAnnualAverages)

arguments
    p struct;
    sim struct = [];
    bCalcAnnualAverages = false; % Whether to calculate annual averages
end
%
% Get the global parameters if they are not already set:
%
if ~isfield(p,'nameModel')
    p = parametersGlobal(p);
end

ixN = p.idxN;
ixDOC = p.idxDOC;

bSilicate = false;
if isfield(p,'idxSi')
    ixSi = p.idxSi;
    bSilicate = true;
end
ixB = p.idxB:p.n;

%Tbc = [];

disp('Preparing simulation')
%
% Check that files exist:
%
path = fileparts(mfilename('fullpath'));
addpath(strcat(path,'/Transport matrix'));

if ~exist(p.pathBoxes,'file')
    error( sprintf('Error: Cannot find transport matrix file: %s',...
        p.pathBoxes));
end
% ---------------------------------
% Load library:
% ---------------------------------
% if p.bParallel
%     if isempty(gcp('nocreate'))
%         parpool('AttachedFiles',...
%             {'../Fortran/NUMmodel_matlab.so',...
%              '../Fortran/NUMmodel_wrap_colmajor4matlab.h'});
%     end
%     %
%     % Set parameters:
%     %
%     h = gcp('nocreate');
%     poolsize = h.NumWorkers;
%     parfor i=1:poolsize
%         loadNUMmodelLibrary();
%         %calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(p.mAdult)), p.mAdult);
%         calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly',int32(10));
%     end
% else
%     loadNUMmodelLibrary();
%     %calllib(loadNUMmodelLibrary(), 'f_setupgeneric', int32(length(p.mAdult)), p.mAdult);
%     calllib(loadNUMmodelLibrary(), 'f_setupgeneralistsonly',int32(10));
% end
% ---------------------------------------
% Initialize run:
% ---------------------------------------
simtime = p.tEnd/p.dtTransport; %simulation time in half days
% Load grid data:
%load(p.pathGrid);
%load(p.pathConfigData);
disp('load boxes')
load(p.pathBoxes, 'nb', 'Ybox', 'Zbox');
load(p.pathGrid,'dznom','nz')
% sim.dznom=dznom;

% Preparing timestepping
Ix = speye(nb,nb);
month = 0;
mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];

%
% Initial conditions:
%
disp('initial conditions')
if (nargin==2)
    disp('Starting from previous simulation.');
    u(:,ixN) = gridToMatrix(squeeze(double(sim.N(:,:,:,end))),[],p.pathBoxes, p.pathGrid);
    u(:, ixDOC) = gridToMatrix(squeeze(double(sim.DOC(:,:,:,end))),[],p.pathBoxes, p.pathGrid);
    if bSilicate
        u(:, ixSi) = gridToMatrix(squeeze(double(sim.Si(:,:,:,end))),[],p.pathBoxes, p.pathGrid);
    end
    for i = 1:p.n -p.idxB+1
        u(:, ixB(i)) = gridToMatrix(squeeze(double(squeeze(sim.B(:,:,:,i,end)))),[],p.pathBoxes, p.pathGrid);
    end
else
    if exist(strcat(p.pathN0),'file')
        load(p.pathN0, 'N_WOAtoMITgcmEcco');
        
        u(:, ixN) = gridToMatrix(squeeze(N_WOAtoMITgcmEcco(:,:,:,1)), [], p.pathBoxes, p.pathGrid);
    else
        u(:, ixN) = 150*ones(nb,1);
    end
    u(:, ixDOC) = zeros(nb,1) + p.u0(ixDOC);
    if bSilicate
        u(:, ixSi) = zeros(nb,1) + p.u0(ixSi);
    end
    u(:, ixB) = ones(nb,1)*p.u0(ixB);
end
u(:,ixDOC)=60*12.011;
Nwoa=gridToMatrix(squeeze(N_WOAtoMITgcmEcco(:,:,:,1)), [], p.pathBoxes, p.pathGrid);
% ini=load(fullfile('..','..','results','20pct','global','mean2022_12_05_0','nrRanditer_1_2022_12_05_10_43_rand0.99469_nr1.mat'));
% ini=load(fullfile('..','..','results','20pct','global','mean2022_12_06','nrRanditer_1_2022_12_06_15_03_rand0.24927_nr1.mat'));
% for i = 1:p.n -p.idxB+1
%      u(:, ixB(i)) = gridToMatrix(squeeze(double(squeeze(ini.sim.B(:,:,:,i,end)))),[],p.pathBoxes, p.pathGrid);
% end
%
% Load temperature:
%
disp('load temp')
load(p.pathTemp, 'Tbc');
Tmat = zeros(nb,12);
for i = 1:12
    Tmat(:,i) = gridToMatrix(Tbc(:,:,:,i), [], p.pathBoxes, p.pathGrid);
end
%
% Load Light:
%
disp('load light')
if p.bUse_parday_light
    load 'Transport Matrix/parday';
end
L0 = zeros(nb,365/p.dtTransport); %divide into half days

% tic
for i = 1:(365/p.dtTransport)
    if p.bUse_parday_light
        L0(:,i) = 1e6*parday(:,i)/(24*60*60).*exp(-p.kw*Zbox);
    else
        % Calculate light:
        
        L0(:,i) = p.EinConv*p.PARfrac*daily_insolation(0,Ybox,i*p.dtTransport,1).*exp(-p.kw*Zbox);
        
    end
end
% toc
%
% Calculate sinking matrix:
%
sim(1).dznom=dznom;
[Asink,p] = calcSinkingMatrix(p, sim, nz);
%
% Matrices for saving the solution:
%
disp('crate output matrix')
iSave = 0;
nSave = floor(p.tEnd/p.tSave) + sign(mod(p.tEnd,p.tSave));
nSave = 12;
nSavefile = 1;
sim = load(p.pathGrid,'x','y','z','dznom','bathy');
sim.N = single(zeros(length(sim.x), length(sim.y), length(sim.z),nSave));
if bSilicate
    sim.Si = sim.N;
end
sim.DOC = sim.N;
sim.B = single(zeros(length(sim.x), length(sim.y), length(sim.z), p.n-p.idxB+1, nSave));
sim.L = sim.N;
sim.T = sim.N;
tSave = [];


jDOC=nan(p.n-p.idxB+1,nb);
jFreal=jDOC;
jLreal=jDOC;
jLossPassive=jDOC;
jPOM=jDOC;
mortpred=jDOC;
mortHTL=jDOC;
mort2=jDOC;
mort=jDOC;

sim.jDOC=sim.B;
sim.jFreal=sim.B;
sim.jLreal=sim.B;
sim.jLossPassive=sim.B;
sim.jPOM=sim.B;
sim.mortpred=sim.B;
sim.mortHTL=sim.B;
sim.mort2=sim.B;
sim.mort=sim.B;






%
% Matrices for annual averages:
%
if bCalcAnnualAverages
    sim.ProdGrossAnnual = zeros( nb,1 );
    sim.ProdNetAnnual = zeros( nb,1 );
    sim.ProdHTLAnnual = zeros( nb,1 );
    sim.BpicoAnnualMean = zeros( nb,1 );
    sim.BnanoAnnualMean = zeros( nb,1 );
    sim.BmicroAnnualMean = zeros( nb,1 );
end

% ---------------------------------------
% Run transport matrix simulation
% ---------------------------------------
disp('Starting simulation')
sLibname = loadNUMmodelLibrary();
dtTransport = p.dtTransport;
n = p.n;

tic
for i=1:simtime
    %
    % Test for time to change monthly transport matrix
    %
    if ismember(mod(i*p.dtTransport,365), p.dtTransport+cumsum(mon))
        % Load TM for this month
        load(strcat(p.pathMatrix, sprintf('%02i.mat',month+1)));
        
        Aexp = function_convert_TM_positive(Aexp);
        Aimp = function_convert_TM_positive(Aimp);
        
        % Preparing for timestepping. 43200s.
        load(p.pathGrid,'deltaT')
        Aexp = Ix + (dtTransport*24*60*60)*Aexp;
        Aimp = Aimp^(dtTransport*24*60*60/deltaT);
        
        % Set monthly mean temperature
        T = Tmat(:,month+1);
        Nwoa=gridToMatrix(squeeze(N_WOAtoMITgcmEcco(:,:,:,month+1)), [], p.pathBoxes, p.pathGrid);
        
        month = mod(month + 1, 12);
    end
    
    %
    % Enforce minimum B concentration
    %
    for k = 1:n
        u(u(:,k)<p.umin(k),k) = p.umin(k);
    end
    % restoring condition on N
    u(:, ixN)=u(:, ixN)+(Nwoa-u(:, ixN)).*(p.dtTransport/p.t_restore);
    % calculate new light
    L = L0(:,mod(i,365/dtTransport)+1).*exp(-(p.kc.*sum(u(:,ixB(1:end-p.ixPOM)),2)).*Zbox);
    dt = p.dt;
    
    if ~isempty(gcp('nocreate'))
        parfor k = 1:nb
            %
            % Integrate ODEs:
            %
            u(k,:) = calllib(sLibname, 'f_simulateeuler', ...
                u(k,:), L(k), T(k), dtTransport, dt);
        end
    else
        for k = 1:nb
            u(k,:) = calllib(sLibname, 'f_simulateeuler', ...
                u(k,:),L(k), T(k), dtTransport, dt);
            %u(k,1) = u(k,1) + 0.5*(p.u0(1)-u(k,1))*0.5;
            % If we use ode23:
            %[t, utmp] = ode23(@fDerivLibrary, [0 0.5], u(k,:), [], L(k));
            %u(k,:) = utmp(end,:);
        end
    end
    if ((floor(i*(p.dtTransport/p.tSave)) > floor((i-1)*(p.dtTransport/p.tSave))) || (i==simtime))
        %tic
        for kkk=1:nb
            calllib(loadNUMmodelLibrary(), 'f_calcderivatives',u(kkk,:),L(kkk),T(kkk), 0.0,0*u(kkk,:)');
            zero = zeros(p.n-p.idxB+1,1);
            [~, jDOC(:,kkk), ~, ~, ~, jFreal(:,kkk), ~, ~, ~, ~, jLossPassive(:,kkk), ~, ...
                jLreal(:,kkk), jPOM(:,kkk), mortpred(:,kkk), mortHTL(:,kkk),...
                mort2(:,kkk), mort(:,kkk)] = calllib(loadNUMmodelLibrary(), 'f_getrates', ...
                zero, zero, zero, zero, zero, zero,zero, zero, zero, zero, zero, ...
                zero,zero, zero, zero, zero, zero, zero);
        end
        %toc
    end
    % Transport
    
    u =  Aimp*(Aexp*u);
    % Sinking
    for j = p.idxSinking
        sinktmp=single(matrixToGrid(u(:,j), [], p.pathBoxes, p.pathGrid));
        b=squeeze(reshape(sinktmp,[],1,size(sinktmp,3)))';
        b=squeeze(Asink(j,:,:))*b;
        c=reshape(b',size(sinktmp,1),size(sinktmp,2),[]);
        u(:, j) = gridToMatrix(c, [], p.pathBoxes, p.pathGrid);
    end
    
    %
    % Save timeseries in grid format
    %
    if ((floor(i*(p.dtTransport/p.tSave)) > floor((i-1)*(p.dtTransport/p.tSave))) || (i==simtime))
        %if ((mod(i*p.dtTransport,p.tSave) < mod((i-1)*p.dtTransport,p.tSave)) || (i==simtime))
        fprintf('t = %u days',floor(i/2))
        
        if any(isnan(u))
            warning('NaNs after running current time step');
            keyboard
        end
        
        
        iSave = iSave + 1;
        sim.N(:,:,:,iSave) = single(matrixToGrid(u(:,ixN), [], p.pathBoxes, p.pathGrid));
        sim.DOC(:,:,:,iSave) = single(matrixToGrid(u(:,ixDOC), [], p.pathBoxes, p.pathGrid));
        if bSilicate
            sim.Si(:,:,:,iSave) = single(matrixToGrid(u(:,ixSi), [], p.pathBoxes, p.pathGrid));
        end
        for j = 1:p.n-p.idxB+1
            sim.B(:,:,:,j,iSave) = single(matrixToGrid(u(:,ixB(j)), [], p.pathBoxes, p.pathGrid));
            sim.jDOC(:,:,:,j,iSave) = single(matrixToGrid(jDOC(j,:)', [], p.pathBoxes, p.pathGrid));
            sim.jFreal(:,:,:,j,iSave) = single(matrixToGrid(jFreal(j,:)', [], p.pathBoxes, p.pathGrid));
            sim.jLreal(:,:,:,j,iSave) = single(matrixToGrid(jLreal(j,:)', [], p.pathBoxes, p.pathGrid));
            sim.jLossPassive(:,:,:,j,iSave) = single(matrixToGrid(jLossPassive(j,:)', [], p.pathBoxes, p.pathGrid));
            sim.jPOM(:,:,:,j,iSave) = single(matrixToGrid(jPOM(j,:)', [], p.pathBoxes, p.pathGrid));
            sim.mortpred(:,:,:,j,iSave) = single(matrixToGrid(mortpred(j,:)', [], p.pathBoxes, p.pathGrid));
            sim.mortHTL(:,:,:,j,iSave) = single(matrixToGrid(mortHTL(j,:)', [], p.pathBoxes, p.pathGrid));
            sim.mort2(:,:,:,j,iSave) = single(matrixToGrid(mort2(j,:)', [], p.pathBoxes, p.pathGrid));
            sim.mort(:,:,:,j,iSave) = single(matrixToGrid(mort(j,:)', [], p.pathBoxes, p.pathGrid));
        end
  
        
        sim.L(:,:,:,iSave) = single(matrixToGrid(L, [], p.pathBoxes, p.pathGrid));
        sim.T(:,:,:,iSave) = single(matrixToGrid(T, [], p.pathBoxes, p.pathGrid));
        tSave = [tSave, i/p.dtTransport];
        fprintf('.\n');
        if iSave==12
            sim.t = tSave;
            sim.p = p;
            thisoutname=[p.outname,'_nr',num2str(nSavefile),'.mat'];
            save(fullfile(string(p.outdir),thisoutname),'sim','-v7.3')
            nSavefile=nSavefile+1;
            iSave=0;
            
        end
        
        
    end
    %
    % Update annual averages:
    %
    if bCalcAnnualAverages
        for k = 1:nb
            [ProdGross1, ProdNet1,ProdHTL1,eHTL,Bpico1,Bnano1,Bmicro1] = ...
                getFunctions(u(k,:), L(k), T(k));
            sim.ProdGrossAnnual(k) = sim.ProdGrossAnnual(k) + ProdGross1/(p.tEnd*2);
            sim.ProdNetAnnual(k) = sim.ProdNetAnnual(k) + ProdNet1/(p.tEnd*2);
            sim.ProdHTLAnnual(k) = sim.ProdHTLAnnual(k) + ProdHTL1/(p.tEnd*2);
            sim.BpicoAnnualMean(k) = sim.BpicoAnnualMean(k) + Bpico1/(p.tEnd*2*365);
            sim.BnanoAnnualMean(k) = sim.BnanoAnnualMean(k) + Bnano1/(p.tEnd*2*365);
            sim.BmicroAnnualMean(k) = sim.BmicroAnnualMean(k) + Bmicro1/(p.tEnd*2*365);
        end
    end
    
end
time = toc;
fprintf('Solving time: %2u:%02u:%02u\n', ...
    [floor(time/3600), mod(floor(time/60),60), floor(mod(time,60))]);
% ---------------------------------------
% Put results into sim structure:
% ---------------------------------------
sim.t = tSave; % days where solution was saved
sim.p = p;
sim.Ntot = calcGlobalN(sim);
sim.B(sim.B<0) = 0.;
sim.DOC(sim.DOC<0) = 0.;

if bCalcAnnualAverages
    tmp = single(matrixToGrid(sim.ProdGrossAnnual, [], p.pathBoxes, p.pathGrid));
    sim.ProdGrossAnnual = squeeze(tmp(:,:,1));
    tmp = single(matrixToGrid(sim.ProdNetAnnual, [], p.pathBoxes, p.pathGrid));
    sim.ProdNetAnnual = squeeze(tmp(:,:,1));
    tmp = single(matrixToGrid(sim.ProdHTLAnnual, [], p.pathBoxes, p.pathGrid));
    sim.ProdHTLAnnual = squeeze(tmp(:,:,1));
    tmp = single(matrixToGrid(sim.BpicoAnnualMean, [], p.pathBoxes, p.pathGrid));
    sim.BpicoAnnualMean = squeeze(tmp(:,:,1));
    tmp = single(matrixToGrid(sim.BnanoAnnualMean, [], p.pathBoxes, p.pathGrid));
    sim.BnanoAnnualMean = squeeze(tmp(:,:,1));
    tmp = single(matrixToGrid(sim.BmicroAnnualMean, [], p.pathBoxes, p.pathGrid));
    sim.BmicroAnnualMean = squeeze(tmp(:,:,1));
end
end
