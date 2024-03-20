%
% Simulate a single water column from the global transport matrix.
% Conservation is enforced rather crudely.
%
% Tranport matrices must be downloaded from http://kelvin.earth.ox.ac.uk/spk/Research/TMM/TransportMatrixConfigs/
% (choose MITgcm_ECCO_v4), and be put into the location 'NUMmodel/TMs'
%
% Input:
%  p: parameter structure from parametersGlobal
%  lat, lon: latitude and longitude
%  sim: (optional) simulation to use for initial conditions
%
% Input options:
%  bCalcAnnualAverages: increases the simulation time by a factor 2-3
%  bExtractcolumn: (logical) forces the extraction of a water column from
%            the transport matrix. Only used when the core code is changed.
%  dayFixed: if non-zero then the run is done with fixed conditions at the
%            specified day.
%
% Output:
%  sim: structure with simulation results
%
function sim = simulateWatercolumn(p, lat, lon, sim, options)

arguments
    p struct;
    lat = 0;
    lon = 0;
    sim struct = [];
    options.bExtractcolumn logical = false; % Extract the watercolumn even though a saved one exists
    options.bRecalcLight logical = false; % Recalc the light (different from the extracted watercolumn)
    options.dayFixed double = 0;
end
%
% Get the watercolumn parameters if they are not already set:
%
if ~isfield(p,'nameModel')
    p = parametersWatercolumn(p);
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

%disp('Preparing simulation')
%
% Set path:
%
path = pwd;
addpath(strcat(path,'/Transport matrix'));

% ---------------------------------------
% Find the indices that corresponds to the selected water column:
% ---------------------------------------
if isempty(sim)
    sim = load(p.pathGrid,'x','y','z','dznom','bathy'); % Load grid
end
load(p.pathBoxes, 'nb', 'Ybox', 'Zbox');
idx = calcGlobalWatercolumn(lat, lon, sim); % Find the indices into matrix
xx = matrixToGrid((1:nb)', [], p.pathBoxes, p.pathGrid); % Find the indices into the grid
idxGrid = squeeze(xx(idx.x, idx.y, idx.z));
nGrid = length(idxGrid);

if nGrid==0
    error('Selected latitude and longitude is on land.')
end
%
% Check if a file with the water column exists; if not extract from TM and
% save:
%
sFile = sprintf('Watercolumn_%s_lat%03i_lon%02i.mat',p.TMname,lat,lon);
if exist(sFile,'file') && ~options.bExtractcolumn
    load(sFile);
end

if ~exist('versionTMcolumn')
    versionTMcolumn=0;
end
versionTMcolumnCurrent = 2; % Current version of the water column
if (versionTMcolumn~=versionTMcolumnCurrent) || options.bExtractcolumn  % Extract water column if the loaded one is too old
    versionTMcolumn = versionTMcolumnCurrent;
    fprintf('-> Extracting water column from transport matrix')
    %
    % Check that transport matrix files exist:
    %
    if ~exist(p.pathBoxes,'file')
        error( 'Error: Cannot find transport matrix file: %s', p.pathBoxes);
    end
    % Load TMs
    for month=1:12
        matrix = load(strcat(p.pathMatrix, sprintf('%02i.mat',month)),'Aimp');
        AimpM(month,:,:) = full(function_convert_TM_positive(matrix.Aimp(idxGrid,idxGrid)));

        % Preparing for timestepping. 43200s.
        temp = load(p.pathGrid,'deltaT');
        %AexpM(month,:,:) = Ix(idxGrid,idxGrid) + squeeze((12*60*60)*AexpM(month,:,:));
        AimpM(month,:,:) = squeeze(AimpM(month,:,:))^(p.dtTransport*24*60*60/temp.deltaT);
        fprintf('.');
    end
    %
    % Load temperature:
    %
    load(p.pathTemp, 'Tbc');
    Tmat = zeros(nb,12);
    for i = 1:12
        Tmat(:,i) = gridToMatrix(Tbc(:,:,:,i), [], p.pathBoxes, p.pathGrid);
    end
    Tmat = Tmat(idxGrid,:); % Use only the specific water column
    %
    % Calc Light:
    %
    L0 = zeros(nGrid,730);
    for i = 1:730
        zup = sim.z - 0.5*sim.dznom; % Top of a box
        zup = zup(1:length(idxGrid));
        dz = sim.dznom(1:length(idxGrid));
        Lup = p.EinConv*p.PARfrac*daily_insolation(0,Ybox(idxGrid),i/2,1).*exp(-p.kw*zup);
%         Lup_save(:,i)=Lup;
        L0(:,i) = Lup.*(1-exp(-p.kw*dz))./(p.kw*dz);
    end

    fprintf('\n');
    save(sFile,'AimpM','Tmat','L0','versionTMcolumn');
end
zup = sim.z - 0.5*sim.dznom; % Top of a box
zup = zup(1:length(idxGrid));
dz = sim.dznom(1:length(idxGrid));
if options.bRecalcLight
    %
    % Calc Light:
    %
    L0 = zeros(nGrid,730);
    for i = 1:730
        zup = sim.z - 0.5*sim.dznom; % Top of a box
        zup = zup(1:length(idxGrid));
        dz = sim.dznom(1:length(idxGrid));
        
        Lup = p.EinConv*p.PARfrac*daily_insolation(0,Ybox(idxGrid),i/2,1).*exp(-p.kw*zup);
        L0(:,i) = Lup.*(1-exp(-p.kw*dz))./(p.kw*dz);
        %Lold = p.EinConv*p.PARfrac*daily_insolation(0,Ybox(idxGrid),i/2,1).*exp(-p.kw*Zbox(idxGrid));
    end
end
Lup1=nan(length(idxGrid),730);
for i = 1:730
    Lup1(:,i) = p.EinConv*p.PARfrac*daily_insolation(0,Ybox(idxGrid),i/2,1);
end

%calculate width between midpoints
zdiff=[0;Zbox(idxGrid)];
zdiff=zdiff(2:end)-zdiff(1:end-1);
% calculate lower radius;
rho=0.4*10^(-6);
radii_lower=(3.*p.mLower./(4*pi*rho)).^(1/3);
ixPOM=(radii_lower.*2)>0.7;
% Get sinking matrix:
[Asink,p] = calcSinkingMatrix(p, sim, nGrid);
% ---------------------------------------
% Initialize run:
% ---------------------------------------
simtime = p.tEnd/p.dtTransport; %simulation time in TM time steps

% Preparing timestepping
mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
%
% Initial conditions:
%
if isfield(sim,'B')
    disp('Starting from previous simulation.');
    u(:,ixN) = sim.N(:,end);
    u(:, ixDOC) = sim.DOC(:,end);
    if bSilicate
        u(:, ixSi) = sim.Si(:,end);
    end
    for i = 1:p.n -p.idxB+1
        u(:, ixB(i)) = sim.B(:,i,end);
    end
else
    if exist(strcat(p.pathN0),'file')
        load(p.pathN0, 'N_WOAtoMITgcmEcco');
        Nmat = zeros(nb,12);
        for i = 1:12
        Nmat(:,i) = gridToMatrix(N_WOAtoMITgcmEcco(:,:,:,i), [], p.pathBoxes, p.pathGrid);
        end
        u(:, ixN) = Nmat(:,1);
        
    else
        u(:, ixN) = p.u0(p.idxN)*ones(nb,1);
    end
    u(:, ixDOC) = zeros(nb,1) + 60*12.011;
    if bSilicate
        u(:, ixSi) = zeros(nb,1) + p.u0(ixSi);
    end
    u(:, ixB) = ones(nb,1)*p.u0(ixB);
    u = u(idxGrid,:); % Use only the specific water column
end
p.u0(ixN) = u(nGrid,ixN); % Use the nitrogen concentration in the last grid cell as BC
Nmat = Nmat(idxGrid,:); % Use only the specific water column
% indsæt ny her og lav p.u0(ixN) til den samme værdi
switch p.thesite
    case 'CCE'
        % data from CCE liter
        N_local=[13.8532777354687 39.3491198206705 33.5374458115977 60.5422539914537 131.474303431804 198.404984424086 279.986493988072 398.083396913394 415.746023622101 483.086770054323 540.464964824627 388.048306774518 401.941571187702 420.222182257682 442.597650207337 467.898015928189 495.684544754558 524.933522466526 554.182500178493 554.182500178493 554.182500178493 554.182500178493 554.182500178493]';
        
    case 'HOT'
        % data from HOT + WOA
        N_local=[0.612707230241040 0.574083642479426 0.541687050639254 0.585276022089345 0.697019876895530 1.11506484276745 5.28336070414508 20.7111351484137 64.5418789508053 149.280221034132 301.705695100030 515.906154154203 515.906154154203 515.906154154203 515.906154154203 515.906154154203 515.906154154203 515.906154154203 515.906154154203 515.906154154203]';

end
N_local=repmat(N_local(1:length(idxGrid)),[1 12]);
Nmat=N_local;
u(:,ixN)=Nmat(:,1);
p.u0(ixN) =u(nGrid,ixN); 

Nwoa=Nmat(:,1);

   
%
% Matrices for saving the solution:
%
iSave = 0;
nSave = floor((p.tEnd-p.tSaveFrom)/p.tSave) + sign(mod(p.tEnd,p.tSave));

sim.N = zeros(nGrid,nSave);
if bSilicate
    sim.Si = sim.N;
end
sim.DOC = sim.N;
sim.B = zeros(nGrid, p.n-p.idxB+1, nSave);
sim.L = sim.N;
sim.T = sim.N;
sim.Nloss = zeros(1,nSave);
sim.NlossHTL = sim.Nloss;
sim.jLreal=sim.B;
sim.jDOC=sim.B;
sim.jFreal=sim.B;

sim.jresp=sim.B;
sim.jMax=sim.B;

sim.jLossPassive=sim.B;
sim.jPOM=sim.B;
sim.mortpred=sim.B;
sim.mortHTL=sim.B;
sim.mort2=sim.B;
sim.mort=sim.B;


tSave = [];

% ---------------------------------------
% Run transport matrix simulation
% ---------------------------------------
%disp('Starting simulation')
%tic
for i = 1:simtime
    %
    % Find the iteration time
    %
    if (options.dayFixed ~= 0)
        iTime = options.dayFixed / p.dtTransport;
    else
        iTime = i;
    end
    %
    % Test for time to change monthly temperature
    %
   
    if (ismember(mod(i*p.dtTransport,365), p.dtTransport+cumsum(mon)) || i==1)
        % Set monthly mean temperature
        month = find(p.dtTransport+cumsum(mon) >= mod(i*p.dtTransport,365),1);
        T = Tmat(:,month);
        if p.restore
        Nwoa=Nmat(:,month);
        end
    end
   
    %
    % Enforce minimum concentration
    %
    for k = 1:nGrid
        u(k,u(k,:)<p.umin) = p.umin(u(k,:)<p.umin);
    end
    
    if p.restore
    % restoring condition on N
    u(:,ixN)=u(:,ixN)+(Nwoa-u(:,ixN)).*(p.dtTransport/p.t_restore);
    else
         % Bottom BC for nutrients:
         u(end, p.idxN) = u(end, p.idxN) +  p.dtTransport* ...
             p.DiffBottom/sim.dznom(nGrid)*(p.u0(p.idxN)-u(end,p.idxN));
    end
    %     if i==300
    %         disp('hi')
    %     end
    % calculate new light
    
    ktot=p.kw+p.kpoc.*cumsum(sum(u(:,ixPOM),2).*zdiff);
    Lup=Lup1(:,mod(i,365/p.dtTransport)+1).*exp(-ktot.*zup);
    L = Lup.*(1-exp(-ktot.*dz))./(ktot.*dz);

    
    


    
    
    
    %
    % Run Euler time step for half a day:
    %
    dt = p.dt;
    dtTransport = p.dtTransport;


    for k = 1:nGrid
        u(k,:) = calllib(loadNUMmodelLibrary(), 'f_simulateeuler', ...
            u(k,:),L(k), T(k), dtTransport, dt);
    end

    if any(isnan(u))
        warning('NaNs after running current time step');
        %return
    end
    %
    % Transport
    %
    u =  squeeze(AimpM(month,:,:)) * u; % Vertical diffusion
    % Sinking:
    for j = p.idxSinking
        u(:,j) = squeeze(Asink(j,:,:)) * u(:,j);
    end

    %
    % Save timeseries in grid format
    %
    if (((floor(i*(p.dtTransport/p.tSave)) > floor((i-1)*(p.dtTransport/p.tSave))) && (i*p.dtTransport>p.tSaveFrom)) || (i==simtime))
        iSave = iSave + 1;
        sim.N(:,iSave) = u(:,ixN);
        sim.DOC(:,iSave) = u(:,ixDOC);
        if bSilicate
            sim.Si(:,iSave) = u(:,ixSi);
        end
        for j = 1:p.n-p.idxB+1
            sim.B(:,j,iSave) = u(:,ixB(j));
        end
        sim.L(:,iSave) = L;
        sim.T(:,iSave) = T;
        % Loss to HTL and POM:
        for j = 1:nGrid
            rates = getRates(p,u(j,:),L(j),T(j));
            sim.jLreal(j,:,iSave)=rates.jLreal;
            sim.jDOC(j,:,iSave)=rates.jDOC;
            sim.jFreal(j,:,iSave)=rates.jFreal;
            
            sim.jresp(j,:,iSave)=rates.jR;
            sim.jMax(j,:,iSave)=rates.jMax;
            
            sim.jLossPassive(j,:,iSave)=rates.jLossPassive;
            
            sim.jPOM(j,:,iSave)=rates.jPOM;
            
            sim.mortpred(j,:,iSave)=rates.mortpred;
            sim.mortHTL(j,:,iSave)=rates.mortHTL;
            sim.mort2(j,:,iSave)=rates.mort2;
            sim.mort(j,:,iSave)=rates.mort;

        end

        tSave = [tSave, i*p.dtTransport];
    end

end
%time = toc;
%fprintf('Solving time: %2u:%02u:%02u\n', ...
%    [floor(time/3600), mod(floor(time/60),60), floor(mod(time,60))]);
% ---------------------------------------
% Put results into sim structure:
% ---------------------------------------
sim.t = tSave; % days where solution was saved
sim.p = p;
%sim.Ntot = calcGlobalN(sim);
% sim.B(sim.B<0) = 0.;
% sim.DOC(sim.DOC<0) = 0.;
sim.z = sim.z(1:length(idx.z));
sim.dznom = sim.dznom(1:length(idx.z));
sim.lat = lat;
sim.lon = lon;
sim.Nmat=Nmat;
sim.idx=idx;

% sim.Ntot = (sum(sim.N.*(sim.dznom*ones(1,length(sim.t)))) + ... % gN/m2 in dissolved phase
%     sum(squeeze(sum(sim.B,2)).*(sim.dznom*ones(1,length(sim.t))))/5.68)/1000; % gN/m2 in biomass
% sim.Nprod = p.DiffBottom*(p.u0(p.idxN)-sim.N(end,:))/1000; % Diffusion in from the bottom; gN/m2/day
end
