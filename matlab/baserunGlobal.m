%% define number of random iterations
nrRandIter=1;
%% Choose number of size classes:
n = 10;
%% Choose number of particle size classes:
nPOM = 2;
%% Choose sumilation time
yearsSim=365*20; %365*1; %1095=3 years

%% Decide which parameters to randomize
theParam={'r_star_d','rand_y','rand_rlstar','aF_random','c_passive_random','alphaMax_rand','alphaR_rand','mortHTL_rand',...
    'rand_rhocn','rand_fracHTLN','rand_epsilonL','rand_alphaN','rand_epsilonF','rand_cF','rand_beta','rand_sigma',...
    'rand_remin2','rand_reminF','rand_rhoval','rand_mort2val','rand_remPOM','rand_vel1','rand_vel2'};
theRandoms=theParam;
theRandoms={};

thisdate = datestr(clock,'YYYY_mm_dd');
time = datestr(clock,'YYYY_mm_dd_HH_MM');

%% Define and assemble random parameters
[s, meanval, sigma]= defineRandomParams;

isnormal=false(1,length(theParam));
isnormal([9,10,12,14,15,16,17,19,20,22,23])=true;
if strcmp(theRandoms,'all')
    theRandoms=theParam;
end
% theRandoms={};
for i=1:length(theParam)
    
    meanval_tmp=exp(meanval.(theParam{i}));
%     sigma_tmp=sigma.(theParam{i});
%     
%     if isnormal(i)
%         pd=makedist('normal','mu',exp(meanval_tmp),'sigma',exp(sigma_tmp));
%         pd=truncate(pd,exp(meanval_tmp)-2*exp(sigma_tmp),exp(meanval_tmp)+2*exp(sigma_tmp));
%         
%     else
%         pd=makedist('lognormal','mu',meanval_tmp,'sigma',sigma_tmp);
%         pd=truncate(pd,exp(meanval_tmp-2*sigma_tmp),exp(meanval_tmp+2*sigma_tmp));
%     end
    
    rand_tmp = (meanval_tmp*1.2-meanval_tmp*0.8).*rand(1,nrRandIter)+meanval_tmp*0.8;
%     rand_tmp=random(pd,1,nrRandIter);
    if isequal(theParam{i},'rand_y')
        % calculate alpha_L: alpha_l=(3.*rand_y)./(4*p_const)
        p_const=0.04;
        rand_tmp=(3.*rand_tmp)./(4*p_const);
    end
    % if the parameter is not to be randomized use the mean value instead:
    if ~any(strcmp(theRandoms,theParam{i}))
        rand_tmp(:)=meanval_tmp;
    end
    randParam(:,i)=rand_tmp';
end
%% Define and create output folder
topLevelFolder = fullfile(pwd,'..','..','results');
files = dir(topLevelFolder);
subFolders = files([files.isdir]); % Extract only those that are directories.
% Get only the folder names into a cell array.
subFolderNames = {subFolders(3:end).name};

resultdir=fullfile(pwd,'..','..','results','20pct','global');

if length(theRandoms)==1 %if only one variable paramenter
%     newdirname=fullfile(resultdir,[theRandoms{1},num2str(dirnr)]);
    outdir=fullfile(resultdir,[theRandoms{1},thisdate]);
elseif length(theRandoms)==length(theParam)
%     newdirname= fullfile(resultdir,['all_random',num2str(dirnr)]);
    outdir= fullfile(resultdir,['all_random',thisdate]);
elseif isempty(theRandoms) %if none
%     newdirname=fullfile(resultdir,[theRandoms{1},num2str(dirnr)]);
    outdir=fullfile(resultdir,['mean',thisdate]);
else
    disp('define outdir')
    return
end

if ~isfolder(outdir)
mkdir(string(outdir));
end
disp(['outdir is: ',outdir])

%% Find the index for mortHTL:
r_star_d_idx=find(contains(theParam,'r_star_d'));
rand_y_idx=find(contains(theParam,'rand_y'));
rand_rlstar_idx=find(contains(theParam,'rand_rlstar'));
aF_random_idx=find(contains(theParam,'aF_random'));
c_passive_random_idx=find(contains(theParam,'c_passive_random'));
alphaMax_rand_idx=find(contains(theParam,'alphaMax_rand'));
alphaR_rand_idx=find(contains(theParam,'alphaR_rand'));
mortHTL_rand_idx=find(contains(theParam,'mortHTL_rand'));

rand_rhocn_idx=find(contains(theParam,'rand_rhocn'));
rand_fracHTLN_idx=find(contains(theParam,'rand_fracHTLN'));
rand_epsilonL_idx=find(contains(theParam,'rand_epsilonL'));
rand_alphaN_idx=find(contains(theParam,'rand_alphaN'));
rand_epsilonF_idx=find(contains(theParam,'rand_epsilonF'));
rand_cF_idx=find(contains(theParam,'rand_cF'));
rand_beta_idx=find(contains(theParam,'rand_beta'));
rand_sigma_idx=find(contains(theParam,'rand_sigma'));
rand_remin2_idx=find(contains(theParam,'rand_remin2'));
rand_reminF_idx=find(contains(theParam,'rand_reminF'));
rand_rhoval_idx=find(contains(theParam,'rand_rhoval'));
rand_mort2val_idx=find(contains(theParam,'rand_mort2val'));
rand_remPOM_idx=find(contains(theParam,'rand_remPOM'));
rand_vel1_idx=find(contains(theParam,'rand_vel1'));
rand_vel2_idx=find(contains(theParam,'rand_vel2'));



%% For each random iteration
tic
for itnr=1:nrRandIter
    arandnr=rand(1);
    copyfile ../input/input_orig.nlm ../input/input.nlm
    
    if ispc
        
        system(['C:\cygwin64\bin\sed.exe "s/.*rNstar.*/      rNstar = ',num2str(randParam(itnr,r_star_d_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*alphaL.*/      alphaL = ',num2str(randParam(itnr,rand_y_idx)),'/" ../input/tmp1.nlm > ../input/tmp.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*rLstar.*/      rLstar = ',num2str(randParam(itnr,rand_rlstar_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*alphaF.*/      alphaF = ',num2str(randParam(itnr,aF_random_idx)),'/" ../input/tmp1.nlm > ../input/tmp.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*cLeakage.*/      cLeakage = ',num2str(randParam(itnr,c_passive_random_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*alphaJ.*/      alphaJ = ',num2str(randParam(itnr,alphaMax_rand_idx)),'/" ../input/tmp1.nlm > ../input/tmp.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*cR.*/      cR = ',num2str(randParam(itnr,alphaR_rand_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        
        system(['C:\cygwin64\bin\sed.exe "s/.*rhoCN.*/      rhoCN = ',num2str(randParam(itnr,rand_rhocn_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*fracHTL_to_N.*/      fracHTL_to_N = ',num2str(randParam(itnr,rand_fracHTLN_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*epsilonL.*/      epsilonL = ',num2str(randParam(itnr,rand_epsilonL_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*alphaN.*/      alphaN = ',num2str(randParam(itnr,rand_alphaN_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*epsilonF*/      epsilonF = ',num2str(randParam(itnr,rand_epsilonF_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*cF.*/      cF = ',num2str(randParam(itnr,rand_cF_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*beta.*/      beta = ',num2str(randParam(itnr,rand_beta_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*sigma.*/      sigma = ',num2str(randParam(itnr,rand_sigma_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*remin2.*/      remin2 = ',num2str(randParam(itnr,rand_remin2_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*reminF.*/      reminF = ',num2str(randParam(itnr,rand_reminF_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*rhoval.*/      rhoval = ',num2str(randParam(itnr,rand_rhoval_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*mort2val.*/      mort2val = ',num2str(randParam(itnr,rand_mort2val_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*remPOM.*/      remPOM = ',num2str(randParam(itnr,rand_remPOM_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*vel1.*/      vel1 = ',num2str(randParam(itnr,rand_vel1_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        system(['C:\cygwin64\bin\sed.exe "s/.*vel2.*/      vel2 = ',num2str(randParam(itnr,rand_vel2_idx)),'/" ../input/tmp.nlm > ../input/tmp1.nlm']);
        
        movefile ../input/tmp1.nlm ../input/input.nlm
    elseif isunix
        system(['sed -i.bak "s/.*rNstar.*/      rNstar = ',num2str(randParam(itnr,r_star_d_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*alphaL.*/      alphaL = ',num2str(randParam(itnr,rand_y_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*rLstar.*/      rLstar = ',num2str(randParam(itnr,rand_rlstar_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*alphaF.*/      alphaF = ',num2str(randParam(itnr,aF_random_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*cLeakage.*/      cLeakage = ',num2str(randParam(itnr,c_passive_random_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*alphaJ.*/      alphaJ = ',num2str(randParam(itnr,alphaMax_rand_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*cR.*/      cR = ',num2str(randParam(itnr,alphaR_rand_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        
        system(['sed -i.bak "s/.*rhoCN.*/     rhoCN = ',num2str(randParam(itnr,rand_rhocn_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*fracHTL_to_N.*/     fracHTL_to_N = ',num2str(randParam(itnr,rand_fracHTLN_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*epsilonL.*/      epsilonL = ',num2str(randParam(itnr,rand_epsilonL_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*alphaN.*/      alphaN = ',num2str(randParam(itnr,rand_alphaN_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*epsilonF.*/      epsilonF = ',num2str(randParam(itnr,rand_epsilonF_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*cF.*/      cF = ',num2str(randParam(itnr,rand_cF_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*beta.*/      beta = ',num2str(randParam(itnr,rand_beta_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*sigma.*/      sigma = ',num2str(randParam(itnr,rand_sigma_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*remin2.*/      remin2 = ',num2str(randParam(itnr,rand_remin2_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*reminF.*/      reminF = ',num2str(randParam(itnr,rand_reminF_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*rhoval.*/      rhoval = ',num2str(randParam(itnr,rand_rhoval_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*mort2val.*/      mort2val = ',num2str(randParam(itnr,rand_mort2val_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*remPOM.*/      remPOM = ',num2str(randParam(itnr,rand_remPOM_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*vel1.*/      vel1 = ',num2str(randParam(itnr,rand_vel1_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
        system(['sed -i.bak "s/.*vel2.*/      vel2 = ',num2str(randParam(itnr,rand_vel2_idx)),'/" ../input/input.nlm > ../input/tmp1.nlm']);
%         movefile ../input/tmp1.nlm ../input/input.nlm
    end
    
p = setupGeneralistsPOM(n, nPOM,true); % Use 10 size groups and parallel execution
p = parametersGlobal(p,2); % Use MITgcm_ECCO
p.tEnd = yearsSim;
p.tSave = 365/12; % How often to save results (monthly)
p.s=s;
p.rand_Param=randParam(itnr,:);
p.outname=['nrRanditer_',num2str(itnr),'_',time,'_rand',num2str(arandnr)];
p.outdir=outdir;
setHTL(randParam(itnr,mortHTL_rand_idx), 1/500^1.5, false, false);
p.kc=0.00037701;
p.kw=0.05;
p.t_restore=15/3;


%
% Simulate
%
if exist(strcat(p.pathInit,'.mat'), 'file')
    % Load decent initial conditions
    disp('Loading initial conditions from file');
    
    simOutput = loadGlobal(p);
    simOutput = simulateGlobal(p,simOutput);
else
    simOutput = simulateGlobal(p);%,sim); % Simulate
end
% simOutput.B(simOutput.B<0)=0; % Get rid of negative biomasses
% disp('save output')
% simOutput.s=s;
% simOutput.rand_Param=randParam(itnr,:);
% outname=['nrRanditer_',num2str(itnr),'_',time,'_rand',num2str(arandnr),'.mat'];
% save(fullfile(string(outdir),outname),'simOutput','-v7.3')
end
toc

disp('done')
