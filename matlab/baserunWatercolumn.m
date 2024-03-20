% paramitnr=1
%% Define points to use:
% lat = -5.64390045591733;
% lon = -136.510730968337;

%% Load datapoints
load(fullfile('data_to_compare_to','lon_lat_TN045.mat'));
load(fullfile('data_to_compare_to','lon_lat_TN050.mat'));
load(fullfile('data_to_compare_to','lon_lat_clarke.mat'));
load(fullfile('data_to_compare_to','lon_lat_boyd.mat'));
load(fullfile('data_to_compare_to','lon_lat_tsuda.mat'));
load(fullfile('data_to_compare_to','lon_lat_cermeno.mat'));
load(fullfile('data_to_compare_to','lon_lat_maranonan.mat'));
load(fullfile('data_to_compare_to','lon_lat_tayler.mat'));

lat=[lat_AMT1; lat_AMT2; lat_Boyd; lat_CCE; lat_CERMENO; lat_clarke; lat_CRD; lat_EB; lat_HOT; lat_TN045; lat_TN050; lat_Tsuda];
lon=[lon_AMT1; lon_AMT2; lon_Boyd; lon_CCE; lon_CERMENO; lon_clarke; lon_CRD; lon_EB; lon_HOT; lon_TN045; lon_TN050; lon_Tsuda];

lat(134)=lat(134)+2;
% i=53;
% lon(i)=[];lat(i)=[];
% i=58;
% lon(i)=[];lat(i)=[];
% i=58;
% lon(i)=[];lat(i)=[];
% i=65;
% lon(i)=[];lat(i)=[];
% i=65;
% lon(i)=[];lat(i)=[];
% i=95;
% lon(i)=[];lat(i)=[];
% i=95;
% lon(i)=[];lat(i)=[];
% lon(95)=-12;%move out from the coast
% lat(96)=-63;%move out from the coast
% lon(96)=-72;%move out from the coast

lon2(1)=lon(11);
lat2(1)=lat(11);

lon2(2)=lon(134);
lat2(2)=lat(134);

lon2(3)=lon(74);
lat2(3)=lat(74);

lon2(4)=lon(118);
lat2(4)=lat(118);

lon2(5)=lon(22);
lat2(5)=lat(22);

lat=lat2(2:3); % HOT  CCE
lon=lon2(2:3); %HOT  CCE

% lat=lat2(3); % CCE
% lon=lon2(3); %CCE


%% define number of random iterations
nrRandIter=1000;
%% Choose number of size classes:
n = 10;
%% Choose number of particle size classes:
nPOM = 2;
%% Choose sumilation time
yearsSim=3650; %1095=3 years
yearsSave = 365*8;
%% Decide which parameters to randomize
%   the parameters are: 'r_star_d','rand_y','rand_rlstar','aF_random',
%   'c_passive_random','alphaMax_rand','alphaR_rand', and 'mortHTL_rand'

theParam={'r_star_d','rand_y','rand_rlstar','aF_random','c_passive_random','alphaMax_rand','alphaR_rand','mortHTL_rand',...
    'rand_rhocn','rand_fracHTLN','rand_epsilonL','rand_alphaN','rand_epsilonF','rand_cF','rand_beta','rand_sigma',...
    'rand_remin2','rand_reminF','rand_rhoval','rand_mort2val','rand_remPOM','rand_vel1','rand_vel2'};

% theRandoms={theParam{paramitnr}};
theRandoms=theParam;

% theRandoms={'r_star_d'};


time = datestr(clock,'YYYY_mm_dd_HH_MM');
%% Define random parameters
[s, meanval, sigma]= defineRandomParamsWide;

%% assemble random parameter
% ------------------------------------------------------------------------
isnormal=false(1,length(theParam));
isnormal([9,10,11,12,13,14,15,16,17,19,20,23])=true;
if strcmp(theRandoms,'all')
    theRandoms=theParam;
end
for i=1:length(theParam)
    
    
    meanval_tmp=meanval.(theParam{i});
    sigma_tmp=sigma.(theParam{i});
    
    if isnormal(i)
        pd=makedist('normal','mu',exp(meanval_tmp),'sigma',exp(sigma_tmp));
        pd=truncate(pd,exp(meanval_tmp)-2*exp(sigma_tmp),exp(meanval_tmp)+2*exp(sigma_tmp));
        
    else
        pd=makedist('lognormal','mu',meanval_tmp,'sigma',sigma_tmp);
        pd=truncate(pd,exp(meanval_tmp-2*sigma_tmp),exp(meanval_tmp+2*sigma_tmp));
    end
    
    rand_tmp=random(pd,1,nrRandIter);
    if isequal(theParam{i},'rand_y')
        % calculate alpha_L: alpha_l=(3.*rand_y)./(4*p_const)
        p_const=0.04;
        rand_tmp=(3.*rand_tmp)./(4*p_const);
    end
    % if the parameter is not to be randomized use the mean value instead:
    if ~any(strcmp(theRandoms,theParam{i}))
        rand_tmp(:)=exp(meanval_tmp);
    end
    randParam(:,i)=rand_tmp';
end
meanval_tmp=0.03;
kc = (meanval_tmp*2-meanval_tmp*0).*rand(1,nrRandIter);
randParam(:,size(randParam,2)+1)=kc';
%% Create output directory
% topLevelFolder = fullfile(pwd,'..','..','results');
% files = dir(topLevelFolder);
% subFolders = files([files.isdir]); % Extract only those that are directories.
% % Get only the folder names into a cell array.
% subFolderNames = {subFolders(3:end).name};

resultdir=fullfile(pwd,'..','..','results','wide');
% resultdir=fullfile(pwd,'results');
% dirnr=0;
% doesthisexist=true;
% while doesthisexist==true
%     dirnr=dirnr+1;

if length(theRandoms)==1 %if only one variable paramenter
    newdirname=fullfile(resultdir,[theRandoms{1},time]);
elseif length(theRandoms)==length(theParam)
    newdirname= fullfile(resultdir,['all_random',time]);
elseif isempty(theRandoms) %if none
    newdirname=fullfile(resultdir,['mean',thisdate]);
else
    disp('define outdir')
    return
end
% doesthisexist=isfolder(newdirname);

% end


% end
if ~isfolder(newdirname)
    mkdir(string(newdirname));
end
outdir= newdirname;
disp(['outdir is: ',outdir])


%find the index for mortHTL:
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
% tic
arandnr=rand(1);
for itnr=1:nrRandIter
    tic
    disp(['itnr: ',num2str(itnr),' of ',num2str(nrRandIter)]);
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
    
    p = setupGeneralistsPOM(n,nPOM);
    p = parametersWatercolumn(p);
    p.tEnd = yearsSim;
    p.tSaveFrom = yearsSave;
    
    restoreon=[false true];
    strout=['false' 'true'];
    for restorenr=1:2
        p.restore = restoreon(restorenr);
        p.t_restore=15/3;
%         p.kc=kc(itnr)/(106/16)/12.011;
        p.kpoc=3*10^(-5);%0.05;
        
        setHTL(randParam(itnr,mortHTL_rand_idx), 1/500^1.5, false, false);
        for p_point=2 %1:length(lon)
            %         tic
            simOutput = simulateWatercolumn(p, lat(p_point),lon(p_point));
            
            %% Output data
            simOutput.s=s;
            simOutput.lon=lon(p_point);
            simOutput.lat=lat(p_point);
            simOutput.rand_Param=randParam(itnr,:);
            %
            outname=['p_point_',num2str(p_point),'_restore_',strout(restorenr),'_nrRanditer_',num2str(itnr),'_',time,'_rand',num2str(arandnr),'.mat'];
            save(fullfile(string(outdir),outname),'simOutput','-v7.3')
            %             toc
        end
    end
    toc
end
% toc
disp('done')

