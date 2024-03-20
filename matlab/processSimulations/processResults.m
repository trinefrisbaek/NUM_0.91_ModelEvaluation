try
    %% Calculate statistics and plot;
    mysite='CCE';
    thissetup1='wide';
%     thissetup=fullfile('final_both1'); % de initiale 100000
    thissetup=fullfile('final_CCEtoHOT'); % de 9 gode CCE genkøres og på HOT
    thissetup=fullfile('final_reduced_first'); % reduktion af range baseret på final_both1 for CCE
    thissetup=fullfile('final_reduced'); % yderligere reduktion af range baseret på final_reduced1 for CCE


    % use ratio AC (ACv1) or full AC (ACv2)?
    ACversion='ACv3';%'ACv2';
    AC_fraction=0.1;
    AC_fraction_F=0.1;
    AC_fraction_O=0.1;
    % use nudging?
    wNudge=false;
    % depth of statistic zone
    mindepth=0;
    maxdepth=200;

    Plot_ACpnm=false;
    Plot_HbacHnano=false;
    Plot_TOTpn=false;
    saveResult=true;

    EinConv=4.57;
    PARfrac=0.4;

    % start of statistic time span
    tstart=1; %365*2+1;

    switch mysite
        case 'HOT'
            sitenr = 1;
        case 'CCE'
            sitenr = 2;
    end

    %check if some files are processed already:
    thisout=['crash_',mysite,'_',thissetup,'.mat'];
    if ~isfile(thisout)
        files = dir(fullfile(thissetup1,thissetup));
        files(1:2)=[];
        thefilenames=string({files.name})';
        site=nan(length(thefilenames),1);
        nudging=nan(length(thefilenames),1);
        for i=1:length(thefilenames)
            site(i,1)=double(string(extractBetween(files(i).name,'p_point_','_restore')));
            iii=string(extractBetween(files(i).name,'_restore_','_firstRandNr'));
            if iii=="a"
                nudging(i,1)=1;
            elseif iii=="f"
                nudging(i,1)=0;
            end
        end


        siteidx=find(site==sitenr);
        if wNudge==true
            nudgidx=find(nudging==1);
        else
            nudgidx=find(nudging==0);
        end

        %% Find the idx for the simulations that fits both site and nudging
        [relevantIdx,~]=intersect(siteidx,nudgidx);
        theruns=1:length(relevantIdx);


        %% Load data for this site
        switch mysite
            case 'HOT'
                TayLanFig1a=load(fullfile('..','..','External_data','Taylor and Landry 2018','Taylor and Landry 2018 data for HOT',['HOT_TayLanFig1a_',num2str(maxdepth),'m.mat']));
                TayLanFig7=load(fullfile('..','..','External_data','Taylor and Landry 2018','Taylor and Landry 2018 data for HOT',['HOT_TayLanFig7_',num2str(maxdepth),'m.mat']));
                GasolFig3=load(fullfile('..','..','External_data','Taylor and Landry 2018','Taylor and Landry 2018 data for HOT',['HOT_GasolFig3_',num2str(maxdepth),'m.mat']));
                totPN=load(fullfile('..','..','External_data','Taylor and Landry 2018','Taylor and Landry 2018 data for HOT',['HOT_totPN_',num2str(maxdepth),'m.mat']));
                load(fullfile('..','..','External_data','Taylor and Landry 2018','daily_insolation_HOT.mat'))
                if maxdepth==200
                    nrBins=[0 138 87 138 96 6 3 2 2 1 0 1 0 0];
                elseif maxdepth==30
                    nrBins=[0 2 24 53 33 2 2 1 1 1 0 0 0 0];
                end
            case 'CCE'
                TayLanFig1a=load(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011',['CCE_TayLanFig1a_',num2str(maxdepth),'m.mat']));
                TayLanFig7=load(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011',['CCE_TayLanFig7_',num2str(maxdepth),'m.mat']));
                GasolFig3=load(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011',['CCE_GasolFig3_',num2str(maxdepth),'m.mat']));
                totPN=load(fullfile('..','..','External_data','Taylor and Landry 2018','Landry CEE 2004-2011',['CCE_totPN_',num2str(maxdepth),'m.mat']));
                load(fullfile('..','..','External_data','Taylor and Landry 2018','daily_insolation_CCE.mat'))
                nrBins=[0 108 127 154 141 88 61 31 23 12 7 0 1 0];
        end

        %% Get relevant info on model setup from first file on the list
        infoSetup=load(fullfile(files(relevantIdx(1)).folder,files(relevantIdx(1)).name));
        %number of generalists
        n=infoSetup.simOutput.p.ixEnd(1)-(infoSetup.simOutput.p.ixStart(1)-1);
        % number of random parameters
        nRandParam=length(infoSetup.simOutput.rand_Param);
        %number of model years
        nyears=length(infoSetup.simOutput.t)/365;

        Nmat=infoSetup.simOutput.Nmat(:,1);
        nlength=length(Nmat);
        diff_ref=Nmat-mean(Nmat,'omitnan');
        std_ref=sqrt(sum(diff_ref.^2,'omitnan')./nlength);

        %% Make statistic output matrices
        %Het vs. Autotrophic cf. Gasol fig 3.
        nOverlap_HC=nan(length(relevantIdx),1);
        %     std_HC_test_log=nOverlap_HC;
        %     std_HC_ref_log=nOverlap_HC;
        %     r_HC=nOverlap_HC;
        %     cRMSd_HC=nOverlap_HC;
        %
        %     nOverlap_HCAC=nOverlap_HC;
        %     std_HCAC_test_log=nOverlap_HC;
        %     std_HCAC_ref_log=nOverlap_HC;
        %     r_HCAC=nOverlap_HC;
        %     cRMSd_HCAC=nOverlap_HC;


        nOverlap_ACpico=nOverlap_HC;
        nOverlap_ACnano=nOverlap_HC;
        nOverlap_ACmicro=nOverlap_HC;
        std_ACpico_test_log=nOverlap_HC;
        std_ACnano_test_log=nOverlap_HC;
        std_ACmicro_test_log=nOverlap_HC;
        std_ACpico_ref_log=nOverlap_HC;
        std_ACnano_ref_log=nOverlap_HC;
        std_ACmicro_ref_log=nOverlap_HC;
        r_ACpico=nOverlap_HC;
        r_ACnano=nOverlap_HC;
        r_ACmicro=nOverlap_HC;
        cRMSd_ACpico=nOverlap_HC;
        cRMSd_ACnano=nOverlap_HC;
        cRMSd_ACmicro=nOverlap_HC;
        RMSd_ACpico=nOverlap_HC;
        RMSd_ACnano=nOverlap_HC;
        RMSd_ACmicro=nOverlap_HC;

        randParam=nan(length(relevantIdx),nRandParam);
        AC_pico_mean_bin_save=nan(length(relevantIdx),14);
        AC_nano_mean_bin_save=nan(length(relevantIdx),14);
        AC_micro_mean_bin_save=nan(length(relevantIdx),14);

        HC_Hbac_tot_mean_bin_save=nan(length(relevantIdx),14);
        HC_Hnano_tot_mean_bin_save=nan(length(relevantIdx),14);

        nOverlap_HC_Hbac=nOverlap_HC;
        std_HC_Hbac_test_log=nOverlap_HC;
        std_HC_Hbac_ref_log=nOverlap_HC;
        r_HC_Hbac=nOverlap_HC;
        cRMSd_HC_Hbac=nOverlap_HC;
        RMSd_HC_Hbac=nOverlap_HC;

        nOverlap_HC_Hnano=nOverlap_HC;
        std_HC_Hnano_test_log=nOverlap_HC;
        std_HC_Hnano_ref_log=nOverlap_HC;
        r_HC_Hnano=nOverlap_HC;
        cRMSd_HC_Hnano=nOverlap_HC;
        RMSd_HC_Hnano=nOverlap_HC;

        TOT_pico_mean_bin_save=nan(length(relevantIdx),14);
        TOT_nano_mean_bin_save=nan(length(relevantIdx),14);
        nOverlap_TOTpico=nOverlap_HC;nOverlap_TOTnano=nOverlap_HC; std_TOTpico_test_log=nOverlap_HC;
        std_TOTnano_test_log=nOverlap_HC; std_TOTpico_ref_log=nOverlap_HC; std_TOTnano_ref_log=nOverlap_HC;
        r_TOTpico=nOverlap_HC;  r_TOTnano=nOverlap_HC; cRMSd_TOTpico=nOverlap_HC; cRMSd_TOTnano=nOverlap_HC;
        RMSd_TOTpico=nOverlap_HC; RMSd_TOTnano=nOverlap_HC;

        N=nan(nlength,length(relevantIdx)+1);
        rmsd_N=nan(1,length(relevantIdx)+1);
        cRMSd_N=nan(1,length(relevantIdx)+1);
        std_N=nan(1,length(relevantIdx)+1);
        r_N=nan(1,length(relevantIdx)+1);


        %% Create bins for Autotrophic Carbon
        % same bins as Taylor & Landry
        %     bin_nr=1:14;
        AC_bin(1:14)=0;
        AC_pico_mean_bin=AC_bin;
        AC_pico_mean_bin(:)=NaN;
        AC_nano_mean_bin=AC_pico_mean_bin;
        AC_micro_mean_bin=AC_pico_mean_bin;
        %     AC_total_mean_bin=AC_pico_mean_bin;

        %     HC_phago_pico_mean_bin=AC_pico_mean_bin;
        %     HC_phago_nano_mean_bin=AC_pico_mean_bin;
        %     HC_phago_micro_mean_bin=AC_pico_mean_bin;

        TOT_pico_mean_bin=AC_pico_mean_bin;
        TOT_nano_mean_bin=AC_pico_mean_bin;

        sumB=nan(length(relevantIdx),n);

        AC_bin(2)=5;
        for ii=3:14
            AC_bin(ii)=AC_bin(ii-1)+0.5*AC_bin(ii-1);
        end
        AC_bin(end)=AC_bin(end)*2;
    else
        load(thisout)
    end

    %% DATA LOAD AND CALC STATISTICS
    for i=theruns
        % display for every 50 simulations
        if mod(i, 50) == 0
            disp([num2str(i),' of ',num2str(length(relevantIdx))])
        end
        thisfile=thefilenames(relevantIdx(i));
        thisfile_name(i)=thisfile;
        load(fullfile(files(relevantIdx(i)).folder,files(relevantIdx(i)).name))
        randParam(i,:)=simOutput.rand_Param;
        thedepths=find(simOutput.z>mindepth & simOutput.z<maxdepth);
        mass=simOutput.p.m(3:end);

        N(:,i)=mean(simOutput.N,2);
    rmsd_N(1,i)=sqrt(sum((N(:,i)-Nmat).^2)./nlength);
    diff_test=N(:,i)-mean(N(:,i),'omitnan');
    cRMSd_N(i)=sum((diff_test-diff_ref).^2,'omitnan')./nlength;
    std_N(i)=sqrt(sum(diff_test.^2,'omitnan')./nlength);
    r_N(i)=(sum(diff_test.*diff_ref,'omitnan')./nlength)./(std_N(i).*std_ref);

        %% remove small biomass for comparison
        % everything below 1.5 um has to be removed
        m_min=mass_function(1.5/2);
        m_max=simOutput.p.m(2+n).*sqrt(simOutput.p.m(2+n)/simOutput.p.m(1+n));
        for tt=1:size(simOutput.B,3)
            for dd=1:size(simOutput.B,1)
                [~,f] = calcBiomassRange(squeeze(simOutput.B(dd,1:n,tt)),simOutput.p.m(3:2+n),m_min,m_max);
                simOutput.B(dd,1:n,tt)=simOutput.B(dd,1:n,tt).*f;
            end
        end
        %% Calculate Biomasses (AC and HC)
        jsum=(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jDOC(thedepths,1:n,tstart:end)+simOutput.jFreal(thedepths,1:n,tstart:end));
        switch ACversion
            case 'ACv1'
                % autotroph carbon is a fraction of the total that is the ratio
                % of the rates
                AC=simOutput.B(thedepths,1:n,tstart:end).*...
                    simOutput.jLreal(thedepths,1:n,tstart:end)./jsum;
                HC=simOutput.B(thedepths,1:n,tstart:end).*...
                    (simOutput.jDOC(thedepths,1:n,tstart:end)+simOutput.jFreal(thedepths,1:n,tstart:end))./jsum;
                HC_osmo=simOutput.B(thedepths,1:n,tstart:end).*simOutput.jDOC(thedepths,1:n,tstart:end)./jsum;
                HC_phago=simOutput.B(thedepths,1:n,tstart:end).*simOutput.jFreal(thedepths,1:n,tstart:end)./jsum;
            case 'ACv2'
                jAH_ratio=simOutput.jLreal(thedepths,1:n,tstart:end)./(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jFreal(thedepths,1:n,tstart:end));
                A_cells=jAH_ratio;
                A_cells(jAH_ratio>AC_fraction)=1; %phototroph
                A_cells(A_cells~=1)=0;
                F_cells = double(~A_cells);

                AC=simOutput.B(thedepths,1:n,tstart:end).*A_cells.*(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jFreal(thedepths,1:n,tstart:end))./jsum...
                    +simOutput.B(thedepths,1:n,tstart:end).*F_cells.*(simOutput.jLreal(thedepths,1:n,tstart:end))./jsum;
                HC=simOutput.B(thedepths,1:n,tstart:end).*simOutput.jDOC(thedepths,1:n,tstart:end)./jsum...
                    +simOutput.B(thedepths,1:n,tstart:end).*F_cells.*simOutput.jFreal(thedepths,1:n,tstart:end)./jsum;
                HC_osmo=simOutput.B(thedepths,1:n,tstart:end).*simOutput.jDOC(thedepths,1:n,tstart:end)./jsum;
                HC_phago=simOutput.B(thedepths,1:n,tstart:end).*F_cells.*simOutput.jFreal(thedepths,1:n,tstart:end)./jsum;
            case 'ACv3'
                I0=ins(1,1:2:end);
                % 1pr day for 2 years
                I0=repmat(I0,[1 nyears]).*EinConv.*PARfrac;
                I0_5pct=I0.*0.05;
                deepcells=simOutput.L(thedepths,:)<=I0_5pct;
                deepcells=repmat(deepcells,[1 1 10]);
                deepcells=permute(deepcells,[1 3 2]);

                %definer ratios af AC til Phago (F) og osmo (O)
                jAF_ratio=simOutput.jLreal(thedepths,1:n,tstart:end)./(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jFreal(thedepths,1:n,tstart:end));
                jAO_ratio=simOutput.jLreal(thedepths,1:n,tstart:end)./(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jDOC(thedepths,1:n,tstart:end));
                %Find celler hvor A er større end ratioen:
                AF_cells=zeros(size(jAF_ratio));
                AF_cells(jAF_ratio>AC_fraction_F)=1; %phototroph i stedet for Phago
                AO_cells=zeros(size(jAO_ratio));
                AO_cells(jAO_ratio>AC_fraction_O & deepcells)=1; %phototroph i stedet for osmo i dybe celler
                %Og marker dem som ikke er
                HF_cells = double(~AF_cells);
                HO_cells = double(~AO_cells);
                % Find celler som både er HF og HO
                H_cells = HF_cells==HO_cells & HF_cells==1;
                % Find celler som både er AF og AO
                AFO_cells = AO_cells==AF_cells & AO_cells==1;
                % fjern dem fra de enkelte
                AF_cells(AFO_cells==1)=0;
                AO_cells(AFO_cells==1)=0;
                % Beregn AC som summen af de celler hvor alt er auto, hvor
                % noget er osmo, hvor noget er phago og hvor der er lidt af
                % alle:
                AC=simOutput.B(thedepths,1:n,tstart:end).*AFO_cells...
                    +simOutput.B(thedepths,1:n,tstart:end).*AF_cells.*(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jFreal(thedepths,1:n,tstart:end))./jsum...
                    +simOutput.B(thedepths,1:n,tstart:end).*AO_cells.*(simOutput.jLreal(thedepths,1:n,tstart:end)+simOutput.jDOC(thedepths,1:n,tstart:end))./jsum...
                    +simOutput.B(thedepths,1:n,tstart:end).*H_cells.*simOutput.jLreal(thedepths,1:n,tstart:end)./jsum;
                HC_phago=simOutput.B(thedepths,1:n,tstart:end).*HF_cells.*simOutput.jFreal(thedepths,1:n,tstart:end)./jsum;
                HC_osmo=simOutput.B(thedepths,1:n,tstart:end).*HO_cells.*simOutput.jDOC(thedepths,1:n,tstart:end)./jsum;

                HC=HC_osmo+HC_phago;
        end

        % sum biomass from size classes
        AC_tot=squeeze(sum(AC,2,'omitnan'));
        %         HC_tot=squeeze(sum(HC,2,'omitnan'));

        %% Data for figure Taylor and Landry 2A

        % Divide AC and HC_phago into pico, nano and micro plankton
        AC_pnm=nan(size(AC,1),size(AC,3),3);
        for tt=1:size(AC,3)
            for dd=1:size(AC,1)
                AC_pnm(dd,tt,:) = calcPicoNanoMicro(squeeze(AC(dd,:,tt)),simOutput.p.m(3:2+n));

            end
        end
        AC_pico=squeeze(AC_pnm(:,:,1));
        AC_nano=squeeze(AC_pnm(:,:,2));
        AC_micro=squeeze(AC_pnm(:,:,3));

        for ii=2:14
            AC_pico_mean_bin(ii)=mean(AC_pico(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
            AC_nano_mean_bin(ii)=mean(AC_nano(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
            AC_micro_mean_bin(ii)=mean(AC_micro(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');

        end
        %make negative data 0:
        AC_pico_mean_bin(AC_pico_mean_bin<0)=0;
        AC_nano_mean_bin(AC_nano_mean_bin<0)=0;
        AC_micro_mean_bin(AC_micro_mean_bin<0)=0;

        AC_pico_mean_bin_save(i,:)=AC_pico_mean_bin;
        AC_nano_mean_bin_save(i,:)=AC_nano_mean_bin;
        AC_micro_mean_bin_save(i,:)=AC_micro_mean_bin;

        TayLanFig1a.AC_pico_mean_bin(nrBins<2)=NaN;
        TayLanFig1a.AC_nano_mean_bin(nrBins<2)=NaN;
        TayLanFig1a.AC_micro_mean_bin(nrBins<2)=NaN;

        % calculate how many overlapping size groups there are
        Overlap_ACpico=((~isnan(AC_pico_mean_bin) & AC_pico_mean_bin~=0)==1 & (~isnan(TayLanFig1a.AC_pico_mean_bin) & TayLanFig1a.AC_pico_mean_bin~=0)==1);
        Overlap_ACnano=((~isnan(AC_nano_mean_bin) & AC_nano_mean_bin~=0)==1 & (~isnan(TayLanFig1a.AC_nano_mean_bin) & TayLanFig1a.AC_nano_mean_bin~=0)==1);
        Overlap_ACmicro=((~isnan(AC_micro_mean_bin) & AC_micro_mean_bin~=0)==1 & (~isnan(TayLanFig1a.AC_micro_mean_bin) & TayLanFig1a.AC_micro_mean_bin~=0)==1);

        nOverlap_ACpico(i)=sum(Overlap_ACpico);
        nOverlap_ACnano(i)=sum(Overlap_ACnano);
        nOverlap_ACmicro(i)=sum(Overlap_ACmicro);

        % convert to log scale
        AC_pico_mean_bin_log=log(AC_pico_mean_bin(1:14));
        AC_nano_mean_bin_log=log(AC_nano_mean_bin(1:14));
        AC_micro_mean_bin_log=log(AC_micro_mean_bin(1:14));
        TayLanFig1a.AC_pico_mean_bin_log=log(TayLanFig1a.AC_pico_mean_bin);
        TayLanFig1a.AC_nano_mean_bin_log=log(TayLanFig1a.AC_nano_mean_bin);
        TayLanFig1a.AC_micro_mean_bin_log=log(TayLanFig1a.AC_micro_mean_bin);
        %change 0(-inf) to nan
        AC_pico_mean_bin_log(isinf(AC_pico_mean_bin_log))=NaN;
        AC_nano_mean_bin_log(isinf(AC_nano_mean_bin_log))=NaN;
        AC_micro_mean_bin_log(isinf(AC_micro_mean_bin_log))=NaN;
        TayLanFig1a.AC_pico_mean_bin_log(isinf(TayLanFig1a.AC_pico_mean_bin_log))=NaN;
        TayLanFig1a.AC_nano_mean_bin_log(isinf(TayLanFig1a.AC_nano_mean_bin_log))=NaN;
        TayLanFig1a.AC_micro_mean_bin_log(isinf(TayLanFig1a.AC_micro_mean_bin_log))=NaN;

        %calculate result - mean
        diff_ACpico_test=AC_pico_mean_bin_log(Overlap_ACpico)-mean(AC_pico_mean_bin_log(Overlap_ACpico),'omitnan');
        diff_ACnano_test=AC_nano_mean_bin_log(Overlap_ACnano)-mean(AC_nano_mean_bin_log(Overlap_ACnano),'omitnan');
        diff_ACmicro_test=AC_micro_mean_bin_log(Overlap_ACmicro)-mean(AC_micro_mean_bin_log(Overlap_ACmicro),'omitnan');

        diff_ACpico_ref=TayLanFig1a.AC_pico_mean_bin_log(Overlap_ACpico)-mean(TayLanFig1a.AC_pico_mean_bin_log(Overlap_ACpico),'omitnan');
        diff_ACnano_ref=TayLanFig1a.AC_nano_mean_bin_log(Overlap_ACnano)-mean(TayLanFig1a.AC_nano_mean_bin_log(Overlap_ACnano),'omitnan');
        diff_ACmicro_ref=TayLanFig1a.AC_micro_mean_bin_log(Overlap_ACmicro)-mean(TayLanFig1a.AC_micro_mean_bin_log(Overlap_ACmicro),'omitnan');
        %calculate standard deviation
        std_ACpico_test_log(i)=sqrt(sum(diff_ACpico_test.^2,'omitnan')./nOverlap_ACpico(i));
        std_ACnano_test_log(i)=sqrt(sum(diff_ACnano_test.^2,'omitnan')./nOverlap_ACnano(i));
        std_ACmicro_test_log(i)=sqrt(sum(diff_ACmicro_test.^2,'omitnan')./nOverlap_ACmicro(i));

        std_ACpico_ref_log(i)=sqrt(sum(diff_ACpico_ref.^2,'omitnan')./nOverlap_ACpico(i));
        std_ACnano_ref_log(i)=sqrt(sum(diff_ACnano_ref.^2,'omitnan')./nOverlap_ACnano(i));
        std_ACmicro_ref_log(i)=sqrt(sum(diff_ACmicro_ref.^2,'omitnan')./nOverlap_ACmicro(i));

        %calculate correlation coefficient
        r_ACpico(i)=(sum(diff_ACpico_test.*diff_ACpico_ref,'omitnan')./nOverlap_ACpico(i))./(std_ACpico_test_log(i).*std_ACpico_ref_log(i));
        r_ACnano(i)=(sum(diff_ACnano_test.*diff_ACnano_ref,'omitnan')./nOverlap_ACnano(i))./(std_ACnano_test_log(i).*std_ACnano_ref_log(i));
        r_ACmicro(i)=(sum(diff_ACmicro_test.*diff_ACmicro_ref,'omitnan')./nOverlap_ACmicro(i))./(std_ACmicro_test_log(i).*std_ACmicro_ref_log(i));

        %calculate cRMSdiff:
        cRMSd_ACpico(i)=sum((diff_ACpico_test-diff_ACpico_ref).^2,'omitnan')./nOverlap_ACpico(i);
        cRMSd_ACnano(i)=sum((diff_ACnano_test-diff_ACnano_ref).^2,'omitnan')./nOverlap_ACnano(i);
        cRMSd_ACmicro(i)=sum((diff_ACmicro_test-diff_ACmicro_ref).^2,'omitnan')./nOverlap_ACmicro(i);

        %calculate RMSdiff:
        RMSd_ACpico(i)=sum((AC_pico_mean_bin_log(Overlap_ACpico)-TayLanFig1a.AC_pico_mean_bin_log(Overlap_ACpico)).^2,'omitnan')./nOverlap_ACpico(i);
        RMSd_ACnano(i)=sum((AC_nano_mean_bin_log(Overlap_ACnano)-TayLanFig1a.AC_nano_mean_bin_log(Overlap_ACnano)).^2,'omitnan')./nOverlap_ACnano(i);
        RMSd_ACmicro(i)=sum((AC_micro_mean_bin_log(Overlap_ACmicro)-TayLanFig1a.AC_micro_mean_bin_log(Overlap_ACmicro)).^2,'omitnan')./nOverlap_ACmicro(i);


        if Plot_ACpnm==true
            %% Plot figure Taylor and Landry 2A
            figure('Name','Taylor & Landry 2018')
            AC_x=(AC_bin(1:14-1)+AC_bin(2:14))./2;
            %         AC_x=exp(log(AC_bin(1:14-1))+0.5*(log(AC_bin(3))-log(AC_bin(2))));
            subplot(1,3,1)
            loglog(AC_x,AC_pico_mean_bin(2:14),'-or'); hold on;
            loglog((TayLanFig1a.AC_bin(1:14-1)+TayLanFig1a.AC_bin(2:14))./2,TayLanFig1a.AC_pico_mean_bin(2:14),'-ok'); hold on;
            set(gca,'Xlim',[2 1000])%,'Ylim',[0.1 1000]);
            subtitle AC_{Pico}
            xlabel('AC_{tot} (\mug C/l)');ylabel('A_{pico}biomass (\mug C/l)');legend('Model','Data')
            ylim([10^-1 10^3])
            subplot(1,3,2)
            loglog(AC_x,AC_nano_mean_bin(2:14),':^r'); hold on;
            loglog((TayLanFig1a.AC_bin(1:14-1)+TayLanFig1a.AC_bin(2:14))./2,TayLanFig1a.AC_nano_mean_bin(2:14),':^k'); hold on;
            set(gca,'Xlim',[2 1000])%,'Ylim',[0.1 1000]);
            subtitle AC_{Nano}
            xlabel('AC_{tot} (\mug C/l)');ylabel('A_{nano}biomass (\mug C/l)');legend('Model','Data')
            ylim([10^-1 10^3])
            subplot(1,3,3)
            loglog(AC_x,AC_micro_mean_bin(2:14),'--sqr'); hold on;
            loglog((TayLanFig1a.AC_bin(1:14-1)+TayLanFig1a.AC_bin(2:14))./2,TayLanFig1a.AC_micro_mean_bin(2:14),'--sqk'); hold on;
            set(gca,'Xlim',[2 1000])%,'Ylim',[0.1 1000]);
            subtitle AC_{Micro}
            xlabel('AC_{tot} (\mug C/l)');ylabel('A_{micro}biomass (\mug C/l)');legend('Model','Data')
            ylim([10^-1 10^3])
        end

        %% Data for figure 7
        % divide phagotrophs into pico, nano, and micro
        HC_phago_pnm=nan(size(HC_phago,1),size(HC_phago,3),3);
        for tt=1:size(AC,3)
            for dd=1:size(AC,1)
                HC_phago_pnm(dd,tt,:) = calcPicoNanoMicro(squeeze(HC_phago(dd,:,tt)),simOutput.p.m(3:2+n));
            end
        end
        HC_phago_pico=squeeze(HC_phago_pnm(:,:,1));
        HC_phago_nano=squeeze(HC_phago_pnm(:,:,2));
        %         HC_phago_micro=squeeze(HC_phago_pnm(:,:,3));
        % Define Hbac as osmo + pico phago
        HC_Hbac_tot=squeeze(sum(HC_osmo,2,'omitnan'))+HC_phago_pico;
        % Define Hnano as phago nano
        HC_Hnano_tot=HC_phago_nano;

        % Find the mean biomass within each AC bin
        AC_tot_mean_bin=nan(size(AC_bin));
        HC_Hbac_tot_mean_bin=nan(size(AC_bin));
        HC_Hnano_tot_mean_bin=nan(size(AC_bin));

        for ii=2:14
            AC_tot_mean_bin(ii)=mean(AC_tot(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
            HC_Hbac_tot_mean_bin(ii)=mean(HC_Hbac_tot(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
            HC_Hnano_tot_mean_bin(ii)=mean(HC_Hnano_tot(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
        end
        %make negative data 0:
        HC_Hbac_tot_mean_bin(HC_Hbac_tot_mean_bin<0)=0;
        HC_Hnano_tot_mean_bin(HC_Hnano_tot_mean_bin<0)=0;

        HC_Hbac_tot_mean_bin_save(i,:)=HC_Hbac_tot_mean_bin;
        HC_Hnano_tot_mean_bin_save(i,:)=HC_Hnano_tot_mean_bin;

        % % Do statitics for HC_Hbak vs AC:
        TayLanFig7.HC_hbac_mean_bin(nrBins<2)=NaN;


        %find overlap:
        Overlap_HC_Hbac=((~isnan(HC_Hbac_tot_mean_bin) & HC_Hbac_tot_mean_bin~=0)==1 & (~isnan(TayLanFig7.HC_hbac_mean_bin) & TayLanFig7.HC_hbac_mean_bin~=0)==1);
        nOverlap_HC_Hbac(i)=sum(Overlap_HC_Hbac);

        % convert to log scale
        HC_Hbac_tot_mean_bin_log=log(HC_Hbac_tot_mean_bin(1:14));
        TayLanFig7.HC_hbac_mean_bin_log=log(TayLanFig7.HC_hbac_mean_bin);

        %change 0(-inf) to nan
        HC_Hbac_tot_mean_bin_log(isinf(HC_Hbac_tot_mean_bin_log))=NaN;
        TayLanFig7.HC_hbac_mean_bin_log(isinf(TayLanFig7.HC_hbac_mean_bin_log))=NaN;

        %calculate result - mean
        diff_HC_Hbac_test=HC_Hbac_tot_mean_bin_log(Overlap_HC_Hbac)-mean(HC_Hbac_tot_mean_bin_log(Overlap_HC_Hbac),'omitnan');
        diff_HC_Hbac_ref=TayLanFig7.HC_hbac_mean_bin_log(Overlap_HC_Hbac)-mean(TayLanFig7.HC_hbac_mean_bin_log(Overlap_HC_Hbac),'omitnan');
        %calculate standard deviation
        std_HC_Hbac_test_log(i)=sqrt(sum(diff_HC_Hbac_test.^2,'omitnan')./nOverlap_HC_Hbac(i));
        std_HC_Hbac_ref_log(i)=sqrt(sum(diff_HC_Hbac_ref.^2,'omitnan')./nOverlap_HC_Hbac(i));
        %calculate correlation coefficient
        r_HC_Hbac(i)=(sum(diff_HC_Hbac_test.*diff_HC_Hbac_ref,'omitnan')./nOverlap_HC_Hbac(i))./(std_HC_Hbac_test_log(i).*std_HC_Hbac_ref_log(i));
        %calculate cRMSdiff:
        cRMSd_HC_Hbac(i)=sum((diff_HC_Hbac_test-diff_HC_Hbac_ref).^2,'omitnan')./nOverlap_HC_Hbac(i);
        %calculate RMSdiff:
        RMSd_HC_Hbac(i)=sum((HC_Hbac_tot_mean_bin_log(Overlap_HC_Hbac)-TayLanFig7.HC_hbac_mean_bin_log(Overlap_HC_Hbac)).^2,'omitnan')./nOverlap_HC_Hbac(i);

        % % Do statitics for HC_nano vs AC:
        TayLanFig7.HC_Hnano_mean_bin(nrBins<2)=NaN;

        %find overlap:
        Overlap_HC_Hnano=((~isnan(HC_Hnano_tot_mean_bin) & HC_Hnano_tot_mean_bin~=0)==1 & (~isnan(TayLanFig7.HC_Hnano_mean_bin) & TayLanFig7.HC_Hnano_mean_bin~=0)==1);
        nOverlap_HC_Hnano(i)=sum(Overlap_HC_Hnano);

        % convert to log scale
        HC_Hnano_tot_mean_bin_log=log(HC_Hnano_tot_mean_bin(1:14));
        TayLanFig7.HC_Hnano_mean_bin_log=log(TayLanFig7.HC_Hnano_mean_bin);

        %change 0(-inf) to nan
        HC_Hnano_tot_mean_bin_log(isinf(HC_Hnano_tot_mean_bin_log))=NaN;
        TayLanFig7.HC_Hnano_mean_bin_log(isinf(TayLanFig7.HC_Hnano_mean_bin_log))=NaN;

        %calculate result - mean
        diff_HC_Hnano_test=HC_Hnano_tot_mean_bin_log(Overlap_HC_Hnano)-mean(HC_Hnano_tot_mean_bin_log(Overlap_HC_Hnano),'omitnan');
        diff_HC_Hnano_ref=TayLanFig7.HC_Hnano_mean_bin_log(Overlap_HC_Hnano)-mean(TayLanFig7.HC_Hnano_mean_bin_log(Overlap_HC_Hnano),'omitnan');
        %calculate standard deviation
        std_HC_Hnano_test_log(i)=sqrt(sum(diff_HC_Hnano_test.^2,'omitnan')./nOverlap_HC_Hnano(i));
        std_HC_Hnano_ref_log(i)=sqrt(sum(diff_HC_Hnano_ref.^2,'omitnan')./nOverlap_HC_Hnano(i));
        %calculate correlation coefficient
        r_HC_Hnano(i)=(sum(diff_HC_Hnano_test.*diff_HC_Hnano_ref,'omitnan')./nOverlap_HC_Hnano(i))./(std_HC_Hnano_test_log(i).*std_HC_Hnano_ref_log(i));
        %calculate cRMSdiff:
        cRMSd_HC_Hnano(i)=sum((diff_HC_Hnano_test-diff_HC_Hnano_ref).^2,'omitnan')./nOverlap_HC_Hnano(i);
        %calculate RMSdiff:
        RMSd_HC_Hnano(i)=sum((HC_Hnano_tot_mean_bin_log(Overlap_HC_Hnano)-TayLanFig7.HC_Hnano_mean_bin_log(Overlap_HC_Hnano)).^2,'omitnan')./nOverlap_HC_Hnano(i);

        if Plot_HbacHnano==true
            %% Plot figure fig 7
            figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(2,1,1)
            h=loglog(AC_tot,HC_Hbac_tot,'.r');
            hold on
            h1=loglog((AC_bin(1:14-1)+AC_bin(2:14))./2, HC_Hbac_tot_mean_bin(2:14),'-or','LineWidth',2);
            h2=loglog(TayLanFig7.AC_tot,TayLanFig7.HC_Hbac,'.k');
            h3=loglog((TayLanFig7.AC_bin(1:14-1)+TayLanFig7.AC_bin(2:14))./2, TayLanFig7.HC_hbac_mean_bin(2:14),'-ok','LineWidth',2);
            set(gca,'xscale','log','yscale','log')
            xlabel('AC_{tot} (\mugC/l)')
            ylabel('Hbac-biomass (\mugC/l)')
            h=h(1);h(2)=h1(1);h(3)=h2(1);h(4)=h3(1);
            legend(h([1 2 3 4]),{'model result','model mean','data','data mean'})
            mystr={['STD = ', num2str(std_HC_Hbac_test_log(i)./std_HC_Hbac_ref_log(i))],...
                ['COR = ', num2str(r_HC_Hbac(i))],...
                ['cRMS = ', num2str(sqrt(cRMSd_HC_Hbac(i))./std_HC_Hbac_ref_log(i))],...
                ['RMS = ', num2str(sqrt(RMSd_HC_Hbac(i)))],...
                ['overlapping bins = ',num2str(nOverlap_HC_Hbac(i))]};
            annotation('textbox', [0.14, 0.68, 0.6, 0], 'string',mystr,'LineStyle','none');
            ylim([10^-2 10^3])


            subplot(2,1,2)
            h=loglog(AC_tot,HC_Hnano_tot,'.r');
            hold on
            h1=loglog((AC_bin(1:14-1)+AC_bin(2:14))./2, HC_Hnano_tot_mean_bin(2:14),'-or','LineWidth',2);
            h2=loglog(TayLanFig7.AC_tot,TayLanFig7.HC_Hnano,'.k');
            h3=loglog((TayLanFig7.AC_bin(1:14-1)+TayLanFig7.AC_bin(2:14))./2, TayLanFig7.HC_Hnano_mean_bin(2:14),'-ok','LineWidth',2);
            set(gca,'xscale','log','yscale','log')
            xlabel('AC_{tot} (\mugC/l)')
            ylabel('Hnano-biomass (\mugC/l)')
            h=h(1);h(2)=h1(1);h(3)=h2(1);h(4)=h3(1);
            legend(h([1 2 3 4]),{'model result','model mean','data','data mean'})
            mystr={['STD = ', num2str(std_HC_Hnano_test_log(i)./std_HC_Hnano_ref_log(i))],...
                ['COR = ', num2str(r_HC_Hnano(i))],...
                ['cRMS = ', num2str(sqrt(cRMSd_HC_Hnano(i))./std_HC_Hnano_ref_log(i))],...
                ['RMS = ', num2str(sqrt(RMSd_HC_Hnano(i)))],...
                ['overlapping bins = ',num2str(nOverlap_HC_Hnano(i))]};
            annotation('textbox', [0.14, 0.2, 0.6, 0], 'string',mystr,'LineStyle','none');
            ylim([10^-2 10^3])


        end
        %% calculate sum pico and sum nano
        TotB=simOutput.B(thedepths,1:n,tstart:end);


        TOT_pnm=nan(size(TotB,1),size(TotB,3),3);
        for tt=1:size(TotB,3)
            for dd=1:size(TotB,1)
                TOT_pnm(dd,tt,:) = calcPicoNanoMicro(squeeze(TotB(dd,:,tt)),simOutput.p.m(3:2+n));

            end
        end

        TOT_pico=squeeze(TOT_pnm(:,:,1));
        TOT_nano=squeeze(TOT_pnm(:,:,2));

        for ii=2:14
            TOT_pico_mean_bin(ii)=mean(TOT_pico(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
            TOT_nano_mean_bin(ii)=mean(TOT_nano(AC_tot>AC_bin(ii-1)&AC_tot<=AC_bin(ii)),'omitnan');
        end

        %make negative data 0:
        TOT_pico_mean_bin(TOT_pico_mean_bin<0)=0;
        TOT_nano_mean_bin(TOT_nano_mean_bin<0)=0;

        TOT_pico_mean_bin_save(i,:)=TOT_pico_mean_bin;
        TOT_nano_mean_bin_save(i,:)=TOT_nano_mean_bin;

        totPN.Tot_pico_mean_bin(nrBins<2)=NaN;
        totPN.Tot_nano_mean_bin(nrBins<2)=NaN;

        % calculate how many overlapping size groups there are
        Overlap_TOTpico=((~isnan(TOT_pico_mean_bin) & TOT_pico_mean_bin~=0)==1 & (~isnan(totPN.Tot_pico_mean_bin) & totPN.Tot_pico_mean_bin~=0)==1);
        Overlap_TOTnano=((~isnan(TOT_nano_mean_bin) & TOT_nano_mean_bin~=0)==1 & (~isnan(totPN.Tot_nano_mean_bin) & totPN.Tot_nano_mean_bin~=0)==1);

        nOverlap_TOTpico(i)=sum(Overlap_TOTpico);
        nOverlap_TOTnano(i)=sum(Overlap_TOTnano);

        % convert to log scale
        TOT_pico_mean_bin_log=log(TOT_pico_mean_bin(1:14));
        TOT_nano_mean_bin_log=log(TOT_nano_mean_bin(1:14));

        totPN.Tot_pico_mean_bin_log=log(totPN.Tot_pico_mean_bin);
        totPN.Tot_nano_mean_bin_log=log(totPN.Tot_nano_mean_bin);

        %change 0(-inf) to nan
        TOT_pico_mean_bin_log(isinf(TOT_pico_mean_bin_log))=NaN;
        TOT_nano_mean_bin_log(isinf(TOT_nano_mean_bin_log))=NaN;

        totPN.Tot_pico_mean_bin_log(isinf(totPN.Tot_pico_mean_bin_log))=NaN;
        totPN.Tot_nano_mean_bin_log(isinf(totPN.Tot_nano_mean_bin_log))=NaN;

        %calculate result - mean
        diff_TOTpico_test=TOT_pico_mean_bin_log(Overlap_TOTpico)-mean(TOT_pico_mean_bin_log(Overlap_TOTpico),'omitnan');
        diff_TOTnano_test=TOT_nano_mean_bin_log(Overlap_TOTnano)-mean(TOT_nano_mean_bin_log(Overlap_TOTnano),'omitnan');

        diff_TOTpico_ref=totPN.Tot_pico_mean_bin_log(Overlap_TOTpico)-mean(totPN.Tot_pico_mean_bin_log(Overlap_TOTpico),'omitnan');
        diff_TOTnano_ref=totPN.Tot_nano_mean_bin_log(Overlap_TOTnano)-mean(totPN.Tot_nano_mean_bin_log(Overlap_TOTnano),'omitnan');

        %calculate standard deviation
        std_TOTpico_test_log(i)=sqrt(sum(diff_TOTpico_test.^2,'omitnan')./nOverlap_TOTpico(i));
        std_TOTnano_test_log(i)=sqrt(sum(diff_TOTnano_test.^2,'omitnan')./nOverlap_TOTnano(i));

        std_TOTpico_ref_log(i)=sqrt(sum(diff_TOTpico_ref.^2,'omitnan')./nOverlap_TOTpico(i));
        std_TOTnano_ref_log(i)=sqrt(sum(diff_TOTnano_ref.^2,'omitnan')./nOverlap_TOTnano(i));

        %calculate correlation coefficient
        r_TOTpico(i)=(sum(diff_TOTpico_test.*diff_TOTpico_ref,'omitnan')./nOverlap_TOTpico(i))./(std_TOTpico_test_log(i).*std_TOTpico_ref_log(i));
        r_TOTnano(i)=(sum(diff_TOTnano_test.*diff_TOTnano_ref,'omitnan')./nOverlap_TOTnano(i))./(std_TOTnano_test_log(i).*std_TOTnano_ref_log(i));

        %calculate cRMSdiff:
        cRMSd_TOTpico(i)=sum((diff_TOTpico_test-diff_TOTpico_ref).^2,'omitnan')./nOverlap_TOTpico(i);
        cRMSd_TOTnano(i)=sum((diff_TOTnano_test-diff_TOTnano_ref).^2,'omitnan')./nOverlap_TOTnano(i);

        %calculate RMSdiff:
        RMSd_TOTpico(i)=sum((TOT_pico_mean_bin_log(Overlap_TOTpico)-totPN.Tot_pico_mean_bin_log(Overlap_TOTpico)).^2,'omitnan')./nOverlap_TOTpico(i);
        RMSd_TOTnano(i)=sum((TOT_nano_mean_bin_log(Overlap_TOTnano)-totPN.Tot_nano_mean_bin_log(Overlap_TOTnano)).^2,'omitnan')./nOverlap_TOTnano(i);


        if Plot_TOTpn==true
            %% Plot figure Taylor and Landry 2A
            figure('Name','TOTpn')
            AC_x=(AC_bin(1:14-1)+AC_bin(2:14))./2;
            subplot(1,2,1)
            loglog(AC_x,TOT_pico_mean_bin(2:14),'-or'); hold on;
            loglog((TayLanFig1a.AC_bin(1:14-1)+TayLanFig1a.AC_bin(2:14))./2,totPN.Tot_pico_mean_bin(2:14),'-ok'); hold on;
            set(gca,'Xlim',[2 1000])%,'Ylim',[0.1 1000]);
            subtitle TOT_{Pico}
            xlabel('AC_{tot} (\mug C/l)');ylabel('TOT_{pico}biomass (\mug C/l)');legend('Model','Data')
            ylim([10^-1 10^3])
            subplot(1,2,2)
            loglog(AC_x,TOT_nano_mean_bin(2:14),':^r'); hold on;
            loglog((TayLanFig1a.AC_bin(1:14-1)+TayLanFig1a.AC_bin(2:14))./2,totPN.Tot_nano_mean_bin(2:14),':^k'); hold on;
            set(gca,'Xlim',[2 1000])%,'Ylim',[0.1 1000]);
            subtitle TOT_{Nano}
            xlabel('AC_{tot} (\mug C/l)');ylabel('TOT_{nano}biomass (\mug C/l)');legend('Model','Data')
            ylim([10^-1 10^3])
        end

        dz=repmat(simOutput.dznom(thedepths), [1 n]);
        meanB=squeeze(mean(TotB,3,'omitnan'));
        sumB(i,:)=squeeze(sum((meanB.*dz),1));



    end
    if saveResult==true
        %% Make table with statistics
        stattable=table(relevantIdx,thefilenames(relevantIdx),nOverlap_HC_Hbac,nOverlap_HC_Hnano,r_HC_Hbac,r_HC_Hnano,std_HC_Hbac_test_log,std_HC_Hbac_ref_log,std_HC_Hnano_test_log,std_HC_Hnano_ref_log,cRMSd_HC_Hbac,RMSd_HC_Hbac,cRMSd_HC_Hnano,RMSd_HC_Hnano,...
            nOverlap_ACpico,nOverlap_ACnano,nOverlap_ACmicro,...
            r_ACpico,r_ACnano,r_ACmicro,std_ACpico_test_log,std_ACnano_test_log,std_ACmicro_test_log,std_ACpico_ref_log,std_ACnano_ref_log,std_ACmicro_ref_log,cRMSd_ACpico,cRMSd_ACnano,cRMSd_ACmicro,RMSd_ACpico,RMSd_ACnano,RMSd_ACmicro,...
            nOverlap_TOTpico,nOverlap_TOTnano,r_TOTpico,r_TOTnano,std_TOTpico_test_log,std_TOTnano_test_log,std_TOTpico_ref_log,std_TOTnano_ref_log,cRMSd_TOTpico,cRMSd_TOTnano,RMSd_TOTpico,RMSd_TOTnano,...
            'VariableNames',["filenr","filename","nOverlapHC_Hbac","nOverlapHC_Hnano","r_HC_Hbac","r_HC_Hnano","std_HC_Hbac","std_HC_Hbac_ref","std_HC_Hnano","std_HC_Hnano_ref","cRMSd_HC_Hbac","RMSd_HC_Hbac","cRMSd_HC_Hnano","RMSd_HC_Hnano",...
            "nOverlapACpico","nOverlapACnano","nOverlapACmicro",...
            "r_ACpico","r_ACnano","r_ACmicro","std_ACpico","std_ACnano","std_ACmicro","std_ACpico_ref","std_ACnano_ref","std_ACmicro_ref","cRMSd_ACpico","cRMSd_ACnano","cRMSd_ACmicro","RMSd_ACpico","RMSd_ACnano","RMSd_ACmicro",...
            "nOverlapTOTpico","nOverlapTOTnano","r_TOTpico","r_TOTnano","std_TOTpico","std_TOTnano","std_TOTpico_ref","std_TOTnano_ref","cRMSd_TOTpico","cRMSd_TOTnano","RMSd_TOTpico","RMSd_TOTnano"]);
        %% Save statistics and randParam
        %     tt = datestr(clock,'YYYY_mm_dd_HH_MM');
        %             tt=extractAfter(thissetup,'random');
        tt=thissetup;
        save(['stattable_',mysite,'_nudge_',num2str(wNudge),'_',ACversion,'_dmax',num2str(maxdepth),'_time_',tt,'.mat'],'stattable','randParam','mass','thisfile_name','ACversion','AC_fraction','AC_fraction_F','AC_fraction_O','maxdepth','tstart','thissetup','thissetup1','AC_pico_mean_bin_save','AC_nano_mean_bin_save','AC_micro_mean_bin_save','HC_Hbac_tot_mean_bin_save','HC_Hnano_tot_mean_bin_save','sumB','TOT_pico_mean_bin_save','TOT_nano_mean_bin_save','N','rmsd_N','cRMSd_N','std_N','r_N','Nmat')

    end

    clearvars -except mysites kkk ntfile files thefilenames site nudging oldfilenames a thissetup thissetup1


catch ME
    %% Make table with statistics
%     stattable=table(relevantIdx,thefilenames(relevantIdx),nOverlap_HC_Hbac,nOverlap_HC_Hnano,r_HC_Hbac,r_HC_Hnano,std_HC_Hbac_test_log,std_HC_Hbac_ref_log,std_HC_Hnano_test_log,std_HC_Hnano_ref_log,cRMSd_HC_Hbac,RMSd_HC_Hbac,cRMSd_HC_Hnano,RMSd_HC_Hnano,...
%         nOverlap_ACpico,nOverlap_ACnano,nOverlap_ACmicro,...
%         r_ACpico,r_ACnano,r_ACmicro,std_ACpico_test_log,std_ACnano_test_log,std_ACmicro_test_log,std_ACpico_ref_log,std_ACnano_ref_log,std_ACmicro_ref_log,cRMSd_ACpico,cRMSd_ACnano,cRMSd_ACmicro,RMSd_ACpico,RMSd_ACnano,RMSd_ACmicro,...
%         nOverlap_TOTpico,nOverlap_TOTnano,r_TOTpico,r_TOTnano,std_TOTpico_test_log,std_TOTnano_test_log,std_TOTpico_ref_log,std_TOTnano_ref_log,cRMSd_TOTpico,cRMSd_TOTnano,RMSd_TOTpico,RMSd_TOTnano,...
%         'VariableNames',["filenr","filename","nOverlapHC_Hbac","nOverlapHC_Hnano","r_HC_Hbac","r_HC_Hnano","std_HC_Hbac","std_HC_Hbac_ref","std_HC_Hnano","std_HC_Hnano_ref","cRMSd_HC_Hbac","RMSd_HC_Hbac","cRMSd_HC_Hnano","RMSd_HC_Hnano",...
%         "nOverlapACpico","nOverlapACnano","nOverlapACmicro",...
%         "r_ACpico","r_ACnano","r_ACmicro","std_ACpico","std_ACnano","std_ACmicro","std_ACpico_ref","std_ACnano_ref","std_ACmicro_ref","cRMSd_ACpico","cRMSd_ACnano","cRMSd_ACmicro","RMSd_ACpico","RMSd_ACnano","RMSd_ACmicro",...
%         "nOverlapTOTpico","nOverlapTOTnano","r_TOTpico","r_TOTnano","std_TOTpico","std_TOTnano","std_TOTpico_ref","std_TOTnano_ref","cRMSd_TOTpico","cRMSd_TOTnano","RMSd_TOTpico","RMSd_TOTnano"]);
%     %% Save statistics and randParam
%     %     tt = datestr(clock,'YYYY_mm_dd_HH_MM');
%     %     tt=extractAfter(thissetup,'random');
%     tt=thissetup;
%     iend=length(thisfile_name)+1;
%     stattable(iend:end,:)=[];
%     randParam(iend:end,:)=[];
%     AC_pico_mean_bin_save(iend:end,:)=[];
%     AC_nano_mean_bin_save(iend:end,:)=[];
%     AC_micro_mean_bin_save(iend:end,:)=[];
%     HC_Hbac_tot_mean_bin_save(iend:end,:)=[];
%     HC_Hnano_tot_mean_bin_save(iend:end,:)=[];
%     sumB(iend:end,:)=[];
%     TOT_pico_mean_bin_save(iend:end,:)=[];
%     TOT_nano_mean_bin_save(iend:end,:)=[];
% 
%     save(['stattable_part',num2str(a+1),'_V1_',mysite,'_nudge_',num2str(wNudge),'_',ACversion,'_dmax',num2str(maxdepth),'_time_',tt,'.mat'],'stattable','randParam','mass','thisfile_name','ACversion','AC_fraction','AC_fraction_F','AC_fraction_O','maxdepth','tstart','thissetup','thissetup1','AC_pico_mean_bin_save','AC_nano_mean_bin_save','AC_micro_mean_bin_save','HC_Hbac_tot_mean_bin_save','HC_Hnano_tot_mean_bin_save','sumB','TOT_pico_mean_bin_save','TOT_nano_mean_bin_save')
%     disp(['Code error at i = ',num2str(i)])
    theruns(1:i)=[];
    error_line = ME.stack(end).line;
    save(thisout)
%     save mycode_errored
end



