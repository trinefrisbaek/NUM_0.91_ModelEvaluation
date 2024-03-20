allparameters={'r^*_d','y','r^*_l','aF','c_{passive}','\alpha_{Max}','\alpha_R','mortHTL',...
    '\rho_{C:N}','fracHTL:N','\epsilon_L','\alpha_N','\epsilon_F','cF','\beta','\sigma',...
    'remin2','reminF','rho','mort2','rem_{POM}','vel_1','vel_2','mHTL'};

allparameters_Label={'r^*_D','\alpha_L','r^*_L','aF','c_{passive}','\alpha_{max}','\alpha_R','\mu_{HTL}',...
    '\rho_{C:N}','remin_{HTL}','\epsilon_L','\alpha_D','\epsilon_F','cF','\beta','\sigma',...
    'remin_2','remin_F','\rho','\mu_{v0}','a','v_1','v_2','mHTL'};



%find the names of the .mat files
thematfiles=dir('Sobol_*.mat');
ifilenames={thematfiles.name}';
filenumbers=str2double(extractBetween(ifilenames,'sobolnr_','_time'));

aretheythere=zeros(20000,49);
rmsd_tpico=nan(20000,49);
rmsd_tnano=nan(20000,49);
rmsd_apico=nan(20000,49);
rmsd_anano=nan(20000,49);
rmsd_amicro=nan(20000,49);

cor_tpico=nan(20000,49);
cor_tnano=nan(20000,49);
cor_apico=nan(20000,49);
cor_anano=nan(20000,49);
cor_amicro=nan(20000,49);

for i=0:2
    %decide which file this is
    this_file_idx=find(filenumbers==i);
    load('thisdata_init.mat');
    thisdata.stattable=stattable;
    for j=1:length(this_file_idx)
        thisdata_tmp=load(ifilenames{this_file_idx(j)},'stattable');
        if strcmp(ifilenames{this_file_idx(j)},'XXSobol_stattable_CCE_nudge_0_ACv3_dmax200_sobolnr_1_time_final_Sobols.mat')
            thisdata2=load('Sobol2_stattable_CCE_nudge_0_ACv3_dmax200_sobolnr_1_time_final_Sobols.mat');
            thisdata_tmp.stattable(72040:85946,:)=thisdata2.stattable(72040:85946,:);
        end
        if contains(thisdata_tmp.stattable.filename(1),'2024')
            randnr=str2double(extractBetween(thisdata_tmp.stattable.filename,'thisrand','_nrRanditer'));
            randnr_str_new=string(randnr+10000);
            thisdata_tmp.stattable.filename = replaceBetween(thisdata_tmp.stattable.filename,'thisrand','_nrRanditer',randnr_str_new);
        end
        thisdata.stattable=[thisdata.stattable;thisdata_tmp.stattable];
    end
    clear stattable randnr thisdata_tmp

    %remove first lines
    thisdata.stattable(1,:) = [];

    %     find the numbers of the files:

    if i==0
%         delete dublicates:
        dep_single=[16203 10049 14928 10008 10006 10004 10001];
        del_10593=[14919	14920	14921	14922];
        del_10594=[14923	14924];
        del_10595=[14925	14926];
        del_11595=[15077 15078];
        del_12598=[15239	15240	15241];
        del_13604=[15397	15398	15399];
        del_14602=[15559	15560	15561];
        del_15605=[15721	15722	15723];
        del_16593=[15880	15881	15882];
        del_19601=[16367    16368];

        del_all=[dep_single';del_10593';del_10594';del_10595';del_11595';del_12598';del_13604';del_14602';del_15605';del_16593';del_19601'];
        thisdata.stattable(del_all,:) = [];


        datastruct_0=thisdata.stattable;
        randnr_0=str2double(extractBetween(thisdata.stattable.filename,'thisrand','_nrRanditer'));
        for ii=1:length(randnr_0)
            aretheythere(randnr_0(ii),1)=aretheythere(randnr_0(ii),1)+1;
        end

        for ii=1:length(randnr_0)
            % tag data fra linje nr ii og put den ind i den der er sorteret
            % på randnr
            rmsd_tpico(randnr_0(ii),1)=datastruct_0.RMSd_TOTpico(ii);
            rmsd_tnano(randnr_0(ii),1)=datastruct_0.RMSd_TOTnano(ii);

            rmsd_apico(randnr_0(ii),1)=datastruct_0.RMSd_ACpico(ii);
            rmsd_anano(randnr_0(ii),1)=datastruct_0.RMSd_ACnano(ii);
            rmsd_amicro(randnr_0(ii),1)=datastruct_0.RMSd_ACmicro(ii);

            cor_tpico(randnr_0(ii),1)=datastruct_0.r_TOTpico(ii);
            cor_tnano(randnr_0(ii),1)=datastruct_0.r_TOTnano(ii);

            cor_apico(randnr_0(ii),1)=datastruct_0.r_ACpico(ii);
            cor_anano(randnr_0(ii),1)=datastruct_0.r_ACnano(ii);
            cor_amicro(randnr_0(ii),1)=datastruct_0.r_ACmicro(ii);
        end

    end
    if i==1
        del_10593=[376451 376452	376453	376454	376455	376456	376457	376458	376459	376460	376461	376462	376463	376464	376465	376466	376467	376468	376469	376470	376471	376472	376473	376474	376475	376476	376477	376478	376479	376480	376481	376482	376483	376484	376485	376486	376487	376488	376489	376490	376491	376492	376493	376494	376495	376496	376497	376498	376499	376500	376501	376502	376503	376504];
        del_10594=[376505	376506	376507	376508	376509	376510	376511	376512	376513	376514	376515	376516	376517	376518	376519	376520	376521	376522	376523	376524	376525	376526	376527	376528	376529	376530	376531	376532	376533	376534	376535	376536	376537	376538	376539	376540	376541	376542	376543	376544	376545	376546	376547	376548	376549	376550	376551	376552];
        del_10595=[376553	376554	376555	376556	376557	376558	376559	376560	376561	376562	376563	376564	376565	376566	376567	376568	376569	376570	376571	376572	376573	376574	376575	376576	376577	376578	376579	376580	376581	376582	376583	376584	376585	376586	376587	376588	376589	376590	376591	376592	376593	376594	376595	376596	376597	376598	376599	376600];
        del_11595=[378142 378140  378137 378138];
        del_12598=[379929	379930	379931	379932	379933	379934	379935	379936	379937	379938	379939	379940	379941	379942	379943	379944	379945	379946	379947	379948	379949	379950	379951	379952	379953	379954	379955	379956	379957	379958	379959	379960	379961	379962	379963	379964	379965	379966	379967	379968	379969	379970	379971	379972	379973	379974	379975	379976	379977	379978	379979];
        del_13604=[381628	381629	381630	381631	381632	381633	381634	381635	381636	381637	381638	381639	381640	381641	381642	381643	381644	381645	381646	381647	381648	381649	381650	381651	381652	381653	381654	381655	381656	381657	381658	381659	381660	381661	381662	381663	381664	381665	381666	381667	381668	381669	381670	381671	381672	381673	381674	381675	381676	381677	381678	381679];
        del_14602=[383411	383412	383413	383414	383415	383416	383417	383418	383419	383420	383421	383422	383423	383424	383425	383426	383427	383428	383429	383430	383431	383432	383433	383434	383435	383436	383437	383438	383439	383440	383441	383442	383443	383444	383445	383446	383447	383448	383449	383450	383451	383452];
        del_15605=[385197	385198	385199	385200	385201	385202	385203	385204	385205	385206	385207	385208	385209	385210	385211	385212	385213	385214	385215	385216	385217	385218	385219	385220	385221	385222	385223	385224	385225	385226	385227	385228	385229	385230	385231	385232	385233	385234	385235	385236	385237	385238	385239	385240	385241	385242	385243	385244	385245	385246	385247];
        del_16593=[386976	386977	386978	386979	386980	386981	386982	386983	386984	386985	386986	386987	386988	386989	386990	386991	386992	386993	386994	386995	386996	386997	386998	386999	387000	387001	387002	387003	387004	387005	387006	387007	387008	387009	387010	387011	387012	387013	387014	387015	387016	387017	387018	387019	387020	387021	387022	387023	387024	387025	387026];
        del_19601=[392493	392494	392495	392496	392497	392498	392499	392500	392501	392502	392503	392504	392505	392506	392507	392508	392509	392510	392511	392512	392513	392514	392515	392516	392517	392518	392519	392520	392521	392522	392523	392524	392525	392526	392527	392528	392529	392530	392531	392532	392533	392534	392535	392536	392537	392538	392539	392540];

        del_all=[del_10593';del_10594';del_10595';del_11595';del_12598';del_13604';del_14602';del_15605';del_16593';del_19601'];

        thisdata.stattable(del_all,:) = [];
        datastruct_1=thisdata.stattable;
        randnr_1=str2double(extractBetween(thisdata.stattable.filename,'thisrand','_nrRanditer'));
        tmpstr=extractBetween(thisdata.stattable.filename,'_nrRanditer_','_rand');
        itnr_1=str2double(extractBefore(tmpstr,'_'));
        for ii=1:length(randnr_1)
            aretheythere(randnr_1(ii),itnr_1(ii))=aretheythere(randnr_1(ii),itnr_1(ii))+1;
        end
        for ii=1:length(randnr_1)
            % tag data fra linje nr ii og put den ind i den der er sorteret
            % på randnr
            rmsd_tpico(randnr_1(ii),itnr_1(ii))=datastruct_1.RMSd_TOTpico(ii);
            rmsd_tnano(randnr_1(ii),itnr_1(ii))=datastruct_1.RMSd_TOTnano(ii);

            rmsd_apico(randnr_1(ii),itnr_1(ii))=datastruct_1.RMSd_ACpico(ii);
            rmsd_anano(randnr_1(ii),itnr_1(ii))=datastruct_1.RMSd_ACnano(ii);
            rmsd_amicro(randnr_1(ii),itnr_1(ii))=datastruct_1.RMSd_ACmicro(ii);

            cor_tpico(randnr_1(ii),itnr_1(ii))=datastruct_1.r_TOTpico(ii);
            cor_tnano(randnr_1(ii),itnr_1(ii))=datastruct_1.r_TOTnano(ii);

            cor_apico(randnr_1(ii),itnr_1(ii))=datastruct_1.r_ACpico(ii);
            cor_anano(randnr_1(ii),itnr_1(ii))=datastruct_1.r_ACnano(ii);
            cor_amicro(randnr_1(ii),itnr_1(ii))=datastruct_1.r_ACmicro(ii);
        end
    end
    if i==2
        del_10593=[376190	376191	376192	376193	376194	376195	376196	376197	376198	376199	376200	376201	376202	376203	376204	376205	376206	376207	376208	376209	376210	376211	376212	376213	376214	376215	376216	376217	376218	376219	376220	376221	376222	376223	376224	376225	376226	376227	376228	376229	376230	376231	376232	376233	376234	376235	376236	376237];
        del_10594=[376238	376239	376240	376241	376242	376243	376244	376245	376246	376247	376248	376249	376250	376251	376252	376253	376254	376255	376256	376257	376258	376259	376260	376261	376262	376263	376264	376265	376266	376267	376268	376269	376270	376271	376272	376273	376274	376275	376276	376277	376278	376279	376280	376281	376282	376283	376284	376285];
        del_10595=[376286	376287	376288	376289	376290	376291	376292	376293	376294	376295	376296	376297	376298	376299	376300	376301	376302	376303	376304	376305	376306	376307	376308	376309	376310	376311	376312	376313	376314	376315	376316	376317	376318	376319	376320	376321	376322	376323	376324	376325	376326	376327	376328	376329	376330	376331	376332	376333];
        del_12598=[379720	379721	379722	379723	379724	379725	379726	379727	379728	379729	379730	379731	379732	379733	379734	379735	379736	379737	379738	379739	379740	379741	379742	379743	379744	379745	379746	379747];
        del_13604=[381428	381429	381430	381431	381432	381433	381434	381435	381436	381437	381438	381439	381440	381441	381442	381443	381444	381445	381446	381447	381448	381449	381450	381451	381452	381453	381454	381455	381456	381457	381458	381459	381460];
        del_14602=[383223	383224	383225	383226	383227	383228	383229	383230	383231	383232	383233	383234	383235	383236	383237	383238	383239	383240	383241	383242	383243	383244	383245	383246];
        del_15605=[385010	385011	385012	385013	385014	385015	385016	385017	385018	385019	385020	385021	385022	385023	385024	385025	385026	385027	385028	385029	385030	385031	385032	385033	385034	385035	385036	385037];
        del_16593=[386778	386779	386780	386781	386782	386783	386784	386785	386786	386787	386788	386789	386790	386791	386792	386793	386794	386795	386796	386797	386798	386799	386800	386801	386802	386803	386804	386805	386806	386807	386808	386809	386810	386811	386812	386813];
        del_19601=[392330	392331	392332	392333	392334	392335	392336	392337	392338	392339	392340	392341	392342	392343	392344	392345	392346	392347	392348	392349	392350	392351	392352	392353	392354	392355	392356	392357	392358];
        del_all=[del_10593';del_10594';del_10595';del_12598';del_13604';del_14602';del_15605';del_16593';del_19601'];
        thisdata.stattable(del_all,:) = [];
        datastruct_2=thisdata.stattable;
        randnr_2=str2double(extractBetween(thisdata.stattable.filename,'thisrand','_nrRanditer'));
        tmpstr=extractBetween(thisdata.stattable.filename,'_nrRanditer_','_rand');
        itnr_2=str2double(extractBefore(tmpstr,'_'));
        for ii=1:length(randnr_2)
            aretheythere(randnr_2(ii),itnr_2(ii))=aretheythere(randnr_2(ii),itnr_2(ii))+1;
        end
        for ii=1:length(randnr_2)
            % tag data fra linje nr ii og put den ind i den der er sorteret
            % på randnr
            rmsd_tpico(randnr_2(ii),itnr_2(ii))=datastruct_2.RMSd_TOTpico(ii);
            rmsd_tnano(randnr_2(ii),itnr_2(ii))=datastruct_2.RMSd_TOTnano(ii);

            rmsd_apico(randnr_2(ii),itnr_2(ii))=datastruct_2.RMSd_ACpico(ii);
            rmsd_anano(randnr_2(ii),itnr_2(ii))=datastruct_2.RMSd_ACnano(ii);
            rmsd_amicro(randnr_2(ii),itnr_2(ii))=datastruct_2.RMSd_ACmicro(ii);

            cor_tpico(randnr_2(ii),itnr_2(ii))=datastruct_2.r_TOTpico(ii);
            cor_tnano(randnr_2(ii),itnr_2(ii))=datastruct_2.r_TOTnano(ii);

            cor_apico(randnr_2(ii),itnr_2(ii))=datastruct_2.r_ACpico(ii);
            cor_anano(randnr_2(ii),itnr_2(ii))=datastruct_2.r_ACnano(ii);
            cor_amicro(randnr_2(ii),itnr_2(ii))=datastruct_2.r_ACmicro(ii);
        end
    end
end
% del_all=[del_10593';del_10594';del_10595';del_12598';del_13604';del_14602';del_15605';del_16593';del_19601'];
[i,~]=find(isnan(rmsd_tpico));
rmsd_tpico(unique(i),:)=[];
[i,~]=find(isnan(rmsd_tnano));
rmsd_tnano(unique(i),:)=[];
[i,~]=find(isnan(rmsd_apico));
rmsd_apico(unique(i),:)=[];
[i,~]=find(isnan(rmsd_anano));
rmsd_anano(unique(i),:)=[];
[i,~]=find(isnan(rmsd_amicro));
rmsd_amicro(unique(i),:)=[];

[i,~]=find(isnan(cor_tpico));
cor_tpico(unique(i),:)=[];
[i,~]=find(isnan(cor_tnano));
cor_tnano(unique(i),:)=[];
[i,~]=find(isnan(cor_apico));
cor_apico(unique(i),:)=[];
[i,~]=find(isnan(cor_anano));
cor_anano(unique(i),:)=[];
[i,~]=find(isnan(cor_amicro));
cor_amicro(unique(i),:)=[];

% increase_nr_runs=1:size(rmsd_apico,1);
% [S1_apico_rmsd,St_apico_rmsd]=Dothethings(increase_nr_runs,rmsd_apico);

%% fig1=figure('name','apico');
% fig1=figure('Name','apico','Color','w','units','centimeters','position',[15,5,16,23]);
% fig1.Renderer='Painters';
% tiledlayout(6,4,'TileSpacing','tight','Padding','tight')
% addtoplot(fig1,increase_nr_runs(1:5:length(increase_nr_runs)),St_apico_rmsd(:,1:5:length(increase_nr_runs)),allparameters)

%%
% fig2=figure('name','anano');
% increase_nr_runs=1:size(rmsd_anano,1);
% [S1_anano_rmsd,St_anano_rmsd]=Dothethings(increase_nr_runs,rmsd_anano);
% addtoplot(fig2,increase_nr_runs,St_anano_rmsd,allparameters_Label)
%%
% fig3=figure('name','anmicro');
% rmsd_amicro(10012,:)=[];
% increase_nr_runs=1:size(rmsd_amicro,1);
% [S1_amicro_rmsd,St_amicro_rmsd]=Dothethings(increase_nr_runs,rmsd_amicro);
% addtoplot(fig3,increase_nr_runs,St_amicro_rmsd,allparameters_Label)

%%
% increase_nr_runs=1:size(rmsd_tpico,1);
% [S1_tpico_rmsd,St_tpico_rmsd]=Dothethings(increase_nr_runs,rmsd_tpico);
% 
% increase_nr_runs=1:size(rmsd_tnano,1);
% [S1_tnano_rmsd,St_tnano_rmsd]=Dothethings(increase_nr_runs,rmsd_tnano);
%% cor
increase_nr_runs=size(cor_tpico,1);%1:size(cor_tpico,1);
[S1_tpico_cor,St_tpico_cor]=Dothethings(increase_nr_runs,cor_tpico);
% 
increase_nr_runs=size(cor_tnano,1);%1:size(cor_tnano,1);
[S1_tnano_cor,St_tnano_cor]=Dothethings(increase_nr_runs,cor_tnano);
%%
% fig2=figure('name','anano');
increase_nr_runs=size(cor_apico,1);%1:size(cor_apico,1);
[S1_apico_cor,St_apico_cor]=Dothethings(increase_nr_runs,cor_apico);
%%
% fig2=figure('name','anano');
% tiledlayout(6,4,'TileSpacing','tight','Padding','tight')
% addtoplot(fig2,increase_nr_runs,St_apico_cor,allparameters_Label)
%%

increase_nr_runs=size(cor_anano,1);%1:size(cor_anano,1);
[S1_anano_cor,St_anano_cor]=Dothethings(increase_nr_runs,cor_anano);

increase_nr_runs=size(cor_amicro,1);%1:size(cor_amicro,1);
[S1_amicro_cor,St_amicro_cor]=Dothethings(increase_nr_runs,cor_amicro);

% addlinetoplot(fig7)
% addlinetoplot(fig8)
% addlinetoplot(fig9)
% addlinetoplot(fig10)
% addlinetoplot(fig11)
%%
allparameters2=allparameters_Label;
allparameters2(23)=[];
% [rankS1_tpico_rmsd, idxrankS1_tpico_rmsd] = sort(S1_tpico_rmsd([1:22 24],end));[rankSt_tpico_rmsd, idxrankSt_tpico_rmsd] = sort(St_tpico_rmsd([1:22 24],end));
% [rankS1_tnano_rmsd, idxrankS1_tnano_rmsd] = sort(S1_tnano_rmsd([1:22 24],end));[rankSt_tnano_rmsd, idxrankSt_tnano_rmsd] = sort(St_tnano_rmsd([1:22 24],end));
% [rankS1_apico_rmsd, idxrankS1_apico_rmsd] = sort(S1_apico_rmsd([1:22 24],end));[rankSt_apico_rmsd, idxrankSt_apico_rmsd] = sort(St_apico_rmsd([1:22 24],end));
% [rankS1_anano_rmsd, idxrankS1_anano_rmsd] = sort(S1_anano_rmsd([1:22 24],end));[rankSt_anano_rmsd, idxrankSt_anano_rmsd] = sort(St_anano_rmsd([1:22 24],end));
% [rankS1_amicro_rmsd, idxrankS1_amicro_rmsd] = sort(S1_amicro_rmsd([1:22 24],end));[rankSt_amicro_rmsd, idxrankSt_amicro_rmsd] = sort(St_amicro_rmsd([1:22 24],end));
%%
[rankS1_tpico_cor, idxrankS1_tpico_cor] = sort(S1_tpico_cor([1:22 24],end));[rankSt_tpico_cor, idxrankSt_tpico_cor] = sort(St_tpico_cor([1:22 24],end));
[rankS1_tnano_cor, idxrankS1_tnano_cor] = sort(S1_tnano_cor([1:22 24],end));[rankSt_tnano_cor, idxrankSt_tnano_cor] = sort(St_tnano_cor([1:22 24],end));
[rankS1_apico_cor, idxrankS1_apico_cor] = sort(S1_apico_cor([1:22 24],end));[rankSt_apico_cor, idxrankSt_apico_cor] = sort(St_apico_cor([1:22 24],end));
[rankS1_anano_cor, idxrankS1_anano_cor] = sort(S1_anano_cor([1:22 24],end));[rankSt_anano_cor, idxrankSt_anano_cor] = sort(St_anano_cor([1:22 24],end));
[rankS1_amicro_cor, idxrankS1_amicro_cor] = sort(S1_amicro_cor([1:22 24],end));[rankSt_amicro_cor, idxrankSt_amicro_cor] = sort(St_amicro_cor([1:22 24],end));


%%
n_p=length(allparameters2);
% fig9=figure('Name','Sobol rmsd','Color','w','units','centimeters','position',[15,5,11.2,17]);
% tiledlayout(5,1,'TileSpacing','compact','Padding','compact');nexttile;
% plottile(rankSt_tpico_rmsd,idxrankSt_tpico_rmsd,n_p,allparameters2);title('tpico');nexttile;
% plottile(rankSt_tnano_rmsd,idxrankSt_tnano_rmsd,n_p,allparameters2);title('tnano');nexttile;
% plottile(rankSt_apico_rmsd,idxrankSt_apico_rmsd,n_p,allparameters2);title('apico');nexttile;
% plottile(rankSt_anano_rmsd,idxrankSt_anano_rmsd,n_p,allparameters2);title('anano');nexttile;
% plottile(rankSt_amicro_rmsd,idxrankSt_amicro_rmsd,n_p,allparameters2);title('amicro')

%%
fig4=figure('Name','Sobol cor','Color','w','units','centimeters','position',[15,5,11.2,17]);
tiledlayout(5,1,'TileSpacing','compact','Padding','compact');nexttile;
plottile(rankSt_tpico_cor,idxrankSt_tpico_cor,n_p,allparameters2);title('tpico');nexttile;
plottile(rankSt_tnano_cor,idxrankSt_tnano_cor,n_p,allparameters2);title('tnano');nexttile;
plottile(rankSt_apico_cor,idxrankSt_apico_cor,n_p,allparameters2);title('apico');nexttile;
plottile(rankSt_anano_cor,idxrankSt_anano_cor,n_p,allparameters2);title('anano');nexttile;
plottile(rankSt_amicro_cor,idxrankSt_amicro_cor,n_p,allparameters2)
% %%
% fig5=figure('Name','Sobol 1 rmsd ','Color','w','units','centimeters','position',[15,5,11.2,17]);
% tiledlayout(5,1,'TileSpacing','compact','Padding','compact');nexttile;
% plottile(rankS1_tpico_rmsd,idxrankS1_tpico_rmsd,n_p,allparameters);nexttile;
% plottile(rankS1_tnano_rmsd,idxrankS1_tnano_rmsd,n_p,allparameters);nexttile;
% plottile(rankS1_apico_rmsd,idxrankS1_apico_rmsd,n_p,allparameters);nexttile;
% plottile(rankS1_anano_rmsd,idxrankS1_anano_rmsd,n_p,allparameters);nexttile;
% plottile(rankS1_amicro_rmsd,idxrankS1_amicro_rmsd,n_p,allparameters);

%% FUNCTIONS
function [f0,D]=calcf0_D(y0)
% y0(isnan(y0))=0;
nlength=length(y0);
f0=sum(y0,'omitnan')./nlength;
D=sum(y0.^2,'omitnan')./nlength;
D=D-f0.^2;
end

function [Di,Ditot]=calcDi_Ditot(D_t,y0,y1,y2)
% y0(isnan(y0))=0; y0=y0(theinterval);
% y1(isnan(y1))=0; y1=y1(:,theinterval);
% y2(isnan(y2))=0; y2=y2(:,theinterval);
nlength=length(y0);
Di=ones(size(y1,2),1).*D_t;
Ditot=zeros(size(y1,2),1);
% figure;
for i=1:length(y0)
    for j=1:length(Di)
%         if i>=9970 && j==23
%         disp('now')
%         end
        Di(j)=Di(j)-(y0(i)-y1(i,j)).^2./(2*nlength);
        Ditot(j)=Ditot(j)+(y0(i)-y2(i,j)).^2./(2*nlength);
    end
end
% Di=D_t-(sum((y0-y1).^2,2,'omitnan')./(2*nlength));
% Ditot=(sum((y0-y2).^2,2,'omitnan')./(2*nlength));
end

function [S1,St]=calc_S(Di,Ditot,D)
S1 = Di./D; %first order effect sensitivity indices
St = Ditot./D; % total effect sensitivity indices
end

function plottile(rankS1_tpico_rmsd,idxrankS1_tpico_rmsd,n_p,allparameters)
plot(1:n_p,rankS1_tpico_rmsd,'LineStyle','-','Color','k')
hold on
plot(1:n_p,rankS1_tpico_rmsd,'.','Color','k','MarkerSize',15)
ax=gca;
ax.XTick=1:n_p;
ax.XTickLabel=allparameters(idxrankS1_tpico_rmsd);
xlim([1 n_p])
ylim([0 0.8])
box off
end

function addtoplot(fignr,increase_nr_runs,St_tnano_rmsd,allparameters)
figure(fignr)

hold on
for j=[1:22,24]
    nexttile(j)
%     subplot(2,12,j)
    hold on
    plot(increase_nr_runs,St_tnano_rmsd(j,:),'-','LineWidth',1.5)
    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.TickDir='out';
    ax.XAxis.LineWidth=1.5;
    ax.YAxis.LineWidth=1.5;
        
    ylim([0 1])
    if ismember(j,[1 5 9 13 17 21])
        ylabel('STi')
    else
%         ax.YAxis.TickLabelColor='w';
    end
    if j>19
        xlabel('No. iterations')
    else
%         ax.XAxis.TickLabelColor='w';
    end

    title(allparameters{j}, 'Position',[4000 0.9 0],'HorizontalAlignment','left')
end
end


function addlinetoplot(fignr)
figure(fignr)
hold on
for j=1:24
    subplot(2,12,j)
    hold on
    plot([10000 10000],ylim,'--k')
end
end
function [S1_tpico_rmsd,St_tpico_rmsd]=Dothethings(increase_nr_runs,rmsd_tpico)
St_tpico_rmsd=nan(24,length(increase_nr_runs));
S1_tpico_rmsd=nan(24,length(increase_nr_runs));
for i=1:length(increase_nr_runs)
    theinterval=1:increase_nr_runs(i);
    [~,D_tpico_rmsd]=calcf0_D(rmsd_tpico(theinterval,1));
    [Di_tpico_rmsd,Ditot_tpico_rmsd]=calcDi_Ditot(D_tpico_rmsd,rmsd_tpico(theinterval,1),rmsd_tpico(theinterval,26:49),rmsd_tpico(theinterval,2:25));
    [S1_tpico_rmsd(:,i),St_tpico_rmsd(:,i)]=calc_S(Di_tpico_rmsd,Ditot_tpico_rmsd,D_tpico_rmsd);
end
end

