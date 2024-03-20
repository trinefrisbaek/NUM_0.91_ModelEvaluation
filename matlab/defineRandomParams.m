function [s, meanval, sigma]= defineRandomParams
%-------------------------------------------------------------------------
% create random initialisation seed and save
%-------------------------------------------------------------------------
s=rng('shuffle');% seed random initial random number
% ------------------------------------------------------------------------
% r*d: Diffusive affility cross-over:
% ------------------------------------------------------------------------
% LOGNORMAL DISTRIBUTION:
meanval.r_star_d=log(2); %0.3);     %=-1.2040
sigma.r_star_d=log(1.2214); %=0.2
%r_star_d(:)=exp(meanval);
% ------------------------------------------------------------------------
% Light harvesting: alphaL and r*L
% aL = αLr−1(1 − exp(−r/r∗L))(1 − ν)
%     set mean and sigmma for rand_y, not alpha_l
% ------------------------------------------------------------------------
meanval.rand_y=log(0.25);
sigma.rand_y=log(1.6487);
p_const=0.4;
% rand_y(:)=exp(meanval);
% alpha_l=(3.*rand_y)./(4*p_const);

meanval.rand_rlstar=log(7.5); %
sigma.rand_rlstar=log(1.6487);
%rand_rlstar(:)=exp(meanval);
% ------------------------------------------------------------------------
% phagotrophic clearance rate: aF
% ------------------------------------------------------------------------
meanval.aF_random=log(0.0189);
sigma.aF_random=log(4.8576);
%aF_random(:)=exp(meanval);
% ------------------------------------------------------------------------
% passive losses: cpassive
% ------------------------------------------------------------------------
meanval.c_passive_random=log(0.03);
sigma.c_passive_random=log(1.8221);
%c_passive_random(:)=exp(meanval);
% ------------------------------------------------------------------------
% maximum synthesis rate: αMax
% ------------------------------------------------------------------------
meanval.alphaMax_rand=log(0.4883);
sigma.alphaMax_rand=log(2.1029);
%alphaMax_rand(:)=exp(meanval);
% ------------------------------------------------------------------------
% basal metabolism coef: αR
% ------------------------------------------------------------------------
meanval.alphaR_rand=log(0.1);
sigma.alphaR_rand=log(1.4918);
%alphaR_rand(:)=exp(meanval);
% ------------------------------------------------------------------------
% mortHTL
% ------------------------------------------------------------------------
meanval.mortHTL_rand=log(0.01);
sigma.mortHTL_rand=log(3.5);
% ------------------------------------------------------------------------
% rhoCN ***new***
% ------------------------------------------------------------------------
meanval.rand_rhocn=log(5.68);
sigma.rand_rhocn=log(1.5);
% ------------------------------------------------------------------------
% fracHTL_to_N
% ------------------------------------------------------------------------
meanval.rand_fracHTLN=log(0.5);
sigma.rand_fracHTLN=log(0.2);
% ------------------------------------------------------------------------
% epsilonL
% ------------------------------------------------------------------------
meanval.rand_epsilonL=log(0.8);
sigma.rand_epsilonL=log(1.1);
% ------------------------------------------------------------------------
% alphaN
% ------------------------------------------------------------------------
meanval.rand_alphaN=log(0.972);
sigma.rand_alphaN=log(0.1);
% ------------------------------------------------------------------------
% epsilonF
% ------------------------------------------------------------------------
meanval.rand_epsilonF=log(0.8); 
sigma.rand_epsilonF=log(1.1);
% ------------------------------------------------------------------------
% cF
% ------------------------------------------------------------------------
meanval.rand_cF=log(30);
sigma.rand_cF=log(5);
% ------------------------------------------------------------------------
% beta
% ------------------------------------------------------------------------
meanval.rand_beta=log(500);
sigma.rand_beta=log(100);
% ------------------------------------------------------------------------
% sigma
% ------------------------------------------------------------------------
meanval.rand_sigma=log(1.3);
sigma.rand_sigma=log(0.2);

% ------------------------------------------------------------------------
% remin2
% ------------------------------------------------------------------------
meanval.rand_remin2=log(0.5);
sigma.rand_remin2=log(0.2);
% ------------------------------------------------------------------------
% reminF
% ------------------------------------------------------------------------
meanval.rand_reminF=log(0.1);
sigma.rand_reminF=log(3);
% ------------------------------------------------------------------------
% rhoval
% ------------------------------------------------------------------------
meanval.rand_rhoval=log(0.4); %missing
sigma.rand_rhoval=log(0.15);
% ------------------------------------------------------------------------
% mort2val
% ------------------------------------------------------------------------
meanval.rand_mort2val=log(0.004);
sigma.rand_mort2val=log(0.0015);
% ------------------------------------------------------------------------
% remPOM
% ------------------------------------------------------------------------
meanval.rand_remPOM=log(0.02);
sigma.rand_remPOM=log(3);
% ------------------------------------------------------------------------
% vel1
% ------------------------------------------------------------------------
meanval.rand_vel1=log(8);
sigma.rand_vel1=log(50);
% ------------------------------------------------------------------------
% vel2
% ------------------------------------------------------------------------
meanval.rand_vel2=log(0.5);
sigma.rand_vel2=log(0.2);
% ------------------------------------------------------------------------
end