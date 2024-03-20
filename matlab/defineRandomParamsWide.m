function [s, meanval, sigma]= defineRandomParamsWide
%-------------------------------------------------------------------------
% create random initialisation seed and save
%-------------------------------------------------------------------------
s=rng('shuffle');% seed random initial random number
% ------------------------------------------------------------------------
% r*d: Diffusive affility cross-over: 1
% ------------------------------------------------------------------------
% LOGNORMAL DISTRIBUTION:
meanval.r_star_d=log(2); %0.3);     %=-1.2040
sigma.r_star_d=log(2.5); %=0.2
%r_star_d(:)=exp(meanval);
% ------------------------------------------------------------------------
% Light harvesting: alphaL and r*L
% aL = αLr−1(1 − exp(−r/r∗L))(1 − ν)
%     set mean and sigma for rand_y, not alpha_l 2 3
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
% phagotrophic clearance rate: aF 4
% ------------------------------------------------------------------------
meanval.aF_random=log(0.0189);
sigma.aF_random=log(4.8576);
%aF_random(:)=exp(meanval);
% ------------------------------------------------------------------------
% passive losses: cpassive 5
% ------------------------------------------------------------------------
meanval.c_passive_random=log(0.03);
sigma.c_passive_random=log(1.8221);
%c_passive_random(:)=exp(meanval);
% ------------------------------------------------------------------------
% maximum synthesis rate: αMax 6
% ------------------------------------------------------------------------
meanval.alphaMax_rand=log(0.4883);
sigma.alphaMax_rand=log(2.1029);
%alphaMax_rand(:)=exp(meanval);
% ------------------------------------------------------------------------
% basal metabolism coef: αR 7
% ------------------------------------------------------------------------
meanval.alphaR_rand=log(0.1);
sigma.alphaR_rand=log(1.4918);
%alphaR_rand(:)=exp(meanval);
% ------------------------------------------------------------------------
% mortHTL 8
% ------------------------------------------------------------------------
meanval.mortHTL_rand=log(0.01);
sigma.mortHTL_rand=log(4.5);%3.5);
% ------------------------------------------------------------------------
% rhoCN ***new*** 9
% ------------------------------------------------------------------------
meanval.rand_rhocn=log(5.68);
sigma.rand_rhocn=log(1.5);
% ------------------------------------------------------------------------
% fracHTL_to_N 10 
% ------------------------------------------------------------------------
meanval.rand_fracHTLN=log(0.5);
sigma.rand_fracHTLN=log(0.2);
% ------------------------------------------------------------------------
% epsilonL 11
% ------------------------------------------------------------------------
meanval.rand_epsilonL=log(0.5);
sigma.rand_epsilonL=log(0.2);
% ------------------------------------------------------------------------
% alphaN 12
% ------------------------------------------------------------------------
meanval.rand_alphaN=log(0.972);
sigma.rand_alphaN=log(0.1);
% ------------------------------------------------------------------------
% epsilonF 13
% ------------------------------------------------------------------------
meanval.rand_epsilonF=log(0.5); 
sigma.rand_epsilonF=log(0.2);
% ------------------------------------------------------------------------
% cF 14
% ------------------------------------------------------------------------
meanval.rand_cF=log(30);
sigma.rand_cF=log(5);
% ------------------------------------------------------------------------
% beta 15 prey-predator mass ratio (500)
% ------------------------------------------------------------------------
meanval.rand_beta=log(500);
sigma.rand_beta=log(100);
% ------------------------------------------------------------------------
% sigma 16  prey-predator with (1.3)
% ------------------------------------------------------------------------
meanval.rand_sigma=log(1.3);
sigma.rand_sigma=log(0.2);

% ------------------------------------------------------------------------
% remin2 17
% ------------------------------------------------------------------------
meanval.rand_remin2=log(0.5);
sigma.rand_remin2=log(0.2);
% ------------------------------------------------------------------------
% reminF 18
% ------------------------------------------------------------------------
meanval.rand_reminF=log(0.1);
sigma.rand_reminF=log(3);
% ------------------------------------------------------------------------
% rhoval 19
% ------------------------------------------------------------------------
meanval.rand_rhoval=log(0.4); %missing
sigma.rand_rhoval=log(0.01);
% ------------------------------------------------------------------------
% mort2val 20
% ------------------------------------------------------------------------
meanval.rand_mort2val=log(0.004);
sigma.rand_mort2val=log(0.002);
% ------------------------------------------------------------------------
% remPOM 21
% ------------------------------------------------------------------------
meanval.rand_remPOM=log(0.02);
sigma.rand_remPOM=log(3);
% ------------------------------------------------------------------------
% vel1  22
% ------------------------------------------------------------------------
meanval.rand_vel1=log(8);
sigma.rand_vel1=log(2.5);
% ------------------------------------------------------------------------
% vel2
% ------------------------------------------------------------------------
meanval.rand_vel2=log(0.5);
sigma.rand_vel2=log(0.0001);
% ------------------------------------------------------------------------
end