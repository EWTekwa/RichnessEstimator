function [Richness_apx,Richness_taylor,Richness_taylor_0,Apx_T_detectP_terms] = RichnessEsts_ideal(TransectAbundance,k,Richness_raw)

%RichnessEsts_ideal.m
%Eden Tekwa Feb 8, 2022 - Apr 11, 2022
%function returns Taylor and Apx versions of omega with survivorship bias eliminated for richness means based on
%the spatial TransectAbundance data: rows=transects, columns=species,
%values=individual counts

TransectAbundance=TransectAbundance(:,sum(TransectAbundance,1)>0); %take out empty species columns

%write observational process model
syms nm P
Dix=(1-exp(-nm/P))*P; %local detection probability of species i at sampled site x
Dis=1-(1-Dix)^k; %detection probability of species i across all sampled sites x in community s

%second-order partial derivatives for 2nd order Taylor expansion of overall detection probability of any species in the community
%     d2Dis_dnm2=diff(Dis,'nm',2);
%     d2Dis_dP2=diff(Dis,'P',2);
%     d2Dis_dnmP=diff(diff(Dis,'nm'),'P');

%hard coded results from above for speed:
d2Dis_dnm2 = - (k*exp(-nm/P)*(P*(exp(-nm/P) - 1) + 1)^(k - 1))/P - k*exp(-(2*nm)/P)*(P*(exp(-nm/P) - 1) + 1)^(k - 2)*(k - 1);
d2Dis_dP2 = - k*(P*(exp(-nm/P) - 1) + 1)^(k - 2)*(k - 1)*(exp(-nm/P) + (nm*exp(-nm/P))/P - 1)^2 - (k*nm^2*exp(-nm/P)*(P*(exp(-nm/P) - 1) + 1)^(k - 1))/P^3;
d2Dis_dnmP = (k*nm*exp(-nm/P)*(P*(exp(-nm/P) - 1) + 1)^(k - 1))/P^2 + k*exp(-nm/P)*(P*(exp(-nm/P) - 1) + 1)^(k - 2)*(k - 1)*(exp(-nm/P) + (nm*exp(-nm/P))/P - 1);

numTrans=size(TransectAbundance,1); %get number of transects
Richness_true=sum(sum(TransectAbundance,1)>0); %get raw richness
n_m_true=mean(TransectAbundance,1);
P_true=zeros(1,Richness_true); %array to record occupancy for each species
for species=1:Richness_true
    P_true(species)=sum(TransectAbundance(:,species)>0)/numTrans; %occupancy as number of transects occupied divided by number of transects
end
P_true(P_true==0)=NaN;
mean_P_true=nanmean(P_true);
var_P_true=nanvar(P_true);
mean_n_m_true=mean(mean(TransectAbundance(:,sum(TransectAbundance)>0)));
var_n_m_true=var(mean(TransectAbundance(:,sum(TransectAbundance)>0)));
cov_nm_P_true=nancov(mean(TransectAbundance(:,sum(TransectAbundance)>0)),P_true);
if size(cov_nm_P_true,2)>1
    cov_nm_P_true=cov_nm_P_true(1,2);
else
    cov_nm_P_true=0;
end


%compute correction terms for omega_T
nm=mean_n_m_true;
P=mean_P_true;
Apx_T_detectP_terms=[eval(Dis)
    eval(d2Dis_dnm2)*var_n_m_true/2
    eval(d2Dis_dP2)*var_P_true/2
    eval(d2Dis_dnmP)*cov_nm_P_true];

if sum(Apx_T_detectP_terms)>0.1 %if sum of correction terms is greater than a threshold
    Ds_taylor=sum(Apx_T_detectP_terms); %use full correction
else
    Ds_taylor=Apx_T_detectP_terms(1); %else, use 0th order correction
end
if Ds_taylor>1
    Ds_taylor=1;
end
if Ds_taylor<0.1
    Ds_taylor=0.1;
end
Richness_taylor=Richness_raw/Ds_taylor;

%compute omega_0
Richness_taylor_0=Richness_raw/Apx_T_detectP_terms(1);

%compute omega_apx
nm=n_m_true;
P=P_true;
Ds_means=eval(Dis);
Ds_mean=mean(Ds_means);
Richness_apx=Richness_raw/Ds_mean;
if Richness_raw==0
    Richness_apx=NaN;
    Richness_taylor=NaN;
    Richness_taylor_0=NaN;
end