function [Richness_raw,Chao1,GP,Chao2,ACE,S_aj2,S_ij2,Richness_apx,Richness_taylor,Richness_taylor_0,Apx_detectP_terms,meanStates] = RichnessEsts(TransectAbundance)

%RichnessEsts.m
%Eden Tekwa Feb 8, 2022 - Apr 11, 2022
%function returns Chao1 and Taylor2 Apx for richness means based on
%the spatial TransectAbundance data: rows=transects, columns=species,
%values=individual counts

TransectAbundance=TransectAbundance(:,sum(TransectAbundance,1)>0); %take out empty species columns

%write observational process model
syms nm k P xi
xi=(1-exp(-nm/P));
Dix=xi*P;
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
Richness_raw=sum(sum(TransectAbundance,1)>0); %get raw richness
Ds=zeros(Richness_raw,1); %store observed species' observation probabilities
P_detected=zeros(1,Richness_raw); %array to record occupancy for each species
for species=1:Richness_raw

    P_detected(species)=sum(TransectAbundance(:,species)>0)/numTrans; %occupancy as number of transects occupied divided by number of transects
end
P_detected(P_detected==0)=NaN;
%compute on full dataset:
mean_P_detected=nanmean(P_detected);
var_P_detected=nanvar(P_detected);
n_m_detected=mean(TransectAbundance,1);
if isnan(n_m_detected)
    n_m_detected=[];
end
mean_n_m_detected=mean(n_m_detected);
var_n_m_detected=var(n_m_detected);
cov_nm_P_detected=nancov(n_m_detected,P_detected);
if size(cov_nm_P_detected,2)>1
    cov_nm_P_detected=cov_nm_P_detected(1,2);
else
    cov_nm_P_detected=0;
end

%compute Chao1 estimate
TransectAbundance=ceil(TransectAbundance); %convert data to discrete count
f1=sum(sum(TransectAbundance,1)==1); %number of singleton species
f2=sum(sum(TransectAbundance,1)==2); %number of doubleton species
Chao1=Richness_raw+f1*(f1-1)/(2*(f2+1)); %Chao1 richness estimator

%compute GP (Gamma-Poisson mixture) estimate (Chiu 2023, peerJ)
f3=sum(sum(TransectAbundance,1)==3); %number of tripleton species
if f3==0
    f3c=1;
else
    f3c=f3;
end
if f1==0
    f1c=1;
else
    f1c=f1;
end
A=2-(2*f2^2)/(3*f1c*f3c);
if f2>0
    f0Chao1=f1c^2/(2*f2);
else
    f0Chao1=f1c*(f1c-1)/2;
end
if A<1
    GP=Richness_raw+f0Chao1*max(.5,A);
else
    GP=Richness_raw+f0Chao1;
end


%compute Chao2 estimate
q1=sum(sum(TransectAbundance>0)==1); %number of species occurring in one transect only
q2=sum(sum(TransectAbundance>0)==2); %number of species occurring in two transect only
m=sum(sum(TransectAbundance>0));
Chao2=Richness_raw+((m-1)/m)*q1*(q1-1)/(2*(q2+1)); %Chao2 richness estimator

%compute abundance-based coverage estimator (ACE)
S_rare=sum(sum(TransectAbundance,1)<=10); %number of rare species (<=10 individuals)
S_abund=sum(sum(TransectAbundance,1)>10); %number of rare species (<=10 individuals)
n_rare=sum(TransectAbundance(:,sum(TransectAbundance,1)<=10),'all'); %total number of individuals in rare species
C_ACE=1-f1/n_rare; %sample coverage estimate
wsumfa=0;
for a=1:10
    wsumfa=wsumfa+a*(a-1)*sum(sum(TransectAbundance,1)==a);
end
V2=max(((S_rare/C_ACE)*wsumfa/(n_rare*(n_rare-1))-1),0); %coefficient of variation
if C_ACE>0
    ACE=S_abund+S_rare/C_ACE+(f1/C_ACE)*V2;
else
    ACE=Chao1;
end

%compute jackknife abundance estimator
S_aj2=Richness_raw+2*f1-f2;

%compute jackknife incidence estimator
if m==1
    S_ij2=Richness_raw;
else
    S_ij2=Richness_raw+(q1*(2*m-3)/m-q2*((m-2)^2)/(m*(m-1)));
end

%compute correction terms for omega_T
k=numTrans;
nm=mean_n_m_detected;
P=mean_P_detected;
Apx_detectP_terms=([eval(Dis)
    eval(d2Dis_dnm2)*var_n_m_detected/2
    eval(d2Dis_dP2)*var_P_detected/2
    eval(d2Dis_dnmP)*cov_nm_P_detected]);

if sum(Apx_detectP_terms)>0.1 %if sum of correction terms is greater than a threshold
    Ds_apx=sum(Apx_detectP_terms); %use full correction
else
    Ds_apx=Apx_detectP_terms(1); %else, use 0th order correction
end
if Ds_apx>1
    Ds_apx=1;
end
if Ds_apx<0.1
    Ds_apx=0.1;
end
Richness_taylor=Richness_raw/Ds_apx;
Richness_taylor_0=Richness_raw/Apx_detectP_terms(1);
meanStates=[mean_n_m_detected,mean_P_detected,var_n_m_detected,var_P_detected,cov_nm_P_detected]';

%compute omega_o
nm=n_m_detected;
P=P_detected;
Ds_means=eval(Dis);
Ds_mean=mean(Ds_means);
Richness_apx=Richness_raw/Ds_mean;

if Richness_raw==0
    Richness_raw=NaN;
    Richness_taylor=NaN;
    Richness_apx=NaN;
    Chao1=NaN;
    Chao2=NaN;
    ACE=NaN;
    S_aj2=NaN;
    S_ij2=NaN;
end