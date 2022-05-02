function [Richness_raw,Chao1,Chao2,ACE,S_aj2,S_ij2,Richness_apx,expectedRichness_raw,expectedChao1,expectedChao2,expectedACE,expectedS_aj2,expectedS_ij2,expectedRichness_apx] = bootRichnessEsts_all(TransectAbundance,numBoot)
%bootRichnessEsts.m
%Ed Tekwa Feb 8, 2022
%function returns Chao1 and Taylor2 Apx for richness means and confidence intervals based on
%the spatial TransectAbundance data: rows=transects, columns=species,
%values=individual counts

%store bootstrapped estimates for the 7 estimators:
expectedRichness_raw=zeros(numBoot,1); %raw
expectedChao1=zeros(numBoot,1); %Chao1
expectedChao2=zeros(numBoot,1); %Chao2
expectedACE=zeros(numBoot,1); %ACE
expectedS_aj2=zeros(numBoot,1); %Jackknife (abundance)
expectedS_ij2=zeros(numBoot,1); %Jackknife (incidence)
expectedRichness_apx=zeros(numBoot,1);

[Richness_raw,Chao1,Chao2,ACE,S_aj2,S_ij2,Richness_apx,~] = RichnessEstsCov(TransectAbundance); %get point estimates from original dataset

maxSampleRichness=length(TransectAbundance);
for resample=1:numBoot
    sampleSetApx=TransectAbundance(:,randi(maxSampleRichness,1,maxSampleRichness)); %resample species with replacement
    for j=1:size(sampleSetApx,1) %for each species
        sampleSetApx(:,j)=sampleSetApx(randi(size(sampleSetApx,1),size(sampleSetApx,1),1),j); %resample with replacement transects for each species
    end
    [Richness_raw_boot,Chao1_boot,Chao2_boot,ACE_boot,S_aj2_boot,S_ij2_boot,Richness_apx_boot,~] = RichnessEstsCov(sampleSetApx);
    expectedRichness_raw(resample)=Richness_raw_boot;
    expectedChao1(resample)=Chao1_boot;
    expectedChao2(resample)=Chao2_boot;
    expectedACE(resample)=ACE_boot;
    expectedS_aj2(resample)=S_aj2_boot;
    expectedS_ij2(resample)=S_ij2_boot;
    expectedRichness_apx(resample)=Richness_apx_boot;
end