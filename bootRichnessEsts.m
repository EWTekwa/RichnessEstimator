function [Richness_raw,Chao1,Richness_apx,expectedRichness_raw,expectedChao1,expectedRichness_apx] = bootRichnessEsts(TransectAbundance,numBoot)
%bootRichnessEsts.m
%Ed Tekwa Feb 8, 2022
%function returns Chao1 and Taylor2 Apx for richness means and confidence intervals based on
%the spatial TransectAbundance data: rows=transects, columns=species,
%values=individual counts

expectedRichness_raw=zeros(numBoot,1);
expectedChao1=zeros(numBoot,1);
expectedRichness_apx=zeros(numBoot,1);

[Richness_raw,Chao1,~,~,~,~,Richness_apx,~] = RichnessEstsCov(TransectAbundance); %get point estimates from original dataset

maxSampleRichness=length(TransectAbundance);
for resample=1:numBoot
    sampleSetApx=TransectAbundance(:,randi(maxSampleRichness,1,maxSampleRichness)); %resample species with replacement
    for j=1:size(sampleSetApx,1) %for each species
        sampleSetApx(:,j)=sampleSetApx(randi(size(sampleSetApx,1),size(sampleSetApx,1),1),j); %resample with replacement transects for each species
    end
    [Richness_raw_boot,Chao1_boot,~,~,~,~,Richness_apx_boot,~] = RichnessEstsCov(sampleSetApx);
    expectedRichness_raw(resample)=Richness_raw_boot;
    expectedChao1(resample)=Chao1_boot;
    expectedRichness_apx(resample)=Richness_apx_boot;
end