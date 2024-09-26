 
addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));

%% In this code script I am determining the grassland grazing capacity
load('grassland_env.mat')
load('livestockdensity.mat')
load('AGB.mat')
load('aggregate_tas.mat')
load('niche.mat')% variable name: landuse_coupled1 (in run4)

%% read population density data
[pop, Rpop] = readgeoraster('gpw_v4_population_density_rev11_2015_2pt5_min.tif')
Rpop = georefcells([-90,90],[-180,180],size(pop));
[resizedpop,resizedRpop] = georesize(pop,Rpop,1800/4320,"bilinear");
resizedpop(resizedpop<0)=0
imagesc(resizedpop)
colorbar
%% livestock data-cattle
[cattle, Rcattle] = readgeoraster('/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/livestockDensity/cattle/Glb_Cattle_CC2006_AD.tif');

[allresizecattle,resizedRcattle] = georesize(cattle,Rcattle,1/12,"bilinear");
allresizecattle(allresizecattle<0)=0

resizecattle = allresizecattle; % First, make allresizecattle a copy of resizecattle
resizecattle(resizecattle>250) = 0;
resizecattle(resizedpop >20) = 0; % Then, set values to 0 where grassland_env is not 0
% cattlenum=sum(sum(resizecattle.*worldarea));
% cattlenumall=sum(sum(allresizecattle.*worldarea));
imagesc(allresizecattle)
colorbar

imagesc(resizecattle)
colorbar
R=georefcells([-90,90],[-180,180],size(resizecattle));

%% livestock-sheep
[sheep, Rsheep] = readgeoraster('/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/livestockDensity/sheep/Glb_SHAD_2006.tif');
[allresizesheep,resizedRsheep] = georesize(sheep,Rsheep,1/12,"bilinear");
allresizesheep(allresizesheep<0) =0;

% keeping only the grassland sheep
resizesheep=allresizesheep;
resizesheep(resizedpop >20)=0
% sheepsnum=sum(sum(resizesheep.*worldarea));
% sheepsnumall=sum(sum(allresizesheep.*worldarea));

%% livestock-goats
[goats, Rgoats] = readgeoraster('/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/livestockDensity/goats/Glb_GTAD_2006.tif');
Rgoats=Rsheep;
[allresizegoats,resizedRgoats] = georesize(goats,Rgoats,1/12,"bilinear");
allresizegoats(allresizegoats<0) =0;
allresizegoats(allresizegoats>250) = 250;
% keeping only the grassland goats
resizegoats=allresizegoats;
resizegoats(resizedpop >20)=0;
% goatsnum=sum(sum(resizegoats.*worldarea));
% goatsnumall=sum(sum(allresizegoats.*worldarea));

%% livestock unit
livestock=resizecattle+0.2*resizegoats+0.2*resizesheep

%% The livestock data processing is now complete. 
% save the livestock data
save('./grazingniche/matdata/livestockdensity.mat','resizecattle','resizedRcattle','resizegoats','resizedRgoats','resizesheep','resizedRsheep','livestock')
load('livestockdensity.mat') %variables: 'resizecattle','resizedRcattle','resizegoats','resizedRgoats','resizesheep','resizedRsheep')

%% the grassland above ground biomass data
%% AGB calculation
%% calculate the fANNP by using mean annual temperature to determine what fraction of 
% biomass is aboveground
load('aggregate_tas.mat')
aggregate_tas=aggregate_data;
fANPP=0.171+0.0129*aggregate_tas
imagesc(fANPP)
colorbar
%% tree cover data 
% This data is 30arcseconds, around 1km, so takes some
% time to load. Do not run this section of code once the resizing is done.
%data=readgeoraster('gm_ve_v1.tif')
% xlimits = [-90,90];
% ylimits = [-180,180]
% dataR = georefcells(xlimits,ylimits,size(data));
% [treecover,treecoverR] = georesize(data,dataR,1/12,"bilinear");
% save('treecover.mat')
load('treecover.mat')
treecoverR=georefcells([-90,90],[-180,180],size(treecover));
treecover(treecover>100) = 0;

imagesc(treecover)
colorbar
treecover = double(treecover);
treecovermultiplier = 1 ./ exp(4.45521 .* 0.01 .* treecover);
imagesc(treecovermultiplier)
colorbar
%% npp data
npp = csvread('MOD17A3H_Y_NPP_2022-01-01_gs_3600x1800.CSV', 0, 0);
%  reading starts from the second row (1 indicates the second row as MATLAB is 0-based indexing) and the first column (0 indicates the first column).
npp(npp>4000)=0
imagesc(npp)
colorbar

AGB=((npp.*fANPP)/0.5).*treecovermultiplier
AGB(AGB>1000)=0
imagesc(AGB)
colorbar
save('./grazingniche/matdata/AGB.mat','AGB')

load('livestockdensity.mat') %variables: 'resizecattle','resizedRcattle','resizegoats','resizedRgoats','resizesheep','resizedRsheep')
load('grassland_env.mat')
load('AGB.mat')
cattlegrass=resizecattle*455*0.03*365*0.001;
goatsgrass=resizegoats*91*0.03*365*0.001;
sheepgrass=resizesheep*91*0.03*365*0.001;

livestockgrass=cattlegrass+goatsgrass+sheepgrass;

livestockgrass(grassland_env==0)=0

livestockpercentage=(livestockgrass)./AGB

livestockpercentage(grassland_env<3)=0

livestockpercentage(livestockpercentage>1) = 1;
livestockpercentage(livestockpercentage<0) = 0;

save('./grazingniche/matdata/livestockpercentage.mat','livestockpercentage')


%% how many is the exceeded grassland capacity outside of the grazing niche
% get a logic layer of overgraze
overgraze_logic = (livestockpercentage >= 0.65 & livestockpercentage <= 1);
niche_logic=(landuse_coupled1>0);

% I now want to get a 1800*3600 matrix of where the value of overgraze_logic is 1 where as the niche_logic value is 0
overgraze_outsideniche_logic = overgraze_logic & ~niche_logic;
imagesc(overgraze_outsideniche_logic)

count1=sum(sum(overgraze_outsideniche_logic));
count2=sum(sum(niche_logic));
percentage_outniche=count1./count2

overgraze_insideniche_logic = overgraze_logic & niche_logic;
imagesc(overgraze_insideniche_logic)



save('./grazingniche/matdata/overgraze_niche.mat','overgraze_insideniche_logic','overgraze_outsideniche_logic','overgraze_logic')


%% grassland that are outside of the niche?

outniche = (grassland_env > 3) & (landuse_coupled1 == 0); % This creates a logical matrix
outniche = double(outniche); % Converts logicals to doubles: 1 for true, 0 for false
outniche=outniche.*grassland_env
imagesc(outniche)
colorbar

%%
save('./grazingniche/matdata/livestock2015.mat','livestockpercentage','outniche','livestock','resizecattle','allresizecattle', 'resizesheep','allresizesheep','resizegoats','allresizegoats','resizedpop');     
%% Here I am dealing with AGBC and BGBC data

load('AGBCBGBC.mat')
%%

[SOC0_30_, Rsoc_] = readgeoraster('SOCS_0_30cm_year_2010AD_10km.tif')
SOC0_30_(SOC0_30_<0)=0

[SOC0_200_, Rsoc_] = readgeoraster('SOCS_0_200cm_year_2010AD_10km.tif')
SOC0_200_(SOC0_200_<0)=0

[SOC0_30,Rsoc] = georesize(SOC0_30_,Rsoc_,3600/4320,"bilinear");
[SOC0_200,Rsoc] = georesize(SOC0_200_,Rsoc_,3600/4320,"bilinear");

imagesc(SOC0_200)
caxis([0 550]);
colorbar

save('./grazingniche/matdata/SOC.mat','SOC0_30','SOC0_200','Rsoc');     


%% carbon loss due to grazing
load('AGBCBGBC.mat')
load('SOC.mat')
load("worldarea.mat")
load('overgraze_niche.mat')
SOC0_200=double(SOC0_200);
AGBC=double(AGBC);
BGBC=double(BGBC);

carbonloss_SOC_map=(overgraze_insideniche_logic.*SOC0_200.*0.1.*grassland_env./100.*worldarea.*100)+((overgraze_outsideniche_logic.*SOC0_200.*0.2.*grassland_env./100.*worldarea.*100));
carbonloss_AGBC_map=overgraze_logic.*AGBC.*grassland_env./100.*worldarea.*100;
carbonloss_BGBC_map=overgraze_logic.*BGBC.*grassland_env./100.*worldarea.*100;

carbonloss_map=carbonloss_SOC_map+carbonloss_AGBC_map+carbonloss_BGBC_map

sum(sum(carbonloss_map))
sum(sum(carbonloss_SOC_map))
sum(sum(carbonloss_AGBC_map))
sum(sum(carbonloss_BGBC_map))


save('./grazingniche/matdata/carbonloss_map.mat','carbonloss_map','carbonloss_SOC','carbonloss_AGBC','carbonloss_BGBC');     

