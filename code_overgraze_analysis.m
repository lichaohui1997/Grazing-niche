addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));

%% import data
load('grassland_env.mat')
load('livestockdensity.mat') 
load('AGB.mat')
load('niche.mat') 
load('futureniche.mat')
load('futureland.mat')
load('livestock2015.mat')
load('AGBCBGBC.mat')
load('SOC.mat')

%% now for continental analysis
x=ncinfo('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc')
lon_nc = ncread('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc', 'lon');
lat_nc = ncread('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc', 'lat');
[lon_nc, lat_nc] = meshgrid(lon_nc, lat_nc);

data_world = ncread('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc', 'm_world');
data_AUS = ncread('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc', 'm_AUS');
data_CHN = ncread('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc', 'm_CHN');

Rmask=georefcells([-90,90],[-180,180],size(data_world));

[resizeworld,resizedRworld] = georesize(data_world,Rmask,5,"bilinear");
[resizeAUS,resizedRAUS] = georesize(data_AUS,Rmask,5,"bilinear");
[resizeCHN,resizedRworld] = georesize(data_CHN,Rmask,5,"bilinear");


% Define path to your file
file_path = '/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc';

% countries = {'m_AFG', 'm_ALB', 'm_DZA', 'm_AND', 'm_AGO', 'm_ATG', 'm_ARG', 'm_ARM', ...
%              'm_AUS', 'm_AUT', 'm_AZE', 'm_BHS', 'm_BHR', 'm_BGD', 'm_BRB', 'm_BLR', ...
%              'm_BEL', 'm_BLZ', 'm_BEN', 'm_BTN', 'm_BOL', 'm_BIH', 'm_BWA', 'm_BRA', ...
%              'm_BRN', 'm_BGR', 'm_BFA', 'm_BDI', 'm_KHM', 'm_CMR', 'm_CAN', 'm_CPV', ...
%              'm_CSID', 'm_CYM', 'm_CAF', 'm_TCD', 'm_CHL', 'm_CHN', 'm_COL', 'm_COM', ...
%              'm_COG', 'm_CRI', 'm_HRV', 'm_CUB', 'm_CYP', 'm_CZE', 'm_CIV', 'm_PRK', ...
%              'm_COD', 'm_DNK', 'm_DJI', 'm_DMA', 'm_DOM', 'm_ECU', 'm_EGY', 'm_SLV', ...
%              'm_GNQ', 'm_ERI', 'm_EST', 'm_ETH', 'm_FLK', 'm_FRO', 'm_FJI', 'm_FIN', ...
%              'm_FRA', 'm_GUF', 'm_PYF', 'm_ATF', 'm_GAB', 'm_GMB', 'm_GEO', 'm_DEU', ...
%              'm_GHA', 'm_GRC', 'm_GRL', 'm_GRD', 'm_GLP', 'm_GUM', 'm_GTM', 'm_GIN', ...
%              'm_GNB', 'm_GUY', 'm_HTI', 'm_HMD', 'm_HND', 'm_HKG', 'm_HUN', 'm_ISL', ...
%              'm_IND', 'm_IOSID', 'm_IDN', 'm_IRN', 'm_IRQ', 'm_IRL', 'm_IMN', 'm_ISR', ...
%              'm_ITA', 'm_JAM', 'm_JKX', 'm_JPN', 'm_JOR', 'm_KAZ', 'm_KEN', 'm_KIR', ...
%              'm_KWT', 'm_KGZ', 'm_LAO', 'm_LVA', 'm_LBN', 'm_LSO', 'm_LBR', 'm_LBY', ...
%              'm_LTU', 'm_LUX', 'm_MDG', 'm_MWI', 'm_MYS', 'm_MLI', 'm_MLT', 'm_MTQ', ...
%              'm_MRT', 'm_MUS', 'm_MYT', 'm_MEX', 'm_FSM', 'm_MDA', 'm_MNG', 'm_MNE', ...
%              'm_MAR', 'm_MOZ', 'm_MMR', 'm_NAM', 'm_NPL', 'm_NLD', 'm_ANT', 'm_NCL', ...
%              'm_NZL', 'm_NIC', 'm_NER', 'm_NGA', 'm_NIU', 'm_NOR', 'm_OMN', 'm_PSID', ...
%              'm_PAK', 'm_PLW', 'm_PSE', 'm_PAN', 'm_PNG', 'm_PRY', 'm_PER', 'm_PHL', ...
%              'm_POL', 'm_PRT', 'm_PRI', 'm_QAT', 'm_KOR', 'm_ROU', 'm_RUS', 'm_RWA', ...
%              'm_REU', 'm_LCA', 'm_SPM', 'm_VCT', 'm_WSM', 'm_STP', 'm_SAU', 'm_SEN', ...
%              'm_SRB', 'm_SLE', 'm_SGP', 'm_SVK', 'm_SVN', 'm_SLB', 'm_SOM', 'm_ZAF', ...
%              'm_SGS', 'm_SSD', 'm_ESP', 'm_LKA', 'm_SDN', 'm_SUR', 'm_SJM', 'm_SWZ', ...
%              'm_SWE', 'm_CHE', 'm_SYR', 'm_TWN', 'm_TJK', 'm_THA', 'm_MKD', 'm_TLS', ...
%              'm_TGO', 'm_TON', 'm_TTO', 'm_TUN', 'm_TUR', 'm_TKM', 'm_GBR', 'm_UGA', ...
%              'm_UKR', 'm_ARE', 'm_TZA', 'm_VIR', 'm_USA', 'm_URY', 'm_UZB', 'm_VUT', ...
%              'm_VEN', 'm_VNM', 'm_ESH', 'm_YEM', 'm_ZMB', 'm_ZWE', 'm_world'};

countries = {x.Variables(3:217).Name};

% Create a structure to hold all country data
country_data = struct();

% Read lat and lon, and create a meshgrid
lon = ncread(file_path, 'lon');
lat = ncread(file_path, 'lat');
[lon, lat] = meshgrid(lon, lat);
Rmask = georefcells([-90,90],[-180,180], size(ncread(file_path, 'm_world')));

% Loop through countries, read and resize each mask
for i = 1:length(countries)
    country_name = countries{i};
    
    % Read country data
    country_data_raw = ncread(file_path, country_name);
    
    % Resize the country data
    [resized_data, ~] = georesize(country_data_raw, Rmask, 5, "bilinear");
    
    % Store in the structure
    country_data.(country_name) = resized_data';
end


Asia={'m_CHN','m_JPN','m_KOR'}
Asiamask=country_data.m_CHN+country_data.m_JPN+country_data.m_KOR
imagesc(Asiamask)
colorbar

% mask for continents
% asia
Asia = {'m_AFG', 'm_ARM', 'm_AZE', 'm_BHR', 'm_BGD', 'm_BTN', 'm_BRN', 'm_KHM', 'm_CHN', 'm_CYM',... 
        'm_GEO', 'm_IND', 'm_IDN', 'm_IRN', 'm_IRQ', 'm_ISR', 'm_JPN', 'm_JOR', 'm_KAZ', 'm_KWT', ...
        'm_KGZ', 'm_LAO', 'm_LBN', 'm_MYS', 'm_MNG', 'm_MMR', 'm_NPL', 'm_PRK', 'm_OMN', ...
        'm_PAK', 'm_PHL', 'm_QAT', 'm_RUS', 'm_SAU', 'm_SGP', 'm_KOR', 'm_LKA', 'm_SYR', 'm_TWN', ...
        'm_TJK', 'm_THA', 'm_TUR', 'm_TKM', 'm_ARE', 'm_UZB', 'm_VNM', 'm_YEM'};

Asiamask = country_data.m_AFG;
for i = 2:length(Asia)
    Asiamask = Asiamask + country_data.(Asia{i});
end

imagesc(Asiamask)
colorbar
% for europe
Europe = {'m_ALB', 'm_AND', 'm_AUT', 'm_BEL', 'm_BIH', 'm_BGR', 'm_HRV', 'm_CYP', 'm_CZE', 'm_DNK',... 
          'm_EST', 'm_FIN', 'm_FRA', 'm_DEU', 'm_GRC', 'm_HUN', 'm_ISL', 'm_IRL', 'm_ITA', ...
          'm_LVA', 'm_LTU', 'm_LUX', 'm_MLT', 'm_MDA', 'm_MNE', 'm_NLD', 'm_MKD', ...
          'm_NOR', 'm_POL', 'm_PRT', 'm_ROU', 'm_SRB', 'm_SVK', 'm_SVN', 'm_ESP', ...
          'm_SWE', 'm_CHE', 'm_UKR', 'm_GBR'};

Europemask = country_data.m_ALB;
for i = 2:length(Europe)
    Europemask = Europemask + country_data.(Europe{i});
end

imagesc(Europemask)
colorbar

% now for south america
SouthAmerica = {'m_ARG', 'm_BOL', 'm_BRA', 'm_CHL', 'm_COL', 'm_ECU', 'm_GUY', 'm_PRY', 'm_PER', 'm_SUR', 'm_URY', 'm_VEN'};
SouthAmericamask = country_data.m_ARG;
for i = 2:length(SouthAmerica)
    SouthAmericamask = SouthAmericamask + country_data.(SouthAmerica{i});
end
imagesc(SouthAmericamask)
colorbar

% north america
NorthAmerica = {'m_ATG', 'm_BHS', 'm_BRB', 'm_BEL', 'm_CAN', 'm_CYM', 'm_CRI', 'm_CUB', 'm_DMA', 'm_DOM', 'm_SLV', 'm_GRL', 'm_GRD', 'm_GLP', 'm_GTM', 'm_HND', 'm_JAM', 'm_MEX', 'm_MTQ', 'm_NIC', 'm_PAN', 'm_PRI', 'm_SPM', 'm_LCA', 'm_VCT', 'm_TTO', 'm_USA', 'm_VIR'};
NorthAmericamask = country_data.m_ATG;
for i = 2:length(NorthAmerica)
    NorthAmericamask = NorthAmericamask + country_data.(NorthAmerica{i});
end
imagesc(NorthAmericamask)
colorbar
% africa
Africa = {'m_DZA', 'm_COD','m_SDN', 'm_AGO', 'm_BEN', 'm_BWA', 'm_BFA', 'm_BDI', 'm_CMR', 'm_CPV', 'm_CAF', 'm_TCD', 'm_COM', 'm_COG', 'm_CIV', 'm_DJI', 'm_EGY', 'm_GNQ', 'm_ERI', 'm_ETH', 'm_GAB', 'm_GMB', 'm_GHA', 'm_GIN', 'm_GNB', 'm_KEN', 'm_LSO', 'm_LBR', 'm_LBY', 'm_MDG', 'm_MWI', 'm_MLI', 'm_MRT', 'm_MUS', 'm_MAR', 'm_MOZ', 'm_NAM', 'm_NER', 'm_NGA', 'm_REU', 'm_RWA', 'm_STP', 'm_SEN', 'm_SLE', 'm_SOM', 'm_ZAF', 'm_SSD', 'm_ESH', 'm_SWZ', 'm_TZA', 'm_TGO', 'm_TUN', 'm_UGA', 'm_ZMB', 'm_ZWE'};
Africamask = country_data.m_DZA;
for i = 2:length(Africa)
    Africamask = Africamask + country_data.(Africa{i});
end
imagesc(Africamask)
colorbar

% Oceania
Oceania = {'m_AUS', 'm_FJI', 'm_KIR', 'm_FSM', 'm_NCL', 'm_NZL', 'm_NIU', 'm_PLW', 'm_PNG', 'm_WSM', 'm_SLB', 'm_TON', 'm_VUT'};

Oceaniamask = country_data.m_AUS;
for i = 2:length(Oceania)
    Oceaniamask = Oceaniamask + country_data.(Oceania{i});
end
imagesc(Oceaniamask)
colorbar


Asiamask=double(Asiamask);
Europemask=double(Europemask);
Oceaniamask=double(Oceaniamask);
Africamask=double(Africamask);
NorthAmericamask=double(NorthAmericamask);
SouthAmericamask=double(SouthAmericamask);

save('./grazingniche/matdata/continentmask.mat','Asiamask','Europemask','Oceaniamask','Africamask','NorthAmericamask','SouthAmericamask')

%% calculating area

R = 6371; % Earth's radius in km
delta_lambda = 0.1 * pi / 180; % change in longitude in radians
phi1 = -pi/2; % starting from -90 degrees
phi2 = phi1 + 0.1 * pi / 180; % add 0.1 degree in radians

area = R^2 * delta_lambda * (sin(phi2) - sin(phi1));
disp(['The area of the (1,1) grid cell is approximately ', num2str(area), ' km^2']);
%% calculate all of the area for all the latitude. The result a matrix of 1800*1 values

R = 6371; % Earth's radius in km
delta_lambda = 0.1 * pi / 180; % change in longitude in radians

areas = zeros(1800, 1);

for i = 1:1800
    phi1 = (-pi/2) + (i-1) * 0.1 * pi / 180; % starting latitude for this step
    phi2 = phi1 + 0.1 * pi / 180; % ending latitude for this step

    areas(i) = R^2 * delta_lambda * (sin(phi2) - sin(phi1));
end

disp(areas);

worldarea=repmat(areas,[1,3600])

save('./grazingniche/matdata/worldarea.mat','worldarea');     
%%
earthsurface=sum(sum(worldarea))
% actually very close to the googled value: 509 600 000 square km
%% land surface area
landsurfacearea=sum(sum(worldarea(resizeworld==1)));
Asiaarea=sum(sum(worldarea(Asiamask==1)));
Africaarea=sum(sum(worldarea(Africamask==1)));
SouthAmericAarea=sum(sum(worldarea(SouthAmericamask==1)));
NorthAmericaArea=sum(sum(worldarea(NorthAmericamask==1)));
Oceaniaarea=sum(sum(worldarea(Oceaniamask==1)));
Europearea=sum(sum(worldarea(Europemask==1)));
Landarea=[Asiaarea;Africaarea;SouthAmericAarea;NorthAmericaArea;Oceaniaarea;Europearea];
%% grassland area
grasslandarea=grassland_env.*worldarea./100
grasslandarea_world=sum(sum(grasslandarea))

grasslandarea_Asia=sum(sum(grasslandarea(Asiamask==1)))
grasslandarea_Africa=sum(sum(grasslandarea(Africamask==1)))
grasslandarea_Oceania=sum(sum(grasslandarea(Oceaniamask==1)))
grasslandarea_SouthAmerica=sum(sum(grasslandarea(SouthAmericamask==1)))
grasslandarea_NorthAmerica=sum(sum(grasslandarea(NorthAmericamask==1)))
grasslandarea_Europe=sum(sum(grasslandarea(Europemask==1)))

grasslandarea=[grasslandarea_Asia;grasslandarea_Africa;grasslandarea_SouthAmerica;grasslandarea_NorthAmerica;grasslandarea_Oceania;grasslandarea_Europe]

% how much percent does that take?
grasslandprc=grasslandarea_world/landsurfacearea;
grasslandprc_Asia=grasslandarea_Asia/Asiaarea;
grasslandprc_Oceania=grasslandarea_Oceania/Oceaniaarea;
grasslandprc_NorthAmerica=grasslandarea_NorthAmerica/NorthAmericaArea;
grasslandprc_SouthAmerica=grasslandarea_SouthAmerica/SouthAmericAarea;
grasslandprc_Africa=grasslandarea_Africa/Africaarea;
grasslandprc_Europe=grasslandarea_Europe/Europearea;

%% how much usable grassland area are there in the world?

nicheland_world=landuse_coupled1.*worldarea./100;
nichelandarea_world=sum(sum(nicheland_world));

nichelandarea_world_prct=nichelandarea_world./landsurfacearea

% how much usable grassland are there in each continent?
nichelandarea_Asia=sum(sum(nicheland_world(Asiamask==1)));
nichelandarea_Africa=sum(sum(nicheland_world(Africamask==1)));
nichelandarea_SouthAmerica=sum(sum(nicheland_world(SouthAmericamask==1)));
nichelandarea_NorthAmerica=sum(sum(nicheland_world(NorthAmericamask==1)));
nichelandarea_Oceania=sum(sum(nicheland_world(Oceaniamask==1)));
nichelandarea_Europe=sum(sum(nicheland_world(Europemask==1)));

save('./grazingniche/matdata/nichelandarea.mat','nichelandarea_Asia','nichelandarea_Africa','nichelandarea_SouthAmerica','nichelandarea_NorthAmerica','nichelandarea_Oceania','nichelandarea_Europe');     


nichelandarea=[nichelandarea_Asia;nichelandarea_Africa;nichelandarea_SouthAmerica;nichelandarea_NorthAmerica;nichelandarea_Oceania;nichelandarea_Europe]
nonnichelandarea=grasslandarea-nichelandarea

figure_nichelandarea_value=[nichelandarea,nonnichelandarea];

save('./grazingniche/matdata/figure_nichelandarea_value.mat','figure_nichelandarea_value');     

%% bar figure of niche area 
figure_nichelandarea_caption={'Asia';'Africa';'South America';'North America';'Oceania';'Europe'}

bar(figure_nichelandarea_value);
% Set x-axis labels
set(gca, 'xticklabel', figure_nichelandarea_caption);
% Add a title and labels
title('GN grassland vs Non-GN grassland');
xlabel('Region');
ylabel('Grassland area');
legend('GN area', 'Non-GN area');%% how much of the unusable grassland are overgrazed?

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/nichelandareabar.svg');


% how much grassland are actually "usable" in each continent?

nicheprc_Asia=nichelandarea_Asia/grasslandarea_Asia;
nicheprc_Africa=nichelandarea_Africa/grasslandarea_Africa;
nicheprc_SouthAmerica=nichelandarea_SouthAmerica/grasslandarea_SouthAmerica;
nicheprc_NorthAmerica=nichelandarea_NorthAmerica/grasslandarea_NorthAmerica;
nicheprc_Oceania=nichelandarea_Oceania/grasslandarea_Oceania;
nicheprc_Europe=nichelandarea_Europe/grasslandarea_Europe;
nicheprc_world=nichelandarea_world./grasslandarea_world;


%% how much grassland is overgrazed? 

overgrazed_logic=double(livestockpercentage >= 0.65 & livestockpercentage <= 1)
mediumgrazed_logic=double(livestockpercentage >= 0.20 & livestockpercentage <= 0.65)
lightgrazed_logic=double(livestockpercentage >= 0 & livestockpercentage <= 0.20)

imagesc(overgrazed_logic)
colorbar


% how much of world's grassland is overgrazed
overgraze_world = double(livestockpercentage >= 0.65);
overgraze_world_area=sum(sum(overgraze_world.*grassland_env.*worldarea/100));
percentage_overgraze=overgraze_world_area./grasslandarea_world;
percentage_overgraze_niche=overgraze_world_area./nichelandarea_world;


% how much grassland is medium-grazed? 
mediumgraze_world = double(livestockpercentage >= 0.20 & livestockpercentage < 0.65);
percentage_medium_area=sum(sum(mediumgraze_world.*grassland_env.*worldarea/100));
percentage_mediumgraze=percentage_medium_area./grasslandarea_world;


% how much grassland is lightly-grazed
lightgraze_world = double(livestockpercentage >= 0 & livestockpercentage < 0.01);
percentage_lightgraze_area=sum(sum(lightgraze_world.*grassland_env.*worldarea/100));
percentage_lightgraze=percentage_lightgraze_area./grasslandarea_world;

land_overgraze=percentage_overgraze.*grasslandarea_world;
land_mediumgraze=percentage_mediumgraze.*grasslandarea_world;
land_lightgraze=percentage_lightgraze.*grasslandarea_world;

%% where are the "empty" niches?

imagesc(lightgraze_world)
colorbar

empty_niche=lightgraze_world.*worldarea.*landuse_coupled1/100

imagesc(empty_niche)
colorbar
save('./grazingniche/matdata/empty_niche.mat','empty_niche')


%% how much grassland is overgrazed in each continent?

overgrazed_Asia = zeros(size(Asiamask)); % Initialize the matrix with zeros, same size as Asiamask
overgrazed_Asia(Asiamask == 1) = overgrazed_logic(Asiamask == 1);

overgrazed_Europe = zeros(size(Europemask)); % Initialize the matrix with zeros, same size as Europemask
overgrazed_Europe(Europemask == 1) = overgrazed_logic(Europemask == 1);

overgrazed_NorthAmerica = zeros(size(NorthAmericamask)); % Initialize the matrix with zeros, same size as NorthAmericamask
overgrazed_NorthAmerica(NorthAmericamask == 1) = overgrazed_logic(NorthAmericamask == 1);

overgrazed_SouthAmerica = zeros(size(SouthAmericamask)); % Initialize the matrix with zeros, same size as SouthAmericamask
overgrazed_SouthAmerica(SouthAmericamask == 1) = overgrazed_logic(SouthAmericamask == 1);

overgrazed_Africa = zeros(size(Africamask)); % Initialize the matrix with zeros, same size as Africamask
overgrazed_Africa(Africamask == 1) = overgrazed_logic(Africamask == 1);

overgrazed_Oceania = zeros(size(Oceaniamask)); % Initialize the matrix with zeros, same size as Oceaniamask
overgrazed_Oceania(Oceaniamask == 1) = overgrazed_logic(Oceaniamask == 1);

imagesc(overgrazed_SouthAmerica)

% what is the overgrazed grassland area in each continent?
overgrazedarea_Asia=sum(sum(overgrazed_Asia.*grassland_env.*worldarea./100))
overgrazedarea_Europe=sum(sum(overgrazed_Europe.*grassland_env.*worldarea./100))
overgrazedarea_Oceania=sum(sum(overgrazed_Oceania.*grassland_env.*worldarea./100))
overgrazedarea_NorthAmerica=sum(sum(overgrazed_NorthAmerica.*grassland_env.*worldarea./100))
overgrazedarea_SouthAmerica=sum(sum(overgrazed_SouthAmerica.*grassland_env.*worldarea./100))
overgrazedarea_Africa=sum(sum(overgrazed_Africa.*grassland_env.*worldarea./100))

% how much of each continents' grassland are overgrazed?
overgrazedprc_Asia=overgrazedarea_Asia/grasslandarea_Asia
overgrazedprc_Europe=overgrazedarea_Europe/grasslandarea_Europe
overgrazedprc_Oceania=overgrazedarea_Oceania/grasslandarea_Oceania
overgrazedprc_NorthAmerica=overgrazedarea_NorthAmerica/grasslandarea_NorthAmerica
overgrazedprc_SouthAmerica=overgrazedarea_SouthAmerica/grasslandarea_SouthAmerica
overgrazedprc_Africa=overgrazedarea_Africa/grasslandarea_Africa

overgrazedarea=[overgrazedarea_Asia;overgrazedarea_Africa;overgrazedarea_NorthAmerica;overgrazedarea_SouthAmerica;overgrazedarea_Oceania;overgrazedarea_Europe]

overgrazedarea_world=sum(overgrazedarea)

%% how much of the world's usable and unusable grassland are overgrazed?

% overgrazed area within the niche
overgrazed_niche = grassland_env .* (outniche == 0 & overgraze_world == 1);
imagesc(overgrazed_niche)
overgrazed_niche_world=overgrazed_niche.*worldarea./100;
overgrazed_niche_area_world=sum(sum(overgrazed_niche_world));
overgrazed_niche_prc_world=overgrazed_niche_area_world./overgraze_world_area;
overgrazed_niche_prc_niche=overgrazed_niche_area_world./nichelandarea_world;

% overgrazed area outside the niche
outniche_area=sum(sum(outniche.*worldarea/100));
overgrazed_outniche = grassland_env .* (outniche ~= 0 & overgraze_world == 1);
imagesc(overgrazed_outniche)
overgrazed_outniche_world=overgrazed_outniche.*worldarea./100
overgrazed_outniche_area_world=sum(sum(overgrazed_outniche_world))
overgrazed_outniche_prc_world=overgrazed_outniche_area_world./overgraze_world_area;

figure_overgrazeGN_value=[]


%% graph of these
figure1 = figure('WindowState','fullscreen');
sgtitle('overgrazed niche grassland', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[-180,180],size(overgrazed_niche));
geoshow(flipud(overgrazed_niche), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/overgrazed_niche.svg');

figure2 = figure('WindowState','fullscreen');
sgtitle('overgrazed grassland outside niche', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[-180,180],size(overgrazed_outniche));
geoshow(flipud(overgrazed_outniche), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/overgrazed_outniche.svg');

%% how much of each continent's usable grassland are overgrazed?
niche_overgrazed_Asia=overgrazed_niche.*Asiamask;
niche_overgrazed_Europe=overgrazed_niche.*Europemask;
niche_overgrazed_SouthAmerica=overgrazed_niche.*SouthAmericamask;
niche_overgrazed_NorthAmerica=overgrazed_niche.*NorthAmericamask;
niche_overgrazed_Oceania=overgrazed_niche.*Oceaniamask;
niche_overgrazed_Africa=overgrazed_niche.*Africamask;

niche_overgrazedarea_Asia=sum(sum(niche_overgrazed_Asia.*worldarea./100))
niche_overgrazedarea_Europe=sum(sum(niche_overgrazed_Europe.*worldarea./100))
niche_overgrazedarea_SouthAmerica=sum(sum(niche_overgrazed_SouthAmerica.*worldarea./100))
niche_overgrazedarea_NorthAmerica=sum(sum(niche_overgrazed_NorthAmerica.*worldarea./100))
niche_overgrazedarea_Oceania=sum(sum(niche_overgrazed_Oceania.*worldarea./100))
niche_overgrazedarea_Africa=sum(sum(niche_overgrazed_Africa.*worldarea./100))

niche_overgrazedprc_Asia=niche_overgrazedarea_Asia./nichelandarea_Asia;
niche_overgrazedprc_Europe=niche_overgrazedarea_Europe./nichelandarea_Europe;
niche_overgrazedprc_SouthAmerica=niche_overgrazedarea_SouthAmerica./nichelandarea_SouthAmerica;
niche_overgrazedprc_NorthAmerica=niche_overgrazedarea_NorthAmerica./nichelandarea_NorthAmerica;
niche_overgrazedprc_Oceania=niche_overgrazedarea_Oceania./nichelandarea_Oceania;
niche_overgrazedprc_Africa=niche_overgrazedarea_Africa./nichelandarea_Africa;


niche_overgrazedarea=[niche_overgrazedarea_Asia;niche_overgrazedarea_Africa;niche_overgrazedarea_NorthAmerica;niche_overgrazedarea_SouthAmerica;niche_overgrazedarea_Oceania;niche_overgrazedarea_Europe]
nonniche_overgrazedarea=overgrazedarea-niche_overgrazedarea

figure_overgrazeGN_value=[niche_overgrazedarea,nonniche_overgrazedarea]

%% save the data
save('./grazingniche/matdata/figure_overgrazeGN_value.mat','figure_overgrazeGN_value')
%% bar chart of overgraze in and out of niche
figure_overgrazedGN_caption={'Asia';'Africa';'South America';'North America';'Oceania';'Europe'}
% Create the bar chart
bar(figure_overgrazeGN_value);
% Set x-axis labels
set(gca, 'xticklabel', figure_overgrazedGN_caption);
% Add a title and labels
title('Overgrazed Area in GN vs. Outside GN');
xlabel('Region');
ylabel('Overgrazed Area');
legend('Overgrazed Area in GN', 'Overgrazed Area outside GN');%% how much of the unusable grassland are overgrazed?
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/overgrazeGNbar.svg');


%% present grassland livestock distribution has concentrated in what percentage of the niche?

imagesc(livestock)
colorbar

% how many grassland cattle in the world?
cattlenum=sum(sum(resizecattle.*worldarea));
cattlenumall=sum(sum(allresizecattle.*worldarea));

goatsnum=sum(sum(resizegoats.*worldarea));
goatsnumall=sum(sum(allresizegoats.*worldarea));

sheepsnum=sum(sum(resizesheep.*worldarea));
sheepsnumall=sum(sum(allresizesheep.*worldarea));

%% how many are these animals inside the niche and how many are outside the niche
% I want the variable cattle_niche to be equal the variable resizecattle, only if landuse_coupled1>0
cattle_niche = zeros(size(resizecattle)); % Initialize cattle_niche with the same size as resizecattle
cattle_niche(landuse_coupled1 > 0) = resizecattle(landuse_coupled1 > 0);
cattle_niche_num=sum(sum(cattle_niche.*worldarea));

imagesc(landuse_coupled1)
imagesc(cattle_niche)
colorbar

cattle_outniche = zeros(size(resizecattle)); % Initialize cattle_niche with the same size as resizecattle
cattle_outniche(landuse_coupled1 == 0) = resizecattle(landuse_coupled1 == 0);
cattle_outniche_num=sum(sum(cattle_outniche.*worldarea));

sheep_niche = zeros(size(resizesheep)); % Initialize sheep_niche with the same size as resizesheep
sheep_niche(landuse_coupled1 > 0) = resizesheep(landuse_coupled1 > 0);
sheep_niche_num=sum(sum(sheep_niche.*worldarea));

sheep_outniche = zeros(size(resizesheep)); % Initialize sheep_niche with the same size as resizesheep
sheep_outniche(landuse_coupled1 == 0) = resizesheep(landuse_coupled1 == 0);
sheep_outniche_num=sum(sum(sheep_outniche.*worldarea));

goats_niche = zeros(size(resizegoats)); % Initialize goats_niche with the same size as resizegoats
goats_niche(landuse_coupled1 > 0) = resizegoats(landuse_coupled1 > 0);
goats_niche_num=sum(sum(goats_niche.*worldarea));

goats_outniche = zeros(size(resizegoats)); % Initialize goats_niche with the same size as resizegoats
goats_outniche(landuse_coupled1 == 0) = resizegoats(landuse_coupled1 == 0);
goats_outniche_num=sum(sum(goats_outniche.*worldarea));


figure_livestockniche_value=[cattle_niche_num,cattle_outniche_num;sheep_niche_num,sheep_outniche_num;goats_niche_num,goats_outniche_num];
figure_livestockniche_caption={'cattle';'sheep';'goat'}

bar(figure_livestockniche_value);
% Set x-axis labels
set(gca, 'xticklabel', figure_livestockniche_caption);
% Add a title and labels
title('Livestock distribution in GN vs. Outside GN');
xlabel('Type of livestock');
ylabel('Number of Livestock');
legend('Inside GN', 'Outside GN');%% how much of the unusable grassland are overgrazed?
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/livestocknichebar.svg');

