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

imagesc(grassland_env)
imagesc(futureland_rcp26)
imagesc(futureniche_struct.rcp26)
grassland_env = grassland_env(:, [ceil(end/2+1):end, 1:floor(end/2)]);
futureland_rcp26 = futureland_rcp26(:, [ceil(end/2+1):end, 1:floor(end/2)]);
futureland_rcp45 = futureland_rcp45(:, [ceil(end/2+1):end, 1:floor(end/2)]);
futureland_rcp60 = futureland_rcp60(:, [ceil(end/2+1):end, 1:floor(end/2)]);
futureland_rcp85 = futureland_rcp85(:, [ceil(end/2+1):end, 1:floor(end/2)]);

%% pixel to area calculation [function]

mapsize=[160,320]

R = 6371; % Earth's radius in km
latitude_diff = 180 / mapsize(1); % in degrees
longitude_diff = 360 / mapsize(2); % in degrees, but note that it remains constant around the globe

% Convert longitude difference to radians
longitude_diff = deg2rad(longitude_diff);

% Initialize the areas matrix
areas = zeros(mapsize(1), 1);

% Compute area for each latitude band
for i = 1:mapsize(1)
    % Calculate the latitude bounds for the current band
    latitude1 = 90 - (i-1) * latitude_diff; % Start from the top and move downwards
    latitude2 = 90 - i * latitude_diff;
    
    % Convert to radians
    latitude1 = deg2rad(latitude1);
    latitude2 = deg2rad(latitude2);
    
    % Compute the area and store in the matrix
    areas(i) = R^2 * longitude_diff * (sin(latitude1) - sin(latitude2));
end

% Display the areas matrix
disp(areas);

worldarea=repmat(areas,[1,mapsize(2)])

sum(sum(worldarea))


%% now we also need to convert our continent masks into the same matrix size.
% the result is mask_futures
% Define path to your file
x=ncinfo('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc')
file_path = '/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc';
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
    [resized_data, ~] = georesize(country_data_raw, Rmask, 320/720, "bilinear");
    
    % Store in the structure
    country_data.(country_name) = resized_data';
end


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

Asiamask = Asiamask(:, [ceil(end/2+1):end, 1:floor(end/2)]);

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

Europemask = Europemask(:, [ceil(end/2+1):end, 1:floor(end/2)]);


% now for south america
SouthAmerica = {'m_ARG', 'm_BOL', 'm_BRA', 'm_CHL', 'm_COL', 'm_ECU', 'm_GUY', 'm_PRY', 'm_PER', 'm_SUR', 'm_URY', 'm_VEN'};
SouthAmericamask = country_data.m_ARG;
for i = 2:length(SouthAmerica)
    SouthAmericamask = SouthAmericamask + country_data.(SouthAmerica{i});
end
imagesc(SouthAmericamask)
colorbar

SouthAmericamask = SouthAmericamask(:, [ceil(end/2+1):end, 1:floor(end/2)]);

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

Africamask = Africamask(:, [ceil(end/2+1):end, 1:floor(end/2)]);

% Oceania
Oceania = {'m_AUS', 'm_FJI', 'm_KIR', 'm_FSM', 'm_NCL', 'm_NZL', 'm_NIU', 'm_PLW', 'm_PNG', 'm_WSM', 'm_SLB', 'm_TON', 'm_VUT'};

Oceaniamask = country_data.m_AUS;
for i = 2:length(Oceania)
    Oceaniamask = Oceaniamask + country_data.(Oceania{i});
end
imagesc(Oceaniamask)
colorbar

Oceaniamask = Oceaniamask(:, [ceil(end/2+1):end, 1:floor(end/2)]);

% convert the masks to double so that it's easier to manipulate them later
% then not only can they be used as logical layers we can directly multiply
% them with matricies, so we don't have to resort to too many conditions
% because the 0 is able to cancle any variable we want

Asiamask=double(Asiamask);
Europemask=double(Europemask);
Oceaniamask=double(Oceaniamask);
Africamask=double(Africamask);
NorthAmericamask=double(NorthAmericamask);
SouthAmericamask=double(SouthAmericamask);


%% land surface area
landsurfacearea=sum(sum(worldarea(country_data.m_world==1)));
Asiaarea=sum(sum(worldarea(Asiamask==1)));
Africaarea=sum(sum(worldarea(Africamask==1)));
SouthAmericAarea=sum(sum(worldarea(SouthAmericamask==1)));
NorthAmericaArea=sum(sum(worldarea(NorthAmericamask==1)));
Oceaniaarea=sum(sum(worldarea(Oceaniamask==1)));
Europearea=sum(sum(worldarea(Europemask==1)));

%% grassland area
% how much grassland are there in the world?
grasslandarea_rcp26=futureland_rcp26.*worldarea
grasslandarea_world_rcp26=sum(sum(grasslandarea_rcp26))

grasslandarea_rcp85=futureland_rcp85.*worldarea
grasslandarea_world_rcp85=sum(sum(grasslandarea_rcp85))

grasslandarea_rcp45 = futureland_rcp45 .* worldarea;
grasslandarea_world_rcp45 = sum(sum(grasslandarea_rcp45));

grasslandarea_rcp60 = futureland_rcp60 .* worldarea;
grasslandarea_world_rcp60 = sum(sum(grasslandarea_rcp60));

imagesc(grasslandarea_rcp26)

%how much grassland are there in each continent?
grasslandarea_Asia_rcp26=sum(sum(grasslandarea_rcp26(Asiamask==1)))
grasslandarea_Africa_rcp26=sum(sum(grasslandarea_rcp26(Africamask==1)))
grasslandarea_Oceania_rcp26=sum(sum(grasslandarea_rcp26(Oceaniamask==1)))
grasslandarea_SouthAmerica_rcp26=sum(sum(grasslandarea_rcp26(SouthAmericamask==1)))
grasslandarea_NorthAmerica_rcp26=sum(sum(grasslandarea_rcp26(NorthAmericamask==1)))
grasslandarea_Europe_rcp26=sum(sum(grasslandarea_rcp26(Europemask==1)))

grasslandarea_Asia_rcp85=sum(sum(grasslandarea_rcp85(Asiamask==1)))
grasslandarea_Africa_rcp85=sum(sum(grasslandarea_rcp85(Africamask==1)))
grasslandarea_Oceania_rcp85=sum(sum(grasslandarea_rcp85(Oceaniamask==1)))
grasslandarea_SouthAmerica_rcp85=sum(sum(grasslandarea_rcp85(SouthAmericamask==1)))
grasslandarea_NorthAmerica_rcp85=sum(sum(grasslandarea_rcp85(NorthAmericamask==1)))
grasslandarea_Europe_rcp85=sum(sum(grasslandarea_rcp85(Europemask==1)))

grasslandarea_Asia_rcp45 = sum(sum(grasslandarea_rcp45(Asiamask == 1)));
grasslandarea_Africa_rcp45 = sum(sum(grasslandarea_rcp45(Africamask == 1)));
grasslandarea_Oceania_rcp45 = sum(sum(grasslandarea_rcp45(Oceaniamask == 1)));
grasslandarea_SouthAmerica_rcp45 = sum(sum(grasslandarea_rcp45(SouthAmericamask == 1)));
grasslandarea_NorthAmerica_rcp45 = sum(sum(grasslandarea_rcp45(NorthAmericamask == 1)));
grasslandarea_Europe_rcp45 = sum(sum(grasslandarea_rcp45(Europemask == 1)));

grasslandarea_Asia_rcp60 = sum(sum(grasslandarea_rcp60(Asiamask == 1)));
grasslandarea_Africa_rcp60 = sum(sum(grasslandarea_rcp60(Africamask == 1)));
grasslandarea_Oceania_rcp60 = sum(sum(grasslandarea_rcp60(Oceaniamask == 1)));
grasslandarea_SouthAmerica_rcp60 = sum(sum(grasslandarea_rcp60(SouthAmericamask == 1)));
grasslandarea_NorthAmerica_rcp60 = sum(sum(grasslandarea_rcp60(NorthAmericamask == 1)));
grasslandarea_Europe_rcp60 = sum(sum(grasslandarea_rcp60(Europemask == 1)));

% how much percent does that take?
grasslandprc_rcp26=grasslandarea_world_rcp26/landsurfacearea;



%% how much "usable" grassland area are there in the world?

nicheland_world_rcp26=futureniche_struct.rcp26.*worldarea;
nichelandarea_world_rcp26=sum(sum(nicheland_world_rcp26));

imagesc(futureniche_struct.rcp26)
colorbar

nicheland_world_rcp85=futureniche_struct.rcp85.*worldarea;
nichelandarea_world_rcp85=sum(sum(nicheland_world_rcp85));

nicheland_world_rcp45=futureniche_struct.rcp45.*worldarea;
nichelandarea_world_rcp45=sum(sum(nicheland_world_rcp45));

nicheland_world_rcp60=futureniche_struct.rcp60.*worldarea;
nichelandarea_world_rcp60=sum(sum(nicheland_world_rcp60));

% how big is the area of usable grassland are there in each continent?
nichelandarea_Asia_rcp26=sum(sum(nicheland_world_rcp26(Asiamask==1)));
nichelandarea_Africa_rcp26=sum(sum(nicheland_world_rcp26(Africamask==1)));
nichelandarea_SouthAmerica_rcp26=sum(sum(nicheland_world_rcp26(SouthAmericamask==1)));
nichelandarea_NorthAmerica_rcp26=sum(sum(nicheland_world_rcp26(NorthAmericamask==1)));
nichelandarea_Oceania_rcp26=sum(sum(nicheland_world_rcp26(Oceaniamask==1)));
nichelandarea_Europe_rcp26=sum(sum(nicheland_world_rcp26(Europemask==1)));

nichelandarea_Asia_rcp85=sum(sum(nicheland_world_rcp85(Asiamask==1)));
nichelandarea_Africa_rcp85=sum(sum(nicheland_world_rcp85(Africamask==1)));
nichelandarea_SouthAmerica_rcp85=sum(sum(nicheland_world_rcp85(SouthAmericamask==1)));
nichelandarea_NorthAmerica_rcp85=sum(sum(nicheland_world_rcp85(NorthAmericamask==1)));
nichelandarea_Oceania_rcp85=sum(sum(nicheland_world_rcp85(Oceaniamask==1)));
nichelandarea_Europe_rcp85=sum(sum(nicheland_world_rcp85(Europemask==1)));

nichelandarea_Asia_rcp45 = sum(sum(nicheland_world_rcp45(Asiamask == 1)));
nichelandarea_Africa_rcp45 = sum(sum(nicheland_world_rcp45(Africamask == 1)));
nichelandarea_SouthAmerica_rcp45 = sum(sum(nicheland_world_rcp45(SouthAmericamask == 1)));
nichelandarea_NorthAmerica_rcp45 = sum(sum(nicheland_world_rcp45(NorthAmericamask == 1)));
nichelandarea_Oceania_rcp45 = sum(sum(nicheland_world_rcp45(Oceaniamask == 1)));
nichelandarea_Europe_rcp45 = sum(sum(nicheland_world_rcp45(Europemask == 1)));

nichelandarea_Asia_rcp60 = sum(sum(nicheland_world_rcp60(Asiamask == 1)));
nichelandarea_Africa_rcp60 = sum(sum(nicheland_world_rcp60(Africamask == 1)));
nichelandarea_SouthAmerica_rcp60 = sum(sum(nicheland_world_rcp60(SouthAmericamask == 1)));
nichelandarea_NorthAmerica_rcp60 = sum(sum(nicheland_world_rcp60(NorthAmericamask == 1)));
nichelandarea_Oceania_rcp60 = sum(sum(nicheland_world_rcp60(Oceaniamask == 1)));
nichelandarea_Europe_rcp60 = sum(sum(nicheland_world_rcp60(Europemask == 1)));

% what percentage grassland are usable in each continent?

nicheprc_Asia_rcp26=nichelandarea_Asia_rcp26/grasslandarea_Asia_rcp26;
nicheprc_Africa_rcp26=nichelandarea_Africa_rcp26/grasslandarea_Africa_rcp26;
nicheprc_SouthAmerica_rcp26=nichelandarea_SouthAmerica_rcp26/grasslandarea_SouthAmerica_rcp26;
nicheprc_NorthAmerica_rcp26=nichelandarea_NorthAmerica_rcp26/grasslandarea_NorthAmerica_rcp26;
nicheprc_Oceania_rcp26=nichelandarea_Oceania_rcp26/grasslandarea_Oceania_rcp26;
nicheprc_Europe_rcp26=nichelandarea_Europe_rcp26/grasslandarea_Europe_rcp26;
nicheprc_world_rcp26=nichelandarea_world_rcp26./grasslandarea_world_rcp26;

nicheprc_Asia_rcp85=nichelandarea_Asia_rcp85/grasslandarea_Asia_rcp85;
nicheprc_Africa_rcp85=nichelandarea_Africa_rcp85/grasslandarea_Africa_rcp85;
nicheprc_SouthAmerica_rcp85=nichelandarea_SouthAmerica_rcp85/grasslandarea_SouthAmerica_rcp85;
nicheprc_NorthAmerica_rcp85=nichelandarea_NorthAmerica_rcp85/grasslandarea_NorthAmerica_rcp85;
nicheprc_Oceania_rcp85=nichelandarea_Oceania_rcp85/grasslandarea_Oceania_rcp85;
nicheprc_Europe_rcp85=nichelandarea_Europe_rcp85/grasslandarea_Europe_rcp85;
nicheprc_world_rcp85=nichelandarea_world_rcp85./grasslandarea_world_rcp85;

nicheprc_Asia_rcp45 = nichelandarea_Asia_rcp45 / grasslandarea_Asia_rcp45;
nicheprc_Africa_rcp45 = nichelandarea_Africa_rcp45 / grasslandarea_Africa_rcp45;
nicheprc_SouthAmerica_rcp45 = nichelandarea_SouthAmerica_rcp45 / grasslandarea_SouthAmerica_rcp45;
nicheprc_NorthAmerica_rcp45 = nichelandarea_NorthAmerica_rcp45 / grasslandarea_NorthAmerica_rcp45;
nicheprc_Oceania_rcp45 = nichelandarea_Oceania_rcp45 / grasslandarea_Oceania_rcp45;
nicheprc_Europe_rcp45 = nichelandarea_Europe_rcp45 / grasslandarea_Europe_rcp45;
nicheprc_world_rcp45 = nichelandarea_world_rcp45 / grasslandarea_world_rcp45;

nicheprc_Asia_rcp60 = nichelandarea_Asia_rcp60 / grasslandarea_Asia_rcp60;
nicheprc_Africa_rcp60 = nichelandarea_Africa_rcp60 / grasslandarea_Africa_rcp60;
nicheprc_SouthAmerica_rcp60 = nichelandarea_SouthAmerica_rcp60 / grasslandarea_SouthAmerica_rcp60;
nicheprc_NorthAmerica_rcp60 = nichelandarea_NorthAmerica_rcp60 / grasslandarea_NorthAmerica_rcp60;
nicheprc_Oceania_rcp60 = nichelandarea_Oceania_rcp60 / grasslandarea_Oceania_rcp60;
nicheprc_Europe_rcp60 = nichelandarea_Europe_rcp60 / grasslandarea_Europe_rcp60;
nicheprc_world_rcp60 = nichelandarea_world_rcp60 / grasslandarea_world_rcp60;

%% Plus present climate
load('nichelandarea.mat')
figure_rcp_value=[
    nichelandarea_Asia,nichelandarea_Asia_rcp26,nichelandarea_Asia_rcp45,nichelandarea_Asia_rcp60,nichelandarea_Asia_rcp85;
    nichelandarea_Africa,nichelandarea_Africa_rcp26,nichelandarea_Africa_rcp45,nichelandarea_Africa_rcp60,nichelandarea_Africa_rcp85;
    nichelandarea_SouthAmerica,nichelandarea_SouthAmerica_rcp26,nichelandarea_SouthAmerica_rcp45,nichelandarea_SouthAmerica_rcp60,nichelandarea_SouthAmerica_rcp85;
    nichelandarea_NorthAmerica,nichelandarea_NorthAmerica_rcp26,nichelandarea_NorthAmerica_rcp45,nichelandarea_NorthAmerica_rcp60,nichelandarea_NorthAmerica_rcp85;
    nichelandarea_Oceania,nichelandarea_Oceania_rcp26,nichelandarea_Oceania_rcp45,nichelandarea_Oceania_rcp60,nichelandarea_Oceania_rcp85;
    nichelandarea_Europe,nichelandarea_Europe_rcp26,nichelandarea_Europe_rcp45,nichelandarea_Europe_rcp60,nichelandarea_Europe_rcp85]

%% save the data
save('./grazingniche/matdata/figure_rcp_value.mat','figure_rcp_value');     

%% bar
figure_rcp_caption={'Asia';'Africa';'South America';'North America';'Oceania';'Europe'}

bar(figure_rcp_value)
% Set x-axis labels
set(gca, 'xticklabel', figure_rcp_caption);
% Add a title and labels
title('Future GN area of continents');
xlabel('Region');
ylabel('Grassland area');
legend('present','rcp2.6', 'rcp4.5','rcp6.0','rcp8.5');
%% turnover rate

Unifying the resolution, unit, also unifying the coordinates 
% Unifying coordination system
halfSize = size(landuse_coupled1, 2) / 2;
landuse_coupled1_resize1 = circshift(landuse_coupled1, [0, halfSize]);
% Unifying unit (into percentage of grassland per pixel, 0-1)
landuse_coupled1_resize2=landuse_coupled1_resize1/100
% Unifying resolution
landuse_coupled1_resize3 = imresize(landuse_coupled1_resize2, [160,320], 'bilinear');
%% Test if they are unified
imagesc(landuse_coupled1_resize3) %美国在右边，值0-1，160*320
colorbar

imagesc(futureniche_struct.rcp26) %美国在右边，值0-1，160*320
colorbar

%% Turn percentage of area per pixel into area per pixel

compare_futureniche_rcp26=futureniche_struct.rcp26.*worldarea;
compare_futureniche_rcp45=futureniche_struct.rcp45.*worldarea;
compare_futureniche_rcp60=futureniche_struct.rcp60.*worldarea;
compare_futureniche_rcp85=futureniche_struct.rcp85.*worldarea;

compare_niche=landuse_coupled1_resize3.*worldarea

%% Compare
compare_turnover_rcp26=compare_futureniche_rcp26-compare_niche
compare_turnover_rcp45=compare_futureniche_rcp45-compare_niche
compare_turnover_rcp60=compare_futureniche_rcp60-compare_niche
compare_turnover_rcp85=compare_futureniche_rcp85-compare_niche

imagesc(compare_niche) 
colorbar

imagesc(compare_futureniche_rcp26)
colorbar

sum(sum(compare_niche))
sum(sum(compare_futureniche_rcp26))
sum(sum(compare_futureniche_rcp85))

%% Calculate turnover rate for each continent

turnover_total_rcp26=sum(sum(compare_turnover_rcp26));
turnover_increase_rcp26=compare_turnover_rcp26;
turnover_increase_rcp26(turnover_increase_rcp26<0)=0
turnover_increase_value_rcp26=sum(sum(turnover_increase_rcp26))

% How much will be lost?
turnover_decrease_rcp26=compare_turnover_rcp26;
turnover_decrease_rcp26(turnover_decrease_rcp26>0)=0
turnover_decrease_value_rcp26=sum(sum(turnover_decrease_rcp26))

% How much in each continent?
turnover_africa_rcp26=sum(sum(compare_turnover_rcp26(Africamask==1)));
turnover_increase_africa_rcp26=sum(sum(turnover_increase_rcp26((Africamask==1))))
turnover_decrease_africa_rcp26=sum(sum(turnover_decrease_rcp26((Africamask==1))))

%% all continents and all rcps
% Define RCPs and Continents
rcps = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
continents = {'Asia', 'Europe', 'Oceania', 'Africa', 'NorthAmerica', 'SouthAmerica'};

% Initialize structures to store results
turnover_total = struct();
turnover_increase_value = struct();
turnover_decrease_value = struct();
turnover_continents = struct();

% Loop over each RCP scenario
for i = 1:length(rcps)
    rcp = rcps{i};
    compare_turnover = eval(['compare_turnover_' rcp]); % Dynamically load the variable

    % Calculate total turnover for the current RCP
    turnover_total.(rcp) = sum(sum(compare_turnover));

    % Calculate turnover increase and decrease
    turnover_increase = compare_turnover;
    turnover_increase(turnover_increase < 0) = 0;
    turnover_increase_value.(rcp) = sum(sum(turnover_increase));

    turnover_decrease = compare_turnover;
    turnover_decrease(turnover_decrease > 0) = 0;
    turnover_decrease_value.(rcp) = sum(sum(turnover_decrease));

    % Calculate turnover for each continent
    for j = 1:length(continents)
        continent = continents{j};
        mask = eval([continent 'mask']); % Dynamically use the mask variable

        turnover_continents.(continent).(rcp).total = sum(sum(compare_turnover(mask == 1)));
        turnover_continents.(continent).(rcp).increase = sum(sum(turnover_increase(mask == 1)));
        turnover_continents.(continent).(rcp).decrease = sum(sum(turnover_decrease(mask == 1)));
    end
end

% Display results (optional)
disp(turnover_total);
disp(turnover_increase_value);
disp(turnover_decrease_value);
disp(turnover_continents);

%% Creat bar chart
%% all rcps

bar([turnover_continents.Africa.rcp26.increase,turnover_continents.Africa.rcp26.decrease,turnover_continents.Africa.rcp45.increase,turnover_continents.Africa.rcp45.decrease,turnover_continents.Africa.rcp60.increase,turnover_continents.Africa.rcp60.decrease,turnover_continents.Africa.rcp85.increase,turnover_continents.Africa.rcp85.decrease])
bar([turnover_continents.Asia.rcp26.increase,turnover_continents.Asia.rcp26.decrease,turnover_continents.Asia.rcp45.increase,turnover_continents.Asia.rcp45.decrease,turnover_continents.Asia.rcp60.increase,turnover_continents.Asia.rcp60.decrease,turnover_continents.Asia.rcp85.increase,turnover_continents.Asia.rcp85.decrease])
bar([turnover_continents.SouthAmerica.rcp26.increase,turnover_continents.SouthAmerica.rcp26.decrease,turnover_continents.SouthAmerica.rcp45.increase,turnover_continents.SouthAmerica.rcp45.decrease,turnover_continents.SouthAmerica.rcp60.increase,turnover_continents.SouthAmerica.rcp60.decrease,turnover_continents.SouthAmerica.rcp85.increase,turnover_continents.SouthAmerica.rcp85.decrease])
bar([turnover_continents.NorthAmerica.rcp26.increase,turnover_continents.NorthAmerica.rcp26.decrease,turnover_continents.NorthAmerica.rcp45.increase,turnover_continents.NorthAmerica.rcp45.decrease,turnover_continents.NorthAmerica.rcp60.increase,turnover_continents.NorthAmerica.rcp60.decrease,turnover_continents.NorthAmerica.rcp85.increase,turnover_continents.NorthAmerica.rcp85.decrease])
bar([turnover_continents.Europe.rcp26.increase,turnover_continents.Europe.rcp26.decrease,turnover_continents.Europe.rcp45.increase,turnover_continents.Europe.rcp45.decrease,turnover_continents.Europe.rcp60.increase,turnover_continents.Europe.rcp60.decrease,turnover_continents.Europe.rcp85.increase,turnover_continents.Europe.rcp85.decrease])
bar([turnover_continents.Oceania.rcp26.increase,turnover_continents.Oceania.rcp26.decrease,turnover_continents.Oceania.rcp45.increase,turnover_continents.Oceania.rcp45.decrease,turnover_continents.Oceania.rcp60.increase,turnover_continents.Oceania.rcp60.decrease,turnover_continents.Oceania.rcp85.increase,turnover_continents.Oceania.rcp85.decrease])

%% Only rcp85

% Set up the tiled layout for the subplots
figure;
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); % 2x3 layout

% Define the colors for gain and loss
colors = [0 1 0; 1 0 0]; % Green for gain, Red for loss

% Plot each continent in a separate subplot
nexttile
bar([turnover_continents.Africa.rcp85.increase, turnover_continents.Africa.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
title("Africa")

nexttile
bar([turnover_continents.Asia.rcp85.increase, turnover_continents.Asia.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
title("Asia")

nexttile
bar([turnover_continents.SouthAmerica.rcp85.increase, turnover_continents.SouthAmerica.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
title("South America")

nexttile
bar([turnover_continents.NorthAmerica.rcp85.increase, turnover_continents.NorthAmerica.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
title("North America")

nexttile
bar([turnover_continents.Oceania.rcp85.increase, turnover_continents.Oceania.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
title("Oceania")

nexttile
bar([turnover_continents.Europe.rcp85.increase, turnover_continents.Europe.rcp85.decrease], 'FaceColor', 'flat', 'CData', colors)
title("Europe")


saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/turnover_rcp85_allcontinents.svg');

%% Create maps
%% Custom colormap just for grazing turnover rate (central white, two spectrum red and green.)
num_colors = 256;  % Total number of colors for smooth transition

color_points = [ 
    0.6, 0.0, 0.0;  % Dark Red
    1.0, 1.0, 1.0;  % White
    0.1, 0.8, 0.3;  % Dark Green
];

% Interpolate to create a smooth transition colormap
reds_to_white = [linspace(color_points(1,1), color_points(2,1), num_colors/2)', ...
                 linspace(color_points(1,2), color_points(2,2), num_colors/2)', ...
                 linspace(color_points(1,3), color_points(2,3), num_colors/2)'];

white_to_greens = [linspace(color_points(2,1), color_points(3,1), num_colors/2)', ...
                   linspace(color_points(2,2), color_points(3,2), num_colors/2)', ...
                   linspace(color_points(2,3), color_points(3,3), num_colors/2)'];

% Combine the two gradients into a single colormap
custom_colormap = [reds_to_white; white_to_greens];

%% Four rcps one plot
% Define the RCP scenarios
rcps = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
titles = {'Turnover for RCP 2.6', 'Turnover for RCP 4.5', 'Turnover for RCP 6.0', 'Turnover for RCP 8.5'};

% Create a figure for the subplots
figure;
sgtitle('Turnover Across Different RCP Scenarios', 'FontSize', 16);

% Loop over each RCP scenario to create subplots
for i = 1:length(rcps)
    rcp = rcps{i};
    title_text = titles{i};
    
    % Create a subplot for each RCP scenario
    subplot(2, 2, i);
    
    % Set up the world map for the current subplot
    ax = worldmap('world');
    setm(ax, 'FFaceColor', [1 1 1]);  % Set the map's background color
    
    % Apply the custom colormap
    set(gcf, 'Colormap', custom_colormap);
    
    % Set color axis limits to ensure zero aligns with white
    caxis([-14000 14000]);  % Adjust according to your data's range
    
    % Load and display the data for the current RCP
    compare_turnover = eval(['compare_turnover_' rcp]); % Dynamically load the variable
    R = georefcells([-90, 90], [0, 360], size(compare_turnover));
    geoshow(flipud(compare_turnover), R, 'DisplayType', 'texturemap');
    
    % Load and plot coastlines
    load coastlines;
    plotm(coastlat, coastlon, 'Color', 'black');
    
    % Add a title to each subplot
    title(title_text, 'FontSize', 12);
    
    % Adjust subplot for tight and compact style
end

% Adjust the figure's size for a tight, compact display
set(gcf, 'Position', [100, 100, 1000, 600]);

% Display a colorbar that applies to all subplots
colorbar('Position', [0.93, 0.1, 0.02, 0.8]);  % Custom position for a single colorbar
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/turnover_allrcp.svg');
