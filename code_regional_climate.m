
ncinfo('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-e_gn_2100-2300.nc')
lon_nc = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-e_gn_2100-2300.nc', 'lon');
lat_nc = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-e_gn_2100-2300.nc', 'lat');
[lon_nc, lat_nc] = meshgrid(lon_nc, lat_nc);

range_rcp60_ = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-GCAM-ssp460-2-1-f_gn_2015-2100.nc', 'range');
range_rcp45_ = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MESSAGE-ssp245-2-1-f_gn_2015-2100.nc', 'range');
range_rcp26_ = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-f_gn_2015-2100.nc', 'range');
range_rcp85_ = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MAGPIE-ssp585-2-1-f_gn_2015-2100.nc', 'range');


range_rcp26=range_rcp26_(:,:,86)'
range_rcp85=range_rcp85_(:,:,86)'
range_rcp45=range_rcp45_(:,:,86)'
range_rcp60=range_rcp60_(:,:,86)'


range_rcp26(range_rcp26>1)=0
range_rcp85(range_rcp85>1)=0
range_rcp45(range_rcp45>1)=0
range_rcp60(range_rcp60>1)=0

R=georefcells([-90,90],[-180,180],size(range_rcp26));
[futureland_rcp26,resizeRrange] = georesize(range_rcp26,R,160/720,"bilinear");
[futureland_rcp85,resizeRrange] = georesize(range_rcp85,R,160/720,"bilinear");
[futureland_rcp45,resizeRrange] = georesize(range_rcp45,R,160/720,"bilinear");
[futureland_rcp60,resizeRrange] = georesize(range_rcp60,R,160/720,"bilinear");

save('./grazingniche/matdata/futureland.mat','futureland_rcp26','futureland_rcp85','futureland_rcp45','futureland_rcp60')


%% MRI
dir_path = '/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/rcp/MRI';
scenarios = {'rcp26', 'rcp60','rcp45','rcp85'};
variables = {'pr', 'tas', 'hurs', 'sfcWind'};


for s = 1:length(scenarios)
    scenario = scenarios{s};
    for v = 1:length(variables)
        var = variables{v};

        filename = [var '_Amon_MRI-CGCM3_' scenario '_r1i1p1_200601-210012.nc'];

        full_path = fullfile(dir_path, filename);

        temp_data = ncread(full_path, var, [1 1 1129], [320 160 12]);

        if strcmp(var, 'tas') % 年均温度
            temp_data = rot90(mean(temp_data(:,:,end-11:end), 3)) - 273.15;
        elseif strcmp(var, 'pr') % 降水
            temp_data = rot90(sum(temp_data(:,:,end-11:end), 3)) * 86400 * 30;
        else
            temp_data = rot90(mean(temp_data(:,:,end-11:end), 3));
        end

        % 创建特定的变量名
        variable_name = ['data_' var '_' scenario 'MRI'];
        eval([variable_name ' = temp_data;']);
    end
end
%% MIROC CHEM的数据 
dir_path = '/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/rcp/MIROCCHEM';
scenarios = {'rcp26', 'rcp60','rcp45','rcp85'};
variables = {'pr', 'tas', 'hurs', 'sfcWind'};


for s = 1:length(scenarios)
    scenario = scenarios{s};
    for v = 1:length(variables)
        var = variables{v};

        % 文件名
        filename = [var '_Amon_MIROC-ESM-CHEM_' scenario '_r1i1p1_200601-210012.nc'];

        % 全路径
        full_path = fullfile(dir_path, filename);

        % 读取数据（仅最后12层）
        temp_data = ncread(full_path, var, [1 1 1128], [128 64 12]);

        % 根据变量处理数据
        if strcmp(var, 'tas') % 年均温度
            temp_data = rot90(mean(temp_data(:,:,end-11:end), 3)) - 273.15;
        elseif strcmp(var, 'pr') % 降水
            temp_data = rot90(sum(temp_data(:,:,end-11:end), 3)) * 86400 * 30;
        else
            temp_data = rot90(mean(temp_data(:,:,end-11:end), 3));
        end

        % 创建特定的变量名
        variable_name = ['data_' var '_' scenario 'CHEM'];
        eval([variable_name ' = temp_data;']);
    end
end

%% MIROC ESM的数据 
dir_path = '/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/rcp/MIROC';
scenarios = {'rcp26', 'rcp60','rcp45','rcp85'};
variables = {'pr', 'tas', 'hurs', 'sfcWind'};


for s = 1:length(scenarios)
    scenario = scenarios{s};
    for v = 1:length(variables)
        var = variables{v};

        % 文件名
        filename = [var '_Amon_MIROC-ESM_' scenario '_r1i1p1_200601-210012.nc'];

        % 全路径
        full_path = fullfile(dir_path, filename);

        % 读取数据（仅最后12层）
        temp_data = ncread(full_path, var, [1 1 1128], [128 64 12]);

        % 根据变量处理数据
        if strcmp(var, 'tas') % 年均温度
            temp_data = rot90(mean(temp_data(:,:,end-11:end), 3)) - 273.15;
        elseif strcmp(var, 'pr') % 降水
            temp_data = rot90(sum(temp_data(:,:,end-11:end), 3)) * 86400 * 30;
        else
            temp_data = rot90(mean(temp_data(:,:,end-11:end), 3));
        end

        % 创建特定的变量名
        variable_name = ['data_' var '_' scenario 'ESM'];
        eval([variable_name ' = temp_data;']);
    end
end
%% resizing the CHEM and ESM data to be the same as MRI data
% example code:
R = georefcells([-90,90],[-180,180],size(data_hurs_rcp26CHEM));


% Define scenarios, variables, and model types
scenarios = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
variables = {'pr', 'tas', 'hurs', 'sfcWind'};
model_types = {'CHEM', 'ESM'};

% Loop over each scenario, variable, and model type
for i = 1:length(scenarios)
    for j = 1:length(variables)
        for k = 1:length(model_types)
            % Construct the variable name for the data
            data_var_name = sprintf('data_%s_%s%s', variables{j}, scenarios{i}, model_types{k});

            % Check if the data variable exists in the workspace
            if evalin('base', ['exist(''', data_var_name, ''', ''var'')'])
                % Command to resize the data
                resize_cmd = sprintf('[%s, ~] = georesize(%s, R, 160/64, "bilinear");', data_var_name, data_var_name);

                % Evaluate the command in the base workspace
                evalin('base', resize_cmd);
            else
                fprintf('Data variable %s not found in workspace.\n', data_var_name);
            end
        end
    end
end
%% taking the average of the data

scenarios = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
variables = {'pr', 'tas', 'hurs', 'sfcWind'};

% Loop over each scenario and variable
for i = 1:length(scenarios)
    for j = 1:length(variables)
        % Construct the variable names for the different models
        data_var_MRI = sprintf('data_%s_%sMRI', variables{j}, scenarios{i});
        data_var_CHEM = sprintf('data_%s_%sCHEM', variables{j}, scenarios{i});
        data_var_ESM = sprintf('data_%s_%sESM', variables{j}, scenarios{i});
        data_var_avg = sprintf('data_%s_%s', variables{j}, scenarios{i});

        % Check if the data variables exist in the workspace
        if evalin('base', ['exist(''', data_var_MRI, ''', ''var'')']) && ...
           evalin('base', ['exist(''', data_var_CHEM, ''', ''var'')']) && ...
           evalin('base', ['exist(''', data_var_ESM, ''', ''var'')'])
            % Command to calculate the average
            avg_cmd = sprintf('%s = (%s + %s + %s) ./ 3;', data_var_avg, data_var_MRI, data_var_CHEM, data_var_ESM);

            % Evaluate the command in the base workspace
            evalin('base', avg_cmd);
        else
            fprintf('One or more data variables for %s %s not found in workspace.\n', variables{j}, scenarios{i});
        end
    end
end


var_titles = {'Precipitation', 'Mean Annual Temperature', 'Near Surface Wind Speed','Relative Humidity'};

% 输出路径
output_path = '/Users/lichaohui/Desktop/calculation/grazingniche/figures/';
data=data_hurs_rcp26;
% 定义xy轴坐标的刻度
y_ticks = [1, size(data, 1)/2, size(data, 1)];
y_tick_labels = {'90°N', '0°', '90°S'};

x_ticks = [1, size(data, 2)/2, size(data, 2)];
x_tick_labels = {'0°', '180°E', '180°W'};

for s = 1:length(scenarios)
    scenario = scenarios{s};

    % 创建一个新图形窗口
    figure('WindowState','fullscreen');

    % 创建整体标题
    sgtitle(['Scenario: ' scenario], 'FontSize', 16);

    for v = 1:length(variables)
        var = variables{v};

        % 获取特定的变量名
        variable_name = ['data_' var '_' scenario];

        % 从工作空间获取数据
        data = eval(variable_name);

        %for the precipitation data, I want the values to only show for 100-2000 the values over 2000 show the same as 2000
        % If the current variable is Precipitation, clip the values over 2000 to be 2000
if contains(var, 'pr')
    data(data > 2000) = 2000;
end

        % 创建子图
        subplot(2, 2, v);

        % 显示图像
% 显示图像
R = georefcells([-90,90],[0,360],size(data));
ax = worldmap('world');      
setm(ax, 'FFaceColor');    

if contains(var, 'pr')
    custom_colormap = [addcolorplus(338)];
elseif contains(var, 'tas')
    custom_colormap = [addcolorplus(275)];
elseif contains(var, 'hurs')
    custom_colormap = [addcolorplus(288)];
elseif contains(var, 'sfcWind')
    custom_colormap = [addcolorplus(306)];
end
colormap(gca, custom_colormap);  % Set the colormap only for the current subplot

colorbar;
geoshow(flipud(data), R, 'DisplayType', 'texturemap');   
load coastlines;
plotm(coastlat, coastlon, 'Color', 'black'); 


        % 创建子图标题
        title(var_titles{v}, 'FontSize', 14);

        % 设置xy轴刻度
        set(gca, 'YTick', y_ticks, 'YTickLabel', y_tick_labels);
        set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);

    end

    % 保存图像
    saveas(gcf, fullfile(output_path, ['figure2100_' scenario '.png']));
end

%% regional masks
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
    [resized_data, ~] = georesize(country_data_raw, Rmask, 320/720, 160/360,"bilinear");
    
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


Oceaniamask = Oceaniamask(:, [ceil(end/2+1):end, 1:floor(end/2)]);
imagesc(Oceaniamask)
colorbar
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

%% 决定niche for Africa 并graph niche
% Define the list of scenarios
scenarios = {'rcp26', 'rcp60','rcp45','rcp85'};

% Initialize a structure to store the niche maps
futureniche_struct_Africa = struct();

% Loop through each scenario
for s = 1:length(scenarios)
    % Get the scenario
    scenario = scenarios{s};

    % Get the data for the scenario
    data_pr = eval(['data_pr_' scenario]);
    data_tas = eval(['data_tas_' scenario]);
    data_sfcWind = eval(['data_sfcWind_' scenario]);
    data_hurs = eval(['data_hurs_' scenario]);
    data_landuse = eval(['futureland_' scenario]);    

data_pr_Africa = data_pr .* Africamask;
data_tas_Africa = data_tas .* Africamask;
data_hurs_Africa = data_hurs .* Africamask;
data_sfcWind_Africa = data_sfcWind .* Africamask;
%     imagesc(data_landuse) %发现数据不一致
%     imagesc(data_tas)
%     imagesc(data_landuse_new)


% 调整数据
    halfSize = size(data_landuse, 2) / 2;
   data_landuse_new = circshift(data_landuse, [0, halfSize]);
 
        % Define the conditions
    cond_pr = (data_pr >= 278 & data_pr <= 3317);
    cond_tas = (data_tas >= 16 & data_tas <= 30);
    cond_sfcWind = (data_sfcWind >= 1.35 & data_sfcWind <= 2.75);
    cond_hurs = (data_hurs >= 42.67 & data_hurs <= 90.3);
    cond_landuse = (data_landuse_new >= 0.01);


    % Combine the conditions to create the niche map
    all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;


    futureniche = data_landuse_new .* all_conditions;

    % Store the niche map in the structure
    futureniche_struct_Africa.(scenario) = futureniche;

    imagesc(futureniche_struct_Africa.rcp26)
    imagesc(data_landuse_new)
    imagesc(data_landuse)
    imagesc(all_conditions)

    % Display the niche map
subplot(2, 2, s);
R=georefcells([-90,90],[0,360],size(futureniche));    
sgtitle('Future livestock niche for Africa', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1;addcolorplus(332)];
ax = worldmap([-34, 37], [-17, 51]);
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
    title(['Niche map for ' scenario]);
    sgtitle('Figure niche of pastureland' )
end

set(gcf, 'Position',  [584,430,966,571])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/futureniche_Africa.svg');

save('./grazingniche/matdata/futureniche.mat','futureniche_struct_Africa')
%% Africa 不用 for loop循环的
%% 决定niche for Africa 并graph niche

subplot(2,2,1)  
data_pr = data_pr_rcp26;
    data_tas = data_tas_rcp26;
    data_sfcWind = data_sfcWind_rcp26;
    data_hurs = data_hurs_rcp26;
    data_landuse = futureland_rcp26;    

data_pr_Africa = data_pr .* Africamask;
data_tas_Africa = data_tas .* Africamask;
data_hurs_Africa = data_hurs .* Africamask;
data_sfcWind_Africa = data_sfcWind .* Africamask;

    cond_pr = (data_pr >= 278 & data_pr <= 3317);
    cond_tas = (data_tas >= 16 & data_tas <= 30);
    cond_sfcWind = (data_sfcWind >= 1.3 & data_sfcWind <= 2.75);
    cond_hurs = (data_hurs >= 42 & data_hurs <= 90);
    cond_landuse = (data_landuse_new >= 0.01);

all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;
futureniche = data_landuse_new .* all_conditions;

R = georefcells([-90, 90], [0, 360], size(futureniche));
% Defining a custom colormap
custom_colormap = [1 1 1; addcolorplus(332)];
% Preparing the figure
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'grassland coverage/pixel');
ax = worldmap([-34, 37], [-17, 51]);; % Adjust these values if needed for better coverage or focus
setm(ax, 'FFaceColor', [1 1 1]);
% Optionally, uncomment the next line to set the colormap if needed
% colormap(ax, custom_colormap);

% Displaying the futureniche data over Africa
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');

% Loading and plotting coastlines for reference
load coastlines
plotm(coastlat, coastlon, 'Color', 'black');

title(['rcp2.6']);

subplot(2,2,2)
    data_pr = data_pr_rcp45;
    data_tas = data_tas_rcp45;
    data_sfcWind = data_sfcWind_rcp45;
    data_hurs = data_hurs_rcp45;
    data_landuse = futureland_rcp45;    

data_pr_Africa = data_pr .* Africamask;
data_tas_Africa = data_tas .* Africamask;
data_hurs_Africa = data_hurs .* Africamask;
data_sfcWind_Africa = data_sfcWind .* Africamask;

    cond_pr = (data_pr >= 278 & data_pr <= 3317);
    cond_tas = (data_tas >= 16 & data_tas <= 30);
    cond_sfcWind = (data_sfcWind >= 1.3 & data_sfcWind <= 2.75);
    cond_hurs = (data_hurs >= 42 & data_hurs <= 90);
    cond_landuse = (data_landuse_new >= 0.01);

all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;
futureniche = data_landuse_new .* all_conditions;

R = georefcells([-90, 90], [0, 360], size(futureniche));
% Defining a custom colormap
custom_colormap = [1 1 1; addcolorplus(332)];
% Preparing the figure
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'grassland coverage/pixel');
ax = worldmap([-34, 37], [-17, 51]);; % Adjust these values if needed for better coverage or focus
setm(ax, 'FFaceColor', [1 1 1]);
% Optionally, uncomment the next line to set the colormap if needed
% colormap(ax, custom_colormap);

% Displaying the futureniche data over Africa
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');

% Loading and plotting coastlines for reference
load coastlines
plotm(coastlat, coastlon, 'Color', 'black');

title(['rcp4.5']);

subplot(2,2,3)
    data_pr = data_pr_rcp60;
    data_tas = data_tas_rcp60;
    data_sfcWind = data_sfcWind_rcp60;
    data_hurs = data_hurs_rcp60;
    data_landuse = futureland_rcp60;    

data_pr_Africa = data_pr .* Africamask;
data_tas_Africa = data_tas .* Africamask;
data_hurs_Africa = data_hurs .* Africamask;
data_sfcWind_Africa = data_sfcWind .* Africamask;

     cond_pr = (data_pr >= 278 & data_pr <= 3317);
    cond_tas = (data_tas >= 16 & data_tas <= 30);
    cond_sfcWind = (data_sfcWind >= 1.3 & data_sfcWind <= 2.75);
    cond_hurs = (data_hurs >= 42 & data_hurs <= 90);
    cond_landuse = (data_landuse_new >= 0.01);

all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;
futureniche = data_landuse_new .* all_conditions;

R = georefcells([-90, 90], [0, 360], size(futureniche));
% Defining a custom colormap
custom_colormap = [1 1 1; addcolorplus(332)];
% Preparing the figure
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'grassland coverage/pixel');
ax = worldmap([-34, 37], [-17, 51]);; % Adjust these values if needed for better coverage or focus
setm(ax, 'FFaceColor', [1 1 1]);
% Optionally, uncomment the next line to set the colormap if needed
% colormap(ax, custom_colormap);

% Displaying the futureniche data over Africa
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');

% Loading and plotting coastlines for reference
load coastlines
plotm(coastlat, coastlon, 'Color', 'black');
title(['rcp6.0']);


subplot(2,2,4)
    data_pr = data_pr_rcp85;
    data_tas = data_tas_rcp85;
    data_sfcWind = data_sfcWind_rcp85;
    data_hurs = data_hurs_rcp85;
    data_landuse = futureland_rcp85;    

data_pr_Africa = data_pr .* Africamask;
data_tas_Africa = data_tas .* Africamask;
data_hurs_Africa = data_hurs .* Africamask;
data_sfcWind_Africa = data_sfcWind .* Africamask;

     cond_pr = (data_pr >= 278 & data_pr <= 3317);
    cond_tas = (data_tas >= 16 & data_tas <= 30);
    cond_sfcWind = (data_sfcWind >= 1.3 & data_sfcWind <= 2.75);
    cond_hurs = (data_hurs >= 42 & data_hurs <= 90);
    cond_landuse = (data_landuse_new >= 0.01);

all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;
futureniche = data_landuse_new .* all_conditions;

R = georefcells([-90, 90], [0, 360], size(futureniche));
% Defining a custom colormap
custom_colormap = [1 1 1; addcolorplus(332)];
% Preparing the figure
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'grassland coverage/pixel');
ax = worldmap([-34, 37], [-17, 51]);; % Adjust these values if needed for better coverage or focus
setm(ax, 'FFaceColor', [1 1 1]);
% Optionally, uncomment the next line to set the colormap if needed
% colormap(ax, custom_colormap);

% Displaying the futureniche data over Africa
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');

% Loading and plotting coastlines for reference
load coastlines
plotm(coastlat, coastlon, 'Color', 'black');
title(['rcp8.5']);

set(gcf, 'Position',  [584,430,966,571])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/futureniche_Africa.svg');


%% 决定niche for Asia 并graph niche

subplot(2,2,1)  
data_pr = data_pr_rcp26;
    data_tas = data_tas_rcp26;
    data_sfcWind = data_sfcWind_rcp26;
    data_hurs = data_hurs_rcp26;
    data_landuse = futureland_rcp26;    

data_pr_Asia = data_pr .* Asiamask;
data_tas_Asia = data_tas .* Asiamask;
data_hurs_Asia = data_hurs .* Asiamask;
data_sfcWind_Asia = data_sfcWind .* Asiamask;

    cond_pr = (data_pr >= 252 & data_pr <= 3098);
    cond_tas = (data_tas >= 2.75 & data_tas <= 21.23);
    cond_sfcWind = (data_sfcWind >= 2.16 & data_sfcWind <= 5.05);
    cond_hurs = (data_hurs >= 38.15 & data_hurs <= 83.10);
    cond_landuse = (data_landuse_new >= 0.01);

all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;
futureniche = data_landuse_new .* all_conditions;

R = georefcells([-90, 90], [0, 360], size(futureniche));
% Defining a custom colormap
custom_colormap = [1 1 1; addcolorplus(332)];
% Preparing the figure
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'grassland coverage/pixel');
ax = worldmap([1.5 81], [26 180]); % Adjust these values if needed for better coverage or focus
setm(ax, 'FFaceColor', [1 1 1]);
% Optionally, uncomment the next line to set the colormap if needed
% colormap(ax, custom_colormap);

% Displaying the futureniche data over Asia
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');

% Loading and plotting coastlines for reference
load coastlines
plotm(coastlat, coastlon, 'Color', 'black');

title(['rcp2.6']);

subplot(2,2,2)
    data_pr = data_pr_rcp45;
    data_tas = data_tas_rcp45;
    data_sfcWind = data_sfcWind_rcp45;
    data_hurs = data_hurs_rcp45;
    data_landuse = futureland_rcp45;    

data_pr_Asia = data_pr .* Asiamask;
data_tas_Asia = data_tas .* Asiamask;
data_hurs_Asia = data_hurs .* Asiamask;
data_sfcWind_Asia = data_sfcWind .* Asiamask;

    cond_pr = (data_pr >= 252 & data_pr <= 3098);
    cond_tas = (data_tas >= 2.75 & data_tas <= 21.23);
    cond_sfcWind = (data_sfcWind >= 2.16 & data_sfcWind <= 5.05);
    cond_hurs = (data_hurs >= 38.15 & data_hurs <= 83.10);
    cond_landuse = (data_landuse_new >= 0.01);

all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;
futureniche = data_landuse_new .* all_conditions;

R = georefcells([-90, 90], [0, 360], size(futureniche));
% Defining a custom colormap
custom_colormap = [1 1 1; addcolorplus(332)];
% Preparing the figure
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'grassland coverage/pixel');
ax = worldmap([1.5 81], [26 180]); % Adjust these values if needed for better coverage or focus
setm(ax, 'FFaceColor', [1 1 1]);
% Optionally, uncomment the next line to set the colormap if needed
% colormap(ax, custom_colormap);

% Displaying the futureniche data over Asia
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');

% Loading and plotting coastlines for reference
load coastlines
plotm(coastlat, coastlon, 'Color', 'black');

title(['rcp4.5']);

subplot(2,2,3)
    data_pr = data_pr_rcp60;
    data_tas = data_tas_rcp60;
    data_sfcWind = data_sfcWind_rcp60;
    data_hurs = data_hurs_rcp60;
    data_landuse = futureland_rcp60;    

data_pr_Asia = data_pr .* Asiamask;
data_tas_Asia = data_tas .* Asiamask;
data_hurs_Asia = data_hurs .* Asiamask;
data_sfcWind_Asia = data_sfcWind .* Asiamask;

    cond_pr = (data_pr >= 252 & data_pr <= 3098);
    cond_tas = (data_tas >= 2.75 & data_tas <= 21.23);
    cond_sfcWind = (data_sfcWind >= 2.16 & data_sfcWind <= 5.05);
    cond_hurs = (data_hurs >= 38.15 & data_hurs <= 83.10);
    cond_landuse = (data_landuse_new >= 0.01);

all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;
futureniche = data_landuse_new .* all_conditions;

R = georefcells([-90, 90], [0, 360], size(futureniche));
% Defining a custom colormap
custom_colormap = [1 1 1; addcolorplus(332)];
% Preparing the figure
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'grassland coverage/pixel');
ax = worldmap([1.5 81], [26 180]); % Adjust these values if needed for better coverage or focus
setm(ax, 'FFaceColor', [1 1 1]);
% Optionally, uncomment the next line to set the colormap if needed
% colormap(ax, custom_colormap);

% Displaying the futureniche data over Asia
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');

% Loading and plotting coastlines for reference
load coastlines
plotm(coastlat, coastlon, 'Color', 'black');
title(['rcp6.0']);


subplot(2,2,4)
    data_pr = data_pr_rcp85;
    data_tas = data_tas_rcp85;
    data_sfcWind = data_sfcWind_rcp85;
    data_hurs = data_hurs_rcp85;
    data_landuse = futureland_rcp85;    

data_pr_Asia = data_pr .* Asiamask;
data_tas_Asia = data_tas .* Asiamask;
data_hurs_Asia = data_hurs .* Asiamask;
data_sfcWind_Asia = data_sfcWind .* Asiamask;

    cond_pr = (data_pr >= 252 & data_pr <= 3098);
    cond_tas = (data_tas >= 2.75 & data_tas <= 21.23);
    cond_sfcWind = (data_sfcWind >= 2.16 & data_sfcWind <= 5.05);
    cond_hurs = (data_hurs >= 38.15 & data_hurs <= 83.10);
    cond_landuse = (data_landuse_new >= 0.01);

all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;
futureniche = data_landuse_new .* all_conditions;

R = georefcells([-90, 90], [0, 360], size(futureniche));
% Defining a custom colormap
custom_colormap = [1 1 1; addcolorplus(332)];
% Preparing the figure
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'grassland coverage/pixel');
ax = worldmap([1.5 81], [26 180]); % Adjust these values if needed for better coverage or focus
setm(ax, 'FFaceColor', [1 1 1]);
% Optionally, uncomment the next line to set the colormap if needed
% colormap(ax, custom_colormap);

% Displaying the futureniche data over Asia
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');

% Loading and plotting coastlines for reference
load coastlines
plotm(coastlat, coastlon, 'Color', 'black');
title(['rcp8.5']);

set(gcf, 'Position',  [584,430,966,571])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/futureniche_Asia.svg');



%% 决定niche for Europe 并graph niche

subplot(2,2,1)  
data_pr = data_pr_rcp26;
    data_tas = data_tas_rcp26;
    data_sfcWind = data_sfcWind_rcp26;
    data_hurs = data_hurs_rcp26;
    data_landuse = futureland_rcp26;    

data_pr_Europe = data_pr .* Europemask;
data_tas_Europe = data_tas .* Europemask;
data_hurs_Europe = data_hurs .* Europemask;
data_sfcWind_Europe = data_sfcWind .* Europemask;

    cond_pr = (data_pr >= 488 & data_pr <= 2235);
    cond_tas = (data_tas >= 5 & data_tas <= 15.63);
    cond_sfcWind = (data_sfcWind >= 2.2 & data_sfcWind <= 5.24);
    cond_hurs = (data_hurs >= 67 & data_hurs <= 88);
    cond_landuse = (data_landuse_new >= 0.01);

all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;
futureniche = data_landuse_new .* all_conditions;

R = georefcells([-90, 90], [0, 360], size(futureniche));
% Defining a custom colormap
custom_colormap = [1 1 1; addcolorplus(332)];
% Preparing the figure
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'grassland coverage/pixel');
ax = worldmap([35 70], [-25 40]); % Adjust these values if needed for better coverage or focus
setm(ax, 'FFaceColor', [1 1 1]);
% Optionally, uncomment the next line to set the colormap if needed
% colormap(ax, custom_colormap);

% Displaying the futureniche data over Europe
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');

% Loading and plotting coastlines for reference
load coastlines
plotm(coastlat, coastlon, 'Color', 'black');

title(['rcp2.6']);

subplot(2,2,2)
    data_pr = data_pr_rcp45;
    data_tas = data_tas_rcp45;
    data_sfcWind = data_sfcWind_rcp45;
    data_hurs = data_hurs_rcp45;
    data_landuse = futureland_rcp45;    

data_pr_Europe = data_pr .* Europemask;
data_tas_Europe = data_tas .* Europemask;
data_hurs_Europe = data_hurs .* Europemask;
data_sfcWind_Europe = data_sfcWind .* Europemask;

    cond_pr = (data_pr >= 488 & data_pr <= 2235);
    cond_tas = (data_tas >= 5 & data_tas <= 15.63);
    cond_sfcWind = (data_sfcWind >= 2.2 & data_sfcWind <= 5.24);
    cond_hurs = (data_hurs >= 67 & data_hurs <= 88);
    cond_landuse = (data_landuse_new >= 0.01);

all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;
futureniche = data_landuse_new .* all_conditions;

R = georefcells([-90, 90], [0, 360], size(futureniche));
% Defining a custom colormap
custom_colormap = [1 1 1; addcolorplus(332)];
% Preparing the figure
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'grassland coverage/pixel');
ax = worldmap([35 70], [-25 40]); % Adjust these values if needed for better coverage or focus
setm(ax, 'FFaceColor', [1 1 1]);
% Optionally, uncomment the next line to set the colormap if needed
% colormap(ax, custom_colormap);

% Displaying the futureniche data over Europe
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');

% Loading and plotting coastlines for reference
load coastlines
plotm(coastlat, coastlon, 'Color', 'black');

title(['rcp4.5']);

subplot(2,2,3)
    data_pr = data_pr_rcp60;
    data_tas = data_tas_rcp60;
    data_sfcWind = data_sfcWind_rcp60;
    data_hurs = data_hurs_rcp60;
    data_landuse = futureland_rcp60;    

data_pr_Europe = data_pr .* Europemask;
data_tas_Europe = data_tas .* Europemask;
data_hurs_Europe = data_hurs .* Europemask;
data_sfcWind_Europe = data_sfcWind .* Europemask;

    cond_pr = (data_pr >= 488 & data_pr <= 2235);
    cond_tas = (data_tas >= 5 & data_tas <= 15.63);
    cond_sfcWind = (data_sfcWind >= 2.2 & data_sfcWind <= 5.24);
    cond_hurs = (data_hurs >= 67 & data_hurs <= 88);
    cond_landuse = (data_landuse_new >= 0.01);

all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;
futureniche = data_landuse_new .* all_conditions;

R = georefcells([-90, 90], [0, 360], size(futureniche));
% Defining a custom colormap
custom_colormap = [1 1 1; addcolorplus(332)];
% Preparing the figure
set(gcf, 'Colormap', custom_colormap);
ax = worldmap([35 70], [-25 40]); % Adjust these values if needed for better coverage or focus
setm(ax, 'FFaceColor', [1 1 1]);
hc=colorbar;
title(hc,'grassland coverage/pixel');
% Optionally, uncomment the next line to set the colormap if needed
% colormap(ax, custom_colormap);

% Displaying the futureniche data over Europe
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');

% Loading and plotting coastlines for reference
load coastlines
plotm(coastlat, coastlon, 'Color', 'black');
title(['rcp6.0']);


subplot(2,2,4)
    data_pr = data_pr_rcp85;
    data_tas = data_tas_rcp85;
    data_sfcWind = data_sfcWind_rcp85;
    data_hurs = data_hurs_rcp85;
    data_landuse = futureland_rcp85;    

data_pr_Europe = data_pr .* Europemask;
data_tas_Europe = data_tas .* Europemask;
data_hurs_Europe = data_hurs .* Europemask;
data_sfcWind_Europe = data_sfcWind .* Europemask;

    cond_pr = (data_pr >= 488 & data_pr <= 2235);
    cond_tas = (data_tas >= 5 & data_tas <= 15.63);
    cond_sfcWind = (data_sfcWind >= 2.2 & data_sfcWind <= 5.24);
    cond_hurs = (data_hurs >= 67 & data_hurs <= 88);
    cond_landuse = (data_landuse_new >= 0.01);

all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;
futureniche = data_landuse_new .* all_conditions;

R = georefcells([-90, 90], [0, 360], size(futureniche));
% Defining a custom colormap
custom_colormap = [1 1 1; addcolorplus(332)];
% Preparing the figure
set(gcf, 'Colormap', custom_colormap);
hc=colorbar;
title(hc,'grassland coverage/pixel');
ax = worldmap([35 70], [-25 40]); % Adjust these values if needed for better coverage or focus
setm(ax, 'FFaceColor', [1 1 1]);
% Optionally, uncomment the next line to set the colormap if needed
% colormap(ax, custom_colormap);

% Displaying the futureniche data over Europe
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');

% Loading and plotting coastlines for reference
load coastlines
plotm(coastlat, coastlon, 'Color', 'black');
title(['rcp8.5']);

set(gcf, 'Position',  [584,430,966,571])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/futureniche_Europe.svg');



