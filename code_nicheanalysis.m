 addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));


%% climate data download
dir_path = '/Users/lichaohui/Desktop/calculation/grazingniche/climate/';

variables = {'pr', 'tas', 'sfcWind', 'hurs'};

for v = 1:length(variables)
    var = variables{v};
        for month = 1:12
        url = sprintf('https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly/%s/CHELSA_%s_%02d_2015_V.2.1.tif', var, var, month);
        filename = sprintf('CHELSA_%s_%02d_2015_V.2.1.tif', var, month);
        full_path = fullfile(dir_path, filename);
        websave(full_path, url);
    end
end

%% load and rescale the present climate data 

dir_path = '/Users/lichaohui/Desktop/calculation/grazingniche/climate/';
variables = {'pr', 'tas', 'sfcWind', 'hurs'};
scaling_factors = [100, 10, 1000, 100];
offsets = [0, -273.15, 0, 0];
for i = 1:length(variables)
    var = variables{i};
    scaling_factor = scaling_factors(i);
    offset = offsets(i);
    
    filename = sprintf('CHELSA_%s_01_2015_V.2.1.tif', var);
    full_path = fullfile(dir_path, filename);
    
    data = imread(full_path);
    data = imresize(data, [1800, 3600]); % Resize data
    [nrows, ncols] = size(data);
    
    data_all = zeros([nrows, ncols, 12]);
    data_all(:,:,1) = data / scaling_factor + offset; % Add the first month to the matrix
    
    for month = 2:12
        filename = sprintf('CHELSA_%s_%02d_2015_V.2.1.tif', var, month);
        full_path = fullfile(dir_path, filename);
        
        data = imread(full_path);
        data = imresize(data, [1800, 3600]); % Resize data
        data_all(:, :, month) = data / scaling_factor + offset;
    end
    
    if strcmp(var, 'pr') % If the variable is 'pr', sum up monthly precipitations
        aggregate_data = sum(data_all, 3);
    else % For all other variables, compute the mean
        aggregate_data = mean(data_all, 3);
    end
    
    % Save the variable to a .mat file
    save(['aggregate_' var '.mat'], 'aggregate_data');
end
% Load the data
if you are running the code not for the first time, you only need to

load('aggregate_hurs.mat')
aggregate_hurs=aggregate_data;
load('aggregate_pr.mat')
aggregate_pr=aggregate_data;
load('aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
load('aggregate_tas.mat')
aggregate_tas=aggregate_data;
%% historical climate data
load('final_data.mat') %variables: 'pr_final', 'tas_final', 'sfcWind_final', 'hurs_final','resizedhyde_vector_s');

middletas = tas_final1(resizedhyde_vector_s >= 1,1);
middlepr = pr_final1(resizedhyde_vector_s >= 1,1);
middlesfcWind = sfcWind_final1(resizedhyde_vector_s >= 1,1);
middlehurs = hurs_final1(resizedhyde_vector_s >= 1,1);
%
lower_bound = prctile(middletas, 2.75);
upper_bound = prctile(middletas, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlepr, 2.75);
upper_bound = prctile(middlepr, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlesfcWind, 2.75);
upper_bound = prctile(middlesfcWind, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%
lower_bound = prctile(middlehurs, 2.75);
upper_bound = prctile(middlehurs, 97.25);
fprintf('95%% of the data lies between %.2f and %.2f\n', lower_bound, upper_bound);
%% historical grazing data
fig_path = '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI';

% Create an array of years
years = 1000:100:1800;

% Loop through the array of years
for i = 1:length(years)
    % Load the data file for the current year
    filename = ['pasture', num2str(years(i)), 'AD.asc'];
    data = importdata(filename);
    
    % Extract the grid values and header information
    gridData = data.data;
    headerInfo = data.textdata;
    
    % Extract header information
    ncols = gridData(1);
    nrows = gridData(2);
    xllcorner = gridData(3);
    yllcorner = gridData(4);
    cellsize = gridData(5);
    
    % Create coordinate vectors
    x = xllcorner + (0:nrows-1) * cellsize;
    y = yllcorner + (0:nrows-1) * cellsize;
    
    % Reshape the grid data into a 2D matrix
    gridMatrix = reshape(gridData(7:size(gridData)), [], nrows)';

    % Display the image
    gridMatrix(gridMatrix < 0) = 0;
    imagesc(gridMatrix);
    colormap([1, 1, 1; jet(255)]); % White color at index 1, followed by the jet colormap
    colorbar;
    
    % Add a title
    title(['Historical Grazing Distribution: Year ', num2str(years(i))]);

    % Save the figure to a PNG file
    saveas(gcf, fullfile(fig_path, ['historical_grazing_distribution_' num2str(years(i)) '.png']), 'png');
    
    % Save the gridMatrix to a variable named gridMatrix and the current year
    assignin('base', ['gridMatrix', num2str(years(i))], gridMatrix);
end

R_hyde = georasterref('LatitudeLimits', [min(y), max(y)], 'LongitudeLimits', [min(x), max(x)], ...
        'RasterSize', [nrows, ncols], 'RasterInterpretation', 'cells');

disp(R_hyde);


%% changing the coordinates of the hyde data from -180-180 to 0-360, 
%% so that the layout of the world map is consistent with the historical climatology maps

% Create an array of years
years = 1000:100:1800;

% Loop through the array of years
for i = 1:length(years)
    % Get the matrix for the current year
    gridMatrix = eval(['gridMatrix', num2str(years(i))]);
    
    % Calculate the half size
    halfSize = size(gridMatrix, 2) / 2;
    
    % Shift the matrix
    shiftedMatrix = circshift(gridMatrix, [0, halfSize]);
    
    % Save the shifted matrix to a variable named shiftedMatrix and the current year
    assignin('base', ['shiftedMatrix', num2str(years(i))], shiftedMatrix);
end


imagesc(shiftedMatrix1800);
colormap([1, 1, 1; jet(255)]); % White color at index 1, followed by the jet colormap
colorbar;  
% alright, the graph shows the conversion is successful.


R_hyde = georasterref('LatitudeLimits', [min(y), max(y)], 'LongitudeLimits', [0,360], ...
        'RasterSize', [nrows, ncols], 'RasterInterpretation', 'cells');

disp(R_hyde)

%% Rescale data
numfile=15
resizedpr = cell(numfile,1);
resizedtas = cell(numfile,1);
resizedhurs = cell(numfile,1);
resizedsfcWind = cell(numfile,1);
resizedhyde = cell(numfile,1);

for i = 1:15
    resizedpr{i} = imresize(eval(['mean_pr_', num2str(i)]), [90 180], 'nearest');
    resizedtas{i} = imresize(eval(['mean_tas_', num2str(i)]), [90 180], 'nearest');
    resizedhurs{i} = imresize(eval(['mean_hurs_', num2str(i)]), [90 180], 'nearest');
    resizedsfcWind{i} = imresize(eval(['mean_sfcWind_', num2str(i)]), [90 180], 'nearest');

end

for i = 1:8
    year = 1000 + (i-1) * 100;
    resizedhyde{i} = imresize(eval(['shiftedMatrix', num2str(year)]), [90, 180]);
end

%% present land use data
shrubs=readgeoraster('consensus_full_class_5.tif')
(save('shrubs.mat', 'shrubs'))
herbaceous=readgeoraster('consensus_full_class_6.tif')
(save('herbaceous.mat', 'herbaceous'))
cultivated=readgeoraster('consensus_full_class_7.tif')
(save('cultivated.mat', 'cultivated'))
builtup=readgeoraster('consensus_full_class_9.tif')
(save('builtup.mat', 'builtup'))
barren=readgeoraster('consensus_full_class_11.tif')
(save('barren.mat', 'barren'))
%
load('shrubs.mat')
load('herbaceous.mat')
load('cultivated.mat')
load('builtup.mat')
load('barren.mat')

herbaceous=double(herbaceous);
shrubs=double(shrubs);
cultivated=double(cultivated);
builtup=double(builtup);
barren=double(barren);

gap=zeros(4080,43200)

herbaceous_full=[herbaceous;gap];
save('herbaceous_full.mat', 'herbaceous_full', '-v7.3');
shrubs_full=[shrubs;gap];
save('shrubs_full.mat', 'shrubs_full', '-v7.3');
cultivated_full=[cultivated;gap];
save('cultivated_full.mat', 'cultivated_full', '-v7.3');
builtup_full=[builtup;gap];
save('builtup_full.mat', 'builtup_full', '-v7.3');
barren_full=[barren;gap];
save('barren_full.mat', 'barren_full', '-v7.3');
%% load the land use data
load('herbaceous_full.mat')
load('shrubs_full.mat')
%% resize the loaded grassland maps (herbaceous+shrubs)
R_env = georefcells([-90,90],[-180,180],size(herbaceous_full));
[herbaceous_resized,R_resized] = georesize(herbaceous_full,R_env,1/12,"bilinear");
[shrubs_resized,R_resized] = georesize(shrubs_full,R_env,1/12,"bilinear");

%% obtaining a grassland map called grassland_env
grassland_env=shrubs_resized+herbaceous_resized;
save('./grazingniche/matdata/grassland_env.mat','grassland_env');

%% grassland map
figure1 = figure('WindowState','fullscreen');
sgtitle('Grassland map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(332)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[-180,180],size(grassland_env));
geoshow(flipud(grassland_env), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/grassland.svg');

%% graph niche
figure;
t = tiledlayout(2, 2);

% mean annual temperature
nexttile;
imagesc(aggregate_tas);
colorbar;
caxis([5 29]);
title('Mean Annual Temperature');

% precipitation
nexttile;
imagesc(aggregate_pr);
caxis([324 2805]);
colorbar;
title('Precipitation');

% near surface wind speed
nexttile;
imagesc(aggregate_sfcWind);
colorbar;
caxis([1.5 6.0]);
title('Near Surface Wind Speed');

% relative humidity
nexttile;
imagesc(aggregate_hurs);
colorbar;
caxis([47 87]);
title('Relative Humidity');

% Adjust the spacing and title of the tiled layout
t.TileSpacing = 'compact';
t.Padding = 'compact';
title(t, 'Climate Data');

% Save the figure as a PNG
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/niche_range.png');
%% couple the niche with the land use data

cond1_pr = (aggregate_pr >= 433 & aggregate_pr <= 2368);
cond1_tas = (aggregate_tas >= 4 & aggregate_tas <= 28);
cond1_sfcWind = (aggregate_sfcWind >= 1.4 & aggregate_sfcWind <= 6);
cond_hurs = (aggregate_hurs >= 46 & aggregate_hurs <= 87);
cond1_landuse = (grassland_env>0);

% Combine the conditions to create the niche map
niche1 = cond1_pr & cond1_tas & cond1_sfcWind & cond_hurs & cond1_landuse;

% Create the coupled landuse map
landuse_coupled1 = double(niche1) .* grassland_env;

% producing the graph of grazing niche
R = georefcells([-90,90],[-180,180],size(landuse_coupled1));

figure1 = figure('WindowState','fullscreen');
sgtitle('Niche map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(332)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
geoshow(flipud(landuse_coupled1), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/niche_general.svg');

save('./grazingniche/matdata/niche.mat',"landuse_coupled1")
