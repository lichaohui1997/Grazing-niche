%% In this code I am using two new climate models for historical data

load('pr_data_miroc.mat');
load('sfcWind_data_miroc.mat');
load('tas_data_miroc.mat');
load('hurs_data_miroc.mat');

years_interest = 1000:100:1800;

start_index = (years_interest - 849) * 12 + 1; % Start from January of each year
end_index = start_index + 11; % Up to December of each year

% Initialize the matrices to store results
tas_yearly_miroc = zeros(128, 64, numel(years_interest));
sfcWind_yearly_miroc = zeros(128, 64, numel(years_interest));
hurs_yearly_miroc = zeros(128, 64, numel(years_interest));
pr_yearly_miroc = zeros(128, 64, numel(years_interest));

for i = 1:numel(years_interest)
    tas_yearly_miroc(:,:,i) = mean(tas_data_miroc(:,:,start_index(i):end_index(i)), 3)-273.15;
    sfcWind_yearly_miroc(:,:,i) = mean(sfcWind_data_miroc(:,:,start_index(i):end_index(i)), 3);
    %hurs_yearly_miroc(:,:,i) = mean(hurs_data_miroc(:,:,start_index(i):end_index(i)), 3);
    capped_hurs_data = min(hurs_data_miroc(:,:,start_index(i):end_index(i)), 100);
    hurs_yearly_miroc(:,:,i) = mean(capped_hurs_data, 3);
    pr_yearly_miroc(:,:,i) = sum(pr_data_miroc(:,:,start_index(i):end_index(i)), 3)*86400*30;
end

directory = '/Users/lichaohui/Desktop/grazingniche/rawdata/climatenewnew/MIROC';
% Save the yearly data
if ~exist(directory, 'dir')
    mkdir(directory);
end

directory = '/Users/lichaohui/Desktop/grazingniche/rawdata/climatenewnew/MIROC';

save(fullfile(directory, 'tas_yearly_miroc.mat'), 'tas_yearly_miroc');
save(fullfile(directory, 'sfcWind_yearly_miroc.mat'), 'sfcWind_yearly_miroc');
save(fullfile(directory, 'hurs_yearly_miroc.mat'), 'hurs_yearly_miroc');
save(fullfile(directory, 'pr_yearly_miroc.mat'), 'pr_yearly_miroc');
%% load the above data. for this section of code, you can start from this section
load('tas_yearly_miroc.mat');
load('sfcWind_yearly_miroc.mat');
load('hurs_yearly_miroc.mat');
load('pr_yearly_miroc.mat');
%% resize and rotate
% rotate
pr_yearly_miroc=rot90(pr_yearly_miroc)
tas_yearly_miroc=rot90(tas_yearly_miroc)
sfcWind_yearly_miroc=rot90(sfcWind_yearly_miroc)
hurs_yearly_miroc=rot90(hurs_yearly_miroc)
% resize
resize_pr_miroc=imresize(pr_yearly_miroc,[90,180],'nearest')
resize_tas_miroc=imresize(tas_yearly_miroc,[90,180],'nearest')
resize_sfcWind_miroc=imresize(sfcWind_yearly_miroc,[90,180],'nearest')
resize_hurs_miroc=imresize(hurs_yearly_miroc,[90,180],'nearest')

%%  graph
figure14 = figure('WindowState','fullscreen');

subplot(2,2,1)
sgtitle('Near surface wind speed', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(342)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[0,360],size(resize_sfcWind_miroc(:,:,1)));
geoshow(flipud(resize_sfcWind_miroc(:,:,1)), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

subplot(2,2,2)
sgtitle('Mean annual temperature', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(277)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[0,360],size(resize_sfcWind_miroc(:,:,1)));
geoshow(flipud(resize_tas_miroc(:,:,1)), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

subplot(2,2,3)
sgtitle('Precipitation', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(284)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[0,360],size(resize_sfcWind_miroc(:,:,1)));
geoshow(flipud(resize_pr_miroc(:,:,1)), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

subplot(2,2,4)
sgtitle('Relative humidity', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(309)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[0,360],size(resize_sfcWind_miroc(:,:,1)));
geoshow(flipud(resize_hurs_miroc(:,:,1)), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/history_climate.svg');

%% another climate model
load('tas_yearly_mri.mat');
load('sfcWind_yearly_mri.mat');
load('hurs_yearly_mri.mat');
load('pr_yearly_mri.mat');
%% rotate and resize
% rotate
tas_yearly_mri=rot90(tas_yearly_mri);
pr_yearly_mri=rot90(pr_yearly_mri);
sfcWind_yearly_mri=rot90(sfcWind_yearly_mri);
hurs_yearly_mri=rot90(hurs_yearly_mri);
% resize
resize_pr_mri=imresize(pr_yearly_mri,[90,180],'nearest')
resize_tas_mri=imresize(tas_yearly_mri,[90,180],'nearest')
resize_sfcWind_mri=imresize(sfcWind_yearly_mri,[90,180],'nearest')
resize_hurs_mri=imresize(hurs_yearly_mri,[90,180],'nearest')
%% turn them into row vector 

resize_hurs_miroc_vector = reshape(resize_hurs_miroc(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector
resize_pr_miroc_vector = reshape(resize_pr_miroc(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector
resize_tas_miroc_vector = reshape(resize_tas_miroc(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector
resize_sfcWind_miroc_vector = reshape(resize_sfcWind_miroc(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector

resize_hurs_mri_vector = reshape(resize_hurs_mri(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector
resize_pr_mri_vector = reshape(resize_pr_mri(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector
resize_tas_mri_vector = reshape(resize_tas_mri(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector
resize_sfcWind_mri_vector = reshape(resize_sfcWind_mri(:,:,1:8), [], 1);  % Reshape into a 12960x1 column vector

%% save these for max's analysis for grassland distribution
resize_hurs_mri_max=resize_hurs_mri_vector;
resize_pr_mri_max=resize_pr_mri_vector;
resize_tas_mri_max=resize_tas_mri_vector;
resize_sfcWind_mri_max=resize_sfcWind_mri_vector;

resize_hurs_miroc_max=resize_hurs_miroc_vector;
resize_pr_miroc_max=resize_pr_miroc_vector;
resize_tas_miroc_max=resize_tas_miroc_vector;
resize_sfcWind_miroc_max=resize_sfcWind_miroc_vector;

resize_hurs_giss_max=resize_hurs_giss_vector
resize_pr_giss_max=resize_pr_giss_vector
resize_tas_giss_max=resize_tas_giss_vector
resize_sfcWind_giss_max=resize_sfcWind_giss_vector

%% get rid of all of the non livestock located land climate points
% here I am doing this for my 3 climate models so far
% Delete the rows where nan is, so that I can produce fits
% becasue distribution fitting does not allow for nan values

resize_hurs_mri_vector(resizedhyde_vector<=0) = []
resize_pr_mri_vector(resizedhyde_vector<=0) = []
resize_sfcWind_mri_vector(resizedhyde_vector<=0) = []
resize_tas_mri_vector(resizedhyde_vector<=0) = []

resize_hurs_miroc_vector(resizedhyde_vector<=0) = []
resize_pr_miroc_vector(resizedhyde_vector<=0) = []
resize_sfcWind_miroc_vector(resizedhyde_vector<=0) = []
resize_tas_miroc_vector(resizedhyde_vector<=0) = []

resize_hurs_giss_vector(resizedhyde_vector<=0) = []
resize_pr_giss_vector(resizedhyde_vector<=0) = []
resize_sfcWind_giss_vector(resizedhyde_vector<=0) = []
resize_tas_giss_vector(resizedhyde_vector<=0) = []

%% calculate the final mean of my three models%%
pr=[resize_pr_miroc_vector,resize_pr_mri_vector,resize_pr_giss_vector]
hurs=[resize_hurs_miroc_vector,resize_hurs_mri_vector,resize_hurs_giss_vector]
sfcWind=[resize_sfcWind_miroc_vector,resize_sfcWind_mri_vector,resize_sfcWind_giss_vector]
tas=[resize_tas_miroc_vector,resize_tas_mri_vector,resize_tas_giss_vector]

pr_max=[resize_pr_miroc_max,resize_pr_mri_max,resize_pr_giss_max]
hurs_max=[resize_hurs_miroc_max,resize_hurs_mri_max,resize_hurs_giss_max]
sfcWind_max=[resize_sfcWind_miroc_max,resize_sfcWind_mri_max,resize_sfcWind_giss_max]
tas_max=[resize_tas_miroc_max,resize_tas_mri_max,resize_tas_giss_max]


