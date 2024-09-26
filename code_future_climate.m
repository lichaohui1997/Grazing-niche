addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));


%% using LUH2 futureland data
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


%% future climate data

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
%% taking the average

scenarios = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
variables = {'pr', 'tas', 'hurs', 'sfcWind'};

for i = 1:length(scenarios)
    for j = 1:length(variables)
        data_var_MRI = sprintf('data_%s_%sMRI', variables{j}, scenarios{i});
        data_var_CHEM = sprintf('data_%s_%sCHEM', variables{j}, scenarios{i});
        data_var_ESM = sprintf('data_%s_%sESM', variables{j}, scenarios{i});
        data_var_avg = sprintf('data_%s_%s', variables{j}, scenarios{i});

        if evalin('base', ['exist(''', data_var_MRI, ''', ''var'')']) && ...
           evalin('base', ['exist(''', data_var_CHEM, ''', ''var'')']) && ...
           evalin('base', ['exist(''', data_var_ESM, ''', ''var'')'])
            avg_cmd = sprintf('%s = (%s + %s + %s) ./ 3;', data_var_avg, data_var_MRI, data_var_CHEM, data_var_ESM);

            evalin('base', avg_cmd);
        else
            fprintf('One or more data variables for %s %s not found in workspace.\n', variables{j}, scenarios{i});
        end
    end
end


%% Future Niche
scenarios = {'rcp26', 'rcp60','rcp45','rcp85'};

futureniche_struct = struct();

for s = 1:length(scenarios)
    scenario = scenarios{s};

    data_pr = eval(['data_pr_' scenario]);
    data_tas = eval(['data_tas_' scenario]);
    data_sfcWind = eval(['data_sfcWind_' scenario]);
    data_hurs = eval(['data_hurs_' scenario]);
    data_landuse = eval(['futureland_' scenario]);


    halfSize = size(data_landuse, 2) / 2;
    data_landuse_new = circshift(data_landuse, [0, halfSize]);
 
        % Define the conditions
    cond_pr = (data_pr >= 433 & data_pr <= 2368);
    cond_tas = (data_tas >= 4 & data_tas <= 28);
    cond_sfcWind = (data_sfcWind >= 1.4 & data_sfcWind <= 6);
    cond_hurs = (data_hurs >= 46 & data_hurs <= 87);
    cond_landuse = (data_landuse_new >= 0.01);


    % Combine the conditions to create the niche map
    all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;


    futureniche = data_landuse_new .* all_conditions;

    % Store the niche map in the structure
    futureniche_struct.(scenario) = futureniche;

    imagesc(futureniche_struct.rcp26)
    imagesc(data_landuse_new)
    imagesc(data_landuse)
    imagesc(all_conditions)


