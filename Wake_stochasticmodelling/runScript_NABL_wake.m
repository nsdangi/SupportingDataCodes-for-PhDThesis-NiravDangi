clear
clc
% script to load LES datasets; estiamte Kramers Moyal coefficients and
% generate corresponding stochastic trajectories

% Author : N.S.Dangi ; TU Delft

rng(0)
addpath(genpath("functions"))
addpath(genpath("DATA/NABL"))
% --- Cases ---
cases = ["6p5", "11p5", ...
         "16p5", "20p5"];

uinfhub = str2double(replace(cases, "p", "."));

% --- Slice folders ---
sliceNames = ["yzslice_5"]; % Modelling and analysis at 5D


turbDia       = 284; % IEA 22 MW
turbHub       = [1120,2500,170]; % ALM-LES based

%% --- Settings/ Masks ---
nComponents   = 1; % Anlayse u component 
max_xByD      = 5; % Modelling and analysis at 5D
min_flowTime  = (turbHub(1)+max_xByD*turbDia)./uinfhub; % To disregard one flow through time


xMin = turbHub(1);
xMax = turbHub(1) + max_xByD*turbDia;

yMin = turbHub(2) - 1*turbDia;
yMax = turbHub(2) + 1*turbDia;

zMin = 0;%turbHub(3) - 0.6*turbDia;
zMax = turbHub(3) + 0.6*turbDia;


mask.x = [xMin, xMax];
mask.y = [yMin, yMax];
mask.z = [zMin, zMax];

for i = 1:numel(cases)
    for j = 1%1:numel(sliceNames)

        sliceFolder = fullfile(cases(i), sliceNames(j));

        % Select mask based on plane
        switch upper(extractBetween(sliceNames(j), 1, 2))
            case "XY"
                planeMask = struct('x', mask.x, 'y', mask.y);
            case "XZ"
                planeMask = struct('x', mask.x, 'z', mask.z);
            case "YZ"
                planeMask = struct('y', mask.y, 'z', mask.z);
        end

        % Import data
        [windfield{i,j}, coords{i,j}, outtime{i,j}] = func_import_probeCSV( ...
            sliceFolder, min_flowTime(i), sliceNames(j), planeMask, nComponents);

        % FPE and Langevin settings
        nBins = 15; % 10 
        [D1_all{i,j}, D2_all{i,j}, centers_ALL{i,j}] = ...
        func_estimateDriftDiffusion(outtime{i,j}, ...
            windfield{i,j}, nBins, nComponents);

        [interpD1_all{i,j}, interpD2_all{i,j}] = func_interpolateDriftDiffusion( ...
            D1_all{i,j}, D2_all{i,j}, centers_ALL{i,j});

        n_transients = 3;
        u_syn_ALL{i,j} = func_sim_Langevin( ...
            outtime{i,j}, windfield{i,j}, ...
            interpD1_all{i,j}, interpD2_all{i,j}, ...
            centers_ALL{i,j}, nComponents,n_transients);

    end
end
