function [windfield, coords, outtime] = func_import_probeCSV( ...
    folder, timeOffset, plane, mask, nComponents)
% Import CSV files with format as present in the DATA folder
% Author : N.S.Dangi ; TU Delft

if nargin < 5 || isempty(nComponents)
    nComponents = 1;
end

nComponents = min(nComponents, 3);

% --- Load time and apply time filter ---
time = readmatrix(fullfile(folder, "Timestamps.csv"), ...
    "NumHeaderLines", 1);

keepTime = time >= time(1) + timeOffset;
outtime = time(keepTime);

% --- Load velocity components ---
componentFiles = ["U.csv", "V.csv", "W.csv"];
windfield = cell(1, nComponents);

for k = 1:nComponents
    velocity = readmatrix(fullfile(folder, componentFiles(k)), ...
        "NumHeaderLines", 1);

    windfield{k} = velocity(:,keepTime)';
end

% --- Load coordinates ---
coords = readmatrix(fullfile(folder, "Coords.csv"), ...
    "NumHeaderLines", 1);

x = coords(:,1);
y = coords(:,2);
z = coords(:,3);

% --- Apply spatial mask ---
switch upper(extractBetween(plane, 1, 2))
    case 'YZ'
        keepCoords = ...
            y >= mask.y(1) & y <= mask.y(2) & ...
            z >= mask.z(1) & z <= mask.z(2);

    case 'XZ'
        keepCoords = ...
            x >= mask.x(1) & x <= mask.x(2) & ...
            z >= mask.z(1) & z <= mask.z(2);

    case 'XY'
        keepCoords = ...
            x >= mask.x(1) & x <= mask.x(2) & ...
            y >= mask.y(1) & y <= mask.y(2);
end

coords = coords(keepCoords,:);

for k = 1:nComponents
    windfield{k} = windfield{k}(:,keepCoords);
end

end
