% Define cell IDs
cellIDs = [2108, 2109, 2110, 2111, 2113, 2114, 2115, 2116, 2118, 2119, ...
           2120, 2121, 2122, 2124, 2125, 2126, 2127, 2128, 2129, 2130, ...
           2131, 2132, 2134, 2135, 2136, 2203, 2210, 2215];

% Load data
distanceMatrix = load_distance_matrix('distance_matrix.csv', cellIDs);
trafficDemand = readtable('traffic_demand.csv');

% Traffic parameters
params.freeFlowSpeed = 60;          % km/h  %~~~~~DATA NEEDED~~~~~
params.congestionWaveSpeed = 20;    % km/h  %~~~~~DATA NEEDED~~~~~
params.jamDensity = 150;            % vehicles/km  %~~~~~DATA NEEDED~~~~~
params.maxFlow = 2000;              % vehicles/hour  %~~~~~DATA NEEDED~~~~~
params.timeStep = 60;               % seconds  %~~~~~DATA NEEDED~~~~~
params.criticalDensity = params.maxFlow / params.freeFlowSpeed;

% Simulation parameters
duration = 2.0;  % hours

% Initial conditions (vehicles/km)
initialDensity = zeros(length(cellIDs), 1);
initialDensity(1) = 30.0;  % Cell 2108  %~~~~~DATA NEEDED~~~~~
initialDensity(2) = 25.0;  % Cell 2109  %~~~~~DATA NEEDED~~~~~
initialDensity(3) = 20.0;  % Cell 2110  %~~~~~DATA NEEDED~~~~~

%% ===== INITIALIZATION =====

fprintf('Initializing Cell Transmission Model...\n');

nCells = length(cellIDs);
nSteps = floor(duration * 3600 / params.timeStep);

% Calculate cell properties
cellLengths = calculate_cell_lengths(distanceMatrix);
maxVehicles = cellLengths * params.jamDensity;
criticalVehicles = cellLengths * params.criticalDensity;

% Initialize state
currentVehicles = initialDensity .* cellLengths;
currentDensity = initialDensity;

% Initialize history tracking
densityHistory = zeros(nSteps, nCells);
flowHistory = cell(nSteps, 1);
timeHistory = zeros(nSteps, 1);

%% ===== SIMULATION LOOP =====

fprintf('Running simulation for %.1f hours...\n', duration);

for step = 1:nSteps
    currentTime = (step - 1) * params.timeStep / 3600.0;  % hours
    
    % Add traffic demand
    currentVehicles = add_demand(currentVehicles, maxVehicles, cellIDs, ...
                                  trafficDemand, currentTime, params.timeStep);
    
    % Compute flows between cells
    flowMatrix = compute_flows(currentVehicles, cellLengths, maxVehicles, ...
                               criticalVehicles, distanceMatrix, params);
    
    % Update state
    currentVehicles = update_state(currentVehicles, flowMatrix, maxVehicles);
    currentDensity = currentVehicles ./ cellLengths;
    
    % Record history
    densityHistory(step, :) = currentDensity';
    flowHistory{step} = flowMatrix;
    timeHistory(step) = currentTime;
    
    % Print progress
    if mod(step, max(1, floor(nSteps/10))) == 0
        fprintf('Progress: %.1f%%\n', 100 * step / nSteps);
    end
end

fprintf('Simulation complete!\n');

%% ===== ANALYSIS & VISUALIZATION =====

fprintf('\nGenerating results...\n');

% Plot results
cellsToPlot = [2108, 2109, 2110, 2111, 2113];
plot_results(densityHistory, timeHistory, cellIDs, cellsToPlot, cellLengths);

% Export results
export_results(densityHistory, timeHistory, cellIDs, cellLengths, 'ctm_results.csv');

% Calculate summary statistics
totalVehicles = densityHistory * cellLengths;
[maxVehicles_total, maxIdx] = max(totalVehicles);
peakTime = timeHistory(maxIdx);

fprintf('\n=== Summary Statistics ===\n');
fprintf('Peak network occupancy: %.0f vehicles at time %.2f hours\n', ...
        maxVehicles_total, peakTime);

% Average density per cell
avgDensities = mean(densityHistory, 1);
for i = 1:min(5, nCells)
    fprintf('Cell %d average density: %.2f vehicles/km\n', ...
            cellIDs(i), avgDensities(i));
end

fprintf('\nResults saved!\n');