% CustomStaticOptimization with Sensitivity Analysis
% Analyzing gastrocnemius medialis and soleus fiber power sensitivity
% Using ±1 SD variations in muscle morphology parameters
% Based on Rajagopal et al. 2016 data

close all; clear all; clc; beep off;

%% Import the OpenSim libraries.
import org.opensim.modeling.*;

%% Define sensitivity analysis parameters based on Rajagopal 2016
% Muscles to analyze
targetMuscles = {'gasmed_r', 'gasmed_l', 'soleus_r', 'soleus_l'};

% Mean and SD values from Rajagopal et al. 2016
rajagopalParams = struct();

% Soleus
rajagopalParams.soleus.Fmax_mean = 6195;      % N
rajagopalParams.soleus.Fmax_sd = 1606;        % N
rajagopalParams.soleus.Lopt_mean = 0.044;     % m (4.4 cm)
rajagopalParams.soleus.Lopt_sd = 0.010;       % m (1.0 cm)
rajagopalParams.soleus.Lts_mean = 0.277;      % m (27.7 cm)
rajagopalParams.soleus.Lts_sd = 0.010;        % m (1.0 cm)
rajagopalParams.soleus.penn_mean = 21.9;      % degrees
rajagopalParams.soleus.penn_sd = 8.0;         % degrees

% Gastrocnemius medialis
rajagopalParams.gasmed.Fmax_mean = 3116;      % N
rajagopalParams.gasmed.Fmax_sd = 727;         % N
rajagopalParams.gasmed.Lopt_mean = 0.051;     % m (5.1 cm)
rajagopalParams.gasmed.Lopt_sd = 0.010;       % m (1.0 cm)
rajagopalParams.gasmed.Lts_mean = 0.399;      % m (39.9 cm)
rajagopalParams.gasmed.Lts_sd = 0.011;        % m (1.1 cm)
rajagopalParams.gasmed.penn_mean = 9.5;       % degrees
rajagopalParams.gasmed.penn_sd = 4.3;         % degrees

% Define test conditions
conditions = {
    'Baseline', ...
    '+1SD_Morphology', ...
    '-1SD_Morphology', ...
    '+1SD_TendonCompliance', ...
    '-1SD_TendonCompliance'
};

% Create valid field names for struct storage (remove special characters)
conditionFieldNames = {
    'Baseline', ...
    'Plus1SD_Morphology', ...
    'Minus1SD_Morphology', ...
    'Plus1SD_TendonCompliance', ...
    'Minus1SD_TendonCompliance'
};

%% Load model and get initial state.
model = Model('subject_walk_adjusted.osim');
state = model.initSystem();

% Store baseline muscle parameters from the model
muscles = model.updMuscles();
baselineParams = struct();
for i = 1:muscles.getSize()
    muscle = Millard2012EquilibriumMuscle.safeDownCast(muscles.get(i-1));
    muscleName = char(muscle.getName());
    baselineParams.(muscleName).Fmax = muscle.get_max_isometric_force();
    baselineParams.(muscleName).Lopt = muscle.get_optimal_fiber_length();
    baselineParams.(muscleName).Lts = muscle.get_tendon_slack_length();
    baselineParams.(muscleName).pennation = muscle.get_pennation_angle_at_optimal();
end

%% Use OpenSim tools to generate data for optimization.
if ~exist('analyze_Kinematics_q.sto', 'file')
    fprintf('Running kinematics analysis...\n\n');
    analyze = AnalyzeTool(model);
    analyze.setName('analyze');
    analyze.setCoordinatesFileName('coordinates.mot');
    analyze.loadStatesFromFile(state);
    analyze.setStartTime(0);
    analyze.setFinalTime(2.37);
    analysisSet = analyze.getAnalysisSet();
    kinematicsAnalysis = Kinematics();
    kinematicsAnalysis.setInDegrees(false);
    analysisSet.cloneAndAppend(kinematicsAnalysis);
    analyze.addAnalysisSetToModel();
    analyze.run();
end

genForcesFile = 'generalized_forces.sto';
if ~exist(genForcesFile, 'file')
    fprintf('Running inverse dynamics...\n\n');
    idtool = InverseDynamicsTool();        
    idtool.setModel(model);
    idtool.setStartTime(0);
    idtool.setEndTime(2.37);
    idtool.setExternalLoadsFileName('grf_walk.xml');
    idtool.setCoordinatesFileName('coordinates.mot');
    idtool.setOutputGenForceFileName(genForcesFile);
    excludedForces = ArrayStr();
    excludedForces.append('muscles');
    idtool.setExcludedForces(excludedForces);
    idtool.run();
end

%% Load data into MATLAB arrays.
lowpassFreq = 6.0;
timeRange = [0.81 1.96];
[coordinates, coordNames, time] = ...
        loadFilterCropArray('analyze_Kinematics_q.sto', lowpassFreq, timeRange);
[speeds, speedNames, ~] = ...
        loadFilterCropArray('analyze_Kinematics_u.sto', lowpassFreq, timeRange);
[genForces, forceNames, ~] = ...
        loadFilterCropArray(genForcesFile, lowpassFreq, timeRange);

% Re-order generalized forces
force2coord = zeros(length(coordNames),1);
for i = 1:length(forceNames)
   forcename = forceNames{i};
   for j = 1:length(coordNames)
       coordname = coordNames{j};
       if contains(forcename, '_moment')
           forcename = forcename(1:end-7);
       elseif contains(forcename, '_force')
           forcename = forcename(1:end-6);
       end
       if strcmp(forcename, coordname)
           force2coord(j) = i;
       end
   end
end
genForces = genForces(:, force2coord);

%% Remove unneeded generalized force data columns.
forcesToRemove = {'beta', 'arm', 'elbow', 'wrist', 'pro_sup', ...
    'mtp','subtalar','pelvis','lumbar'};
colsToRemove = [];
for i = 1:length(coordNames)
    for j = 1:length(forcesToRemove)
        if contains(coordNames{i}, forcesToRemove{j})
            colsToRemove = [colsToRemove i];
        end
    end
end
coordNamesAll = coordNames;
genForces(:, colsToRemove) = [];  
coordNames(colsToRemove) = [];

%% MAIN SENSITIVITY ANALYSIS LOOP
numConditions = length(conditions);

% Storage for results
sensitivityResults = struct();
for m = 1:length(targetMuscles)
    muscleName = targetMuscles{m};
    sensitivityResults.(muscleName).peakPower = zeros(numConditions, 1);
    sensitivityResults.(muscleName).avgPower = zeros(numConditions, 1);
    sensitivityResults.(muscleName).totalWork = zeros(numConditions, 1);
    sensitivityResults.(muscleName).powerTimeSeries = cell(numConditions, 1);
    sensitivityResults.(muscleName).activationTimeSeries = cell(numConditions, 1);
    sensitivityResults.(muscleName).fiberLengthTimeSeries = cell(numConditions, 1);
    sensitivityResults.(muscleName).fiberVelocityTimeSeries = cell(numConditions, 1);
end

% Store applied parameter modifications for reporting
appliedParams = struct();

for condIdx = 1:numConditions
    conditionName = conditions{condIdx};
    
    fprintf('\n========================================\n');
    fprintf('Running condition %d/%d: %s\n', condIdx, numConditions, conditionName);
    fprintf('========================================\n\n');
    
    % Reload model for each condition
    model = Model('subject_walk_adjusted.osim');
    state = model.initSystem();
    muscles = model.updMuscles();
    
    % Modify muscle parameters based on condition
    for i = 1:muscles.getSize()
        muscle = Millard2012EquilibriumMuscle.safeDownCast(muscles.get(i-1));
        muscleName = char(muscle.getName());
        
        % Check if this is a target muscle
        isTargetMuscle = false;
        muscleType = '';
        for m = 1:length(targetMuscles)
            if strcmp(muscleName, targetMuscles{m})
                isTargetMuscle = true;
                if contains(muscleName, 'soleus')
                    muscleType = 'soleus';
                elseif contains(muscleName, 'gasmed')
                    muscleType = 'gasmed';
                end
                break;
            end
        end
        
        if isTargetMuscle
            % Get baseline values
            baseline_Fmax = baselineParams.(muscleName).Fmax;
            baseline_Lopt = baselineParams.(muscleName).Lopt;
            baseline_Lts = baselineParams.(muscleName).Lts;
            baseline_penn = baselineParams.(muscleName).pennation;
            
            % Apply modifications based on condition
            switch conditionName
                case 'Baseline'
                    % No changes
                    new_Fmax = baseline_Fmax;
                    new_Lopt = baseline_Lopt;
                    new_Lts = baseline_Lts;
                    new_penn = baseline_penn;
                    
                case '+1SD_Morphology'
                    % Increase all parameters by +1 SD
                    new_Fmax = baseline_Fmax + rajagopalParams.(muscleType).Fmax_sd;
                    new_Lopt = baseline_Lopt + rajagopalParams.(muscleType).Lopt_sd;
                    new_Lts = baseline_Lts + rajagopalParams.(muscleType).Lts_sd;
                    new_penn = baseline_penn + rajagopalParams.(muscleType).penn_sd * (pi/180); % Convert to radians
                    
                case '-1SD_Morphology'
                    % Decrease all parameters by -1 SD
                    new_Fmax = baseline_Fmax - rajagopalParams.(muscleType).Fmax_sd;
                    new_Lopt = baseline_Lopt - rajagopalParams.(muscleType).Lopt_sd;
                    new_Lts = baseline_Lts - rajagopalParams.(muscleType).Lts_sd;
                    new_penn = baseline_penn - rajagopalParams.(muscleType).penn_sd * (pi/180);
                    
                case '+1SD_TendonCompliance'
                    % Only increase tendon slack length (more compliant)
                    new_Fmax = baseline_Fmax;
                    new_Lopt = baseline_Lopt;
                    new_Lts = baseline_Lts + rajagopalParams.(muscleType).Lts_sd;
                    new_penn = baseline_penn;
                    
                case '-1SD_TendonCompliance'
                    % Only decrease tendon slack length (stiffer)
                    new_Fmax = baseline_Fmax;
                    new_Lopt = baseline_Lopt;
                    new_Lts = baseline_Lts - rajagopalParams.(muscleType).Lts_sd;
                    new_penn = baseline_penn;
            end
            
            % Apply modifications to muscle
            muscle.set_max_isometric_force(new_Fmax);
            muscle.set_optimal_fiber_length(new_Lopt);
            muscle.set_tendon_slack_length(new_Lts);
            muscle.set_pennation_angle_at_optimal(new_penn);
            
            % Store applied parameters (use valid field name)
            if condIdx == 1
                appliedParams.(muscleName) = struct();
            end
            fieldName = conditionFieldNames{condIdx};
            appliedParams.(muscleName).(fieldName).Fmax = new_Fmax;
            appliedParams.(muscleName).(fieldName).Lopt = new_Lopt;
            appliedParams.(muscleName).(fieldName).Lts = new_Lts;
            appliedParams.(muscleName).(fieldName).pennation = new_penn * (180/pi);
            
            fprintf('  %s:\n', muscleName);
            fprintf('    Fmax: %.1f N, Lopt: %.4f m, Lts: %.4f m, Penn: %.1f deg\n', ...
                    new_Fmax, new_Lopt, new_Lts, new_penn*(180/pi));
        end
        
        % Disable muscle dynamics for static optimization
        muscle.set_ignore_tendon_compliance(true);
        muscle.set_ignore_activation_dynamics(true);
    end
    
    state = model.initSystem();
    
    % Store max isometric forces
    Fmax = zeros(muscles.getSize(), 1);
    for i = 1:muscles.getSize()
        muscle = Millard2012EquilibriumMuscle.safeDownCast(muscles.get(i-1));
        Fmax(i) = muscle.get_max_isometric_force();
    end
    
    %% Perform static optimization
    timeInterval = 5;
    N = size(coordinates, 1);
    coords_subset = coordinates(1:timeInterval:N, :);
    speeds_subset = speeds(1:timeInterval:N, :);
    genForces_subset = genForces(1:timeInterval:N, :);
    numTimePoints = size(coords_subset, 1);
    
    % FMINCON setup
    options = optimoptions('fmincon','Display','off', ...
         'TolCon',1e-4,'TolFun',1e-12,'TolX',1e-8,'MaxFunEvals',100000, ...
         'MaxIter',5000,'Algorithm','interior-point');
    
    numCoords = length(coordNames);
    numMuscles = muscles.getSize();
    lb = [zeros(1,numMuscles), -ones(1,numCoords)];
    ub = [ones(1,numMuscles), ones(1,numCoords)];
    x0 = zeros(1,numMuscles+numCoords);
    
    w = [ones(1,numMuscles), 100*ones(1,numCoords)];
    cost = @(x) sum(w.*(x.^2));
    
    % Pre-allocate solution and fiber mechanics arrays
    xsol = zeros(numTimePoints, length(x0));
    fiberForce = zeros(numTimePoints, numMuscles);
    fiberVelocity = zeros(numTimePoints, numMuscles);
    fiberLength = zeros(numTimePoints, numMuscles);
    activations = zeros(numTimePoints, numMuscles);
    
    coords = model.getCoordinateSet();
    
    % Optimization loop
    for i = 1:numTimePoints
        if mod(i, 10) == 0
            fprintf('  Time step %i/%i\n', i, numTimePoints);
        end
        
        % Set kinematics
        for j = 1:length(coordNamesAll)
            coord = coords.get(coordNamesAll{j});
            coord.setValue(state, coords_subset(i,j));
            coord.setSpeedValue(state, speeds_subset(i,j));
        end
        
        % Compute muscle properties
        fl = zeros(1, numMuscles);
        fv = zeros(1, numMuscles);
        fp = zeros(1, numMuscles);
        cosPenn = zeros(1, numMuscles);
        r = zeros(length(coordNames), numMuscles);
        
        model.realizeVelocity(state);
        for j = 1:length(coordNames)
            coord = coords.get(coordNames{j});
            for k = 1:muscles.getSize()
                muscle = Millard2012EquilibriumMuscle.safeDownCast(muscles.get(k-1));
                r(j,k) = muscle.computeMomentArm(state, coord);
                if j == 1
                    fl(k) = muscle.getActiveForceLengthMultiplier(state);            
                    fv(k) = muscle.getForceVelocityMultiplier(state);
                    fp(k) = muscle.getPassiveForceMultiplier(state);
                    cosPenn(k) = muscle.getCosPennationAngle(state);
                    
                    % Store fiber mechanics for power calculation
                    fiberLength(i,k) = muscle.getFiberLength(state);
                    fiberVelocity(i,k) = muscle.getFiberVelocity(state);
                end
            end 
        end
        
        % Create linear equality constraints
        Amusc = bsxfun(@times, r, fl.*fv.*Fmax'.*cosPenn);
        Ares = eye(numCoords);
        Aeq = [Amusc Ares];
        
        M = genForces_subset(i, :)';
        Mpass = r*(Fmax.*fp'.*cosPenn');
        Beq = M - Mpass;
        
        % Solve optimization
        x = fmincon(cost, x0, [], [], Aeq, Beq, lb, ub, [], options);
        
        xsol(i, :) = x;
        x0 = x;
        
        % Extract activations and compute fiber force
        activations(i,:) = x(1:numMuscles);
        fiberForce(i,:) = activations(i,:) .* fl .* fv .* Fmax' + fp .* Fmax';
    end
    
    %% Calculate fiber power for target muscles
    % Power = Force * Velocity
    fiberPower = fiberForce .* fiberVelocity; % Watts
    
    % Get muscle name mapping
    muscleNames = ArrayStr();
    muscles.getNames(muscleNames);
    
    % Calculate time step for work integration
    timeVec = linspace(timeRange(1), timeRange(2), numTimePoints);
    dt = timeVec(2) - timeVec(1);
    
    for m = 1:length(targetMuscles)
        muscleName = targetMuscles{m};
        
        % Find muscle index
        muscleIdx = -1;
        for k = 1:muscles.getSize()
            if strcmp(char(muscleNames.get(k-1)), muscleName)
                muscleIdx = k;
                break;
            end
        end
        
        if muscleIdx > 0
            powerTS = fiberPower(:, muscleIdx);
            sensitivityResults.(muscleName).powerTimeSeries{condIdx} = powerTS;
            sensitivityResults.(muscleName).activationTimeSeries{condIdx} = activations(:, muscleIdx);
            sensitivityResults.(muscleName).fiberLengthTimeSeries{condIdx} = fiberLength(:, muscleIdx);
            sensitivityResults.(muscleName).fiberVelocityTimeSeries{condIdx} = fiberVelocity(:, muscleIdx);
            
            % Calculate metrics
            sensitivityResults.(muscleName).peakPower(condIdx) = max(powerTS);
            sensitivityResults.(muscleName).avgPower(condIdx) = mean(abs(powerTS));
            
            % Calculate total work (integral of power over time)
            % Only positive work (muscle doing work)
            positiveWork = trapz(timeVec, max(powerTS, 0));
            sensitivityResults.(muscleName).totalWork(condIdx) = positiveWork;
        end
    end
end

%% Visualize sensitivity results
pgc = linspace(0, 100, numTimePoints);

% Create comprehensive figure for each muscle
for m = 1:length(targetMuscles)
    muscleName = targetMuscles{m};
    
    figure('Position', [50 + (m-1)*50, 50, 1600, 1000]);
    
    % Power time series
    subplot(2,3,1);
    colors = lines(numConditions);
    for c = 1:numConditions
        plot(pgc, sensitivityResults.(muscleName).powerTimeSeries{c}, ...
             'Color', colors(c,:), 'LineWidth', 2, 'DisplayName', conditions{c});
        hold on;
    end
    xlabel('Gait Cycle (%)');
    ylabel('Fiber Power (W)');
    title('Fiber Power');
    legend('Location', 'best');
    grid on;
    
    % Activation time series
    subplot(2,3,2);
    for c = 1:numConditions
        plot(pgc, sensitivityResults.(muscleName).activationTimeSeries{c}, ...
             'Color', colors(c,:), 'LineWidth', 2, 'DisplayName', conditions{c});
        hold on;
    end
    xlabel('Gait Cycle (%)');
    ylabel('Activation');
    title('Muscle Activation');
    legend('Location', 'best');
    grid on;
    ylim([0 1]);
    
    % Fiber length time series
    subplot(2,3,3);
    for c = 1:numConditions
        plot(pgc, sensitivityResults.(muscleName).fiberLengthTimeSeries{c}, ...
             'Color', colors(c,:), 'LineWidth', 2, 'DisplayName', conditions{c});
        hold on;
    end
    xlabel('Gait Cycle (%)');
    ylabel('Fiber Length (m)');
    title('Fiber Length');
    legend('Location', 'best');
    grid on;
    
    % Fiber velocity time series
    subplot(2,3,4);
    for c = 1:numConditions
        plot(pgc, sensitivityResults.(muscleName).fiberVelocityTimeSeries{c}, ...
             'Color', colors(c,:), 'LineWidth', 2, 'DisplayName', conditions{c});
        hold on;
    end
    xlabel('Gait Cycle (%)');
    ylabel('Fiber Velocity (m/s)');
    title('Fiber Velocity');
    legend('Location', 'best');
    grid on;
    
    % Bar plot of peak power
    subplot(2,3,5);
    bar(sensitivityResults.(muscleName).peakPower);
    set(gca, 'XTickLabel', conditions, 'XTickLabelRotation', 45);
    ylabel('Peak Power (W)');
    title('Peak Power Comparison');
    grid on;
    
    % Bar plot of total positive work
    subplot(2,3,6);
    bar(sensitivityResults.(muscleName).totalWork);
    set(gca, 'XTickLabel', conditions, 'XTickLabelRotation', 45);
    ylabel('Total Positive Work (J)');
    title('Total Positive Work');
    grid on;
    
    sgtitle(sprintf('%s Sensitivity Analysis', muscleName), 'Interpreter', 'none');
end

% Create summary comparison figure
figure('Position', [100, 100, 1400, 800]);

% Peak power comparison
subplot(2,3,1);
peakPowerMatrix = zeros(length(targetMuscles), numConditions);
for m = 1:length(targetMuscles)
    peakPowerMatrix(m,:) = sensitivityResults.(targetMuscles{m}).peakPower';
end
bar(peakPowerMatrix');
legend(targetMuscles, 'Interpreter', 'none', 'Location', 'best');
set(gca, 'XTickLabel', conditions, 'XTickLabelRotation', 45);
ylabel('Peak Power (W)');
title('Peak Power Comparison');
grid on;

% Average power comparison
subplot(2,3,2);
avgPowerMatrix = zeros(length(targetMuscles), numConditions);
for m = 1:length(targetMuscles)
    avgPowerMatrix(m,:) = sensitivityResults.(targetMuscles{m}).avgPower';
end
bar(avgPowerMatrix');
legend(targetMuscles, 'Interpreter', 'none', 'Location', 'best');
set(gca, 'XTickLabel', conditions, 'XTickLabelRotation', 45);
ylabel('Average Absolute Power (W)');
title('Average Power Comparison');
grid on;

% Total work comparison
subplot(2,3,3);
workMatrix = zeros(length(targetMuscles), numConditions);
for m = 1:length(targetMuscles)
    workMatrix(m,:) = sensitivityResults.(targetMuscles{m}).totalWork';
end
bar(workMatrix');
legend(targetMuscles, 'Interpreter', 'none', 'Location', 'best');
set(gca, 'XTickLabel', conditions, 'XTickLabelRotation', 45);
ylabel('Total Positive Work (J)');
title('Total Work Comparison');
grid on;

% Percent change from baseline - Peak Power
subplot(2,3,4);
peakPowerPctChange = zeros(length(targetMuscles), numConditions-1);
for m = 1:length(targetMuscles)
    baseline = sensitivityResults.(targetMuscles{m}).peakPower(1);
    peakPowerPctChange(m,:) = ((sensitivityResults.(targetMuscles{m}).peakPower(2:end) - baseline) / baseline) * 100;
end
bar(peakPowerPctChange');
legend(targetMuscles, 'Interpreter', 'none', 'Location', 'best');
set(gca, 'XTickLabel', conditions(2:end), 'XTickLabelRotation', 45);
ylabel('Change from Baseline (%)');
title('Peak Power % Change');
grid on;
yline(0, 'k--');

% Percent change from baseline - Total Work
subplot(2,3,5);
workPctChange = zeros(length(targetMuscles), numConditions-1);
for m = 1:length(targetMuscles)
    baseline = sensitivityResults.(targetMuscles{m}).totalWork(1);
    workPctChange(m,:) = ((sensitivityResults.(targetMuscles{m}).totalWork(2:end) - baseline) / baseline) * 100;
end
bar(workPctChange');
legend(targetMuscles, 'Interpreter', 'none', 'Location', 'best');
set(gca, 'XTickLabel', conditions(2:end), 'XTickLabelRotation', 45);
ylabel('Change from Baseline (%)');
title('Total Work % Change');
grid on;
yline(0, 'k--');

sgtitle('Plantarflexor Sensitivity Analysis Summary');

%% Print detailed results table
fprintf('\n========================================\n');
fprintf('SENSITIVITY ANALYSIS RESULTS\n');
fprintf('========================================\n\n');

for m = 1:length(targetMuscles)
    muscleName = targetMuscles{m};
    fprintf('\n%s:\n', muscleName);
    fprintf('%-25s %12s %12s %12s\n', 'Condition', 'Peak Power', 'Avg Power', 'Total Work');
    fprintf('%-25s %12s %12s %12s\n', '', '(W)', '(W)', '(J)');
    fprintf('-----------------------------------------------------------\n');
    
    for c = 1:numConditions
        fprintf('%-25s %12.2f %12.2f %12.2f\n', ...
                conditions{c}, ...
                sensitivityResults.(muscleName).peakPower(c), ...
                sensitivityResults.(muscleName).avgPower(c), ...
                sensitivityResults.(muscleName).totalWork(c));
    end
    
    fprintf('\nPercent Change from Baseline:\n');
    fprintf('%-25s %12s %12s %12s\n', 'Condition', 'Peak Power', 'Avg Power', 'Total Work');
    fprintf('%-25s %12s %12s %12s\n', '', '(%)', '(%)', '(%)');
    fprintf('-----------------------------------------------------------\n');
    
    baseline_peak = sensitivityResults.(muscleName).peakPower(1);
    baseline_avg = sensitivityResults.(muscleName).avgPower(1);
    baseline_work = sensitivityResults.(muscleName).totalWork(1);
    
    for c = 2:numConditions
        pct_peak = ((sensitivityResults.(muscleName).peakPower(c) - baseline_peak) / baseline_peak) * 100;
        pct_avg = ((sensitivityResults.(muscleName).avgPower(c) - baseline_avg) / baseline_avg) * 100;
        pct_work = ((sensitivityResults.(muscleName).totalWork(c) - baseline_work) / baseline_work) * 100;
        
        fprintf('%-25s %+11.1f%% %+11.1f%% %+11.1f%%\n', ...
                conditions{c}, pct_peak, pct_avg, pct_work);
    end
end

%% Print applied parameter modifications
fprintf('\n========================================\n');
fprintf('APPLIED PARAMETER MODIFICATIONS\n');
fprintf('========================================\n\n');

for m = 1:length(targetMuscles)
    muscleName = targetMuscles{m};
    fprintf('\n%s:\n', muscleName);
    fprintf('%-25s %12s %12s %12s %12s\n', 'Condition', 'Fmax (N)', 'Lopt (m)', 'Lts (m)', 'Penn (deg)');
    fprintf('-------------------------------------------------------------------------\n');
    
    for c = 1:numConditions
        fieldName = conditionFieldNames{c};
        fprintf('%-25s %12.1f %12.4f %12.4f %12.1f\n', ...
                conditions{c}, ...
                appliedParams.(muscleName).(fieldName).Fmax, ...
                appliedParams.(muscleName).(fieldName).Lopt, ...
                appliedParams.(muscleName).(fieldName).Lts, ...
                appliedParams.(muscleName).(fieldName).pennation);
    end
end

%% Save results
save('sensitivity_results.mat', 'sensitivityResults', 'conditions', ...
     'targetMuscles', 'rajagopalParams', 'appliedParams', 'baselineParams');

fprintf('\n\nResults saved to sensitivity_results.mat\n');
fprintf('Analysis complete!\n');