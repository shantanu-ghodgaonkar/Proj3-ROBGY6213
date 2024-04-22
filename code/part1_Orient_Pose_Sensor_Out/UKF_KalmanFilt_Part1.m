clear; % Clear variables
addpath('../data')
datasetNum = 1; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime,proj2Data] = init(datasetNum);

Z = sampledVicon(1:6,:);
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = 0.1*eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %Just for saving state his.
prevTime = 0; %last time step in real time
pos = proj2Data.position;
pose = proj2Data.angle;
for i = 1:length(sampledTime)
% MY IMPLEMENTATION START -------------------------------------------------
    if sampledData(i).is_ready == 1
        dt = sampledTime(i) - prevTime; % Calculate time interval dt
        prevTime = sampledTime(i); % Update the previous time variable
        % Perform the prediction step
        [covarEst,uEst] = pred_step(uPrev,covarPrev,sampledData(i).omg,sampledData(i).acc,dt);
        % Perform the update step
        [uCurr,covar_curr] = upd_step([pos(i,:), pose(i,:)]',covarEst,uEst);
        % Store updated state for plotting
        savedStates(:, i) = uCurr;
        % Update previous values
        uPrev = uCurr;
        covarPrev = covar_curr;
    end
% MY IMPLEMENTATION END ---------------------------------------------------
end
plotData(savedStates, sampledTime, sampledVicon, 1, datasetNum);