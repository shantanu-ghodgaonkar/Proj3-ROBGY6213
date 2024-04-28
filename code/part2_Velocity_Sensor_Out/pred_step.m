 function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
    %% Parameter Definition
    % uPrev - is the mean of the prev state
    %covarPrev - covar of the prev state
    %angVel - angular velocity input at the time step
    %acc - acceleration at the timestep
    %dt - difference in time 
% MY IMPLEMENTATION START -------------------------------------------------
    % Set nDash = n + nq = 15 + 12 = 27
    nDash = 27;
    % Set UKF Parameters
    alpha = 0.001;
    k = 1;
    beta = 2;
    % Compute lambda dash
    lambdaDash = (alpha^2) * (nDash + k) - nDash;
    % Create the mean augmented matrix
    uAugPrev = [uPrev; zeros(12,1)];
    % Set Q as 0.02
    Q = 0.004;
    % Calculate Qd
    Qd = dt*(eye(12) * Q);
    % Create the covariance augmented matrix
    PAug = [covarPrev, zeros(15,12); zeros(12,15), Qd];
    % Find the square root of the covariance augmented matrix
    sqrtPAug = chol(PAug);
    % initialise the sigma points matrix if size 27x55
    XAug = zeros(nDash, ((2*nDash)+1));
    % Set the first value of the sigma points
    XAug(:,1) = uAugPrev;
    % Compute all other sigma points
    for i = 2:1:(nDash+1)
        XAug(:,i) = uAugPrev + sqrt(nDash + lambdaDash) * sqrtPAug(:,i-1);
        XAug(:,(i+nDash)) = uAugPrev - sqrt(nDash + lambdaDash) * sqrtPAug(:,i-1);
    end
    % Initialise state vector matrix
    Xt = zeros(15, ((2*nDash)+1));
    % Propagate the sigma points through the non-linear process function
    for i = 1:1:((2*nDash)+1)
        x1 = XAug(1:3, i);
        x2 = XAug(4:6, i);
        x3 = XAug(7:9, i);
        x4 = XAug(10:12, i);
        x5 = XAug(13:15, i);
        ng = XAug(16:18, i);
        na = XAug(19:21, i);
        nbg = XAug(22:24, i);
        nba = XAug(25:27, i);

        % Precompute required sines and cosines to reduce computation time
        sxo = sin(x2(1));
        syo = sin(x2(2));
        szo = sin(x2(3));
        cxo = cos(x2(1));
        cyo = cos(x2(2));
        czo = cos(x2(3));

        % Calculate inverse of Euler Rate Parameterisation (ZYX) for f(x,u,n)
        G_inv=[(czo*syo)/cyo,(syo*szo)/cyo,1;-szo,czo,0;czo/cyo,szo/cyo,0];
        % G_inv = round(G_inv, 4);
        % Calculate Rotation matrix (ZYX) for f(x,u,n)
        R=[cyo*czo,czo*sxo*syo-cxo*szo,sxo*szo+cxo*czo*syo;
            cyo*szo,cxo*czo+sxo*syo*szo,cxo*syo*szo-czo*sxo;
            -syo,cyo*sxo,cxo*cyo];
        % R = round(R, 4);
        f1 = x3;
        f2 = G_inv * (angVel - x4);
        f3 = [0;0;-9.8] + R *(acc - x5 - na); % 9.8 is only for Z axis
        f4 = nbg;
        f5 = nba; 
        % Perform euler integration and find the estimated state vectors
        Xt(:,i) = [x1 + (dt * f1); 
            x2 + (dt * f2) - ng;
            x3 + (dt * f3) - na;
            x4 + (f4);
            x5 + (f5)];
    end
    % Compute the weight array for mean and covariances
    WmDash = [(lambdaDash / (nDash + lambdaDash)), ( ones(1, 2*nDash)*(1 / (2 * (nDash + lambdaDash))) )];
    WcDash = [(( lambdaDash / (nDash + lambdaDash) ) + (1 - (alpha^2) + beta)), ( ones(1, 2*nDash)*(1 / (2 * (nDash + lambdaDash))) )];
    % Initialise the estimated mean state vector
    uEst = zeros(15,1);
    % Find the weighted sum of the propagated estimated state vectors
    for i = 1:1:((2*nDash)+1)
        uEst = uEst + (WmDash(i) * Xt(:,i));
    end
    % Initialise the estimated noise covariance matrix
    covarEst = zeros(15,15);
    % Find the weighted sum of the estimated noise covariance matrix
    for i = 1:1:((2*nDash)+1)
        covarEst = covarEst + (WcDash(i) * ((Xt(:,i) - uEst) * (Xt(:,i) - uEst)'));
    end
% MY IMPLEMENTATION END ---------------------------------------------------
end

