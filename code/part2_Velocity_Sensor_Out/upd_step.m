function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    % R = Noise in System
    R =  0.1 * eye(3);
    % vt =  Measurement Update Model noise
    vt = normrnd(0,(0.0005),[3,1]);
    % n = size of state vector
    n = 15;
    % alpha, k and beta are UKF Parameters
    alpha = 0.001;
    k = 1;
    beta = 2;
    % lambda is the spread in sigma points
    lambda = ((alpha^2)*(n + k)) - n;
    % initialise sigma points matrix of vectors to 0
    Xt = zeros(15, ((2*n)+1));
    % Set sigma point at 0 to be equal to estimated mean
    Xt(:,1) = uEst;
    % Find sqaure root of estimated covariance matrix using Cholesky
    % Decomposition
    sqrtCovarEst = chol(covarEst, "lower");
    % Find all 2n+1 = 31 sigma points
    for i = 1:n
        Xt(:, i+1) = uEst + (sqrt(n+lambda) * sqrtCovarEst(:,i));
        Xt(:, i+n+1) = uEst - (sqrt(n+lambda) * sqrtCovarEst(:, i));
    end
    % Rotation matrix to go from IMU to camera frame
    Rb2c = transpose([0.707 -0.707 0; -0.707 -0.707 0; 0 0 -1]);%[0.7071, -0.7071, 0; -0.7071, -0.7071, 0; 0, 0, -1]; % eul2rotm([-pi/4,pi,0])
    % Translation matrix to go from IMU to camera frame
    Tc2b = [0.0283; -0.0283; -0.0300]; % eul2rotm([-pi/4,pi,0]) * [-0.04, 0.0, -0.03]';
    % Skew Symmetric version of Translation matrix
    Skew_Tc2b = [0, -Tc2b(3), Tc2b(2); Tc2b(3), 0, -Tc2b(1); -Tc2b(2), Tc2b(1), 0];
    % Rotation matrix to go from camera to IMU frame
    Rc2b = Rb2c';
    % Rotation matrix to go from IMU to world frame
    Rb2w = eul2rotm(uEst(4:6)'); 
    % Initialise update vector matrix to accomodate all sigma points after
    % being passed through non-linear update model
    Zt = zeros(3, ((2*n)+1));
    % Find all update vectors through iteration
    for i = 1:((2*n)+1)
        % x1 = Xt(1:3, i);
        % x2 = Xt(4:6, i); 
        x3 = Xt(7:9, i);
        % x4 = Xt(10:12, i);
        % x5 = Xt(13:15, i);
        l_x2_x3 = Rb2w \ x3;
        % Non-linear update model
        Zt(:,i) = (Rb2c * l_x2_x3) - (Rb2c * Skew_Tc2b * (Rc2b * z_t(4:6,1))) + vt;
    end
    % Find weights for mean and covariances 
    Wm = zeros(1, ((2*n)+1));
    Wm(1) = lambda / (n + lambda);
    Wc = zeros(1, ((2*n)+1));
    Wc(1) = ((lambda) / (n + lambda)) + (1 - (alpha^2) + beta);
    for i = 2:((2*n)+1)
        Wm(i) = 1 / (2 * (n + lambda));
        Wc(i) = 1 / (2 * (n + lambda));
    end
    % Initialise estimated measurement update vector
    Zut = zeros(3,1);
    % Find weighted sum of all update vectors obtained by passing all sigma 
    % points through non-linear update model
    for i = 1:((2*n)+1)
        Zut = Zut + (Wm(i) * Zt(:,i));
    end
    % Initialise covariance and cross-covariance matrices
    Ct = zeros(15,3);
    St = zeros(3,3);
    % Computed weighted sum of covariance and cross-covariance matrices
    for i = 1:((2*n)+1)
        Ct = Ct + (Wc(i) * (Xt(:,i) - uEst) * (Zt(:,i) - Zut)');
        St = St + (Wc(i) * (Zt(:,i) - Zut) * (Zt(:,i) - Zut)');
    end
    St = St + R;

    % Compute Kalman Gain
    K = Ct / St;
    % Compute current mean
    uCurr = uEst + K * (z_t(1:3) - Zut);
    % Compute current covariance
    covar_curr = covarEst - (K * St * K');
end

