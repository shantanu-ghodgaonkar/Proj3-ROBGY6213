 function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
    %% Parameter Definition
    % uPrev - is the mean of the prev state
    %covarPrev - covar of the prev state
    %angVel - angular velocity input at the time step
    %acc - acceleration at the timestep
    %dt - difference in time 
% MY IMPLEMENTATION START -------------------------------------------------
    nDash = 27;
    alpha = 0.001;
    k = 1;
    beta = 2;
    lambdaDash = (alpha^2) * (nDash + k) - nDash;
    uAugPrev = [uPrev; zeros(12,1)];
    % Set Q as 0.01
    Q = 0.01;
    % Calculate Qd
    Qd = dt*(eye(12) * Q);
    PAug = [covarPrev, zeros(15,12); zeros(12,15), Qd];
    sqrtPAug = chol(PAug);
    XAug = zeros(nDash, ((2*nDash)+1));
    XAug(:,1) = uAugPrev;
    
    for i = 2:1:(nDash+1)
        XAug(:,i) = uAugPrev + sqrt(nDash + lambdaDash) * sqrtPAug(:,i-1);
        XAug(:,(i+nDash)) = uAugPrev - sqrt(nDash + lambdaDash) * sqrtPAug(:,i-1);
    end
    % disp("DEBUG POINT");
    Xt = zeros(15, ((2*nDash)+1));
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
        Xt(:,i) = [x1 + (dt * f1); 
            x2 + (dt * f2) - (dt * G_inv * ng); 
            x3 + (dt * f3) - (dt * R * na); 
            x4 + (f4); ... x4 + (dt * f4); 
            x5 + (f5)]; ... x5 + (dt * f5)];
    end
    % disp("DEBUG POINT");

    WmDash = [(lambdaDash / (nDash + lambdaDash)), ( ones(1, 2*nDash)*(1 / (2 * (nDash + lambdaDash))) )];
    WcDash = [(( lambdaDash / (nDash + lambdaDash) ) + (1 - (alpha^2) + beta)), ( ones(1, 2*nDash)*(1 / (2 * (nDash + lambdaDash))) )];
    uEst = zeros(15,1);
    for i = 1:1:((2*nDash)+1)
        uEst = uEst + (WmDash(i) * Xt(:,i));
    end
    covarEst = zeros(15,15);
    for i = 1:1:((2*nDash)+1)
        covarEst = covarEst + (WcDash(i) * ((Xt(:,i) - uEst) * (Xt(:,i) - uEst)'));
    end
    % disp("DEBUG POINT");
% MY IMPLEMENTATION END ---------------------------------------------------
end

