clc
clear all

electrode_global;

global wk vk R1 i C xkpredhat Pk Yk w1 w2 w3 b out C_est res_matrix_hybrid

% Loading saved noise from noise_generator file

load('noise_generator.mat', 'vk');
load('noise_generator.mat', 'wk');

% Sampling time in minutes (6sec)

sample_T = 0.1;

% step length (used while updating NN parameters)

alpha = 0.1;

N_samples = 1000; % number of time instants

% Parameter values

F = 96487;         % C/mol
R = 8.314;           % J/mol K
T = 298.15;        % K
phi_eq1 = 0.42;    % V
phi_eq2 = 0.303;   % V
rho = 3.4;         % g/cm^3
W = 92.7;          % g/mol
V = 1*10^(-5);      % cm
i_app = 1*10^(-5);  % A/cm^2
i01 = 1*10^(-4);    % A/cm^2
i02 = 1*10^(-8);    % A/cm^2


% finding steady state solution

yss = fsolve(@(X) electrode_dae(X), [0.35024 0.4071], optimoptions('fsolve','Display','iter'));

% Standard Deviation of noise in the differential states of the system 

state_sigma = [sqrt(10^(-5))];      %  (1% of steady state values)
meas_sigma = [0.01];

% Measurement noise covariance matrix 

R1 = diag([0.0001]);

% State noise covariance matrix
Q1 = diag([state_sigma(1)^2]);
Q = Q1;



% Measurement bias 

b = 0.1; 

% Measurement matrix

C = [0 1];          % only y2 is measured (C matrix for plant simulation)

C_est = [0 1 0 0 0]; % only y2 is measured (C matrix for DAE-EKF-NN simulation)

x_ss = [0.35024 0.4071]';    % Initial state used in the simulation study
Xk = x_ss; 


% Defining Initial weights

w1 = 0.2;
w2 = 0.4;
w3 = 0.6;

% Augmented noise matrix

 % Q_aug = [Q1 0 0 0 0; 0 Q1 0 0 0; 0 0 w1 0 0; 0 0 0 w2 0; 0 0 0 0 w3];

Q_aug = [Q1 0 0 0 0; 0 Q1 0 0 0; 0 0 10 0 0; 0 0 0 10 0; 0 0 0 0 10];

% Initial state of the estimator

Xkhat  = [0.5322; 0.4254; w1; w2; w3];   % Initial state of the estimator
% X = [0.35024; 0.4071; w1; w2; w3];

X = [0.35024; 0.4071]; % Initial state of the plant

% Neural network for the measurement model

NN = electrode_NN(Xkhat);

% Measurement model function

h = @(Xkhat) C_est*Xkhat + electrode_NN(Xkhat);         % Measurement model for the estimator

Yk(:,1) = C*X + 0*NN + vk(:,1);                % plant measurement model 

% Defining mass matrix for solving the DAE system

M = [1 0; 0 0];   %(2x2)

options = odeset('Mass',M);



% Initial covariance matrix

% Pk = [0.005  0 0 0 0; 0 0.005 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];

%Pk = [0.005  0 0 0 0; 0 0.005 0 0 0; 0 0 0.005 0 0; 0 0 0 0.005 0; 0 0 0 0 0.005];   
% Since initial error covariance is higher

Pk = [10^-1  0 0 0 0; 0 10^-1 0 0 0; 0 0 10^-1 0 0; 0 0 0 10^-1 0; 0 0 0 0 10^-1];
Pk_ekf = [0.005  0; 0 0.005];

% For plotting the 1st point (initial guess)
res_matrix_hybrid(1,:) = [1 X(1) X(2) Xkhat(1) Xkhat(2) Yk(1,1) 0 0 0];

for i = 2:1:500
    
    % ------Plant calculations-------
    func = @(t,Xk) electrode_dae(Xk);
    [t,Xt] = ode15s(func, [0,15], Xk, options);
    Xk = Xt(length(t),:);       % Generate x(k+1)
    Xk(1,1) = Xk(1,1) + wk(1,i);

    X = [Xk(1,1); Xk(1,2)]; 
    % X = [Xk(1,1); Xk(1,2); w1; w2; w3];  % In plant simulation neural
    % network model should not come into picture (NN should come only in
    % the estimator design)

    if i <=200

    Yk(:,i) = C*X   + vk(1,i);  % No bias in the estimates

    else

    Yk(:,i) = C*X + 0.5*b + vk(1,i);  %  bias in the estimates    

    end

     
    % ------Soft sensor calculations-----
    
    % Extended Kalman Filter 
    % Step - 1: Prediction step
    xkhat = [Xkhat(1,1), Xkhat(2,1)];
    [t,Xj] = ode15s(func, [0,15], xkhat, options);
    xkpredhat = Xj(end,:)';
    
    % Generating linear perturbation model at the current steady state 
    if i <= 100
    
    Z_vec = xkhat;
    [A1_mat] = electrode_Jacob(Z_vec);   
    phy1 = expm(A1_mat*sample_T);                      
    Pk_ekf = phy1 * Pk_ekf * phy1' + Q;      
    
    Vk1 = R1 + C*Pk_ekf*C';
    
    % Calculating Kalman Gain (Lk)
    Lk1 = Pk_ekf*C'*inv(Vk1); 
    
    % Step - 2: Update step
    ek1 = Yk(:,i)-C*xkpredhat;           % error matrix
    Xkhat = xkpredhat + Lk1*ek1;          % updating mean and covariance in this step
    Pk_ekf = (eye(2,2)-Lk1*C)*Pk_ekf;

    % Simulation result
    res_matrix_hybrid(i,:) = [i Xk(1) Xk(2) Xkhat(1) Xkhat(2) Yk(1,i) 0 0 0];
    
    else
    
    Xkhat = [Xkhat(1,1);  Xkhat(2,1); w1; w2; w3];    
    Z_vec = Xkhat;
    [A_mat] = electrode_hybrid_Jacob(Z_vec);               % A matrix  (Jacobian of the function F)
    phy = expm(A_mat*sample_T);                            % discrete model
    Pk = phy * Pk * phy' + Q_aug;     
    


    % Calculating Kalman Gain (Lk)
    
    z0 = [xkpredhat(1,1); xkpredhat(2,1); w1; w2; w3];
    [H] = electrode_nn_Jacob(h, z0);                        % Jacobian of the measurement model
    Vk = R1 + H*Pk*H';
    Lk = Pk*H'*inv(Vk); 
     
    % Step - 2: Update step
    
    X_pred_hat = [xkpredhat(1,1); xkpredhat(2,1); w1; w2; w3];
    ek = Yk(:,i)-(C_est*X_pred_hat + electrode_NN(X_pred_hat));   % error matrix
    Xkhat = X_pred_hat + Lk*ek;                                                 % updating mean and covariance in this step
    Pk = (eye(5,5)-Lk*H)*Pk;

    w1 = Xkhat(3);  % updated weight   (weights are updated through augmentation approach, not through training)
    w2 = Xkhat(4);  % updated weight  (weights are updated through augmentation approach, not through training)
    w3 = Xkhat(5);  % updated weight  (weights are updated through augmentation approach, not through training)


    % % Updating the NN parameters - w1, w2, w3
    % df_by_dw1 = -ek*w3*(1-(out.^2))*Xkhat(1,1);
    % df_by_dw2 = -ek*w3*(1-(out.^2))*Xkhat(2,1);
    % df_by_dw3 = -ek*out;
    % w1 = w1 - (alpha*df_by_dw1);
    % w2 = w2 - (alpha*df_by_dw2);
    % w3 = w3 - (alpha*df_by_dw3);
    
    % Simulation result
    res_matrix_hybrid(i,:) = [i Xk(1) Xk(2) Xkhat(1) Xkhat(2) Yk(1,i) Xkhat(3) Xkhat(4) Xkhat(5)];
    
    end
    
    % save results_hybrid_model
end


% Plotting the states: y1 and y2
SetGraphics 

figure(1), 

subplot(211),
plot(res_matrix_hybrid(:,1),res_matrix_hybrid(:,2),'k',res_matrix_hybrid(:,1),res_matrix_hybrid(:,4),'r'),
grid on
ylabel('y1'), 
legend('True State', 'Estimated State')

subplot(212),
plot(res_matrix_hybrid(:,1),res_matrix_hybrid(:,3),'k',res_matrix_hybrid(:,1),res_matrix_hybrid(:,6),'b',res_matrix_hybrid(:,1),res_matrix_hybrid(:,5),'r'),
grid on
ylabel('y2'), 
xlabel('Sampling Instant'),
legend('True (Actual) state','Measurement', 'Estimated State')


figure(2),

plot(res_matrix_hybrid(:,1),res_matrix_hybrid(:,3),'k',res_matrix_hybrid(:,1),res_matrix_hybrid(:,6),'b',res_matrix_hybrid(:,1),res_matrix_hybrid(:,5),'r'),
grid on
ylabel('y2'), 
xlabel('Sampling Instant'),
legend('True (Actual) state','Measurement', 'Estimated State')


% mean of error matrix of y1
mean(res_matrix_hybrid(:,2)-res_matrix_hybrid(:,4));

% mean of error matrix of y2
mean(res_matrix_hybrid(:,3)-res_matrix_hybrid(:,5));

% SSE -  sum of squared of estimation error  for y1
SSE1  = (res_matrix_hybrid(:,2)- res_matrix_hybrid(:,4))'*(res_matrix_hybrid(:,2)-res_matrix_hybrid(:,4))

% SSE -  sum of squared of estimation error  for y2
SSE2  = (res_matrix_hybrid(:,3)-res_matrix_hybrid(:,5))'*(res_matrix_hybrid(:,3)-res_matrix_hybrid(:,5))

% RMSE -  y1
rmse_y1 = sqrt(mean((res_matrix_hybrid(:,2)- res_matrix_hybrid(:,4)).^2));

% RMSE -  y2
rmse_y2 = sqrt(mean((res_matrix_hybrid(:,3)- res_matrix_hybrid(:,5)).^2));


% Estimation Error Plots for y1 and y2
figure(3),
subplot(211),
plot(res_matrix_hybrid(:,1),res_matrix_hybrid(:,2)-res_matrix_hybrid(:,4),'m'),
grid on
ylabel('Estimation error - y1'), 
xlabel('Sampling Instant')

subplot(212),

plot(res_matrix_hybrid(:,1),res_matrix_hybrid(:,3)-res_matrix_hybrid(:,5),'m'),
grid on
ylabel('Estimation error - y2'), 
xlabel('Sampling Instant')
