
% --------------------------------------------------
% File for RANDOM NOISE GENERATION 
%---------------------------------------------------


% Standard Deviation of noise in the differential states of the system 
state_sigma = [sqrt(10^(-5))];      %  (1% of steady state values)
meas_sigma = [0.01];

% Noise generation in states 
for i = 1:1:1000 
    wk(1,i) = state_sigma(1)*randn;
end

% Noise generation in measurements
for i = 1:1:1000 
    vk(1,i) = meas_sigma(1)*randn;
end

save noise_generator