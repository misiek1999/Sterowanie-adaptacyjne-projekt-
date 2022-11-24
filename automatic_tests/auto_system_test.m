%% create plot of intergral observer for diffrent T and Beta
% load initialization script
run('../initialization.m');

%% create vector with T and noise power compramision
T_list = [0.05 0.1 0.2 0.5 1 2 5 10 20 50];

noise_power_list= [0 1e-7 1e-6 1e-5 1e-4 1e-3];


%% run simulation to colect data
%% iterate for each T
rms_list = zeros([length(T_list) length(noise_power_list)]);
itr_T = 1;  % create T iterator 
itr_noise = 1;  % create noise iterator 
for T = T_list
    % calculate new time steps
    time_stamps = 0 : time_step :T;
    time_samples_count = length(time_stamps);
    % find best G1 and G2 matrix
    % initialize P matrix's
    Pi1 = [];
    Pi2 = [];
    M = {}; % cell array
    W = {}; 
    % find all Mi matrix's and save it to M array
    for i = 1:n
        % calculate W matrix
        Wi = [A,     betas(i)*B*B.'
             C.'*C, -A.'    ];
        W{i} = Wi;
        % find Grama matrix's
        Mi_integral_fun = @(tau) expm(-A.'*(T-tau))*C.'*C*getPhi11(tau, i);
        Mi = integral(@(tau)Mi_integral_fun(tau), 0, T, 'ArrayValued', true);
        M{i} = Mi;
    end
    itr = 1;
    for t = time_stamps
        G1_vec(:,:, itr) = getG1(t);
        G2_vec(:,:, itr) = getG2(t); 
        itr = itr + 1;
    end
    % check diffrent noise power 
    itr_noise = 1;  % reset noise iterator
    for noise_power= noise_power_list
        collected_data = sim('model_sim_2020b.slx');    % run simulation to collect data
        rms_list(itr_T, itr_noise) = collected_data.error_sum(end); % save collected data to array
        itr_noise = itr_noise + 1;%incraese noise iterator
    end
    
    itr_T = itr_T + 1;  % increase T iterator
end

%% plot collected data
figure();
h = gca;
surf(T_list, noise_power_list, rms_list.');
grid on;
colorbar;
title('Porównanie jakości dla różnych parametrów obserwatora i obiektu')
xlabel('Długość okna obserwatora [s]');
ylabel('Wsółczynnik zaszumienia sygnału')
zlabel('Całka z kwadratu błedu')
set(h,'zscale','log')
set(h,'xscale','log')
set(h,'yscale','log')