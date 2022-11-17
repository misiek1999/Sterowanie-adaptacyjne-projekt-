function state = fcn(observation_data, control_data, time_stamps,G1_vec, G2_vec, g1_m, g2_m)
%% calculate vectors
G1xY = [];
G2xU = [];
itr = 1;
for i = time_stamps
    G1xY = [G1xY G1_vec(:, itr) * observation_data(itr)];
    G2xU = [G2xU G2_vec(:, itr) * control_data(itr)];
    itr = itr + 1;
end

disp('siema')

%% make itegration of both vector
state_obs = trapz(time_stamps, G1xY, 2);
control_obs = trapz(time_stamps, G2xU, 2);

state = state_obs + control_obs;
