%% load system from other directory
path_to_model = 'model/load_model.m';
run(path_to_model);

x0 = [0 0];

% determine the row of the matrix
n = length(A);
global A B C n W M 
global T time_step time_stamps time_samples_count
global G1_vec G2_vec
%% Define integral observer parameters
betas = [1 1];
T = 5;
time_step = 0.01;
time_stamps = 0 : time_step :T;
time_samples_count = length(time_stamps);

%% find best G1 and G2 matrix
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

%%  Create vector with G1 and G2
G1_vec = [];
G2_vec = [];
[g1_n g1_m] = size(getG1(0));
[g2_n g2_m] = size(getG2(0));
itr = 1;
for t = time_stamps
    G1_vec = [G1_vec getG1(t)];
    G2_vec = [G2_vec getG2(t)]; 
    itr = itr + 1;
end
