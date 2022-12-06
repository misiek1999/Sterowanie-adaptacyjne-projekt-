%% load system from other directory
path_to_model = 'model/load_model.m';
run(path_to_model);
%% add this folder to matlab path 
path_to_main_folder = pwd;
addpath(path_to_main_folder);
%% Argument initialization
x0 = [0 0 0];

% determine the row of the matrix
n = length(A);
n_C = size(C);
n_C = n_C(1);
n_B = size(B);
n_B = n_B(2);

global A B C n W M 
global T time_step time_stamps time_samples_count
global G1_vec G2_vec
%% Define integral observer parameters
betas = [1 1];
betas = [1 1 1];
T = 2;
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
    G1_vec(:,:, itr) = getG1(t);
    G2_vec(:,:, itr) = getG2(t); 
    itr = itr + 1;
end

%% plot G1 G2 matrices

% G1
figure
for i = 1:n
    for j = 1:n_C
        subplot(n, n_C, j + (i-1)*n_C)
        plot(time_stamps, squeeze(G1_vec(i,j,:)))
        grid on;
        title("G1["+num2str(i)+','+num2str(j)+']')
    end
end

% G2
figure
for i = 1:n
    for j = 1:n_B
        subplot(n, n_B, j + (i-1)*n_B)
        plot(time_stamps, squeeze(G2_vec(i,j,:)))
        grid on;
        title("G2["+num2str(i)+','+num2str(j)+']')
    end
end