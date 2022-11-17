%% load system from other directory
path_to_model = 'model/load_model.m';
run(path_to_model);
% determine the row of the matrix
n = length(A);
global A B C n W M 
%% Define integral observer parameters
betas = [1 1];
T = 5;


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


% % transpose Pi1 and Pi2
% Pi1 = Pi1.';
% Pi2 = Pi2.';
% 
% %% calculate G matrix's
% G1 = P1 * C.';
% G2 = P2 * B;

