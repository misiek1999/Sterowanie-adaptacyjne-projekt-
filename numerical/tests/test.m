addpath("../");

% % parameters of dynamic damper
% m = 10; % mass of object 
% c = 1;  % stiffness of object
% k = 5;  % damping of object
% A = [0 1
%     -k/m -c/m];
% 
% B = [0;1];
% 
% C = [1 0]; % we can meansure only mass position

% A = [-1 0 0
%     1 -1 0
%     0 1 -1];
% B = [1;0;0];
% C = [1 1 1];

% A = [0 1
%     0 0];
% B = [0
%     1];
% C = [2 0];

% A = [0 1
%     0 0];
% B = [0 1
%     1 0];
% C = [1 0
%      0 1
%      1 1];

A = [-0.05 1 0
    0  -0.9 1
    0 0 -0.1];
B = [2
    1
    1];
C = [0.5 0.5 0
    0 0 2];

%%
time_step = 0.01;
startpoint = [1 0.5 0.5 0.5];  % T0, beta1_0, beta2_0
ploton = false;
searchiter = 150;

% [optiG1, optiG2, time_samples, optibetas, optiT, x, fval, exitflag, output] = getoptiG1G2(A, B, C, time_step, startpoint, ploton, searchiter);

test_betas = [0.1 0.5 1 2 10];
no_betas = length(test_betas);
test_Ts = [0.5 1 2 10];
no_Ts = length(test_Ts);

output_noise_powers = [0 1e-3 1e-2];
no_output_noise_powers = length(output_noise_powers);
input_noise_powers = [0 1e-3 1e-2];
no_input_noise_powers = length(input_noise_powers);

n = length(A);
outputsize = size(C);
outputsize = outputsize(1);
inputsize  = size(B);
inputsize  = inputsize(2);
x0 = ones(1, n);

i=0;
error_table = zeros(length(test_betas), length(test_Ts), length(output_noise_powers), length(input_noise_powers));
for test_beta = test_betas
    i=i+1
    j=0;
    for test_T = test_Ts
        j=j+1;
        time_samples = 0:time_step:test_T;
        time_samples_count = length(time_samples);
        [G1, G2] = getG1G2(A, B, C, test_beta*ones(1,n), time_samples);
        errors = zeros(length(output_noise_powers), length(input_noise_powers));
        k=0;
        for output_noise_power = output_noise_powers
            k=k+1;
            l=0;
            for input_noise_power = input_noise_powers
                l=l+1;
                out    = sim("test_model.slx");
                simtime = out.tout;
                error  = out.error.signals.values;
                % dzielone przez długość symulacji minus długość okna w celu normalizacji
                errors(k, l) = trapz(simtime(simtime>test_T), sum(error(simtime>test_T),2).^2) / (100 - test_T);
                states = out.states;
            end
        end
        error_table(i, j, :, :) = errors;
    end
end

%%
markers = ['rx', 'r.', 'ro'];
figure();
h = gca;
hold on;
i=0;
for test_beta = test_betas
    i=i+1
    j=0;
    for test_T = test_Ts
        j=j+1;
        k=0;
        for output_noise_power = output_noise_powers
            k=k+1;
            l=0;
            for input_noise_power = input_noise_powers
                l=l+1;
                plot3(test_beta, test_T, error_table(i, j, k, l), markers(k));
            end
        end
    end
end
xlabel("\beta");
ylabel("T")
grid on;
set(h,'zscale','log')
set(h,'xscale','log')
set(h,'yscale','log')