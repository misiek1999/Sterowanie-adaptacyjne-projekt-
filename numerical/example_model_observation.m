example_params;

if (length(A)-rank(obsv(A,C))) == 0  % sprawdzenie obserwowalności
    sprintf("System obserwowalny")
else
    sprintf("System nie obserwowalny")
    return
end

options = optimset('MaxIter',100,'PlotFcns',@optimplotfval);

params0 = [1, 0.2, 0.2 0.2];  % T0, beta1, beta2, ..., beta_n_B - długość zależy od wymiaru macierzy B
optifhand = @(p)optiG1G2(A, B, C, p);

[x,fval,exitflag,output] = fminsearch(optifhand, params0, options);

T = x(1);
betas = x(2:end);
input_noise_power = 0.01;
output_noise_power = 0.01;

time_step = 0.01;
time_samples = 0:time_step:T;
time_samples_count = length(time_samples);
n = length(A);
outputsize = size(C);
outputsize = outputsize(1);
inputsize  = size(B);
inputsize  = inputsize(2);
x0 = ones(1, n);

[G1, G2] = getG1G2(A, B, C, betas, time_samples);
