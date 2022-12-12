% przykład obserwacji rzeczywistego obiektu 
% silnik DC obsługujący serwomechanizm
% obserwowane stany:
% 1) prąd silnika (z czujnika w sterowniku mocy)
% 2) prędkość obrotowa (z tachoprądnicy)
% 3) położenie wału (z enkodera)
% wejście jest mierzone w zakresie -1 do 1 co odpowiada napięciom -12, +12 V
% wyjściem jest położenie

% TODO: Dość poważny problem - dla dużych wartości własnych A algorytm 
% wykrzacza się przy liczeniu macierzy Gramma M
% na tamtym etapie liczy się e do macierzy, co może dać Inf a potem przy 
% całkowaniu NaN

load example_measurements.mat  % zawiera też linearyzację A, B, C, D

if (length(A)-rank(obsv(A,C))) == 0  % sprawdzenie obserwowalności
    sprintf("System obserwowalny")
else
    sprintf("System nie obserwowalny")
    return
end

% wyznaczenie optymalnych parametrów obserwatora dla 

options = optimset('MaxIter',100,'PlotFcns',@optimplotfval);

params0 = [1, 0.2, 0.2, 0.2];  % T0, beta1, beta2, ..., beta_n_B - długość zależy od wymiaru macierzy B
optifhand = @(p)optiG1G2(A, B, C, p);

[x,fval,exitflag,output] = fminsearch(optifhand, params0, options);

T = x(1);
betas = x(2:end);

% time_step = 0.01;
time_step = time(2) - time(1);
time_samples = 0:time_step:T;
time_samples_count = length(time_samples);
n = length(A);
n_C = size(C);
n_C = n_C(1);
n_B = size(B);
n_B = n_B(2);
x0 = zeros(1, n);

[G1, G2] = getG1G2(A, B, C, betas, time_samples);
plotG1G2(G1,G2,time_samples,A,B,C)

% symulacja obserwatora

inzeros = zeros(time_samples_count, n_C);
outzeros = zeros(time_samples_count, n_B);
extended_output = cat(n_C, inzeros, y);
extended_input  = cat(n_B, outzeros, u);

observator_states = [];
for i = 1:length(time)
    observation_data = extended_output(i:i+time_samples_count,:);
    control_data     = extended_input(i:i+time_samples_count,:);
    observator_state = observatorFnc(observation_data, control_data, time_samples, G1, G2);
    observator_states = [observator_states; observator_state];
end
