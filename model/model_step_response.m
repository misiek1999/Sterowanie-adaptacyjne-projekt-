%% init model 
load_model;

%% craete system model from state space matrix
sys = ss(A,B,C, 0);
step_response = step(sys);

plot(step_response); grid on;
xlabel('Time')
ylabel('Mass position')
title('Model response for step')