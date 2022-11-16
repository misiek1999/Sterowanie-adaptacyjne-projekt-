%% this script create linear model of dynamic damper
% parameters of dynamic damper
m = 10; % mass of object 
c = 1;  % stiffness of object
k = 5;  % damping of object

%% state space matrix
A = [0 1
    -k/m -c/m];

B = [0;1];

C = [1 0]; % we can meansure only mass position


