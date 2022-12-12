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
% A = [-1 0 0
%     1 -1 0
%     0 1 -1];
% B = [1;0;0];
% C = [1 1 1];
% 
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

% %% real tank model 
% A  =[-0.0220329569615680 0 0
%     0.0225287188113119 -0.0263853892550530 0
%     0 0.0182410642165413 -0.0169457838246552];
% B = [0.0113378684807256; 0;0];
% C = [0 0 1];