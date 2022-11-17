clear all

% run("..\model\load_model.m");  % run in 'symbolic' folder as workspace
run("model_params.m");

sizeA = size(A);
system_size = sizeA(2);

syms t tau T 
% syms beta [1 system_size]
beta = [1 1];

%%

%%%%%%%%%%%%
% M0 = int(expm(A'*tau)*C'*C*expm(A*tau), tau, 0, T);  % M0

% G1 = expm(A*T)*inv(M0)*expm(A'*t)*C';  % G1(t)
% G2 = int(subs(G1, t, s)*C*expm(A*(s-t)), s, 0, t);  % G2(t)

%%%%%%%%%%%%%%
% P1 = [];
% P2 = [];
% for i = 1:system_size
%     W = [A,    beta(i)*B*B'
%          C'*C, -A']
%     Phi = simplify(expm(W * t),'steps',1000)
%     Phi11 = Phi(1:system_size,1:system_size)
%     Phi21 = Phi(system_size+1:2*system_size,1:system_size)
%     int_internal = simplify(expm(-A'*(T-tau)*C'*C*subs(Phi11,t,tau)))
%     Mi = int(int_internal, tau, 0, T)
%     ei = zeros([system_size, 1])
%     ei(i) = 1
%     iMi = inv(Mi);
%     p1 = Phi11 * iMi  * ei
%     p2 = Phi21 * iMi * ei
%     P1 = [P1, p1]
%     P2 = [P2, p2]
% end
% 
% G1 = P1 * C';
% G2 = P2 * B;

%%%%%%%%%%%%%
% MI = expm(-A'*t) * int(expm(A'*tau)*C'*C*Phi, tau, 0, T);
% 
% Phi11 = ;
% Phi21 = ;

% p

%%%%%%%%%%%%%
W = [A,    B*B'
    C'*C, -A']
Phi = expm(W*t)
Phi11 = Phi(1:system_size,1:system_size)
Phi21 = Phi(system_size+1:2*system_size,1:system_size)
M = int(subs(Phi11.',t,tau)*C.'*C*expm(A*tau),tau,0,T) * expm(-A*T)
M = simplify(M,'steps',1000)
iM = inv(M)
G1 = iM * Phi11.' * C.'
G2 = iM * Phi21 * B