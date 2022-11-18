clear all;

% run("..\model\load_model.m");  % uruchom skrypt z folderu 'symbolic'
run("model_params.m");

sizeA = size(A);
system_size = sizeA(2);  % liczba zmiennych stanu

syms t tau T real positive

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wzór ogólny (około 10 minut liczenia dla najogólniejszego przypadku)
syms beta [1 system_size] real positive
beta(1) = 1;
beta(2) = 1;
syms J
P1 = [];
P2 = [];
for i = 1:system_size
    W = [A,    beta(i)*B*B.'
         C.'*C, -A.']
    Phi = simplify(expm(W * t),'steps',1000)
    Phi11 = Phi(1:system_size,1:system_size)
    Phi21 = Phi(system_size+1:2*system_size,1:system_size)
    int_internal = simplify(expm(-A.'*(T-tau))*C.'*C*subs(Phi11,t,tau), 'steps', 1000)
    Mi = int(int_internal, tau, 0, T)
    Mi = simplify(Mi, 'steps', 1000)
    ei = zeros([system_size, 1])
    ei(i) = 1
    iMi = inv(Mi);
    iMi = simplify(iMi, 'steps', 1000)
    p1 = Phi11 * iMi * ei
    p1 = simplify(p1, 'steps', 1000)
    p2 = Phi21 * iMi * ei
    p1 = simplify(p2, 'steps', 1000)
    P1 = [P1, p1]  % to jest tak na prawdę transpozycja P1
    P2 = [P2, p2]
end

G1 = P1' * C.';  % transponuję zmienną P1, ale w teorii jest to macierz P1 nietransponowana
G2 = P2' * B;

G1 = simplify(G1,'steps',1000);
G2 = simplify(G2,'steps',1000);

% przykladowe przebiegi składowych macierzy G1 i G2

Treal = 0.5;
figure(1)
figure(2)
for i = 1:system_size
    G1tempf = matlabFunction(subs(G1(i), T, Treal));
    G2tempf = matlabFunction(subs(G2(i), T, Treal));
    figure(1)
    subplot(system_size,1,i)
    plot(0:0.01:Treal,G1tempf(0:0.01:Treal))
    title(["G1"+num2str(i)])
    figure(2)
    subplot(system_size,1,i)
    plot(0:0.01:Treal,G2tempf(0:0.01:Treal))
    title(["G2"+num2str(i)])
end


% dodatkowy FOR do liczenia wskaźnika jakości, bardzo wydłuża czas obliczeń
for i = 1:system_size
    J = G1(i,1).^2 + beta(i)*G2(i,1).^2
end
J = simplify(J,'steps',1000)
J = int(J, t, 0, T)
J = simplify(J,'steps',1000)
Jf = matlabFunction(J)

figure(3)
plot(0:0.01:10, Jf(0:0.01:10))
title("Funkcja kosztu w zależności od T dla \beta_1="+num2str(double(beta(1)))+...
    " \beta_2="+num2str(double(beta(2))))
xlabel("T")
ylabel("Wartość funkcji kosztu")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  przykład strona 318
% beta = [1 1];
% W = [A,    B*B'
%     C'*C, -A']
% Phi = expm(W*t)
% Phi = simplify(Phi,'steps',1000)
% Phi11 = Phi(1:system_size,1:system_size)
% Phi21 = Phi(system_size+1:2*system_size,1:system_size)
% M = int(subs(Phi11.',t,tau)*C.'*C*expm(A*tau),tau,0,T) * expm(-A*T)
% M = simplify(M,'steps',1000)
% iM = inv(M)
% G1 = iM * Phi11.' * C.'
% G2 = iM * Phi21.' * B
% 
% J = G1(1,1)^2 + G1(2,1)^2;
% J = J + G2(1,1)^2 + G2(2,1)^2;
% J = int(J, t, 0, T)
% Jsimp = simplify(J,'steps',1000)