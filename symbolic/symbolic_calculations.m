clear all;

% run("..\model\load_model.m");  % uruchom skrypt z folderu 'symbolic'
run("model_params.m");

sizeA = size(A);
system_size = sizeA(2);  % liczba zmiennych stanu
sizeC = size(C);
output_size = sizeC(1);  % liczba wyjsc
sizeB = size(B);
input_size = sizeB(2);  % liczba wejsc

syms t tau T real positive
% assume(t>=0);
% assume(tau>=0);
% assume(T>0);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wzór ogólny (około 10 minut liczenia dla najogólniejszego przypadku)
syms beta [1 system_size] real
for i=1:system_size
    assume(beta(i)>0);
end
% beta(1) = 0;
% beta(2) = 0;
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
    p2 = simplify(p2, 'steps', 1000)
    P1 = [P1, p1]  % to jest tak na prawdę transpozycja P1
    P2 = [P2, p2]
end

G1 = P1' * C.'  % transponuję zmienną P1, ale w teorii jest to macierz P1 nietransponowana
G2 = P2' * B

G1 = simplify(G1,'steps',1000)
G2 = simplify(G2,'steps',1000)

% przykladowe przebiegi składowych macierzy G1 i G2
%%
Treal = 0.5;
betasreal = [0.5 0.5];

figure(1)
for i = 1:system_size
    for j = 1:output_size
        G1tempf = subs(G1(i,j), T, Treal);
        for k = 1:system_size
            G1tempf = subs(G1tempf, beta(k), betasreal(k));
        end
        G1tempf = matlabFunction(G1tempf);
        subplot(system_size,output_size,j+(i-1)*output_size)
        plot(0:0.01:Treal,G1tempf(0:0.01:Treal))
        grid on;
        title(["G1 "+num2str(i)+num2str(j)])
    end
end

figure(2)
for i = 1:system_size
    for j = 1:input_size
        G2tempf = subs(G2(i,j), T, Treal);
        for k = 1:system_size
            G2tempf = subs(G2tempf, beta(k), betasreal(k));
        end
        G2tempf = matlabFunction(G2tempf);
        subplot(system_size,output_size,j+(i-1)*output_size)
        plot(0:0.01:Treal,G2tempf(0:0.01:Treal))
        grid on;
        title(["G2 "+num2str(i)+num2str(j)])
    end
end

%%
% dodatkowy FOR do liczenia wskaźnika jakości, bardzo wydłuża czas obliczeń
% syms betashat [1 system_size] real positive
for i = 1:system_size
    J = sum(G1(i,:).^2) + 1*sum(G2(i,:).^2)  % koszt liczymy jakby na wejściu były zakłócenia
end
J = simplify(J,'steps',1000)
J = int(J, t, 0, T)
% J = simplify(J,'steps',1000)

%%
Jfs    = [];
Trange = 0:0.1:10;

for betasfortest = 0.01:0.01:1.01
    Jf = J;
    for k = 1:system_size  % liczba parametrów beta
        Jf = subs(J, beta(k), betasfortest);
    end
    Jf = matlabFunction(Jf);
    Jfs = [Jfs; Jf(Trange)];
end
%%
sizeJfs = size(Jfs);
figure(4)
X = meshgrid(Trange);
X = X(1:101,:);
Y = meshgrid(0.0:0.01:1.0)';
Y = Y(1:101,:);
mesh(X, Y, Jfs)
grid on;
xlabel("T")
ylabel("$\hat\beta(1)=\hat\beta(2)=...=Y$", 'Interpreter', 'latex')
zlabel("Wartość funkcji kosztu")
zlim([0, 5])
caxis([0 3])
colorbar
