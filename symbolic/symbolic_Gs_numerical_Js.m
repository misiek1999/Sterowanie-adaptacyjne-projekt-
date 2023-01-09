% clear all;

% run("..\model\load_model.m");  % uruchom skrypt z folderu 'symbolic'
% run("model_params.m");
A=A-mod(A,0.01); B=B-mod(B,10);
if (length(A)-rank(obsv(A,C))) == 0  % sprawdzenie obserwowalności
    sprintf("System obserwowalny")
else
    sprintf("System nie obserwowalny")
    return
end

time_step = 0.01;

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
syms beta [1 system_size] real positive
for i=1:system_size
    assume(beta(i)>0);
end

% zakomentowanie daje pełną analizę dla różnych beta (dłuzej się liczy)
for i = 1:system_size
    beta(i) = 0;
end

syms J
P1 = [];
P2 = [];
for i = 1:system_size
    W = [A,    beta(i)*B*B.'
         C.'*C, -A.']
    Phi = simplify(expm(W * t),'steps',10)
%     Phi = expm(W * t)
    Phi11 = Phi(1:system_size,1:system_size)
    Phi21 = Phi(system_size+1:2*system_size,1:system_size)
    int_internal = simplify(expm(-A.'*(T-tau))*C.'*C*subs(Phi11,t,tau), 'steps', 10)
%     int_internal = (-A.'*(T-tau))*C.'*C*subs(Phi11,t,tau)
    Mi = int(int_internal, tau, 0, T)
    Mi = simplify(Mi, 'steps', 10)
    ei = zeros([system_size, 1])
    ei(i) = 1
    iMi = inv(Mi);
    iMi = simplify(iMi, 'steps', 10)
    p1 = Phi11 * iMi * ei
    p1 = simplify(p1, 'steps', 10)
    p2 = Phi21 * iMi * ei
    p2 = simplify(p2, 'steps', 10)
    P1 = [P1, p1]  % to jest tak na prawdę transpozycja P1
    P2 = [P2, p2]
end

G1 = P1' * C.'  % transponuję zmienną P1, ale w teorii jest to macierz P1 nietransponowana
G2 = P2' * B

G1 = simplify(G1,'steps',10)
G2 = simplify(G2,'steps',10)

% przykladowe przebiegi składowych macierzy G1 i G2
%%
Treal = 2;
betasreal = ones(1,system_size);

figure
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
        title("G1["+num2str(i)+','+num2str(j)+']')
    end
end

figure
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
        title("G2["+num2str(i)+','+num2str(j)+']')
    end
end

%% numerycznie obliczony wskaźnik jakości J
Tslen = 10;  % liczba analizowanych długości okien
bslen = 20;  % liczba analizowanych wartości współczynników beta
Ts    = logspace(-1,log10(10),Tslen);  % analizowane długości okna
betas = logspace(-2,log10(1),bslen);  % analizowane współczynniki beta DLA MACIERZY G1 i G2
J = sum(G1(1,:).^2) + 1*sum(G2(1,:).^2)
for i = 2:system_size
    J = J + sum(G1(i,:).^2) + 1*sum(G2(i,:).^2)  % koszt liczymy jakby na wejściu były zakłócenia (tak samo ważne jak na wyjściu) + dla każdego pojedynczego wejścia/wyjścia współczynniki są takie same
end

try
    double(beta);
    syms betafnc [1 system_size] real positive
    Jsymfun(T,betafnc,t) = J;
catch
    Jsymfun(T,beta,t) = J;
end

Jfun = matlabFunction(Jsymfun);
Jvalues = zeros(length(Ts), length(betas));
i = 0;
for Ttest = Ts  % długość okna
    i = i + 1;
    j = 0;
    for betatest = betas  % bety, które przyjmujemy do obliczenia G1 i G2 (przyjmujemy jednakowe wartości dla każdej z nich)
        j = j + 1;
        % przekazywanie różnej liczby argumentów do funkcji kosztów (w zależności od wymiaru systemu)
        args = {Ttest};
        for k = 2:(1+system_size)
            args(k) = {betatest};
        end
        Jnumerical = integral(@(tau)Jfun(args{:}, tau), 0, Ttest);
        Jvalues(i, j) = Jnumerical;
    end
end

%%

figure
[X, Y] = meshgrid(Ts, betas);
mesh(X, Y, (real(Jvalues))')
grid on;
xlabel("T")
ylabel("$\beta(1)=\beta(2)=...=Y$", 'Interpreter', 'latex')
zlabel("Wartość funkcji kosztu")
zlim([0, 5])
caxis([0 3])
colorbar