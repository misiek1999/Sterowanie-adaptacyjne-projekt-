example_params;

n = length(A);
n_C = size(C);
n_C = n_C(1);
n_B = size(B);
n_B = n_B(2);

betashat  = [1 1 1];
betashat2 = [0 0 0];

Js = [];
J2s = [];
Ts = 0.1:0.5:10;
betasrange = [0, 5.^(linspace(-3,0,9))];
for T = Ts
    T
    time_samples = 0:0.01:T;
    for beta = betasrange 
        betas = beta * ones(1,n);  
        [G1, G2] = getG1G2(A, B, C, betas, time_samples);
        J = calculateJ(G1, G2, betashat, time_samples);  % gdy uwzględniamy zakłócenia na wejściu
        J2 = calculateJ(G1, G2, betashat2, time_samples);  % gdy nie uwzględniamy zakłóceń na wejściu
        Js = [Js J];
        J2s = [J2s J2];
    end
end

%% 
% wyrysowanie wyników
Jvalues = reshape(Js, length(betasrange), length(Ts));
figure
[X, Y] = meshgrid(Ts, betasrange);
mesh(X, Y, Jvalues, 'FaceColor', 'r', 'EdgeColor', 'k')
grid on;
xlabel("T")
ylabel("$\beta(1)=\beta(2)=...=Y$", 'Interpreter', 'latex')
zlabel("Wartość funkcji kosztu")
% legend("$\hat\beta(1)=\hat\beta(2)=...=1$", "$\hat\beta(1)=\hat\beta(2)=...=0$", 'Interpreter', 'latex')
zlim([0, 5])
% caxis([0 3])
% colorbar

%% 
% jeśli dodatkowo bierzemy pod uwagę koszt obliczeń 
Tcost = 1 + Ts;
Tcost = bsxfun(@times,Jvalues,reshape(Tcost, 1, []));
figure
[X, Y] = meshgrid(Ts, betasrange);
mesh(X, Y, Tcost, 'FaceColor', 'r', 'EdgeColor', 'k')
grid on;
xlabel("T")
ylabel("$\beta(1)=\beta(2)=...=Y$", 'Interpreter', 'latex')
zlabel("Wartość funkcji kosztu")
% legend("$\hat\beta(1)=\hat\beta(2)=...=1$", "$\hat\beta(1)=\hat\beta(2)=...=0$", 'Interpreter', 'latex')
zlim([1, 10])

%%
% a tu jeszcze inny przypadek, bierzemy funkcje kosztu dla beta=1 i beta=0
% oraz bierzemy pod uwagę czas oblieczeń, w takim wypadku zawsze powinniśmy
% dostać jakieś minimum lokalne
J2values = reshape(J2s, length(betasrange), length(Ts));
Tcost2 = 1 + (n^2).^(Ts);
Tcost2 = bsxfun(@times,Jvalues.*J2values,reshape(Tcost2, 1, []));
mincost = min(min(Tcost2));

figure
[X, Y] = meshgrid(Ts, betasrange);
bestT = (X(Tcost2==mincost))
bestbeta = (Y(Tcost2==mincost))
mesh(X, Y, Tcost2, 'FaceColor', 'r', 'EdgeColor', 'k')
grid on;
xlabel("T")
ylabel("$\beta(1)=\beta(2)=...=Y$", 'Interpreter', 'latex')
zlabel("Wartość funkcji kosztu")
% legend("$\hat\beta(1)=\hat\beta(2)=...=1$", "$\hat\beta(1)=\hat\beta(2)=...=0$", 'Interpreter', 'latex')
zlim([mincost, (mincost+10)*1.5])