global n time_stamps G1_vec G2_vec T 
plotGs = 0;  % nie rysuj macierzy G1, G2

path_to_model = 'model/load_model.m';
run(path_to_model);

n = length(A);
n_C = size(C);
n_C = n_C(1);
n_B = size(B);
n_B = n_B(2);

Js = [];
J2s = [];
Ts = 0.1:0.5:10;
betasrange = [0, 5.^(linspace(-3,0,9))];
for T = Ts
    T
    for beta = betasrange 
        betas = beta * ones(1,n);
        J  = 0;  % gdy uwzględniamy zakłócenia na wejściu
        J2 = 0;  % gdy nie uwzględniamy zakłóceń na wejściu
        initialization;
        for i = 1:n
            for j = 1:n_C
                J = J + trapz(time_stamps, (squeeze(G1_vec(i,j,:))).^2);
                J2 = J2 + trapz(time_stamps, (squeeze(G1_vec(i,j,:))).^2);
            end
            for k = 1:n_B
                J = J + trapz(time_stamps, (squeeze(G2_vec(i,k,:))).^2);
            end
        end
        J2s = [J2s J2];
        Js = [Js J];
    end
end

%%
Jvalues = reshape(Js, length(betasrange), length(Ts));
J2values = reshape(J2s, length(betasrange), length(Ts));
figure
[X, Y] = meshgrid(Ts, betasrange);
mesh(X, Y, Jvalues, 'FaceColor', 'r', 'EdgeColor', 'k')
hold on;
mesh(X, Y, J2values, 'FaceColor', 'b', 'EdgeColor', 'k')
grid on;
xlabel("T")
ylabel("$\beta(1)=\beta(2)=...=Y$", 'Interpreter', 'latex')
zlabel("Wartość funkcji kosztu")
legend("$\hat\beta(1)=\hat\beta(2)=...=1$", "$\hat\beta(1)=\hat\beta(2)=...=0$", 'Interpreter', 'latex')
zlim([0, 5])
% caxis([0 3])
% colorbar