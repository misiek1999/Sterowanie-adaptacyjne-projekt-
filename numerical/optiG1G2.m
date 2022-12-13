function y = optiG1G2(A, B, C, params)
    % funkcja kosztu do optymalizacji parametrów macierzy G1 i G2
    % A - macierz przestrzeni stanów A
    % B - macierz przestrzeni stanów B
    % C - macierz przestrzeni stanów C
    % params - parametry optymalizacji w postaci [T, beta1, beta2, ..., beta_n_B]
    
    
    Tend = params(1);
    betas = params(2:end);
    
    n = length(A);
    n_C = size(C);
    n_C = n_C(1);
    n_B = size(B);
    n_B = n_B(2);
    
    betashat  = ones(1,n);
    betashat2 = zeros(1,n);
    
    time_samples = 0:0.01:Tend;
    
    Tcost = 1 + (n^2).^(Tend);
    
    [G1, G2] = getG1G2(A, B, C, betas, time_samples);
    J = calculateJ(G1, G2, betashat, time_samples);  % gdy uwzględniamy zakłócenia na wejściu
    J2 = calculateJ(G1, G2, betashat2, time_samples);  % gdy nie uwzględniamy zakłóceń na wejściu
    
    y = J * J2 * Tcost;
    
end