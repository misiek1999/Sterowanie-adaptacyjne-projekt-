function [G1, G2, time_samples, betas, T, x, fval, exitflag, output] = getoptiG1G2(A, B, C, time_step, startpoint, ploton, searchiter)
    % funkcja zwracająca optymalne macierze G1 G2 w sensie zdefiniowanego
    % przez nas wskaźnika
    % A - macierz stanu A
    % B - macierz stanu B
    % C - macierz stanu C
    % time_step  - okres próbkowania, co tyle wyliczane są wartości
    %              macierzy G1 G2
    % start_point - punkt startowy algorytmu optymalizacji w formie
    %               [T, beta1, beta2, ..., beta_n_B]
    % ploton - czy pokazywać przebieg optymalizacji i rysować macierze G1 G2
    % searchiter - maksymalna liczba iteracji algorytmu optymalizacji

    if ploton
        options = optimset('MaxIter',searchiter,'PlotFcns',@optimplotfval);
    else
        options = optimset('MaxIter',searchiter);
    end

    params0 = startpoint;  % T0, beta1, beta2, ..., beta_n_B - długość zależy od wymiaru macierzy B
    optifhand = @(p)optiG1G2(A, B, C, p);

    [x,fval,exitflag,output] = fminsearch(optifhand, params0, options);

    T = x(1);
    betas = x(2:end);

    time_samples = 0:time_step:T;
    
    [G1, G2] = getG1G2(A, B, C, betas, time_samples);
    
    if ploton
        plotG1G2(G1, G2, time_samples);
    end
    
end