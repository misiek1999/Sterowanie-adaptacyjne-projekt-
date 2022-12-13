function plotG1G2(G1_vec, G2_vec, time_samples)
    % funkcja rysująca przebiegi macierzy G1 G2
    % G1 - macierz G1
    % G1 - macierz G2
    % time_samples - próbki, dla których zostały wyliczone macierze G1 G2

    s = size(G1_vec);
    n = s(1);
    n_C = s(2);
    s = size(G2_vec);
    n_B = s(2);

    % G1
    figure
    for i = 1:n
        for j = 1:n_C
            subplot(n, n_C, j + (i-1)*n_C)
            plot(time_samples, squeeze(G1_vec(i,j,:)))
            grid on;
            title("G1["+num2str(i)+','+num2str(j)+']')
        end
    end

    % G2
    figure
    for i = 1:n
        for j = 1:n_B
            subplot(n, n_B, j + (i-1)*n_B)
            plot(time_samples, squeeze(G2_vec(i,j,:)))
            grid on;
            title("G2["+num2str(i)+','+num2str(j)+']')
        end
    end
    
end