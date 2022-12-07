function [J] = calculateJ(G1, G2, betashat, time_samples)
    % obliczenie wskaźnika jakości J dla danych macierzy G1, G2 i wag beta
    % G1 - macierz 3 wymiarowa obserwatora wyjścia
    % G1 - macierz 3 wymiarowa obserwatora wejścia
    % betashat - wagi wskaźnika jakości
    % time_samples - próbki czasu dla których wyliczone są macierze G1 i G2

    s = size(G1);
    n = s(1);
    n_C = s(2);
    s = size(G2);
    n_B = s(2);

    J = 0;
    for i = 1:n
        for j = 1:n_C
            J = J + trapz(time_samples, (squeeze(G1(i,j,:))).^2);
        end
        for k = 1:n_B
            J = J + betashat(k) * trapz(time_samples, (squeeze(G2(i,k,:))).^2);
        end
    end
    
end