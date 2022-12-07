function [phi11] = getPhi11(t, i, W, n)
    % t - próbka czasu
    % i - numer próbki
    % W - macierz W
    % n - rozmiar obserwowanego systemu
    phis = expm(W(:,:,i)*t);
    phi11 = phis(1:n, 1:n);
end

