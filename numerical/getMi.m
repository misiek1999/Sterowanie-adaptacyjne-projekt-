function Mi = getMi(A, C, T, i, W, n)
    % funkcja pomocnicza do liczenia macierzy Gramma M
    % A - macierz przestrzeni stanów A
    % C - macierz przestrzeni stanów C
    % T - długość okna
    % i -numer iteracji
    % W - macierz W
    % n - liczba stanów
    Mi_integral_fun = @(tau) expm(-A.'*(T-tau))*C.'*C*getPhi11(tau, i, W, n);
	Mi = integral(@(tau)Mi_integral_fun(tau), 0, T, 'ArrayValued', true);
end