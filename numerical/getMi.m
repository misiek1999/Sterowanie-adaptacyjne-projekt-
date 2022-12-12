function Mi = getMi(A, C, T, i, W, n)
    Mi_integral_fun = @(tau) expm(-A.'*(T-tau))*C.'*C*getPhi11(tau, i, W, n);
	Mi = integral(@(tau)Mi_integral_fun(tau), 0, T, 'ArrayValued', true);
end