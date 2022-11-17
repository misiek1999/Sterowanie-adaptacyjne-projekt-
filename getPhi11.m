function [phi11] = getPhi11(t, i)
    global W n
    phis = expm(W{i}*t);
    phi11 = phis(1:n, 1:n);
end

