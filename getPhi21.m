function [phi21] = getPhi21(t, i)
    global W n
    phis = expm(W{i}*t);
    phi21 = phis(n+1:end, 1:n);
end

