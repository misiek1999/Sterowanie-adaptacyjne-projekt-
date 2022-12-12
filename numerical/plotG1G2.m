function plotG1G2(G1_vec, G2_vec, time_samples, A, B, C)

    n = length(A);
    n_C = size(C);
    n_C = n_C(1);
    n_B = size(B);
    n_B = n_B(2);

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