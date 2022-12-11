function [G1, G2] = getG1G2(A, B, C, betas, time_samples)
    % funkcja wyliczająca numerycznie macierze G1 G2 dla danego systemu
    % A - macierz stanu A
    % B - macierz stanu B
    % C - macierz stanu C
    % betas - współczynniki wagowe przy projektowaniu macierzy G1 i G2
    % time_samples - próbki czasu dla których wyliczane są macierze G1 i G2
    % G1 - macierz 3 wymiarowa obserwatora wyjścia
    % G1 - macierz 3 wymiarowa obserwatora wejścia
    
    n = length(A);
    if (n-rank(obsv(A,C))) == 0  % sprawdzenie obserwowalności
        % sprintf("System obserwowalny")
    else
        % sprintf("System nie obserwowalny")
        return
    end
    
    % A, B, C i wektora czasu time_samples (zawarta w nim jest długość okna)
    n_C = size(C);
    n_C = n_C(1);
    n_B = size(B);
    n_B = n_B(2);
    T = time_samples(end);
    samples_no = length(time_samples);
    
    M = zeros(n,n,n);
    W = zeros(2*n,2*n,n);
    % wyliczenie wszystkich macierzy Mi i zapisanie ich w M
    for i = 1:n
        % calculate W matrix
        Wi = [A,     betas(i)*B*B.'
             C.'*C, -A.'    ];
        W(:,:,i) = Wi;
        % find Grama matrix's
        Mi_integral_fun = @(tau) expm(-A.'*(T-tau))*C.'*C*getPhi11(tau, i, W, n);
        Mi = integral(@(tau)Mi_integral_fun(tau), 0, T, 'ArrayValued', true);
        M(:,:,i) = Mi;
    end
    
    % utworzenie macierzy 3 wymiarowych G1, G2 (trzeci wymiar odpowiada czasowi)
    G1_vec = zeros(n,n_C,samples_no);
    G2_vec = zeros(n,n_B,samples_no);
    itr = 1;
    for t = time_samples
        P1 = [];
        P2 = [];
        for i = 1:n
            ei = zeros([n 1]);
            ei(i) = 1;
            
            phis = expm(W(:,:,i)*t);
            phi11 = phis(1:n, 1:n);
            phi21 = phis(n+1:end, 1:n);

            p1 = phi11 / M(:,:,i) * ei;
            p2 = phi21 / M(:,:,i) * ei;
        
            P1 = [P1 p1];
            P2 = [P2 p2];
        end

        G1_vec(:,:, itr) = P1.'*C.';
        G2_vec(:,:, itr) = P2.'*B; 
        itr = itr + 1;
    end
    
    G1 = G1_vec;
    G2 = G2_vec;
    
end