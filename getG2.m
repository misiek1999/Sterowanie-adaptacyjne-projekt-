function [G1] = getG2(t)
    % function to calcluate new G1 matrix
    global A B C n W M 
    
    P2 = zeros(n);
    for i = 1:n
        ei = zeros([n 1]);
        ei(i) = 1;

        p = getPhi21(t, i) /M{i}*ei;
        
        P2(:,i) = p;
    end

    G1 = P2.'*B;
end