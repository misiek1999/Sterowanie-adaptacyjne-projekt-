function [G1] = getG1(t)
    % function to calcluate new G1 matrix
    global A B C n W M 
    
    P1 = zeros(n);
    for i = 1:n
        ei = zeros([n 1]);
        ei(i) = 1;

        p = getPhi11(t, i) /M{i}*ei;
        
        P1(:,i) = p;
    end

    G1 = P1.'*C.';
end