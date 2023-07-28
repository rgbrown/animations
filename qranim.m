%QRANIM: A simple visualisation of the QR factorisation

n = 15;
A = randn(n);
Q = eye(n);

for j = 1:n
    for i = j+1:n
        % Givens rotation 
        G = givens(A(j, j), A(i, j));
        A([j, i], :) = G*A([j, i], :);
        Q([j, i], :) = G'*Q([j, i], :);
        A(abs(A) < 1e-15) = 0;
        clc()
        disp('Q = ')
        disp(Q)
        disp('R = ')
        disp(A)
        pause(0.15)
    end
end
