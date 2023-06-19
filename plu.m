% A simple animation that factorises a matrix by the PLU decomposition

n = 4
A = rand(n);
A0 = A;

p = 1:n;
s = max(abs(A'))';

for k = 1:n-1
	[~, idx] = max(abs(A(p(k:end), k)) ./ s(p(k:end)));
	j = idx + k - 1;
	p([j, k]) = p([k, j]);
	for i = k + 1:n
		z = A(p(i), k) / A(p(k), k);
		A(p(i), k) = z;
		for j = k+1:n
			A(p(i), j) = A(p(i), j) - z*A(p(k), j);
		end
	end
end
L = tril(A, -1) + eye(n);
U = triu(A);
L*U
A0


