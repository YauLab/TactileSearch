function Y = shuffleVector(X,N)
% Shuffle the vector X for N times, and concatenate the each shuffle to
% generate a new vector Y.

[n,m] = size(X);
if n > 1 && m > 1
    error('X is not a Vector.')
end

if nargin < 2
    N = 1;
end

Y = zeros(1,N*length(X));
for k = 1:N
    [~,index] = sort(rand(size(X)));
    Y((k-1)*length(X)+1:k*length(X)) = X(index); 
end