function Y = OuterProduct(X)
% OUTERPRODUCT row-by-row outer product of matrix X.

N = size(X,1); p = size(X,2);
Y = zeros(N, p, p);
for ii=1:N
  Y(ii,:,:) = X(ii,:)' * X(ii,:);
end

end

