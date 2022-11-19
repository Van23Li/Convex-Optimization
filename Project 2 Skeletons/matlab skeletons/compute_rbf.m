function [K] = compute_rbf(X1,X2)
%COMPUTE_RBF compute rbf matrix between X1 and X2
K = [];
for idx = 1:size(X2,1)
    K = [K, exp(-sum((X1 - X2(idx,:)).^2, 2)/(2*1.5))];
end
end

