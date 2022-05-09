function [X,X_mean] = preprocess(X)

% Center the data
X_mean = mean(X,2);
X = (X - X_mean);

% Whiten the data
covariance = cov(X');
[V, D] = eig(covariance);
d_inv = sqrt(1./diag(D));
X = V*diag(d_inv,0)*V'*X;

end