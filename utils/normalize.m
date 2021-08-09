function x = normalize(X)
% will normalize any matrix by its abs max value
x = X./max(abs(X(:)));
