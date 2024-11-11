function ans = is_in_set(A,x,tol)

% auxiliary function that checks (using tolerance tol) if x is already a
% member of set A

if nargin<3
    tol = 1e-6;
end

m = size(A,1);
if m==0
    ans = 0;
else
    y = zeros(m,1);
    for k = 1:m
        tmp = A(k,:)-x;
        y(k,1) = tmp*tmp'; % norm(A(k,:)-x);
    end
    ans = min(sqrt(y))<=tol;
end
