function ans = is_in_set(A,x,tol)

% auxiliary function that checks (using tolerance tol) if x is already a
% member of set A

class_t = superiorfloat(A,x);

if nargin<3
    tol = 1e-6;
end

m = size(A,1);
if m==0
    ans = 0;
else
    if strcmp(class_t,'mp')
        ans = min(vecnorm(double(A-x),2,2)) <= tol*(1+norm(x));
    else
        ans = min(vecnorm(A-x,2,2)) <= tol*(1+norm(x));
    end
end
