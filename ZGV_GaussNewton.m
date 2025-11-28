function [lambda,mu,x,y,flag,err,step,init] = ZGV_GaussNewton(A,B,C,lambda,mu,x,y,opts)

% [lambda,mu,x,y,flag,err] = ZGV_GaussNewton(A,B,C,lambda,mu,x,y,opts) 
% returns ZGV (2D) point for a twoparameter pencil (A + lambda*B + mu*C)x = 0 
% using the Gauss-Newton method.
% 
% Output: 
%    - a point (lambda,mu) and unit vectors x and y such that:
%      a) (A + lambda*B + mu*C)*x = 0  
%      b) y'*(A + lambda*B + mu*C) = 0
%      c) y'*B*x = 0
%    - flag: convergence (1) or not (0)
%    - err: norms of the residual
%    - step: number of steps
%    - init: 0 - initial vectors were given, 1 - last singular vectors were used, 2 -
%      a combination of two last singular vectors was used
%
% Input: 
%    - A,B,C: square matrices
%    - lambda, mu: initial approximation for the ZGV point
%    - x, y: initial approximation for the right and left eigenvector, if
%      not provided or empty then singular vectors corresponding to the 
%      smallest singular values of (A+lambda*B+mu*C) are used
%   - opts : options
%
% Options in opts:
%   - maxsteps: (15) maximal number of steps 
%   - tol: (1e2*eps) tolerance (relative to the norm of matrices)
%   - maxdif: (1e-1), maximal relative change of (lambda,mu) in one iteration
%   - show: (0), set to 1 to display values and residuals

% This is Algorithm 4 in B. Plestenjak: On properties and numerical 
% computation of critical points of eigencurves of bivariate matrix pencils

% Bor Plestenjak 2024, revised in 2025

narginchk(5, 8);
if nargin < 6, x = [];    end
if nargin < 7, y = [];    end
if nargin < 8, opts = []; end

if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A,B,C);
end

if isfield(opts,'showiter'), showiter = opts.showiter; else, showiter = 0;                       end
if isfield(opts,'tol'),      tol      = opts.tol;      else, tol      = 1e2*eps(class_t);        end
if isfield(opts,'maxsteps'), maxsteps = opts.maxsteps; else, maxsteps = 15;                      end
if isfield(opts,'maxdif'),   maxdif   = opts.maxdif;   else, maxdif   = 1e-1;                    end

init = 0;

if isempty(x) || isempty(y)
    if (size(A,1)>500) && (~strcmp(class_t,'mp'))
        [U,S,V] = svds(A + lambda*B + mu*C,5,'smallest');
    else
        [U,S,V] = svd(full(A + lambda*B + mu*C));
    end
    s = diag(S);
    if length(s)<3
        s = [0;s];
    end 
    % heuristic criteria when to use a combination of two last singular vectors
    % or just the last singular vector
    if (s(end-1)<s(end-2)*1e-2 || (s(end-1)<s(end)*10) && (s(end-1)<s(end)*1e4))
        if showiter
            disp('Taking a combination of the two last singular vectors')
        end
        init = 2;
        if isempty(x), x = V(:,end-1:end)*randn(2,1); x = x/norm(x); end
        if isempty(y) 
            % we select y so that y'*B*x = 0
            y = (-U(:,end-1)'*B*x)*U(:,end) + (U(:,end)'*B*x)*U(:,end-1); y = y/norm(y); 
        end
    else
        if showiter
            disp('Taking the last singular vector')
        end
        init = 1;
        if isempty(x), x = V(:,end); end
        if isempty(y), y = U(:,end); end
    end
end

% we use relative tolerance
reltol = tol * norm(A + lambda*B +mu*C,'fro');
n = size(A,1);
step = 0;
flag = 0;
err = [];
% random vectors used for the normalization
a = randn(n,1,class_t) + 1i*randn(n,1,class_t); a = a/norm(a);
b = randn(n,1,class_t) + 1i*randn(n,1,class_t); b = b/norm(b);
% we normalize initial vectors
x = x/(a'*x);
y = y/(b'*y);
bc = conj(b);
dif_eig = 0;
while (step < maxsteps) 
    step = step + 1;
    y = conj(y); 
    M = A + lambda*B + mu*C;
    J = [M                  zeros(n,class_t)   B*x   C*x;
         zeros(n,class_t)   M.'                B.'*y C.'*y;
         y.'*B              x.'*B.'            0     0;
         a'                 zeros(1,n,class_t) 0     0;
         zeros(1,n,class_t) bc'                0     0];
    res = [M*x; M.'*y; y.'*B*x; (a'*x-1); (bc'*y-1)];
    err(step) = norm(res);
    if showiter
        % we compute the smallest singular value of the Jacobian
        if (size(A,1)>500) && (~strcmp(class_t,'mp'))
            s = svds(J,5,'smallest');
        else
            s = svd(J);
        end
        zgv = y.'*C*x;
        fprintf('step %3d | err: %8.3e | lambda (%14.7e,%14.7e) | mu (%14.7e,%14.7e) |dif_eig: %8.3e |s_min: %8.3e |zgv: %8.3e\n',step,err(step),real(lambda),imag(lambda),real(mu),imag(mu),dif_eig,s(end),zgv);
    end
    if dif_eig>maxdif
        % if the relative change of (lambda,mu) is too large, this
        % indicates accidental convergence  
        flag = 0;
        y = conj(y);
        break
    end
    if err(step)<reltol 
        % convergence
        flag = 1;
        y = conj(y);
        break
    end
    % turn off the warning when J is rank deficient
    warning off
    upd = -J\res;
    warning on
    x = x + upd(1:n);
    y = y + upd(n+1:2*n);
    dif_eig = norm(upd(2*n+1:2*n+2))/(1+norm([lambda mu]));
    lambda = lambda + upd(2*n+1);
    mu = mu + upd(2*n+2);
    y = conj(y);
end    

