function [lambda, mu] = critical_points_project(A,B,C,opts)

% [lambda, mu] = critical_points_project(A,B,C,opts) finds critical points (lambda,mu) 
% for the eigenvalue problem (A + lambda*B + mu*C)x = 0, such that 
% a) 2D points: there exist nonzero vectors x and y such that 
%    (A + lambda*B + mu*C)*x=0, y'*(A + lambda*B + mu*C)=0, y'*B*x=0
% b) ZGV points: 2D points where mu is a simple eigenvalue 
%
% Keep matrices of the problem small, i.e., not larger than 40 x 40.
%
% Input:
%   - A,B,C : square matrices n x n
%   - opts : options
%
% Options in opts:
%   - goal ('2D'): set to 'ZGV' to compute only ZGV points
%   - cplx_rnd (1): set to 0 to use random projections (faster for real matrices but less reliable)
%   - gensize (0): set to 1 to always extract the generic number n(n-1) of finite eigenvalues
%   - sep (1e4*sqrt(eps)): threshold for the separation of regular eigenvalues
%   - sepinf (1e4*eps): threshold for the separation of finite eigenvalues
%   - membtol (1e-6): treshold for preventing to output multiple instances of critical points
%   - YCXtol (1e-6): treshold for the condition |y'*C*x| > YCXtol*norm(C) for the ZGV point
%   - GStol  (1e-6): treshold for the condition sigma(n-1) > GStol*sigma(1) for geometric simplicity
%   - GN_refine (0): set to 1 to apply Gauss-Newton refinement of computed points
%   - show (0): display values used for the identificaton (1) or no (0)
%   - mprefine (0): refine final results in higher precision (requires Advanpix MCT)
%   - other options for twopareig
%
% Output:
%   - lambda, mu: vectors with critical points

% This is Algorithm 3 in B. Plestenjak: On properties and numerical 
% computation of critical points of eigencurves of bivariate matrix pencils

% Bor Plestenjak 2024, revised in 2025

if nargin<4, opts=[]; end

if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A,B,C);
end

if isfield(opts,'goal'),      goal = opts.goal;           else, goal = '2D';       end
if isfield(opts,'cplx_rnd'),  cplx_rnd = opts.cplx_rnd;   else, cplx_rnd = 1;      end
if isfield(opts,'gensize'),   gensize = opts.gensize;     else, gensize = 0;       end
if isfield(opts,'show'),      show = opts.show;           else, show = 0;          end
if isfield(opts,'sep'),       sep      = opts.sep;        else, sep      = sqrt(1e4*eps(class_t));      end
if isfield(opts,'sepinf'),    sepinf   = opts.sepinf;     else, sepinf   = 1e4*eps(class_t);        end
if isfield(opts,'GN_refine'), GN_refine  = opts.GN_refine; else, GN_refine   = 0;      end
if isfield(opts,'membtol'),   membtol = opts.membtol;     else, membtol = 1e-6;    end
if isfield(opts,'YBXtol'),    YBXtol = opts.YBXtol;       else, YBXtol = 1e-6;  end
if isfield(opts,'YCXtol'),    YCXtol = opts.YCXtol;       else, YCXtol = 1e-6;     end
if isfield(opts,'GStol'),     GStol = opts.GStol;         else, GStol = 1e-6;      end
if isfield(opts,'multtol'),   multtol = opts.multtol;     else, multtol = 1e-4; end
if isfield(opts,'use_eigs'),  use_eigs = opts.use_eigs;   else, use_eigs = 0;      end
if isfield(opts,'target'),    target = opts.target;       else, target = 0;        end
if isfield(opts,'neigs'),     neigs = opts.neigs;         else, neigs = 10;        end
if isfield(opts,'opts_eigs'), opts_eigs = opts.opts_eigs; else, opts_eigs = [];    end
if isfield(opts,'mprefine'),  mprefine = opts.mprefine;   else, mprefine = 0;      end
if isfield(opts,'opts_mp'),   opts_mp = opts.opts_mp;     else, opts_mp = [];        end

% Make sure all inputs are of the same numeric type.
if ~isa(A,class_t), A = numeric_t(A,class_t); end
if ~isa(B,class_t), B = numeric_t(B,class_t); end
if ~isa(C,class_t), C = numeric_t(C,class_t); end

lambda = numeric_t([],class_t); 
mu = numeric_t([],class_t); 
n = size(A,2);

if cplx_rnd
    [U,~] = qr(randn(2*n,class_t)+1i*randn(2*n,class_t),0); % random complex orthonormal columns
    [V,~] = qr(randn(2*n,class_t)+1i*randn(2*n,class_t),0); % random complex orthonormal columns
else
    [U,~] = qr(randn(2*n,class_t),0); % random ncols orthonormal columns
    [V,~] = qr(randn(2*n,class_t),0); % random ncols orthonormal columns
end

Z = zeros(n,class_t);
A2 = [A Z; B A];
B2 = [B Z; Z B];
C2 = [C Z; Z C];

AA = U'*A2*V;
BB = U'*B2*V;
CC = U'*C2*V;

AP = AA(1:end-1,1:end-1);
BP = BB(1:end-1,1:end-1);
CP = CC(1:end-1,1:end-1);

% we solve the projected problem, but since Delta0 can be singular and
% Delta1 is nonsingular, we use substitution of variables

if use_eigs
    [lambda,mu,x1,x2] = twopareigs_ks(A,-B,-C,AP,-BP,-CP,neigs,opts);
    m = length(lambda);
    for j = 1:m
        y1(:,j) = min_sing_vec((A+lambda(j)*B+mu(j)*C)');
        y2(:,j) = min_sing_vec((AP+lambda(j)*BP+mu(j)*CP)');
    end
else
    [tmplambda,tmpmu,x1,x2,y1,y2] = twopareig(B,-A,-C,BP,-AP,-CP,opts);
    lambda = 1./tmplambda;
    mu = tmpmu./tmplambda;
end
solution = numeric_t([],class_t); 

m = length(lambda);

nrmAA = norm(AA,'fro');
nrmBB = norm(BB,'fro');
nrmCC = norm(CC,'fro');
a = zeros(m,1);
b = zeros(m,1);
c = zeros(m,1);
condB = zeros(m,1);
condC = zeros(m,1);
tip = zeros(m,1);
for j = 1:m
    tmpABCnorm = nrmAA + abs(lambda(j))*nrmBB + abs(mu(j))*nrmCC;
    ABCnorm = norm([tmpABCnorm 0; nrmBB tmpABCnorm],'fro');	
    a(j) = (AA(end,1:end-1) + lambda(j)*BB(end,1:end-1) + mu(j)*CC(end,1:end-1))*x2(:,j);
    b(j) = y2(:,j)'*(AA(1:end-1,end) + lambda(j)*BB(1:end-1,end) + mu(j)*CC(1:end-1,end));
    a(j) = abs(a(j))/ABCnorm;
    b(j) = abs(b(j))/ABCnorm;
    c(j) = (y1(:,j)'*B*x1(:,j))*(y2(:,j)'*CP*x2(:,j)) - (y1(:,j)'*C*x1(:,j))*(y2(:,j)'*BP*x2(:,j));
    c(j) = abs(c(j))/sqrt(1+abs(lambda(j))^2);
end
if gensize
    mx = sort(max(a,b));
    sep = sqrt(mx(n*n)*mx(n*n+1));
end
normB = norm(B);
normC = norm(C);
for j = 1:m
    condB(j) = abs(y1(:,j)'*B*x1(:,j))/(norm(y1(:,j))*norm(x1(:,j))*normB);
    condC(j) = abs(y1(:,j)'*C*x1(:,j))/(norm(y1(:,j))*norm(x1(:,j))*normC);
    if (a(j)<sep) && (b(j)<sep)
        if c(j)>sepinf
            tip(j) = 1;
        else
            tip(j) = 2;
        end
    else
        tip(j) = 3;
    end
    if (tip(j) == 1) 
        if GN_refine
            % optional refinement with Gauss-Newton method
            [ll,mm,x,y,flag,err,step] = ZGV_GaussNewton(A,B,C,lambda(j),mu(j),x1(:,j),y1(:,j),opts);
            if err(end)<err(1)
                if flag && mprefine
                    [ll,mm,x,y,flagmp,errmp,stepmp,initmp] = ZGV_GaussNewton(mp(A),mp(B),mp(C),mp(ll),mp(mm),mp(x),mp(y),opts_mp);
                    ll = double(ll);
                    mm = double(mm);
                    x = double(x);
                    y = double(y);
                end
                lambda(j) = ll;
                mu(j) = mm;
                y1(:,j) = y;
                x1(:,j) = x;
                condB(j) = abs(y1(:,j)'*B*x1(:,j))/(norm(y1(:,j))*norm(x1(:,j))*normB);
                condC(j) = abs(y1(:,j)'*C*x1(:,j))/(norm(y1(:,j))*norm(x1(:,j))*normC);
            end
        end
        if ~is_in_set(solution,[lambda(j) mu(j)],membtol)
            if strcmp(goal,'ZGV')
                % if we are computing ZGV points, we first check if |y'*B*x|=0 and |y'*C*x| is nonzero
                if (condB(j)<YBXtol) && (condC(j)>YCXtol)
                    % if this is true, we additionally check if mu is geometrically simple
                    s = svd(A + lambda(j)*B + mu(j)*C);
                    if s(n-1) > GStol*s(1)
                        tip(j) = 0;
                        solution = [solution; [lambda(j) mu(j)]];
                    end
                end
            else
                % if we are computing 2D points, we check that either |y'*B*x|=0 or 
                % geometric multiplicity is at least two
                if (condB(j) < YBXtol)  
                    tip(j) = 0;
                    solution = [solution; [lambda(j) mu(j)]];
                else
                    % we estimate multiplicity from the number of close eigenvalues
                    dist = vecnorm([lambda mu]-[lambda(j) mu(j)],2,2);
                    pos = dist < multtol*(1+norm([lambda(j) mu(j)]));
                    alg_mult = sum(pos);
                    if alg_mult>1
                        s = svd(A + lambda(j)*B + mu(j)*C);
                        if s(n-1) <= GStol*s(1)
                            tip(j) = -1;
                            solution = [solution; [lambda(j) mu(j)]];
                        end
                    end
                end
            end
        end
    end
end

if show
    mab = max(a,b);
    [tmp, ord] = sort(tip);
    sub1 = find(tip(ord)<3);
    [tmp, ordsub1] = sort(-c(ord(sub1)));
    sub2 = find(tip(ord)==3);
    [tmp, ordsub2] = sort(mab(ord(sub2)));
    tmp = ord(sub1);
    ord(sub1) = ord(sub1(ordsub1));
    tmp = ord(sub2);
    ord(sub2) = ord(sub2(ordsub2));

    [tmp, ord] = sort(real(lambda));
    if show == 2
        disp(' ')
        disp(' i   | p(i)|     re(lambda)     im(lambda) |     re(mu)          im(mu)    |     s(eig) |   ||V''*x|| |   ||U''*y|| |   |y''*B*x| |   |y''*C*x| | tip')
        disp('-------------------------------------------------------------------------------------------------------------------------------------------------------')
        for j = 1:m
            fprintf('%4d |%4d | %14.6e %14.6e | %14.6e %14.6e | %10.2e | %10.2e | %10.2e | %10.2e | %10.2e | %1d \n',j, ord(j), ...
                real(lambda(ord(j))),imag(lambda(ord(j))),real(mu(ord(j))),imag(mu(ord(j))),abs(c(ord(j))),abs(a(ord(j))),abs(b(ord(j))),abs(condB(ord(j))),abs(condC(ord(j))),tip(ord(j)));
        end
    end
    sub0 = find(tip==0);
    sub1 = find(tip==1);
    sub3 = find(tip==3);
    fprintf('max(a,b)   in group 0: %10.2e,    sep: %10.2e, min(a,b)   in group 1: %10.2e,   min(a,b) in group 3: %10.2e\n',...
        max(mab(sub0)),sep,min(mab(sub1)),min(mab(sub3)))
    fprintf('max(condB) in group 0: %10.2e, YBXtol: %10.2e, min(condB) in group 1: %10.2e, min(condB) in group 3: %10.2e\n',...
        max(condB(sub0)),YBXtol,min(condB(sub1)),min(condB(sub3)))
end

if ~isempty(solution)
    lambda = solution(:,1);
    mu = solution(:,2);
else
    lambda = [];
    mu =  [];
end
       


