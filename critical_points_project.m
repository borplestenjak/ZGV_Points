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
%   - cplx_rnd (1): set to 0 to use random perturbations or projections
%   - gensize (0): set to 1 to always extract the generic number of finite eigenvalues 
%   - sep (1e2*sqrt(eps)): threshold for the separation of regular eigenvalues
%   - sepinf (1e4*eps): threshold for the separation of finite eigenvalues
%   - multtol (1e-4): treshold for identifying the multiplicity of mu as an eigenvalus of (A+lambda*B,-C)
%   - refine (0): set to 1 to apply Gauss-Newton refinement of computed points
%   - show (0): display values used for the identificaton (1) or no (0)
%   - other options for twopareig
%
% Output:
%   - lambda, mu: vectors with critical points
%   - mult: multiplicities of mu as an eigenvalue of A+lambda*B+mu*C for a fixed lambda

% This is Algorithm 3 in B. Plestenjak: On properties and numerical 
% computation of critical points of eigencurves of bivariate matrix pencils

% Bor Plestenjak 2024

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
if isfield(opts,'sep'),       sep      = opts.sep;        else, sep      = sqrt(1e2*eps(class_t));      end
if isfield(opts,'sepinf'),    sepinf   = opts.sepinf;     else, sepinf   = 1e4*eps(class_t);        end
if isfield(opts,'refine'),    refine   = opts.refine;     else, refine   = 1;      end
if isfield(opts,'multtol'),   multtol = opts.multtol;     else, multtol = 1e-4;    end
if isfield(opts,'use_eigs'),  use_eigs = opts.use_eigs;   else, use_eigs = 0;          end
if isfield(opts,'target'),    target = opts.target;       else, target = 0;            end
if isfield(opts,'neigs'),     neigs = opts.neigs;         else, neigs = 10;            end
if isfield(opts,'opts_eigs'), opts_eigs = opts.opts_eigs; else, opts_eigs = [];        end

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
    % [mu,lambda,x1,x2] = twopareigs(A,-C,-B,AP,-CP,-BP,neig,opts);
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
if show
    disp(' ')
    disp(' i   | p(i)|     re(lambda)     im(lambda) |     re(mu)          im(mu)    |     s(eig) |   ||V''*x|| |   ||U''*y|| | tip')
    disp('-------------------------------------------------------------------------------------------------------------------------')
end

nrmAA = norm(AA,'fro');
nrmBB = norm(BB,'fro');
nrmCC = norm(CC,'fro');
a = zeros(m,1);
b = zeros(m,1);
c = zeros(m,1);
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
    % diff1 = [mx(n*n) mx(n*n+1)]
end
for j = 1:m
    if a(j)<sep && b(j)<sep
        if c(j)>sepinf
            tip(j) = 1;
        else
            tip(j) = 2;
        end
    else
        tip(j) = 3;
    end
    if tip(j) == 1
        if strcmp(goal,'ZGV')
            M = A + lambda(j)*B;
            [X,D,Y] = eig(M,-C);
            mu_pos = diag(D);
            pos = abs(mu_pos-mu(j))<multtol*(1+abs(mu(j)));
            alg_mult = sum(pos);
            if alg_mult == 1
                if refine
                    [ll,mm,x,y,flag,err,step] = ZGV_GaussNewton(A,B,C,lambda(j),mu(j),x1(:,j),y1(:,j));
                    lambda(j) = ll;
                    mu(j) = mm;
                end
                solution = [solution; [lambda(j) mu(j)]];
            end
        else
            if refine
                [ll,mm,x,y,flag,err,step] = ZGV_GaussNewton(A,B,C,lambda(j),mu(j),x1(:,j),y1(:,j));
                lambda(j) = ll;
                mu(j) = mm;
            end
            solution = [solution; [lambda(j) mu(j)]];
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

    for j = 1:m
        fprintf('%4d |%4d | %14.6e %14.6e | %14.6e %14.6e | %10.2e | %10.2e | %10.2e | %1d \n',j, ord(j), ...
                real(lambda(ord(j))),imag(lambda(ord(j))),real(mu(ord(j))),imag(mu(ord(j))),abs(c(ord(j))),abs(a(ord(j))),abs(b(ord(j))),tip(ord(j)));
    end
    % diff = [max(mab(ord(sub1))) min(mab(ord(sub2)))]
end

if ~isempty(solution)
    lambda = solution(:,1);
    mu = solution(:,2);
else
    lambda = [];
    mu =  [];
end
       


