function [lambda,mu,cand_lambda,cand_mu] = critical_points_MFRD(A,B,C,opts)

% [lambda, mu] = critical_points_MFRD(A,B,C) finds critical points (lambda,mu) 
% for the eigenvalue problem (A + lambda*B + mu*C)x = 0, such that 
% a) 2D points: there exist nonzero vectors x and y such that 
%    (A + lambda*B + mu*C)*x=0, y'*(A + lambda*B + mu*C)=0, y'*B*x=0
% b) ZGV points: 2D points where mu is a simple eigenvalue and y'*C*x ~=0
%
% Input:
%   - A,B,C : square matrices n x n
%   - opts : options
%
% Options in opts:
%   - goal ('2D'): set to 'ZGV' to compute only ZGV points
%   - delta (1e-2+1e-2i): perturbation parameter for the MFRD
%   - membtol (1e-6): treshold for preventing to output multiple instances of critical points
%   - multtol (1e-4): treshold for identifying the multiplicity of mu as an eigenvalus of (A+lambda*B,-C)
%   - show (0): display values used for the identificaton (1) or no (0)
%   - mprefine (0): refine final results in multiple precision (requires Advanpix MCT)
%   - showplot (0): plots convergence graph for Gauss-Newton method
%   - other options for twopareig (e.g., some problem require singular=1)
%     and ZGV_GaussNewton
%
% Output:
%   - lambda, mu: vectors with critical points
%   - lambda, mu: vectors with candidates from the MFRD

% This method is from B. Plestenjak: Numerical methods for the computation 
% of critical points of eigencurves of bivariate matrix pencils

% MFRD method is based on: E. Jarlebring, S. Kvaal, and W. Michiels: Computing 
% all pairs (λ, μ) such that λ is a double eigenvalue of A + μB. 
% SIAM J. Matrix Anal. Appl., 32(3):902–927, 2011.

% Bor Plestenjak 2024

if nargin<4, opts=[]; end

if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A,B,C);
end

if isfield(opts,'goal'),      goal = opts.goal;           else, goal = '2D';           end
if isfield(opts,'delta'),     delta = opts.delta;         else, delta = 1e-2 + 1e-2i;  end
if isfield(opts,'membtol'),   membtol = opts.membtol;     else, membtol = 1e-6;        end
if isfield(opts,'multtol'),   multtol = opts.multtol;     else, multtol = 1e-4;        end
if isfield(opts,'mprefine'),  mprefine = opts.mprefine;   else, mprefine = 0;          end
if isfield(opts,'showplot'),  showplot = opts.showplot;   else, showplot = 0;          end
if isfield(opts,'use_eigs'),  use_eigs = opts.use_eigs;   else, use_eigs = 0;          end
if isfield(opts,'target'),    target = opts.target;       else, target = 0;            end
if isfield(opts,'neigs'),     neigs = opts.neigs;         else, neigs = 10;            end
if isfield(opts,'opts_eigs'), opts_eigs = opts.opts_eigs; else, opts_eigs = [];        end
if isfield(opts,'use_jd'),    use_jd = opts.use_jd;       else, use_jd = 0;            end
if isfield(opts,'opts_jd'),   opts_jd = opts.opts_jd;     else, opts_jd = [];          end


% Make sure all inputs are of the same numeric type.
if ~isa(A,class_t), A = numeric_t(A,class_t); end
if ~isa(B,class_t), B = numeric_t(B,class_t); end
if ~isa(C,class_t), C = numeric_t(C,class_t); end

lambda = numeric_t([],class_t); 
mu = numeric_t([],class_t); 
n = size(A,2);

% we get initial candidates from the MFRD 
% linearization for the second equation of the 2EP
A2 = A;
B2 = (1+delta)*B;
C2 = C;

if use_eigs
    [cand_mu, cand_lambda] = twopareigs(A+target*B,-C,-B,A2+target*B2,-C2,-B2,neigs,opts_eigs);
    cand_lambda = cand_lambda + target;
elseif use_jd
    [cand_lambda, cand_mu] = twopareigs_jd(A,-B,-C,A2,-B2,-C2,neigs,opts_jd);
else
    [cand_lambda, cand_mu] = twopareig(A,-B,-C,A2,-B2,-C2,opts);
end

m = length(cand_lambda);
solution = numeric_t([],class_t); 

if showplot
    figure
end

% we test and refine each candidates with a Gauss-Newton method
for j = 1:m
    [ll,mm,x,y,flag,err,step,init] = ZGV_GaussNewton(A,B,C,cand_lambda(j),cand_mu(j),[],[],opts);
    if mprefine
        [ll,mm,x,y,flag,err,step,init] = ZGV_GaussNewton(mp(A),mp(B),mp(C),mp(ll),mp(mm),mp(x),mp(y));
        ll = double(ll);
        mm = double(mm);
    end
    dif = double(norm([ll mm]- [cand_lambda(j) cand_mu(j)]));
    if flag && ~is_in_set(solution,[ll mm],membtol)  
        if showplot && dif<0.5
            if init==1
                semilogy(0:step-1,err,'r.-','MarkerSize',40)
            else
                semilogy(0:step-1,err,'b.:','MarkerSize',40)
            end
            hold on
        end
        if strcmp(goal,'ZGV')
            M = A + ll*B;
            [X,D,Y] = eig(M,-C);
            mu_pos = diag(D);
            pos = abs(mu_pos-mm)<multtol*(1+abs(mm));
            alg_mult = sum(pos);
            if alg_mult == 1
                solution = [solution; [ll mm]];
            end
        else
            solution = [solution; [ll mm]];
        end
    end
end

if showplot
    hold off
end


if ~isempty(solution)
    lambda = solution(:,1);
    mu = solution(:,2);
else
    lambda = [];
    mu =  [];
end

end

