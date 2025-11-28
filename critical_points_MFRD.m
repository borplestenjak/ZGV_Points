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
%   - delta (1e-3): perturbation parameter for the MFRD
%   - membtol (1e-6): treshold for preventing to output multiple instances of critical points
%   - YCXtol (1e-6): treshold for the condition |y'*C*x| > YCXtol*norm(C) for the ZGV point
%   - GStol  (1e-6): treshold for the condition sigma(n-1) > GStol*sigma(1) for geometric simplicity
%   - show (0): display values used for the identificaton (1,2) or no (0)
%   - showplot (0): plots convergence graph for Gauss-Newton method
%   - other options for twopareig (e.g., some problems require singular=1)
%     and ZGV_GaussNewton
%   - use_jd (0): set of 1 to use twopareigs_jd (JAcobi-Davidson)
%   - opts_jd : options for Jacobi-Davidson method (see twopareigs_jd)
%   - use_eigs (0): set of 1 to use twopareigs (Krylov-Schur)
%   - target: target for twopareigs
%   - opts_eigs : options for Sylvester-Arnoldi method (see twopareigs)
%   - neigs: number of critical points for twopareigs and twopareigs_jd
%   - mprefine (0): refine final results in higher precision (requires Advanpix MCT)
%   - opts_mp : options for GN refinement in higher precision
%
% Output:
%   - lambda, mu: vectors with critical points
%   - cand_lambda, cand_mu: vectors with candidates from the MFRD

% This method is from B. Plestenjak: Numerical methods for the computation 
% of critical points of eigencurves of bivariate matrix pencils

% MFRD method is based on: E. Jarlebring, S. Kvaal, and W. Michiels: Computing 
% all pairs (λ, μ) such that λ is a double eigenvalue of A + μB. 
% SIAM J. Matrix Anal. Appl., 32(3):902–927, 2011.

% Bor Plestenjak 2024, revised in 2025

if nargin<4, opts=[]; end

if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A,B,C);
end

if isfield(opts,'goal'),      goal = opts.goal;           else, goal = '2D';       end
if isfield(opts,'delta'),     delta = opts.delta;         else, delta = 1e-3;      end
if isfield(opts,'membtol'),   membtol = opts.membtol;     else, membtol = 1e-6;    end
if isfield(opts,'YCXtol'),    YCXtol = opts.YCXtol;       else, YCXtol = 1e-6;     end
if isfield(opts,'GStol'),     GStol = opts.GStol;         else, GStol = 1e-6;      end
if isfield(opts,'showplot'),  showplot = opts.showplot;   else, showplot = 0;      end
if isfield(opts,'use_eigs'),  use_eigs = opts.use_eigs;   else, use_eigs = 0;      end
if isfield(opts,'target'),    target = opts.target;       else, target = 0;        end
if isfield(opts,'neigs'),     neigs = opts.neigs;         else, neigs = 10;        end
if isfield(opts,'opts_eigs'), opts_eigs = opts.opts_eigs; else, opts_eigs = [];    end
if isfield(opts,'mprefine'),  mprefine = opts.mprefine;   else, mprefine = 0;      end
if isfield(opts,'opts_mp'),   opts_mp = opts.opts_mp;     else, opts_mp = [];      end
if isfield(opts,'use_jd'),    use_jd = opts.use_jd;       else, use_jd = 0;        end
if isfield(opts,'opts_jd'),   opts_jd = opts.opts_jd;     else, opts_jd = [];      end
if ~isfield(opts,'svmult'),   opts.svmult = delta;        end

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
    [cand_lambda, cand_mu,XR,YR,XL,YL] = twopareig(A,-B,-C,A2,-B2,-C2,opts);
end

m = length(cand_lambda);
solution = numeric_t([],class_t); 

if showplot
    figure
end

normC = norm(C);
normB = norm(B);

% we test and refine each candidate with a Gauss-Newton method using up to
% two different normalization vectors
for j = 1:m
     if ~use_eigs && ~use_jd && abs(XL(:,j)'*B*XR(:,j))/normB<1e-2
        [ll,mm,x,y,flag,err,step,init] = ZGV_GaussNewton(A,B,C,cand_lambda(j),cand_mu(j),XR(:,j),XL(:,j),opts);
     else
        [ll,mm,x,y,flag,err,step,init] = ZGV_GaussNewton(A,B,C,cand_lambda(j),cand_mu(j),[],[],opts);
    end
    % one more try if there is no convergence due to unlucky choice of
    % normalizing vectors in the GN method
    if ~flag
        [ll,mm,x,y,flag,err,step,init] = ZGV_GaussNewton(A,B,C,cand_lambda(j),cand_mu(j),[],[],opts);
    end
    if flag && mprefine
        [ll,mm,x,y,flagmp,errmp,stepmp,initmp] = ZGV_GaussNewton(mp(A),mp(B),mp(C),mp(ll),mp(mm),mp(x),mp(y),opts_mp);
        ll = double(ll);
        mm = double(mm);
        x = double(x);
        y = double(y);
    end
    dif = double(norm([ll mm]- [cand_lambda(j) cand_mu(j)]))/double(1+norm([cand_lambda(j) cand_mu(j)]));
    if flag && ~is_in_set(solution,[ll mm],membtol) && (dif<0.2)
        if showplot 
            if init<=1 
                semilogy(0:step-1,err,'r.-','MarkerSize',40)
            else
                semilogy(0:step-1,err,'kd--','MarkerSize',14,'MarkerFaceColor','w')
            end
            hold on
        end
        if strcmp(goal,'ZGV')
            % if we are computing ZGV points, we first check if |y'*C*x| is
            % nonzero and mu is geometrically simple
            condC = abs(y'*C*x)/(norm(y)*norm(x));
            if (condC>YCXtol*normC)
                % if this is true, we additionally check if mu is geometrically simple
                M = A + ll*B + mm*C;
                s = svd(M);
                if s(n-1) > GStol*s(1)
                    solution = [solution; [ll mm]];
                end
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

