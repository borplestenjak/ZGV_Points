function [lambda, mu, mult] = critical_points(A,B,C,opts)

% [lambda, mu, mult] = critical_points(A,B,C,opts) finds critical points (lambda,mu) 
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
%   - sc_steps (0): steps of staircase algorithm for an optional initial
%     reduction of the delta matrices (when it works it might speed up the computation)
%   - membtol (1e-6): treshold for preventing to output multiple instances of critical points
%   - YBXtol (1e-6): treshold for checking the condition y'*B*x=0 for the 2D point
%   - multtol (1e-4): treshold for identifying the multiplicity of mu as an eigenvalus of (A+lambda*B,-C)
%   - other options for singgep or staircase_step_cr_np
%   - inviter (1): use inverse iteration for eigenvectors or slow svd (0)
%   - refine (1): Newton refinement steps to improve the accuracy of simple
%     eigenvalues of a regular MEP - requires computation of eigenvectors 
%   - refineeps (eps): relative tolerance for Newton refinement
%
% Output:
%   - lambda, mu: vectors with critical points
%   - mult: multiplicities of mu as an eigenvalue of A+lambda*B+mu*C for a fixed lambda

% This is Algorithm 1 in B. Plestenjak: On properties and numerical 
% computation of critical points of eigencurves of bivariate matrix pencils

% Bor Plestenjak 2024

if nargin<4, opts=[]; end

if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A,B,C);
end

if isfield(opts,'goal'),      goal = opts.goal;         else,  goal = '2D';    end
if isfield(opts,'sc_steps'),  sc_steps = opts.sc_steps; else,  sc_steps = 0;   end
if isfield(opts,'membtol'),   membtol = opts.membtol;   else,  membtol = 1e-6; end
if isfield(opts,'YBXtol'),    YBXtol = opts.YBXtol;     else,  YBXtol = 1e-6;  end
if isfield(opts,'multtol'),   multtol = opts.multtol;   else,  multtol = 1e-4; end

% Make sure all inputs are of the same numeric type.
if ~isa(A,class_t), A = numeric_t(A,class_t); end
if ~isa(B,class_t), B = numeric_t(B,class_t); end
if ~isa(C,class_t), C = numeric_t(C,class_t); end

lambda = numeric_t([],class_t); 
mu = numeric_t([],class_t); 
n = size(A,2);
Z = zeros(n,class_t);

% linearization for the second equation
A2 = [A Z;  B A];
B2 = [B Z;  Z B];
C2 = [C Z;  Z C];

[Delta0,Delta1,Delta2] = twopar_delta(A,-B,-C,A2,-B2,-C2);

% optional initial staircase reduction of singular pencils
if sc_steps>0
    Delta = {Delta0,Delta1,Delta2};
    for k=1:sc_steps
        Delta = staircase_step_cr_np(Delta,opts);
    end
    Delta0 = Delta{1};  
    Delta1 = Delta{2};  
    Delta2 = Delta{3};
end

% We solve a singular two-parameter eigenvalue problem by first applying
% singgep to the pencil (Delta1,Delta0) and then by inserting computed
% lambda's in the GEP and checking the 2D criteria for the computed mu's. 
lambda_pos = singgep(Delta1,Delta0,opts);
m = length(lambda_pos);
normB = norm(B,'fro');
solution = [];
multiplicities = [];
for k = 1:m
    M = A + lambda_pos(k)*B;
    [X,D,Y] = eig(M,-C);
    mu_pos = diag(D);
    for j = 1:length(mu_pos)
        if isfinite(mu_pos(j))
            if ~is_in_set(solution,[lambda_pos(k) mu_pos(j)],membtol)
                % we estimate multiplicity from the number of close eigenvalues
                pos = abs(mu_pos-mu_pos(j))<multtol*(1+abs(mu_pos(j)));
                alg_mult = sum(pos);
                % condition y'*B*x = 0
                if abs(Y(:,j)'*B*X(:,j))/(norm(Y(:,j))*norm(X(:,j)))<YBXtol*normB 
                    if strcmp(goal,'2D') || ((strcmp(goal,'ZGV') && alg_mult==1)) 
                        solution = [solution; [lambda_pos(k) mu_pos(j)]];
                        multiplicities = [multiplicities; alg_mult];
                    end
                else
                    if strcmp(goal,'2D') && (alg_mult>1)
                        % we estimate geometric multiplicity from the rank of
                        % the matrix composed of eigenvectors
                        subX = X(:,find(pos==1));
                        if rank(subX)>1
                            % geometric multiplicity of mu_0 is larger than 1, this is a 2D point   
                            solution = [solution; [lambda_pos(k) mu_pos(j)]];
                            multiplicities = [multiplicities; alg_mult];
                        end
                    end
                end
            end
        end
    end
end

if ~isempty(solution)
    lambda = solution(:,1);
    mu = solution(:,2);
else
    lambda = [];
    mu =  [];
end
