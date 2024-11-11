function p = matchrows_short(A,B,mu1,mu2,p1,p2)

% MATCHROWS(A,B)  returns (subset of) perturbation p that minimizes norm(A-B(p,:))

% This is the assigmnent problem that could be solved by the Hungarian
% algorithm. Here we use a random approach using two random projections to
% a line. We repeat the method until we get the same ordering from both
%
% Optional arguments:
% mu1 and mu2 are initial random vectors (since random generator is slow,
% this way we can this use the same projection for all examples)
% p1 and p2 are initial vectors with zeros and ones

% Bor Plestenjak 2024

class_t = superiorfloat(A,B);
[nA,d] = size(A);
[nB,d] = size(A);

if nargin<3
    mu1 = randn(d,1,class_t)+1i*randn(d,1,class_t);
end
if nargin<4
    mu2 = randn(d,1,class_t)+1i*randn(d,1,class_t);
end
if nargin<5
    p1 = zeros(nA,1,class_t);
end
if nargin<6
    p2 = ones(nA,1,class_t);
end

runs = 0;
while norm(p1-p2)>0 && runs<5
    runs = runs + 1;
    if runs>1
        % if we did not get the same order for both projections we repeat the
        % test with two new projections
        mu1 = randn(d,1,class_t)+1i*randn(d,1,class_t);
        mu2 = randn(d,1,class_t)+1i*randn(d,1,class_t);
    end
    pA1 = A*mu1;
    pB1 = B*mu1;
    pA2 = A*mu2;
    pB2 = B*mu2;
    for row = 1:nA
        % we go through rows in a loop and assign each row to the nearest
        % nonassigned row from B
        [tmp,pos1] = min(abs(pB1-pA1(row)));
        [tmp,pos2] = min(abs(pB2-pA2(row)));
        p1(row) = pos1;
        pB1(pos1) = Inf;
        p2(row) = pos2;
        pB2(pos2) = Inf;
    end
end

if norm(p1-p2)>0
    disp('Permutations are not equal, computing row norms')
    PA = A;
    PB = B;
    p = zeros(nA,1);
    for row = 1:nA
        % we go through rows in a loop and assign each row to the nearest
        % nonassigned row from B
        [tmp,pos1] = min(vecnorm(PB-PA(row,:),2,2));
        p(row) = pos1;
        pB(pos1) = Inf;
    end
else
    p = p1;
end


end