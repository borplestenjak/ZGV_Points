% Example 7.5 (right hand side Figure 7.5) from B. Plestenjak: On properties 
% and numerical computation of critical points of eigencurves of bivariate 
% matrix pencils

% Matrix A0 is from Example 5 in M. A. Freitag and A. Spence. A Newton-based 
% method for the calculation of the distance to instability. Linear Algebra 
% Appl., 435(12):3189–3205, 2011.

% Bor Plestenjak 2024

% -----------------------------------------------------------------------
% next setting increases chances of reproducible results 
% comment it if speed is more important then reproducibility
maxNumCompThreads(1); % could 
rng(1,'twister')
% -----------------------------------------------------------------------

A0 = diag([-0.4 + 6i, -0.1+1i, -1-3i, -5+1i]) + diag([1 1 1],1) + diag([1 1 1],-1);
n0 = size(A0,1);

A = [zeros(n0) A0; A0' zeros(n0)];
B = -[zeros(n0) 1i*eye(n0); -1i*eye(n0) zeros(n0)];
C = -eye(2*n0);

opts = [];
opts.GN_refine = 1;
[lambda0,mu0] = critical_points(A,B,C,opts);

indr = find(abs(imag(lambda0))<1e-6);
lambdar = real(lambda0(indr));
mur = real(mu0(indr));

pos1 = find(mur>0);
lambdarr = lambdar(pos1);
murr = mur(pos1);
[mu0,k] = min(murr);
lambda0 = lambdarr(k);
dti = mu0

close all

PlotSettings

% first figure
fig1 = figure;
n = 2*n0;
pts = 3000;
lam = linspace(-3,6,pts);
H = [];
for k = 1:pts
    tmp = sort(eig(A+lam(k)*B,-C)); 
    H(:,k) = tmp; 
end
pl = kron(ones(2*n0,1),lam);
scatter(pl(:),real(H(:)),10,abs(imag(H(:))),'filled')
hold on

% find point that are not optimal points
g = zeros(length(lambdar),1);
for k = 1:length(lambdar)
    g(k) = is_in_set([lambda0 mu0],[lambdar(k) mur(k)]);
end
ind2 = find(g==0);

plot(lambdar(ind2),mur(ind2),'.r','MarkerSize',40)
plot(lambda0,mu0,'pb','MarkerSize',16,'MarkerFaceColor','b')
hold off

axis([-3 6 -6.5 6.5])
xlabel('\lambda')
ylabel('\mu')


% second figure
fig2 = copyobj(fig1, groot);
axis([0 1.8 4.5 6.5])
xlabel('\lambda')
ylabel('\mu')

% third figure
fig3 = copyobj(fig1, groot);
axis([-1.5 4 -0.25 1.75])
xlabel('\lambda')
ylabel('\mu')
