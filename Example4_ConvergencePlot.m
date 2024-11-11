% Example 7.2 from B. Plestenjak: On properties and numerical computation 
% of critical points of eigencurves of bivariate matrix pencils
%
% This example requires Advanpix MCT

% Bor Plestenjak 2024

A = [1 2 3 0; 2 0 1 0; 3 1 1 0; 0 0 0 -3];
B = [1 0 1 0; 0 1 1 0; 1 1 0 0; 0 0 0 -3];
C = [2 1 0 0; 1 3 0 0; 0 0 1 0; 0 0 0 1];
n = 4;

mp.Digits(100);
rng(1);
opts = [];
opts.delta = 1e-2;
opts.tol = 1e-85;
opts.showplot = 1;

PlotSettings

[lambda, mu, cand_lambda, cand_mu] = critical_points_MFRD(mp(A),mp(B),mp(C),opts);
legend('ZGV points','','','','intersections','Location','SouthWest')
xlabel('iteration')
ylabel('error')