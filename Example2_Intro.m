% Example 2.9 from B. Plestenjak: On properties and numerical computation 
% of critical points of eigencurves of bivariate matrix pencils

% Bor Plestenjak 2024

A = [0 1; 0 0];
B = [1 0; 0 1];
C = [1 0; 0 2];
n = size(A,1);

% computation of all critical points
[lambda,mu] = critical_points(A, B, C);
crit_points = [lambda mu]

opts = [];
opts.goal = 'ZGV';
[lambdaZGV,muZGV] = critical_points(A, B, C, opts);
ZGV_points = [lambdaZGV muZGV]

indr = find(abs(imag(lambda))<1e-6);
lambdar = real(lambda(indr));
mur = real(mu(indr));

indr = find(abs(imag(lambdaZGV))<1e-6);
lambdarZGV = real(lambdaZGV(indr));
murZGV = real(muZGV(indr));

PlotSettings

pts = 1000;
lam = linspace(-1,1,pts);
H = [];
for k = 1:pts
   H(:,k) = eig(A+lam(k)*B,-C); 
end
pl = kron(ones(n,1),lam);
scatter(pl(:),real(H(:)),15,abs(imag(H(:))),'filled')
hold on
plot(lambdar,mur,'.k','MarkerSize',40)
plot(lambdarZGV,murZGV,'.r','MarkerSize',40)
hold off
axis([-1 1 -1 1])