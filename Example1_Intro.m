% Example 2.8 from B. Plestenjak: On properties and numerical computation 
% of critical points of eigencurves of bivariate matrix pencils

% Bor Plestenjak 2024

A = [1 2 3 0; 2 0 1 0; 3 1 1 0; 0 0 0 -3];
B = [1 0 1 0; 0 1 1 0; 1 1 0 0; 0 0 0 -3];
C = [2 1 0 0; 1 3 0 0; 0 0 1 0; 0 0 0 1];
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
lam = linspace(-4,4,pts);
H = [];
for k = 1:pts
   H(:,k) = eig(A+lam(k)*B,-C); 
end
pl = kron(ones(n,1),lam);
scatter(pl(:),real(H(:)),15,abs(imag(H(:))),'filled')
hold on

% find 2D point that are not ZGV points
g = zeros(length(lambdar),1);
for k = 1:length(lambdar)
    g(k) = is_in_set([lambdarZGV murZGV],[lambdar(k) mur(k)]);
end
ind2 = find(g==0);
plot(lambdarZGV,murZGV,'.r','MarkerSize',40)
plot(lambdar(ind2),mur(ind2),'ok','MarkerSize',14,'MarkerFaceColor','w')
hold off
axis([-3 1 -4 4])
xlabel('\lambda')
ylabel('\mu')