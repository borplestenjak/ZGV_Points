% Example 7.1 from B. Plestenjak: On properties and numerical computation 
% of critical points of eigencurves of bivariate matrix pencils

% Bor Plestenjak 2024

A = [ 3  0;  0  0];
B = [ 0  1; -1 -1];
C = [-2 -2;  2  0];
n = 2;

opts = [];
rng(1) % so that results are reproducible

% first option - we compute ZGV points using Algorithm 1 
[lambda1,mu1] = critical_points(A,B,C);
sol1 = [lambda1 mu1]

indr = find(abs(imag(lambda1))<1e-6);
lambdar = real(lambda1(indr));
mur = real(mu1(indr));

% second option - we compute ZGV points using Algorithm 2
opts2 = [];
opts2.show = 1;
[lambda2,mu2] = critical_points_project(A,B,C,opts2);
sol2 = [lambda2 mu2]

% third option - we compute ZGV points using the MRFD and Gauss-Newton method
opts3 = [];
opts3.delta = 1e-2;
[lambda3,mu3,cand_lambda3,cand_mu3] = critical_points_MFRD(A,B,C,opts3);
candidates3 = [cand_lambda3 cand_mu3]
sol3 = [lambda3 mu3]

PlotSettings

figure
pts = 1000;
lam = linspace(0,4,pts);
H = [];
for k = 1:pts
   H(:,k) = eig(A+lam(k)*B,-C); 
end
pl = kron(ones(n,1),lam);
scatter(pl(:),real(H(:)),15,0*abs(imag(H(:))),'filled')
hold on
lam2 = linspace(-0.5,1.5,pts);
Hy = [];
for k = 1:pts
   Hy(:,k) = eig(A+lam2(k)*C,-B); 
end
pl2 = kron(ones(n,1),lam2);
scatter(real(Hy(:)),pl2(:),15,0*abs(imag(Hy(:))),'filled')
plot(lambdar,mur,'.r','MarkerSize',40)
hold off
axis equal
axis([-0.5 4.5 -1 2])
