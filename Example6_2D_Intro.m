% Example 7.4 from B. Plestenjak: On properties and numerical computation 
% of critical points of eigencurves of bivariate matrix pencils
%
% This example is taken from Lu, Su, Bai: SIMAX 45 (2024) 1455-1486

A = [2 0 1;0 0 1;1 1 0];
C = -eye(3);
B = -[1 0 1;0 1 1;1 1 0];
n = size(A,1);

% computation of all critical points
rng(1)
opts.membtol = 1e-15; % otherwise the point (1,0) is computed only once
[lambda,mu] = critical_points(A, B, C, opts);
crit_points = [lambda mu]

indr = find(abs(imag(lambda))<1e-6);
lambdar = real(lambda(indr));
mur = real(mu(indr));

PlotSettings

pts = 1000;
lam = linspace(0,2,pts);
H = [];
for k = 1:pts
   H(:,k) = eig(A+lam(k)*B,-C); 
end
pl = kron(ones(n,1),lam);
scatter(pl(:),real(H(:)),15,abs(imag(H(:))),'filled')
hold on
plot(lambdar,mur,'.r','MarkerSize',40)
hold off
axis([0 2 -3 3])

