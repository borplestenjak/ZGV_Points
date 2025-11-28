% Example 7.8 from B. Plestenjak: On properties and numerical computation 
% of critical points of eigencurves of bivariate matrix pencils

% Bor Plestenjak 2024

% -----------------------------------------------------------------------
% next setting increases chances of reproducible results 
% comment it if speed is more important then reproducibility
maxNumCompThreads(1); % could 
rng(1,'twister')
% -----------------------------------------------------------------------

% we set up parameters for bde2mep so that we can get a discretizaton by
% Chebyshev collocation for the
% y''(x) - 2*lambda*cos(2x)*y(x)+mu*y(x) = 0
% and b.c. y'(0)=y'(pi/2)=0

a = 0;
b = pi/2;
p = 1;
q = 0;
r = 0;
s = @(x) -2*cos(2*x);
t = 1;
bc = [1 0; 1 0];

n = 30;
opts = [];

[z,A,B,C,G,k,rr] = bde2mep(a,b,p,q,r,s,t,bc,n,opts);

[lambda,mu] = critical_points(A,B,C);

indr = find(abs(imag(lambda))<1e-6);
lambdar = real(lambda(indr));
mur = real(mu(indr));

n = 50;
[z,A,B,C,G,k,rr] = bde2mep(a,b,p,q,r,s,t,bc,n,opts);

n = size(A,1);
pts = 2000;
lam = linspace(-80,80,pts);
H = [];
for k = 1:pts
    tmp = sort(eig(A+lam(k)*B,-C)); 
    H(:,k) = tmp; 
end
pl = kron(ones(n,1),lam);
scatter(pl(:),real(H(:)),10,abs(imag(H(:))),'filled')
hold on
plot(lambdar,mur,'.r','MarkerSize',40)
hold off
axis([-80 80 -10 90])
xlabel('\lambda')
ylabel('\mu')

% we can iteratively refine the solution using finer discretization
target = [11 17];
n = 25;
[z,A,B,C,G,k,rr] = bde2mep(a,b,p,q,r,s,t,bc,n,opts);
[lambda5a,mu5a,x,y,flag,err] = ZGV_GaussNewton(A,B,C,target(1),target(2));
sol1a = [lambda5a mu5a]
n = 50;
[z,A,B,C,G,k,rr] = bde2mep(a,b,p,q,r,s,t,bc,n,opts);
[lambda5b,mu5b,x,y,flag,err] = ZGV_GaussNewton(A,B,C,lambda5a,mu5a);
sol1b = [lambda5b mu5b]
n = 100;
[z,A,B,C,G,k,rr] = bde2mep(a,b,p,q,r,s,t,bc,n,opts);
[lambda5c,mu5c,x,y,flag,err] = ZGV_GaussNewton(A,B,C,lambda5b,mu5b);
sol1c = [lambda5c mu5c]

target = [31 42];
n = 25;
[z,A,B,C,G,k,rr] = bde2mep(a,b,p,q,r,s,t,bc,n,opts);
[lambda5a,mu5a,x,y,flag,err] = ZGV_GaussNewton(A,B,C,target(1),target(2));
sol2 = [lambda5a mu5a]
n = 50;
[z,A,B,C,G,k,rr] = bde2mep(a,b,p,q,r,s,t,bc,n,opts);
[lambda5b,mu5b,x,y,flag,err] = ZGV_GaussNewton(A,B,C,lambda5a,mu5a);
sol2b = [lambda5b mu5b]
n = 100;
[z,A,B,C,G,k,rr] = bde2mep(a,b,p,q,r,s,t,bc,n,opts);
[lambda5c,mu5c,x,y,flag,err] = ZGV_GaussNewton(A,B,C,lambda5b,mu5b);
sol2c = [lambda5c mu5c]

target = [60 78];
n = 25;
[z,A,B,C,G,k,rr] = bde2mep(a,b,p,q,r,s,t,bc,n,opts);
[lambda5a,mu5a,x,y,flag,err] = ZGV_GaussNewton(A,B,C,target(1),target(2));
sol3 = [lambda5a mu5a]
n = 50;
[z,A,B,C,G,k,rr] = bde2mep(a,b,p,q,r,s,t,bc,n,opts);
[lambda5b,mu5b,x,y,flag,err] = ZGV_GaussNewton(A,B,C,lambda5a,mu5a);
sol3b = [lambda5b mu5b]
n = 100;
[z,A,B,C,G,k,rr] = bde2mep(a,b,p,q,r,s,t,bc,n,opts);
[lambda5c,mu5c,x,y,flag,err] = ZGV_GaussNewton(A,B,C,lambda5b,mu5b);
sol3c = [lambda5c mu5c]
