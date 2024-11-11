% Example 7.3 from B. Plestenjak: On properties and numerical computation 
% of critical points of eigencurves of bivariate matrix pencils
%
% The timings and error were obtained by running Example5_Figures for 
% n = 5:5:45

% Bor Plestenjak 2024

n = 5:5:45;
t1 = [6.02E-03 8.91E-02 6.36E-01 2.85E+00 1.25E+01 3.99E+01 1.14E+02 2.14E+02 5.24E+02];
t2 = [8.70E-03 5.87E-02 2.10E-01 6.53E-01 1.70E+00 4.40E+00 1.09E+01 2.31E+01 4.55E+01];
t3 = [1.70E-02 7.37E-02 1.91E-01 4.98E-01 1.19E+00 2.20E+00 3.88E+00 6.64E+00 1.16E+01];

e1 = [6.70E-13 1.50E-12 1.10E-11 4.30E-11 7.40E-11 1.60E-10 1.30E-10 1.20E-10 6.10E-10];
e2 = [2.00E-12 9.00E-11 3.20E-09 3.70E-08 5.30E-08 2.30E-07 2.00E-07 1.30E-07 2.70E-06];
e3 = [3.80E-15 1.30E-13 1.40E-13 2.90E-13 3.90E-13 4.90E-13 6.20E-13 7.50E-13 7.20E-13];

PlotSettings

semilogy(n,t1,'r*-','MarkerSize',10)
hold on
semilogy(n,t2,'gs--',n,t3,'bd:')
legend('Alg. 1','Alg. 3','MFRD + Alg. 4','Location','NorthWest')
ylabel('time in seconds')
xlabel('n')

figure
semilogy(n,e1,'r*-','MarkerSize',10)
hold on
semilogy(n,e2,'gs--',n,e3,'bd:')
hold off
legend('Alg. 1','Alg. 3','MFRD + Alg. 4','Location','NorthWest')
ylabel('error')
xlabel('n')
axis([0 50 1e-15 1e-1])
