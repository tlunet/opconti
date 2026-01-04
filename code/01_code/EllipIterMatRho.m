clc; clear; close all

nu = 1e-1; %penalization parameters
N = 2^5; h = 1/(N+1); %space grid
a1 = 1; a2 = 1; na = a1+a2+1; %nb overlap points
n1 = N/2+a1; n2 = N/2+a2; %subdomains interior points 

%define matrix
e = ones(N-1, 1);
A = spdiags([-e/h^2 2/h^2*e -e/h^2], [-1 0 1], N-1, N-1);
H = nu*speye(N-1, N-1) + (A*A')^(-1);

%restriction matrices
R1 = [speye(n1) sparse(n1, N-1-n1)];
R2 = [sparse(n2, N-1-n2) speye(n2)];

%each block matrix for u1
H1ff = R1*H*R1'; H1cc = H(n1+1:end, n1+1:end);
H1fc = H(1:n1, n1+1:end); H1cf = H(n1+1:end, 1:n1);
[H1ff H1fc; H1cf H1cc] - H %check if zero

%each block matrix for u2
H2ff = R2*H*R2'; H2cc = H(1:N-1-n2, 1:N-1-n2);
H2fc = H(N-n2:end, 1:N-1-n2); H2cf = H(1:N-1-n2, N-n2:end);
[H2cc H2cf; H2fc H2ff] - H %check if zero

%partition
R1t = R1; R1t(n1-na+1:n1, n1-na+1:n1) = 1/2*speye(na);
R2t = R2; R2t(1:na, N-n2:N-n2+na-1) = 1/2*speye(na);
R1t'*R1 + R2t'*R2 %check if identity

%two Schur complements
S1 = H1ff - H1fc*(H1cc\H1cf);
S2 = H2ff - H2fc*(H2cc\H2cf);

%iterative matrix for u
al = 10;
rho = R1t'*(speye(n1) - al*S1)*R1 - R2t'*(speye(n2) - al*S2)*R2;
max(abs(eig(full(rho))))
