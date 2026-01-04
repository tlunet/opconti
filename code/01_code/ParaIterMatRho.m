clc; clear; close all

nu = 1e0; %penalization parameters
L = 1; Nx = 2^5; x = linspace(0, L, Nx+1); hx = x(2) - x(1); %space grid
T = 1; Nt = 2^5; t = linspace(0, T, Nt+1); ht = t(2) - t(1); %time grid
a1 = 1; a2 = 2; na = a1+a2; %nb overlap points in time
n1 = Nt/2+a1; n2 = Nt/2+a2; %subdomains points in time 

%space matrix
e = ones(Nx-1, 1);
A = spdiags([-e/hx^2 2/hx^2*e -e/hx^2], [-1 0 1], Nx-1, Nx-1);
%state matrix
At = kron(speye(Nt), speye(Nx-1)+ht*A/2) ...
    - kron(spdiags(e, -1, Nt, Nt), speye(Nx-1)-ht*A/2);
H = nu*speye(Nt*(Nx-1), Nt*(Nx-1)) + (At*At')^(-1);

%restriction matrices
R1 = [speye(n1*(Nx-1)) sparse(n1*(Nx-1), (Nt-n1)*(Nx-1))];
R2 = [sparse(n2*(Nx-1), (Nt-n2)*(Nx-1)) speye(n2*(Nx-1))];

%each block matrix for u1
H1ff = R1*H*R1'; H1cc = H(n1*(Nx-1)+1:end, n1*(Nx-1)+1:end);
H1fc = H(1:n1*(Nx-1), n1*(Nx-1)+1:end); 
H1cf = H(n1*(Nx-1)+1:end, 1:n1*(Nx-1));
[H1ff H1fc; H1cf H1cc] - H %check if zero

%each block matrix for u2
H2ff = R2*H*R2'; H2cc = H(1:(Nt-n2)*(Nx-1), 1:(Nt-n2)*(Nx-1));
H2fc = H((Nt-n2)*(Nx-1)+1:end, 1:(Nt-n2)*(Nx-1)); 
H2cf = H(1:(Nt-n2)*(Nx-1), (Nt-n2)*(Nx-1)+1:end);
[H2cc H2cf; H2fc H2ff] - H %check if zero

%partition
R1t = R1; 
R1t((n1-na)*(Nx-1)+1:n1*(Nx-1), (n1-na)*(Nx-1)+1:n1*(Nx-1)) ...
    = 1/2*speye(na*(Nx-1));
R2t = R2; 
R2t(1:na*(Nx-1), (Nt-n2)*(Nx-1)+1:(Nt-n2+na)*(Nx-1)) ...
    = 1/2*speye(na*(Nx-1));
R1t'*R1 + R2t'*R2 %check if identity

%two Schur complements
S1 = H1ff - H1fc*(H1cc\H1cf);
S2 = H2ff - H2fc*(H2cc\H2cf);

%iterative matrix for u
al = 0.19;
rho = R1t'*(speye(n1*(Nx-1)) - al*S1)*R1 ...
    - R2t'*(speye(n2*(Nx-1)) - al*S2)*R2;
max(abs(eig(full(rho))))