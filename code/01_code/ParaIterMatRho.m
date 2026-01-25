clc; clear; close all

nu = 1e0; %penalization parameters
L = 1; Nx = 2^5; x = linspace(0, L, Nx+1); hx = x(2) - x(1); %space grid
T = 1; Nt = 2^5; t = linspace(0, T, Nt+1); ht = t(2) - t(1); %time grid
a1 = 2; a2 = 2; na = a1+a2; %nb overlap points in time
n1 = Nt/2+a1; n2 = Nt/2+a2; %subdomains points in time 

%space matrix
e = ones(Nx-1, 1);
A = spdiags([-e/hx^2 2/hx^2*e -e/hx^2], [-1 0 1], Nx-1, Nx-1);
%build space-time matrix with Crank-Nicolson (optimize-then-discretize)
[At11, At13, At21, At22, At33] = BuildCNMatrixOtd(Nt, ht, Nx, A, nu, 0);
H = At33 + At22\(At21*(At11\At13));
% %build space-time matrix with Crank-Nicolson (discretize-then-optimize)
% [At11, At13, At21, At22, At32, At33] = BuildCNMatrixDto(Nt, ht, Nx, A, nu, 0);
% H = At33 + At32*(At22\(At21*(At11\At13)));

%restriction matrices
R1 = [speye(n1*(Nx-1)) sparse(n1*(Nx-1), (Nt-n1)*(Nx-1))];
R2 = [sparse(n2*(Nx-1), (Nt-n2)*(Nx-1)) speye(n2*(Nx-1))];

%each block matrix for u1
H1ff = R1*H*R1'; H1cc = H(n1*(Nx-1)+1:end, n1*(Nx-1)+1:end);
H1fc = H(1:n1*(Nx-1), n1*(Nx-1)+1:end); 
H1cf = H(n1*(Nx-1)+1:end, 1:n1*(Nx-1));
% [H1ff H1fc; H1cf H1cc] - H %check if zero

%each block matrix for u2
H2ff = R2*H*R2'; H2cc = H(1:(Nt-n2)*(Nx-1), 1:(Nt-n2)*(Nx-1));
H2fc = H((Nt-n2)*(Nx-1)+1:end, 1:(Nt-n2)*(Nx-1)); 
H2cf = H(1:(Nt-n2)*(Nx-1), (Nt-n2)*(Nx-1)+1:end);
% [H2cc H2cf; H2fc H2ff] - H %check if zero

%partition
R1t = R1; 
R1t((n1-na)*(Nx-1)+1:n1*(Nx-1), (n1-na)*(Nx-1)+1:n1*(Nx-1)) ...
    = 1/2*speye(na*(Nx-1));
R2t = R2; 
R2t(1:na*(Nx-1), (Nt-n2)*(Nx-1)+1:(Nt-n2+na)*(Nx-1)) ...
    = 1/2*speye(na*(Nx-1));
% for j = 1:na
%     R1t(end-j*(Nx-1)+1:end-(j-1)*(Nx-1), ...
%         (n1-j)*(Nx-1)+1:(n1-j+1)*(Nx-1)) = 1/2*speye(Nx-1);
%     R2t((j-1)*(Nx-1)+1:j*(Nx-1), ...
%         (Nt-n2+j-1)*(Nx-1)+1:(Nt-n2+j)*(Nx-1)) = 1/2*speye(Nx-1);
% end
% R1t'*R1 + R2t'*R2 %check if identity

S1 = H1ff - H1fc*(H1cc\H1cf); %Schur complements for u1ff
S2 = H2ff - H2fc*(H2cc\H2cf); %Schur complements for u2ff

%iterative matrix from gradient descenet 
al = 1; %step size
rho = R1t'*(speye(n1*(Nx-1)) - al*S1)*R1 ...
    + R2t'*(speye(n2*(Nx-1)) - al*S2)*R2;
max(abs(eig(full(rho)))) %spectral radius

%------------------------------%
%   Optimize-then-discretize   %
%------------------------------%
function [At11, At13, At21, At22, At33] = BuildCNMatrixOtd(N, dt, J, A, nu, gam)
e = ones(N, 1);
%construct four matrices
At11 = kron(speye(N), speye(J-1)+dt*A/2) ...
    - kron(spdiags(e, -1, N, N), speye(J-1)-dt*A/2);
At21 = kron(spdiags([e e], [0 1], N, N), -dt*speye(J-1)/2);
At13 = At21';
At22 = At11';
At33 = nu*speye(N*(J-1), N*(J-1));

%adapt the initial and final condition
At11(1:J-1, 1:J-1) = speye(J-1);
At13(1:J-1, 1:J-1) = sparse(J-1, J-1);
At21(end-J+2:end, end-J+2:end) = -gam*speye(J-1);
At22(end-J+2:end, end-J+2:end) = speye(J-1);
end
%------------------------------%
%   Discretize-then-optimize   %
%------------------------------%
function [At11, At13, At21, At22, At32, At33] = BuildCNMatrixDto(N, dt, J, A, nu, gam)
e = ones(N-1, 1);
%construct four matrices
At11 = kron(speye(N-1), speye(J-1)+dt*A/2) ...
    - kron(spdiags(e, -1, N-1, N-1), speye(J-1)-dt*A/2);
At13 = kron(spdiags(ones(N, 2), [0 1], N-1, N), -dt*speye(J-1)/2);
At21 = kron(speye(N-1), -dt*speye(J-1));
At22 = At11';
At32 = kron(spdiags(ones(N, 2), [-1 0], N, N-1), speye(J-1)/2);
At33 = nu*speye(N*(J-1), N*(J-1));

%adapt the final condition
At21(end-J+2:end, end-J+2:end) = At21(end-J+2:end, end-J+2:end)/2 ...
    - gam*speye(J-1);
At32(1:J-1, 1:J-1) = speye(J-1);
At32(end-J+2:end, end-J+2:end) = speye(J-1);
end