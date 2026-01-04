clc; clear; close all

nu = 1e0; %gam = 1e0; %penalization parameters
L = 1; Nx = 2^3; x = linspace(0, L, Nx+1); hx = x(2) - x(1); %space grid
T = 1; Nt = 2^3; t = linspace(0, T, Nt+1); ht = t(2) - t(1); %time grid
gam = (pi^2*(2*T^2+T) + 4*T + 1) / (pi^4*(2*T^2+T) - 4); %special gamma

%laplace matrix
e = ones(Nx-1, 1); 
A = spdiags([-e/hx^2 2/hx^2*e -e/hx^2], [-1 0 1], Nx-1, Nx-1);

Niter = 10; %nb of iterations


%Discretize-then-optimize (Dto)
AtDto = BuildCNMatrixDto(Nt, ht, Nx, A, nu, gam);
FDto = BuildCNRhsDto(Nt, ht, Nx, yExact(x(2:end-1), t(1)), ...
    yTarget(x(2:end-1), t(2:end), nu), gam);
UDto = AtDto\FDto; UDto = reshape(UDto, Nx-1, 2*Nt);
yDto = UDto(:, 1:Nt); lamDto = UDto(:, Nt+1:end);

% figure
% mesh(t, x, yExact(x, t))
% figure
% mesh(t(2:end), x(2:end-1), yDto)
% figure
% mesh(t, x, lamExact(x, t, nu))
% figure
% mesh(t(2:end), x(2:end-1), lamDto)

norm(ht*(yDto - yExact(x(2:end-1), t(2:end))))
norm(ht*(lamDto - lamExact(x(2:end-1), t(2:end), nu)))

eyDto = [0.0324 0.0087 0.0023 5.8807e-04 1.4898e-04 3.7485e-05];
elamDto = [1.3517 0.7274 0.3776 0.1924 0.0971 0.0488];
k = 3:8; hk = 2.^(-k);
figure
loglog(hk, eyDto, '-o', hk, elamDto, '-s', ...
    hk, hk/hk(1), '-*', hk, (hk/hk(1)).^2, '-+', ...
    'LineWidth', 2, 'MarkerSize', 10);
xlabel('mesh size $h_x = h_t$', 'Interpreter', 'latex');
ylabel('error $l^2$ norm', 'Interpreter', 'latex');
legend({'error in $y$', 'error in $\lambda$', 'order 1', 'order 2'}, ...
    'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', 25); set(gca, 'linewidth', 1.5);


%Optimize-then-discretize (Otd)
AtOtd = BuildCNMatrixOtd(Nt, ht, Nx, A, nu, gam);
FOtd = BuildCNRhsOtd(Nt, ht, Nx, yExact(x(2:end-1), t(1)), ...
    yTarget(x(2:end-1), t, nu), gam);
UOtd = AtOtd\FOtd; UOtd = reshape(UOtd, Nx-1, 2*Nt+2);
yOtd = UOtd(:, 1:Nt+1); lamOtd = UOtd(:, Nt+2:end);

% figure
% mesh(t, x, yExact(x, t))
% figure
% mesh(t, x(2:end-1), yOtd)
% figure
% mesh(t, x, lamExact(x, t, nu))
% figure
% mesh(t, x(2:end-1), lamOtd)

norm(ht*(yOtd - yExact(x(2:end-1), t)))
norm(ht*(lamOtd - lamExact(x(2:end-1), t, nu)))

eyOtd = [0.0259 0.0062 0.0015 3.7515e-04 9.3349e-05 2.3283e-05];
elamOtd = [0.1484 0.0358 0.0089 0.0022 5.5318e-04 1.3827e-04];
k = 3:8; hk = 2.^(-k);
figure
loglog(hk, eyOtd, '-o', hk, elamOtd, '-s', ...
    hk, hk/hk(1), '-*', hk, (hk/hk(1)).^2, '-+', ...
    'LineWidth', 2, 'MarkerSize', 10);
xlabel('mesh size $h_x = h_t$', 'Interpreter', 'latex');
ylabel('error $l^2$ norm', 'Interpreter', 'latex');
legend({'error in $y$', 'error in $\lambda$', 'order 1', 'order 2'}, ...
    'Interpreter', 'latex', 'Location', 'best');
set(gca, 'FontSize', 25); set(gca, 'linewidth', 1.5);

%------------------------------%
%   Test with exact solution   %
%------------------------------%
function g = yExact(x, t)
g = sin(pi*x').*(2*t.^2 + t);
end
function g = lamExact(x, t, nu)
g = -nu*sin(pi*x').*(pi^2*(2*t.^2 + t) + 4*t + 1);
end
function g = yTarget(x, t, nu)
g = nu*sin(pi*x').*((pi^4+1/nu)*(2*t.^2 + t) - 4);
end
%------------------------------%
%   Discretize-then-optimize   %
%------------------------------%
function At = BuildCNMatrixDto(N, dt, J, A, nu, gam)
e = ones(N, 1);
%construct four matrices
At11 = kron(speye(N), speye(J-1)+dt*A/2) ...
    - kron(spdiags(e, -1, N, N), speye(J-1)-dt*A/2);
At21 = kron(speye(N), -dt*speye(J-1));
At12 = kron(spdiags([e e], [-1 1], N, N), dt*speye(J-1)/4/nu) - At21/2/nu;
At22 = At11';

%adapt the initial and final condition
At12(1:J-1, 1:J-1) = 3*dt*speye(J-1)/4/nu;
At12(end-J+2:end, end-J+2:end) = 3*dt*speye(J-1)/4/nu;
At21(end-J+2:end, end-J+2:end) = At21(end-J+2:end, end-J+2:end)/2 ...
    - gam*speye(J-1);

%concatenate four matrices
At = [At11, At12; At21, At22];
end
function F = BuildCNRhsDto(N, dt, J, Y0, Yhat, gam)
%construct F1 to solve Y
F1 = zeros(N*(J-1), 1);
F1(1:J-1) = Y0;

%construct F2 to solve Lambda
Yhat = Yhat(:);
F2 = -dt*Yhat;
F2(end-J+2:end) = F2(end-J+2:end)/2 - gam*Yhat(end-J+2:end);

% concatenate F1 and F2
F = [F1; F2];
end
%------------------------------%
%   Optimize-then-discretize   %
%------------------------------%
function At = BuildCNMatrixOtd(N, dt, J, A, nu, gam)
e = ones(N+1, 1);
%construct four matrices
At11 = kron(speye(N+1), speye(J-1)+dt*A/2) ...
    - kron(spdiags(e, -1, N+1, N+1), speye(J-1)-dt*A/2);
At21 = kron(spdiags([e e], [0 1], N+1, N+1), dt*speye(J-1)/2);
At12 = At21'/nu;
At22 = -At11';

%adapt the initial and final condition
At11(1:J-1, 1:J-1) = speye(J-1);
At12(1:J-1, 1:J-1) = sparse(J-1, J-1);
At21(end-J+2:end, end-J+2:end) = -gam*speye(J-1);
At22(end-J+2:end, end-J+2:end) = speye(J-1);

%concatenate four matrices
At = [At11, At12; At21, At22];
end
function F = BuildCNRhsOtd(N, dt, J, Y0, Yhat, gam)
%construct F1 to solve Y
F1 = zeros((N+1)*(J-1), 1);
F1(1:J-1) = Y0;

%construct F2 to solve Lambda
F2 = dt*(Yhat(:, 1:end-1) + Yhat(:, 2:end))/2;
F2 = [F2(:); -gam*Yhat(:, end)];

% concatenate F1 and F2
F = [F1; F2];
end