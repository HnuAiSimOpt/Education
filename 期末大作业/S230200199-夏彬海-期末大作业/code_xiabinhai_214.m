%---main---
clc
clear;

E = 100;        % Young's modulus in MPa,
r = 50;         % radius of the rod in mm,
I = pi * r^4 / 4;   % area moment of inertia in mm^4,
A = pi * r^2;       % cross-section area in mm^2,

% node information. 1st column is node number,
% 2nd, 3rd and 4th column are X, Y, Z coordinate.
node = [1 0 0 0;
        2 171.010 0 469.846;
        3 457.798 0 879.422;
        4 890.811 0 1129.422;
        5 1390.811 0 1129.422];

% element information. 1st column is element number,
% 2nd and 3rd column are the numbers of the nodes it includes.
ele = [1 1 2;
       2 2 3;
       3 3 4;
       4 4 5];

n_ele = length(ele(:, 1));       % total number of the elements
L = zeros(1, n_ele);         % beam length
ang = zeros(1, n_ele);           % beam horizontal angle

for i = 1 : n_ele
    L(i) = sqrt((node(ele(i, 3), 2) - node(ele(i, 2), 2))^2 ...
         + (node(ele(i, 3), 4) - node(ele(i, 2), 4))^2);
    ang(i) = atand((node(ele(i, 3), 4) - node(ele(i, 2), 4)) ...
         / (node(ele(i, 3), 2) - node(ele(i, 2), 2)));
end

% obtain the global stiffness matrix of the rod
dof = length(node(:, 1)) * 3;      % degree of freedom
F = zeros(dof, 1) * nan;           % global force vector
D = ones(dof, 1) * nan;            % deformation vector in global form
K = zeros(dof);                % global stiffness matrix

for i = 1 : n_ele
    k_ele = transform(E, I, A, L(i), ang(i));
    K = combine(K, k_ele, ele(i, 2), ele(i, 3));
end

% Force boundary conditions, units N, mm.
F5x = 0; F5z = -20; T5 = 0;
F4x = 0; F4z = 0; T4 = 0;
F3x = 0; F3z = 0; T3 = 0;
F2x = 0; F2z = 0; T2 = 0;
F(4 : end) = [F2x F2z T2 F3x F3z T3 F4x F4z T4 F5x F5z T5];

% Displacement boundary conditions, units mm, degree.
D1x = 0; D1z = 0; M1 = 0;
D(1 : 3) = [D1x D1z M1];

% Solve the equation
index = D ~= 0;
D(index) = K(index, index) \ F(index);
F = K * D;

% print the nodal deformation
disp('The nodal deformations are:');
disp(D);

% plot the deformed and undefromed fishing rod
A = 0.001 * node(:, 2);
B = 0.001 * node(:, 4);
A1 = zeros(1, length(A));
B1 = zeros(1, length(B));
for i = 1 : length(D(:, 1))/3
    A1(i) = A(i) + 0.001 * D(3 * i - 2, 1);
    B1(i) = B(i) + 0.001 * D(3 * i - 1, 1);
end
plot(A, B, '-k', A1, B1, '-r');
axis([0, 1.6, 0, 1.2]);
set(gca, 'XTick', [0 : 0.2 : 1.6]);
set(gca, 'YTick', [0 : 0.2 : 1.2]);

function K_ele = transform(E, I, A, L, ang)
% This is to transform local stiffness matrix into global form.
%   E is Young's modulus,
%   I is area moment of inertia,
%   A is cross-section area,
%   L is length,
%   ang is horizontal angle.
%   the size of the beam stiffness matrix is 6*6.

% local stiffness matrix
k1 = A * E / L;
k2 = E * I / L^3;
k_local = [k1, 0, 0, -k1, 0, 0;
           0, 12 * k2, 6 * k2 * L, 0, -12 * k2, 6 * k2 * L;
           0, 6 * k2 * L, 4 * k2 * L^2, 0, -6 * k2 * L, 2 * k2 * L^2;
           -k1, 0, 0, k1, 0, 0;
           0, -12 * k2, -6 * k2 * L, 0, 12 * k2, -6 * k2 * L;
           0, 6 * k2 * L, 2 * k2 * L^2, 0, -6 * k2 * L, 4 * k2 * L^2];

% transformation matrix
c = cosd(ang);
s = sind(ang);
T = [c s 0 0 0 0;
     -s c 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 c s 0;
     0 0 0 -s c 0;
     0 0 0 0 0 1];

% transform local stiffness matrix to global form
K_ele = (T') * k_local * T;

end

function k_t = combine(k_t, k_ele, n1, n2)
% This is to obtain the global stiffness matrix.
%   stiffness matrix of a beam with node i and j are used as input.

d = zeros(1, 6);
d(1) = 3 * n1 - 2;
d(2) = 3 * n1 - 1;
d(3) = 3 * n1;
d(4) = 3 * n2 - 2;
d(5) = 3 * n2 - 1;
d(6) = 3 * n2;
for m = 1 : 6
    for n = 1 : 6
        k_t(d(m), d(n)) = k_t(d(m), d(n)) + k_ele(m, n);
    end
end

end
