clear all;
close all;
clc;

%% SSH model parameters setup
t1=0.5;    % Intra-cell hopping strength
t2=1.0;    % Inter-cell hopping strength
N=40;      % Number of unit cells (for finite chain calculation)
k_points=100; % Number of k-points

%% 1. Calculate band structure for infinite SSH chain (implicit periodic boundary conditions)
k=linspace(-pi, pi, k_points);
E_bands=zeros(2, length(k));

for i=1:length(k)
    H_k=[0, t1 + t2*exp(-1j*k(i)); 
           t1 + t2*exp(1j*k(i)), 0];
    eigenvalues=eig(H_k);
    E_bands(:, i)=sort(real(eigenvalues));
end

% Plot band structure
figure;
plot(k, E_bands(1,:), 'b', 'LineWidth', 1.5);
hold on;
plot(k, E_bands(2,:), 'r', 'LineWidth', 1.5);
xlabel('Wave vector k');
ylabel('Energy');
title('SSH Model Band Structure (PBC)');
grid on;
legend('Lower band', 'Upper band');
set(gca, 'FontSize', 12);

%% 2. Finite SSH chain calculation - Periodic Boundary Conditions (PBC)
H_pbc=zeros(2*N);

% Construct Hamiltonian with periodic boundary conditions
for i=1:N
    % Intra-cell hopping (t1)
    H_pbc(2*i-1, 2*i)=t1;
    H_pbc(2*i, 2*i-1)=t1;
    
    % Inter-cell hopping (t2) - regular terms
    if i < N
        H_pbc(2*i, 2*i+1)=t2;
        H_pbc(2*i+1, 2*i)=t2;
    end
    
    % Periodic boundary terms - connect last cell to first cell
    if i == N
        H_pbc(2*N, 1)=t2;
        H_pbc(1, 2*N)=t2;
    end
end

% Diagonalization
[psi_pbc, E_pbc]=eig(H_pbc);
E_pbc=diag(E_pbc);

% Plot finite chain PBC energy spectrum
figure;
plot(E_pbc, 'o', 'MarkerSize', 8);
xlabel('State index');
ylabel('Energy');
title(['Energy spectrum of finite SSH chain (PBC), N=' num2str(N)]);
grid on;
set(gca, 'FontSize', 12);

%% 3. Finite SSH chain calculation - Open Boundary Conditions (OBC) comparison
H_obc=zeros(2*N);

% Construct Hamiltonian with open boundary conditions
for i=1:N
    % Intra-cell hopping (t1)
    H_obc(2*i-1, 2*i)=t1;
    H_obc(2*i, 2*i-1)=t1;
    
    % Inter-cell hopping (t2) - internal connections only
    if i < N
        H_obc(2*i, 2*i+1)=t2;
        H_obc(2*i+1, 2*i)=t2;
    end
end

% Diagonalization
[psi_obc, E_obc]=eig(H_obc);
E_obc=diag(E_obc);

% Plot finite chain OBC energy spectrum
figure;
plot(E_obc, 'o', 'MarkerSize', 8);
xlabel('State index');
ylabel('Energy');
title(['Energy spectrum of finite SSH chain (OBC), N=' num2str(N)]);
grid on;
set(gca, 'FontSize', 12);

%% 4. Compare PBC and OBC energy spectra
figure;
plot(sort(E_pbc), 'ro', 'MarkerSize', 8, 'DisplayName', 'PBC');
hold on;
plot(sort(E_obc), 'bo', 'MarkerSize', 8, 'DisplayName', 'OBC');
xlabel('State index');
ylabel('Energy');
title(['Comparison of PBC and OBC spectra, N=' num2str(N)]);
legend;
grid on;
set(gca, 'FontSize', 12);

%% 5. Calculate topological invariant - Winding number
d_vector=zeros(2, length(k));
for i=1:length(k)
    dx=t1 + t2*cos(k(i));
    dy=t2*sin(k(i));
    d_vector(:, i)=[dx; dy];
end

% Plot d-vector trajectory in Brillouin zone
figure;
plot(d_vector(1,:), d_vector(2,:), 'LineWidth', 1.5);
hold on;
plot(0, 0, 'ro'); % Origin
xlabel('d_x');
ylabel('d_y');
title('Trajectory of d(k) vector');
grid on;
axis equal;
set(gca, 'FontSize', 12);

% Calculate winding number
winding_number=0;
for i=1:length(k)-1
    theta1=atan2(d_vector(2,i), d_vector(1,i));
    theta2=atan2(d_vector(2,i+1), d_vector(1,i+1));
    delta_theta=theta2 - theta1;
    
    % Handle phase jumps
    if delta_theta > pi
        delta_theta=delta_theta - 2*pi;
    elseif delta_theta < -pi
        delta_theta=delta_theta + 2*pi;
    end
    
    winding_number=winding_number + delta_theta;
end
winding_number=winding_number / (2*pi);

fprintf('Winding number: %.2f\n', winding_number);

% Modified part: using if-else instead of ternary operator
if abs(t1) < abs(t2)
    phase='Non-trivial (winding=1)';
else
    phase='Trivial (winding=0)';
end
fprintf('Topological phase: %s\n', phase);