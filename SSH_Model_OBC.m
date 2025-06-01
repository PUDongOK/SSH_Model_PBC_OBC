clear all;
close all;
clc;

%% SSH Model Parameters Setup
t1=0.5;    % Intra-cell hopping strength
t2=1;      % Inter-cell hopping strength
N=40;      % Number of unit cells (for finite chain calculation)
k_points=100; % Number of k-points

%% Calculate Band Structure for Infinite SSH Chain
k=linspace(-pi, pi, k_points);
E_bands=zeros(2, length(k));

for i=1:length(k)
    H_k=[0, t1 + t2*exp(-1j*k(i)); 
           t1 + t2*exp(1j*k(i)), 0];
    eigenvalues=eig(H_k);
    E_bands(:, i)=sort(real(eigenvalues));
end

% Plot Band Structure
figure;
plot(k, E_bands(1,:), 'b', 'LineWidth', 1.5);
hold on;
plot(k, E_bands(2,:), 'r', 'LineWidth', 1.5);
xlabel('Wave vector k');
ylabel('Energy');
title('SSH Model Band Structure');
grid on;
legend('Lower band', 'Upper band');
set(gca, 'FontSize', 12);

%% Calculate Winding Number
d_vector=zeros(2, length(k));
for i=1:length(k)
    dx=t1 + t2*cos(k(i));
    dy=t2*sin(k(i));
    d_vector(:, i)=[dx; dy];
end

% Plot d-vector Trajectory in Brillouin Zone
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

% Calculate Winding Number
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

%% Finite SSH Chain Calculation (Open Boundary Conditions)
H_finite=zeros(2*N);

% Construct Hamiltonian
for i=1:N
    % Intra-cell hopping (t1)
    if mod(i,2) == 1  % A sublattice to B sublattice
        H_finite(2*i-1, 2*i)=t1;
        H_finite(2*i, 2*i-1)=t1;
    end
    
    % Inter-cell hopping (t2)
    if i < N
        H_finite(2*i, 2*i+1)=t2;
        H_finite(2*i+1, 2*i)=t2;
    end
end

% Diagonalization
[psi, E]=eig(H_finite);
E=diag(E);

% Plot Finite Chain Energy Spectrum
figure;
plot(E, 'o', 'MarkerSize', 8);
xlabel('State index');
ylabel('Energy');
title('Energy spectrum of finite SSH chain');
grid on;
set(gca, 'FontSize', 12);

% Plot Edge State Wavefunctions (if they exist)
if abs(t1) < abs(t2)  % Topologically non-trivial phase
    figure;
    subplot(2,1,1);
    plot(abs(psi(:,1)).^2, 'LineWidth', 1.5);
    title('Left edge state');
    xlabel('Site index');
    ylabel('Probability density');
    grid on;
    
    subplot(2,1,2);
    plot(abs(psi(:,end)).^2, 'LineWidth', 1.5);
    title('Right edge state');
    xlabel('Site index');
    ylabel('Probability density');
    grid on;
    set(gca, 'FontSize', 12);
end