# SSH_Model_PBC_OBC
This MATLAB code simulates the Su-Schrieffer-Heeger (SSH) model, a fundamental topological model in condensed matter physics. 
Below is a detailed explanation of each section:

1. Parameter Setup

t1: Hopping amplitude within a unit cell (weaker bond)
t2: Hopping amplitude between unit cells (stronger bond)
N: Size of finite chain for numerical calculations
k_points: Resolution for momentum space calculations

2. Infinite Chain Band Structure (Periodic Boundary Conditions)
matlab
k = linspace(-pi, pi, k_points);
E_bands=zeros(2, length(k));

for i = 1:length(k)
    H_k = [0, t1 + t2*exp(-1j*k(i)); 
           t1 + t2*exp(1j*k(i)), 0];
    eigenvalues = eig(H_k);
    E_bands(:, i) = sort(real(eigenvalues));
end
Calculates the band structure for an infinite chain by diagonalizing the Hamiltonian in momentum space

The Hamiltonian in k-space is a 2×2 matrix representing the two-site unit cell

Results in two energy bands (valence and conduction bands)

3. Finite Chain with Periodic Boundary Conditions (PBC)
matlab
H_pbc = zeros(2*N);
for i = 1:N
    % Intra-cell hopping
    H_pbc(2*i-1, 2*i) = t1;
    H_pbc(2*i, 2*i-1) = t1;
    
    % Inter-cell hopping
    if i < N
        H_pbc(2*i, 2*i+1) = t2;
        H_pbc(2*i+1, 2*i) = t2;
    end
    
    % Connect last cell to first
    if i == N
        H_pbc(2*N, 1) = t2;
        H_pbc(1, 2*N) = t2;
    end
end
Constructs a finite Hamiltonian with periodic boundary conditions

The chain forms a ring where the last cell connects back to the first

Results in discrete energy levels that approximate the infinite chain

4. Finite Chain with Open Boundary Conditions (OBC)
matlab
H_obc = zeros(2*N);
for i = 1:N
    % Intra-cell hopping
    H_obc(2*i-1, 2*i) = t1;
    H_obc(2*i, 2*i-1) = t1;
    
    % Inter-cell hopping (no connection at ends)
    if i < N
        H_obc(2*i, 2*i+1) = t2;
        H_obc(2*i+1, 2*i) = t2;
    end
end
Constructs a finite Hamiltonian with open ends

The chain has true endpoints with no connection between first and last sites

May show edge states in the topological phase

5. Topological Invariant Calculation (Winding Number)
matlab
d_vector = zeros(2, length(k));
for i = 1:length(k)
    dx = t1 + t2*cos(k(i));
    dy = t2*sin(k(i));
    d_vector(:, i) = [dx; dy];
end

winding_number = 0;
for i = 1:length(k)-1
    theta1 = atan2(d_vector(2,i), d_vector(1,i));
    theta2 = atan2(d_vector(2,i+1), d_vector(1,i+1));
    delta_theta = theta2 - theta1;
    
    % Handle phase jumps
    if delta_theta > pi
        delta_theta = delta_theta - 2*pi;
    elseif delta_theta < -pi
        delta_theta = delta_theta + 2*pi;
    end
    
    winding_number = winding_number + delta_theta;
end
winding_number = winding_number / (2*pi);
Calculates the winding number by tracking the d-vector trajectory

The d-vector components come from rewriting the Hamiltonian as H(k) = d(k)·σ

The winding number counts how many times d(k) wraps around the origin

Integer values (0 or 1) indicate distinct topological phases

Key Physics Concepts Illustrated:
Band Structure: Shows the energy-momentum relationship in periodic systems

Bulk-Boundary Correspondence: Relates bulk topology (winding number) to edge states

Topological Invariants: The winding number remains quantized under continuous deformations

Boundary Condition Effects: PBC vs OBC show different spectral features

The code demonstrates how to:

Construct tight-binding Hamiltonians

Diagonalize systems with different boundary conditions

Calculate topological invariants

Visualize results for physical interpretation

The SSH model serves as a paradigmatic example of a topological insulator in 1D, and this code provides a complete numerical implementation of its key features.


Here's an explanation of the key differences between the Open Boundary Condition (OBC) implementation in this code and what would be done with Periodic Boundary Conditions (PBC), along with their physical significance:

1. Hamiltonian Construction Differences
OBC Implementation (Current Code):

matlab
% Inter-cell hopping (t2)
if i < N  % Only connect interior cells
    H_finite(2*i, 2*i+1) = t2;
    H_finite(2*i+1, 2*i) = t2;
end
PBC Implementation (Comparison):

matlab
% Inter-cell hopping (t2) - regular terms
if i < N
    H_pbc(2*i, 2*i+1) = t2;
    H_pbc(2*i+1, 2*i) = t2;
end
% Additional PBC connection:
if i == N  % Connect last cell back to first
    H_pbc(2*N, 1) = t2;
    H_pbc(1, 2*N) = t2;
end
Physical Significance:

OBC creates a finite chain with true endpoints where electrons cannot hop past the boundary

PBC forms a closed loop where the chain has no endpoints, mimicking infinite systems

This difference crucially affects edge state formation in topological phases

2. Spectral Differences
OBC Features (Shown in Code):

matlab
plot(E, 'o', 'MarkerSize', 8);  % Will show:
% - Possible zero-energy edge states when |t1| < |t2|
% - Bulk states that don't form continuous bands
PBC Features (For Comparison):

matlab
plot(E_pbc, 'o', 'MarkerSize', 8);  % Would show:
% - No zero-energy states (gap remains open)
% - States approximate the infinite band structure
Physical Significance:

OBC spectrum reveals the bulk-boundary correspondence through edge states

PBC spectrum better approximates the infinite system's bulk properties

The winding number calculation (done earlier in the code) predicts OBC edge states

3. Edge State Visualization (Unique to OBC)
Current OBC Code:

matlab
if abs(t1) < abs(t2)  % Topological phase check
    plot(abs(psi(:,1)).^2);  % Shows localized edge state
end
PBC Would Have:

No equivalent section because PBC has no edges

All eigenstates would be extended Bloch waves

Physical Significance:

The edge state plots demonstrate:

Spatial localization at chain ends (exponential decay)

Topological protection (exists despite disorder)

These are the hallmark features of topological insulators

4. Mathematical Formalism Differences
OBC Treatment:

Uses real-space Hamiltonian with missing connections

Diagonalization yields both bulk and edge solutions

Requires explicit finite-size calculations

PBC Treatment:

Can use momentum-space representation (k-space)

Solutions are Bloch waves with well-defined momentum

Bulk properties can be calculated analytically

Key Physical Implications:
Boundary Effects Matter:

OBC reveals phenomena invisible in PBC (edge states)

Real experimental systems are always OBC-like

Topological Detection:

Winding number (calculated from infinite/PBC system) predicts OBC edge states

This bulk-boundary correspondence is fundamental to topology

Size Dependence:

OBC results are more sensitive to system size N

Finite-size gaps appear that would close in PBC

The code intentionally uses OBC to demonstrate these important physical differences that arise when considering real systems with boundaries versus idealized periodic systems.
