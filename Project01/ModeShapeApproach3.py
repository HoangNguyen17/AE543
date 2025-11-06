import numpy as np
import matplotlib.pyplot as plt
# Note: scipy.linalg.eigh is no longer imported

# --- 1. Define Parameters ---
L = 5.0       # m, total length
b_o = 0.70    # m, outer chord
h_o = 0.25    # m, outer height
t = 0.010     # m, wall thickness
E = 69e9      # Pa (N/m^2), Young's modulus
rho = 2700    # kg/m^3, density
m_p = 100.0   # kg, tip payload
n = 5         # number of segments

print(f"--- 1. Parameters Defined ---")
print(f"L={L}m, E={E/1e9} GPa, rho={rho} kg/m^3, m_p={m_p} kg, n={n} segments")

# --- 2. Calculate Beam Properties ---
print("\n--- 2. Calculating Beam Properties ---")

# Inner dimensions
b_i = b_o - 2 * t
h_i = h_o - 2 * t

# Area Moment of Inertia (I) for bending about the chord-wise axis
I = (1/12) * (b_o * h_o**3 - b_i * h_i**3)
EI = E * I

# Cross-sectional Area (A) and Total Wing Mass (m_w)
A = (b_o * h_o) - (b_i * h_i)
m_w = rho * A * L

# Segment properties
delta_L = L / n
m_seg = m_w / n

print(f"Area Moment of Inertia (I): {I:.4e} m^4")
print(f"Flexural Rigidity (EI): {EI:.4e} N*m^2")
print(f"Total Wing Mass (m_w): {m_w:.2f} kg")
print(f"Segment Mass (m_seg): {m_seg:.2f} kg")

# --- 3. Build Stiffness Matrix (K) ---
print("\n--- 3. Building Stiffness Matrix [K] ---")

# We build the Flexibility Matrix (C) first, then K = C^-1
C = np.zeros((n, n))
# Node positions (x=1, 2, 3, 4, 5)
x_nodes = np.arange(1, n + 1) * delta_L

# Use cantilever beam deflection formulas C_ij = deflection at i for force at j
for i in range(n):
    for j in range(n):
        xi = x_nodes[i]
        xj = x_nodes[j]
        
        if xi <= xj:
            # Deflection at x=xi due to load P=1 at a=xj
            C[i, j] = (xi**2 * (3*xj - xi)) / (6 * EI)
        else: # xi > xj
            # Deflection at x=xi due to load P=1 at a=xj
            # (By reciprocity, C[i,j] = C[j,i], this is equivalent)
            C[i, j] = (xj**2 * (3*xi - xj)) / (6 * EI)

# Stiffness Matrix K = C^-1
K = np.linalg.inv(C)

print("Stiffness Matrix K (N/m) (first 5x5 shown):")
print(np.array2string(K, formatter={'float_kind':lambda x: "%.2e" % x}))

# --- 4. Build Mass Matrices (M) ---
print("\n--- 4. Building Mass Matrices [M] ---")

# Case 1: Without Payload (M1)
m_lumped_1 = np.full(n, m_seg)
M1 = np.diag(m_lumped_1)
print("Mass Matrix M1 (Without Payload) (kg):")
print(np.array2string(M1, formatter={'float_kind':lambda x: "%.2f" % x}))

# Case 2: With Payload (M2)
m_lumped_2 = np.full(n, m_seg)
m_lumped_2[-1] += m_p  # Add payload to the last node
M2 = np.diag(m_lumped_2)
print("\nMass Matrix M2 (With Payload) (kg):")
print(np.array2string(M2, formatter={'float_kind':lambda x: "%.2f" % x}))

# --- 5. Solve Eigenvalue Problem ---
print("\n--- 5. Solving Eigenvalue Problems ---")

def solve_and_print(K, M, case_name):
    """Solves the generalized eigenvalue problem (K*phi = lambda*M*phi) 
       using only numpy and prints formatted results."""
    
    # Convert generalized problem K*phi = w^2*M*phi
    # to standard problem K_tilde*u = w^2*u
    
    # 1. Get M^(-1/2)
    # Since M is diagonal, M^(-1/2) is diag(1/sqrt(m_i))
    M_inv_sqrt = np.diag(1.0 / np.sqrt(np.diag(M)))
    
    # 2. Form the modified stiffness matrix K_tilde
    # K_tilde = M^(-1/2) * K * M^(-1/2)
    K_tilde = M_inv_sqrt @ K @ M_inv_sqrt
    
    # 3. Solve the standard eigenvalue problem using numpy.linalg.eigh
    # K_tilde is symmetric, so we can use eigh
    # eigvals are w^2
    # eigvecs_u are the eigenvectors in the {u} coordinate system
    eigvals, eigvecs_u = np.linalg.eigh(K_tilde)
    
    # 4. Sort eigenvalues (and corresponding eigenvectors) from low to high
    idx = eigvals.argsort()
    eigvals = eigvals[idx]
    eigvecs_u_sorted = eigvecs_u[:, idx]
    
    # 5. Transform eigenvectors {u} back to physical {phi}
    # phi = M^(-1/2) * u
    eigvecs = M_inv_sqrt @ eigvecs_u_sorted
    
    # 6. Calculate frequencies
    freqs_rad = np.sqrt(eigvals)      # omega (rad/s)
    freqs_hz = freqs_rad / (2 * np.pi) # f (Hz)
    
    print(f"\n--- Results: {case_name} ---")
    print("Natural Frequencies (Hz):")
    print(np.array2string(freqs_hz, formatter={'float_kind':lambda x: "%.2f" % x}))
    
    print("\nMode Shapes (Eigenvectors, normalized to max component):")
    modes = []
    vec_formatter = {'float_kind':lambda x: "% .3f" % x} # Formatter for vectors

    for i in range(n):
        vec = eigvecs[:, i]
        # Normalize by dividing by the component with the largest absolute value
        normalized_vec = vec / np.max(np.abs(vec))
        modes.append(normalized_vec)
        print(f"  Mode {i+1} (f={freqs_hz[i]:.2f} Hz):")
        print("    " + np.array2string(normalized_vec, formatter=vec_formatter))
        
    return freqs_hz, modes

# Solve for both cases
freqs1, modes1 = solve_and_print(K, M1, "Without Payload (mp=0)")
freqs2, modes2 = solve_and_print(K, M2, "With Payload (mp=100kg)")

# --- 6. Plot Mode Shapes ---
print("\n--- 6. Generating Plots ---")

def plot_modes(x_nodes, modes, freqs, title, filename):
    """Plots the first 5 mode shapes."""
    plt.figure(figsize=(10, 8))
    
    # Add the root (x=0, y=0) for plotting
    x_plot = np.insert(x_nodes, 0, 0) 
    
    for i in range(n):
        # Add the root displacement (0) to the mode shape vector
        mode_shape_plot = np.insert(modes[i], 0, 0)
        plt.plot(x_plot, mode_shape_plot, 'o-', 
                 label=f'Mode {i+1} (f={freqs[i]:.2f} Hz)')
        
    plt.title(title, fontsize=16)
    plt.xlabel('Wing Span (m)', fontsize=12)
    plt.ylabel('Normalized Deflection', fontsize=12)
    plt.legend(loc='best')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.axhline(0, color='black', linewidth=0.5) # Add x-axis
    plt.savefig(filename)
    print(f"Saved plot to {filename}")

# Plot for Case 1
plot_modes(x_nodes, modes1, freqs1, 
           'Mode Shapes - Without Payload (mp=0)', 
           'mode_shapes_no_payload.png')

# Plot for Case 2
plot_modes(x_nodes, modes2, freqs2, 
           'Mode Shapes - With Payload (mp=100kg)', 
           'mode_shapes_with_payload.png')

print("\n--- Analysis Complete ---")