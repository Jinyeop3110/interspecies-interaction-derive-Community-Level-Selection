"""
Generalized Lotka-Volterra (gLV) Simulation for Community Coalescence

This script implements a gLV model to simulate microbial community dynamics
and coalescence events between communities.

Reference: Interspecies interactions derive community-level selection
"""

import numpy as np
from scipy.integrate import solve_ivp


def gLV_dynamics(t, y, interaction_matrix, growth_rates, carrying_capacities):
    """
    Generalized Lotka-Volterra equations.

    dN_i/dt = r_i * N_i * (1 - sum_j(alpha_ij * N_j) / K_i)

    Parameters
    ----------
    t : float
        Time (not used, for ODE solver compatibility)
    y : array
        Species abundances
    interaction_matrix : 2D array
        Interaction coefficients alpha_ij
    growth_rates : array
        Intrinsic growth rates r_i
    carrying_capacities : array
        Carrying capacities K_i

    Returns
    -------
    dydt : array
        Rate of change for each species
    """
    y = np.maximum(y, 0)  # Ensure non-negative abundances
    dydt = growth_rates * y * (1 - interaction_matrix @ y / carrying_capacities)
    return dydt


def simulate_community(y0, t_span, interaction_matrix, growth_rates, carrying_capacities,
                       t_eval=None, extinction_threshold=1e-6):
    """
    Simulate community dynamics using gLV model.

    Parameters
    ----------
    y0 : array
        Initial abundances
    t_span : tuple
        (t_start, t_end)
    interaction_matrix : 2D array
        Interaction coefficients
    growth_rates : array
        Growth rates
    carrying_capacities : array
        Carrying capacities
    t_eval : array, optional
        Time points for output
    extinction_threshold : float
        Below this abundance, species considered extinct

    Returns
    -------
    result : OdeResult
        Solution object with .t (times) and .y (abundances)
    """
    def f(t, y):
        return gLV_dynamics(t, y, interaction_matrix, growth_rates, carrying_capacities)

    result = solve_ivp(f, t_span, y0, method='RK45', t_eval=t_eval,
                       max_step=0.1, dense_output=True)

    # Apply extinction threshold
    result.y[result.y < extinction_threshold] = 0

    return result


def coalescence_event(community_A, community_B, mixing_ratio=0.5):
    """
    Mix two communities in a coalescence event.

    Parameters
    ----------
    community_A : array
        Abundances of community A
    community_B : array
        Abundances of community B
    mixing_ratio : float
        Fraction of community A (0 to 1)

    Returns
    -------
    mixed : array
        Initial condition for coalesced community
    """
    return mixing_ratio * community_A + (1 - mixing_ratio) * community_B


def generate_random_interactions(n_species, interaction_strength=0.5, symmetric=False):
    """
    Generate random interaction matrix.

    Parameters
    ----------
    n_species : int
        Number of species
    interaction_strength : float
        Mean interaction strength (0 to 1)
    symmetric : bool
        If True, make interactions symmetric

    Returns
    -------
    interactions : 2D array
        Interaction matrix with 1s on diagonal
    """
    # Random interactions scaled by strength
    interactions = np.random.uniform(0, 2 * interaction_strength, (n_species, n_species))

    if symmetric:
        interactions = (interactions + interactions.T) / 2

    # Self-interaction = 1 (normalized carrying capacity)
    np.fill_diagonal(interactions, 1.0)

    return interactions


# =============================================================================
# Example usage
# =============================================================================

if __name__ == "__main__":

    # Parameters
    n_species = 10
    interaction_strength = 0.3
    t_assembly = 500    # Time for initial assembly
    t_coalescence = 500 # Time after coalescence

    print(f"Simulating {n_species} species community coalescence")
    print(f"Interaction strength: {interaction_strength}")
    print("-" * 50)

    # Generate shared species pool parameters
    np.random.seed(42)
    interaction_matrix = generate_random_interactions(n_species, interaction_strength)
    growth_rates = np.ones(n_species)
    carrying_capacities = np.ones(n_species)

    # Initial conditions: two different communities
    y0_A = np.random.uniform(0.01, 0.1, n_species)
    y0_B = np.random.uniform(0.01, 0.1, n_species)

    # Assemble communities separately
    print("Assembling Community A...")
    result_A = simulate_community(y0_A, (0, t_assembly), interaction_matrix,
                                   growth_rates, carrying_capacities)
    community_A = result_A.y[:, -1]

    print("Assembling Community B...")
    result_B = simulate_community(y0_B, (0, t_assembly), interaction_matrix,
                                   growth_rates, carrying_capacities)
    community_B = result_B.y[:, -1]

    # Coalescence event
    print("\nCoalescence event (50:50 mixing)...")
    y0_mixed = coalescence_event(community_A, community_B, mixing_ratio=0.5)

    result_coalesced = simulate_community(y0_mixed, (0, t_coalescence), interaction_matrix,
                                          growth_rates, carrying_capacities)
    community_final = result_coalesced.y[:, -1]

    # Results summary
    print("\n" + "=" * 50)
    print("RESULTS")
    print("=" * 50)

    richness_A = np.sum(community_A > 1e-6)
    richness_B = np.sum(community_B > 1e-6)
    richness_final = np.sum(community_final > 1e-6)

    print(f"Community A richness: {richness_A}")
    print(f"Community B richness: {richness_B}")
    print(f"Coalesced richness:   {richness_final}")

    # Species retention analysis
    present_A = community_A > 1e-6
    present_B = community_B > 1e-6
    present_final = community_final > 1e-6

    retained_from_A = np.sum(present_A & present_final)
    retained_from_B = np.sum(present_B & present_final)

    print(f"\nSpecies retained from A: {retained_from_A}/{richness_A}")
    print(f"Species retained from B: {retained_from_B}/{richness_B}")
