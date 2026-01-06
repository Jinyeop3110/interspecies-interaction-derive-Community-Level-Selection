# Interspecies Interactions Derive Community-Level Selection

Code associated with the manuscript on community-level selection driven by interspecies interactions in microbial communities.

## Overview

This repository contains Python scripts for simulating microbial community coalescence using generalized Lotka-Volterra (gLV) models.

## Repository Structure

```
├── gLV_simulation.py    # Core gLV simulation and coalescence functions
├── requirements.txt     # Python dependencies
└── README.md
```

## Requirements

```bash
pip install -r requirements.txt
```

## Quick Start

```python
from gLV_simulation import (
    simulate_community,
    coalescence_event,
    generate_random_interactions
)
import numpy as np

# Setup
n_species = 10
interaction_matrix = generate_random_interactions(n_species, interaction_strength=0.3)
growth_rates = np.ones(n_species)
carrying_capacities = np.ones(n_species)

# Simulate two communities
y0_A = np.random.uniform(0.01, 0.1, n_species)
result_A = simulate_community(y0_A, (0, 500), interaction_matrix, growth_rates, carrying_capacities)

# Coalescence
y0_mixed = coalescence_event(result_A.y[:, -1], result_B.y[:, -1], mixing_ratio=0.5)
result_coalesced = simulate_community(y0_mixed, (0, 500), interaction_matrix, growth_rates, carrying_capacities)
```

## Run Example

```bash
python gLV_simulation.py
```

## Model Description

The generalized Lotka-Volterra model:

$$\frac{dN_i}{dt} = r_i N_i \left(1 - \frac{\sum_j \alpha_{ij} N_j}{K_i}\right)$$

Where:
- $N_i$ = abundance of species $i$
- $r_i$ = intrinsic growth rate
- $\alpha_{ij}$ = interaction coefficient (effect of species $j$ on species $i$)
- $K_i$ = carrying capacity

## Data Availability

Raw sequencing data: [TBD - Dryad DOI]

## Citation

[TBD]

## License

MIT License
