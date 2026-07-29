"""
COLORMAP.py - Centralized Color Definitions for Coalescence Analysis

This module provides consistent color schemes across all analysis scripts.
Colors are standardized for nutrient conditions and phase diagram classifications.

NUTRIENT CONDITION COLORS:
--------------------------
- LN (Low Nitrogen): Rich magenta-purple (#A7216A)
- MN (Medium Nitrogen): Dark brown-red (#802000)  
- HN (High Nitrogen): Deep orange-red (#E24912)

PHASE DIAGRAM COLORS:
--------------------
- Dominance: Red (#D32F2F) - One community dominates
- Mixing: Green (#388E3C) - Communities mix/balance
- Restructuring: Purple (#7B1FA2) - New structure emerges

USAGE:
------
# Nutrient condition colors
from COLORMAP import MEDIUM_COLORS, get_medium_color, get_medium_colors

color = get_medium_color('L')  # or 'LN', 'M', 'MN', 'H', 'HN'
colors = get_medium_colors()  # Returns [LN_color, MN_color, HN_color]
ln_color = MEDIUM_COLORS['LN']

# Phase diagram colors  
from COLORMAP import PHASE_DIAGRAM_COLORS, get_phase_diagram_color, get_phase_diagram_colors

color = get_phase_diagram_color('dominance')
colors = get_phase_diagram_colors()  # Returns [dominance, mixing, restructuring]
dom_color = PHASE_DIAGRAM_COLORS['dominance']

All colors provide good contrast and are colorblind-friendly.
"""

# Primary color scheme for nutrient conditions
MEDIUM_COLORS = {
    'LN': '#A7216A',   # Low Nitrogen - Rich magenta-purple
    'MN': '#802000',   # Medium Nitrogen - Dark brown-red
    'HN': '#E24912',   # High Nitrogen - Deep orange-red
    'L': '#A7216A',    # Alias for LN
    'M': '#802000',    # Alias for MN  
    'H': '#E24912',    # Alias for HN
}

# Color list in LN, MN, HN order (for legacy compatibility)
MEDIUM_COLOR_LIST = [
    MEDIUM_COLORS['LN'],  # Index 0: Low Nitrogen
    MEDIUM_COLORS['MN'],  # Index 1: Medium Nitrogen  
    MEDIUM_COLORS['HN'],  # Index 2: High Nitrogen
]

# Phase diagram colors for coalescence outcome classification (lighter/fainter version)
PHASE_DIAGRAM_COLORS = {
    'dominance': '#E57373',      # Light red - One community dominates (was #D32F2F)
    'mixing': '#81C784',         # Light green - Communities mix/balance (was #388E3C)
    'restructuring': '#BA68C8',  # Light purple - New structure emerges (was #7B1FA2)
    'Dominance': '#E57373',      # Alias with capitalization
    'Mixing': '#81C784',         # Alias with capitalization  
    'Restructuring': '#BA68C8',  # Alias with capitalization
}

# Phase diagram color list in order [Dominance, Mixing, Restructuring]
PHASE_DIAGRAM_COLOR_LIST = [
    PHASE_DIAGRAM_COLORS['dominance'],    # Index 0: Dominance
    PHASE_DIAGRAM_COLORS['mixing'],       # Index 1: Mixing
    PHASE_DIAGRAM_COLORS['restructuring'] # Index 2: Restructuring  
]

# Alternative color schemes for special cases
SECONDARY_COLORS = {
    'grey': '#808080',
    'light_grey': '#CCCCCC',
    'dark_grey': '#404040',
    'black': '#000000',
    'white': '#FFFFFF',
}

# Seaborn tab10 palette mapping (for backward compatibility)
TAB10_COLORS = {
    'tab10_0': '#1f77b4',  # blue
    'tab10_1': '#ff7f0e',  # orange  
    'tab10_2': '#2ca02c',  # green
    'tab10_3': '#d62728',  # red
    'tab10_4': '#9467bd',  # purple
    'tab10_5': '#8c564b',  # brown
    'tab10_6': '#e377c2',  # pink
    'tab10_7': '#7f7f7f',  # grey
    'tab10_8': '#bcbd22',  # olive
    'tab10_9': '#17becf',  # cyan
}

def get_medium_color(medium):
    """
    Get color for a specific medium condition.
    
    Parameters:
    -----------
    medium : str
        Medium condition ('L', 'LN', 'M', 'MN', 'H', 'HN')
        
    Returns:
    --------
    str : Hex color code
    
    Examples:
    ---------
    >>> get_medium_color('L')
    '#A7216A'
    >>> get_medium_color('MN') 
    '#802000'
    """
    medium_upper = medium.upper()
    if medium_upper in MEDIUM_COLORS:
        return MEDIUM_COLORS[medium_upper]
    else:
        raise ValueError(f"Unknown medium: {medium}. Use 'L'/'LN', 'M'/'MN', or 'H'/'HN'")

def get_medium_colors():
    """
    Get all medium colors as a list in LN, MN, HN order.
    
    Returns:
    --------
    list : List of hex color codes [LN, MN, HN]
    
    Examples:
    ---------  
    >>> get_medium_colors()
    ['#A7216A', '#802000', '#E24912']
    """
    return MEDIUM_COLOR_LIST.copy()

def get_medium_color_by_index(index):
    """
    Get medium color by index (0=LN, 1=MN, 2=HN).
    
    Parameters:
    -----------
    index : int
        Index (0, 1, or 2)
        
    Returns:
    --------
    str : Hex color code
    
    Examples:
    ---------
    >>> get_medium_color_by_index(0)  # LN
    '#A7216A'
    >>> get_medium_color_by_index(2)  # HN  
    '#E24912'
    """
    if 0 <= index < len(MEDIUM_COLOR_LIST):
        return MEDIUM_COLOR_LIST[index]
    else:
        raise ValueError(f"Index {index} out of range. Use 0 (LN), 1 (MN), or 2 (HN)")

def get_tab10_compatible_colors():
    """
    Get medium colors in tab10-compatible format for legacy scripts.
    
    Returns:
    --------
    list : Colors matching tab10 indices used in existing scripts
    """
    return [
        MEDIUM_COLORS['LN'],  # tab10[2] replacement
        MEDIUM_COLORS['MN'],  # tab10[0] replacement  
        MEDIUM_COLORS['HN'],  # tab10[3] replacement
    ]

# Utility function for color mapping
def map_medium_to_index(medium):
    """
    Map medium string to index for color array access.
    
    Parameters:
    -----------
    medium : str
        Medium condition ('L', 'LN', 'M', 'MN', 'H', 'HN')
        
    Returns:
    --------
    int : Index (0=LN, 1=MN, 2=HN)
    """
    medium_upper = medium.upper()
    mapping = {'L': 0, 'LN': 0, 'M': 1, 'MN': 1, 'H': 2, 'HN': 2}
    if medium_upper in mapping:
        return mapping[medium_upper]
    else:
        raise ValueError(f"Unknown medium: {medium}. Use 'L'/'LN', 'M'/'MN', or 'H'/'HN'")

# Phase diagram color functions
def get_phase_diagram_color(category):
    """
    Get color for a specific phase diagram category.
    
    Parameters:
    -----------
    category : str
        Phase diagram category ('dominance', 'mixing', 'restructuring')
        Case insensitive.
        
    Returns:
    --------
    str : Hex color code
    
    Examples:
    ---------
    >>> get_phase_diagram_color('dominance')
    '#D32F2F'
    >>> get_phase_diagram_color('Mixing')
    '#7B1FA2'
    """
    category_lower = category.lower()
    if category_lower in PHASE_DIAGRAM_COLORS:
        return PHASE_DIAGRAM_COLORS[category_lower]
    elif category in PHASE_DIAGRAM_COLORS:  # Check capitalized version
        return PHASE_DIAGRAM_COLORS[category]
    else:
        raise ValueError(f"Unknown phase diagram category: {category}. Use 'dominance', 'mixing', or 'restructuring'")

def get_phase_diagram_colors():
    """
    Get all phase diagram colors as a list in dominance, mixing, restructuring order.
    
    Returns:
    --------
    list : List of hex color codes [dominance, mixing, restructuring]
    
    Examples:
    ---------  
    >>> get_phase_diagram_colors()
    ['#D32F2F', '#7B1FA2', '#388E3C']
    """
    return PHASE_DIAGRAM_COLOR_LIST.copy()

def get_phase_diagram_color_by_index(index):
    """
    Get phase diagram color by index (0=dominance, 1=mixing, 2=restructuring).
    
    Parameters:
    -----------
    index : int
        Index (0, 1, or 2)
        
    Returns:
    --------
    str : Hex color code
    
    Examples:
    ---------
    >>> get_phase_diagram_color_by_index(0)  # Dominance
    '#D32F2F'
    >>> get_phase_diagram_color_by_index(2)  # Restructuring  
    '#388E3C'
    """
    if 0 <= index < len(PHASE_DIAGRAM_COLOR_LIST):
        return PHASE_DIAGRAM_COLOR_LIST[index]
    else:
        raise ValueError(f"Index {index} out of range. Use 0 (dominance), 1 (mixing), or 2 (restructuring)")

def map_phase_to_index(category):
    """
    Map phase diagram category string to index for color array access.
    
    Parameters:
    -----------
    category : str or int
        Phase diagram category ('dominance', 'mixing', 'restructuring') or class index (0, 1, 2)
        
    Returns:
    --------
    int : Index (0=dominance, 1=mixing, 2=restructuring)
    
    Examples:
    ---------
    >>> map_phase_to_index('mixing')
    1
    >>> map_phase_to_index(0)  # Pass-through for class indices
    0
    """
    # If already an integer (class index), return as-is
    if isinstance(category, int):
        if 0 <= category <= 2:
            return category
        else:
            raise ValueError(f"Class index {category} out of range. Use 0, 1, or 2")
    
    # If string, map to index
    category_lower = category.lower()
    mapping = {
        'dominance': 0, 
        'mixing': 1, 
        'restructuring': 2,
        'dominant': 0,  # Alternative naming
        'mix': 1,       # Alternative naming
        'restructure': 2, # Alternative naming
    }
    
    if category_lower in mapping:
        return mapping[category_lower]
    else:
        raise ValueError(f"Unknown phase diagram category: {category}. Use 'dominance', 'mixing', or 'restructuring'")

# For debugging and verification
if __name__ == "__main__":
    print("COLORMAP.py - Color Definitions")
    print("=" * 40)
    
    print("Medium Colors:")
    for medium, color in MEDIUM_COLORS.items():
        if len(medium) <= 2:  # Only show main keys, not aliases
            print(f"  {medium}: {color}")
    
    print(f"\nMedium Color List: {MEDIUM_COLOR_LIST}")
    
    print("\nPhase Diagram Colors:")
    for phase, color in PHASE_DIAGRAM_COLORS.items():
        if phase.islower():  # Only show lowercase keys, not aliases
            print(f"  {phase}: {color}")
    
    print(f"\nPhase Diagram Color List: {PHASE_DIAGRAM_COLOR_LIST}")
    
    print("\nTesting medium functions:")
    print(f"get_medium_color('L'): {get_medium_color('L')}")
    print(f"get_medium_colors(): {get_medium_colors()}")
    print(f"get_medium_color_by_index(1): {get_medium_color_by_index(1)}")
    print(f"map_medium_to_index('HN'): {map_medium_to_index('HN')}")
    
    print("\nTesting phase diagram functions:")
    print(f"get_phase_diagram_color('dominance'): {get_phase_diagram_color('dominance')}")
    print(f"get_phase_diagram_colors(): {get_phase_diagram_colors()}")
    print(f"get_phase_diagram_color_by_index(1): {get_phase_diagram_color_by_index(1)}")
    print(f"map_phase_to_index('mixing'): {map_phase_to_index('mixing')}")
    print(f"map_phase_to_index(2): {map_phase_to_index(2)}")