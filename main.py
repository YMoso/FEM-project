from files.FEM import initialize_fem
from CONFIG.config import config_setup

def main():
    # Domain Dimensions [m]
    width = 1.0
    height = 1.0

    # Mesh Density (Number of divisions)
    nx = 20
    ny = 20

    # Boundary Conditions (Dirichlet [V], Neumann [V/m])
    boundary_conditions_dirichlet = {"left": 20.0, "right": 10.0}
    boundary_conditions_neumann = {"top": 0.0, "bottom": 0.0}

    # Material Permittivity Values [F/m]
    Layer_Left = 5.0
    Center_Bottom = 2.0
    Center_Middle = 10.0
    Center_Top = 2.0
    Layer_Right = 5.0

    # Visualization Settings
    line_values = {"line_visible": True, "x0": 0.0, "x1": 1.0, "y0": 0.0, "y1": 1.0}

    config = config_setup(width, height, nx, ny, boundary_conditions_dirichlet, boundary_conditions_neumann, line_values, Layer_Left, Center_Bottom, Center_Middle, Center_Top, Layer_Right)

    initialize_fem(config)

if __name__ == "__main__":
    main()