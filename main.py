from files.FEM import initialize_fem
from CONFIG.config import config_setup

def main():
    # Domain Dimensions [m]
    width = 2.0
    height = 2.0

    # Mesh Density (Number of divisions)
    nx = 20
    ny = 20

    # Boundary Conditions (Dirichlet [V], Neumann [V/m])
    boundary_conditions_dirichlet = {"left": 20.0, "right": 10.0}
    boundary_conditions_neumann = {"top": -5, "bottom": 5}

    # Material Permittivity Values [F/m]
    Layer_Left = 5.0
    Center_Bottom = 5.0
    Center_Middle = 5.0
    Center_Top = 5.0
    Layer_Right = 5.0

    # Visualization Settings
    line_values = {"line_visible": False, "x0": 0.0, "x1": 1.0, "y0": width / 3, "y1": width / 3}

    config = config_setup(width, height, nx, ny, boundary_conditions_dirichlet, boundary_conditions_neumann, line_values, Layer_Left, Center_Bottom, Center_Middle, Center_Top, Layer_Right)

    initialize_fem(config)

if __name__ == "__main__":
    main()