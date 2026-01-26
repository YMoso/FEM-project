def config_setup(width, height, nx, ny, boundary_conditions_dirichlet, boundary_conditions_neumann, line_values, Layer_Left, Center_Bottom, Center_Middle, Center_Top, Layer_Right):
    xb = [0.0, 0.33 * width, 0.66 * width, width]  # x-breakpoints
    yb = [0.0, 0.33 * height, 0.66 * height, height] # y-breakpoints
    config = {
        "width": width,
        "height": height,
        "nx": nx,
        "ny": ny,
        "boundary_conditions_dirichlet": boundary_conditions_dirichlet,
        "boundary_conditions_neumann": boundary_conditions_neumann,
        "line_values": line_values,
        "materials": [
            {"name": "Layer_Left",    "permittivity": Layer_Left, "x_range": [xb[0], xb[1]], "y_range": [0.0, height]},
            {"name": "Center_Bottom", "permittivity": Center_Bottom, "x_range": [xb[1], xb[2]], "y_range": [yb[0], yb[1]]},
            {"name": "Center_Middle", "permittivity": Center_Middle, "x_range": [xb[1], xb[2]], "y_range": [yb[1], yb[2]]},
            {"name": "Center_Top",    "permittivity": Center_Top, "x_range": [xb[1], xb[2]], "y_range": [yb[2], yb[3]]},
            {"name": "Layer_Right",   "permittivity": Layer_Right, "x_range": [xb[2], xb[3]], "y_range": [0.0, height]}
        ]
    }
    return config
