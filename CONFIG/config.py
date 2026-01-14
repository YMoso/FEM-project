config = {
    "width": 1.0,
    "height": 1.0,
    "nx": 10,
    "ny": 10,
    "materials": [
        # {
        #     "name": "Square_5_inner",
        #     "permittivity": 1.0,
        #     "x_range": [0.0, 0.],
        #     "y_range": [0.4, 0.6]
        # },
        {
            "name": "Square_4",
            "permittivity": 1000.0,
            "x_range": [0.0, 0.5],
            "y_range": [0.0, 0.5]
        },
        {
            "name": "Square_3",
            "permittivity": 1.0,
            "x_range": [0.0, 0.5],
            "y_range": [0.5, 1.0]
        },
        {
            "name": "Square_2",
            "permittivity": 5.0,
            "x_range": [0.5, 1.0],
            "y_range": [0.0, 0.5]
        },
        {
            "name": "Square_1",
            "permittivity": 200.0,
            "x_range": [0.5, 1.0],
            "y_range": [0.5, 1.0]
        }
    ],

    "boundary_conditions": {
    "left":   {"type": "dirichlet", "value": 300.0},
    "right":  {"type": "dirichlet", "value": 200.0},
    "top":    {"type": "neumann", "value": -30.0},
    "bottom": {"type": "dirichlet", "value": 200.0}
}
}