config = {
    "width": 1.0,
    "height": 1.0,
    "nx": 10,
    "ny": 10,
    "materials": [
        {
            "name": "Square_5_inner",
            "permittivity": 1.0,
            "x_range": [0.4, 0.6],
            "y_range": [0.4, 0.6]
        },
        {
            "name": "Square_4",
            "permittivity": 2.0,
            "x_range": [0.3, 0.7],
            "y_range": [0.3, 0.7]
        },
        {
            "name": "Square_3",
            "permittivity": 3.0,
            "x_range": [0.2, 0.8],
            "y_range": [0.2, 0.8]
        },
        {
            "name": "Square_2",
            "permittivity": 4.0,
            "x_range": [0.1, 0.9],
            "y_range": [0.1, 0.9]
        },
        {
            "name": "Square_1",
            "permittivity": 5.0,
            "x_range": [0.0, 1.0],
            "y_range": [0.0, 1.0]
        }
    ],

    "boundary_conditions": {
    "left":   {"type": "dirichlet", "value": 100.0},
    "right":  {"type": "dirichlet", "value": 80.0},
    "top":    {"type": "dirichlet", "value": 0.0},
    "bottom": {"type": "dirichlet", "value": 0.0}
}
}