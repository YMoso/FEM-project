config = {
    "width": 1.0,
    "height": 1.0,
    "nx": 10,
    "ny": 10,
    "materials": [

        {
            "name": "Square_1",
            "permittivity": 1.0,
            "x_range": [0.0, 0.3],
            "y_range": [0.0, 1.0]
        },
        {
            "name": "Square_2",
            "permittivity": 1.0,
            "x_range": [0.3, 0.7],
            "y_range": [0.0, 0.35]
        },
        {
            "name": "Square_3",
            "permittivity": 1.0,
            "x_range": [0.3, 0.7],
            "y_range": [0.35, 0.7]
        },
        {
            "name": "Square_4",
            "permittivity": 1.0,
            "x_range": [0.3, 0.7],
            "y_range": [0.7, 1.0]
        },
{
            "name": "Square_5",
            "permittivity": 1.0,
            "x_range": [0.7, 1.0],
            "y_range": [0.0, 1.0]
        }
    ],

    "boundary_conditions_dirichlet": {
    "left": 0.0,
    "right":5.0
},
    "boundary_conditions_neumann": {
    "top": 0,
    "bottom":0
}
}