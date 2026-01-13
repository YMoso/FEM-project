class Material:
    def __init__(self, name, permittivity, x_range, y_range):
        self.name = name
        self.permittivity = permittivity
        self.x_min, self.x_max = x_range
        self.y_min, self.y_max = y_range

    def contains(self, x, y):
        return (self.x_min <= x <= self.x_max and self.y_min <= y <= self.y_max)