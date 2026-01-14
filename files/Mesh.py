from files.Node import Node
from files.TriangularElement import TriangularElement


class Mesh:
    def __init__(self, width, height, nx, ny, materials):
        self.width = width
        self.height = height
        self.nx = nx
        self.ny = ny
        self.materials = materials
        self.hx = width / nx
        self.hy = height / ny
        self.nodes = {}
        self.elements = []
        self.corners = {"left": [],"right": [],"top": [],"bottom": []}
        self.edges = {"left": [],"right": [],"top": [],"bottom": []}

        self.node_map = self.create_nodes()
        self.create_elements(self.node_map)

    def create_nodes(self):
        node_map = {}
        node_id = 0
        for j in range(self.ny + 1):
            for i in range(self.nx + 1):
                x = i * self.hx
                y = j * self.hy
                self.nodes[node_id] = Node(node_id, x, y)
                node_map[(i, j)] = node_id
                node_id += 1

        return node_map

    def create_elements(self, node_map):
        element_id = 0
        for j in range(self.ny):
            for i in range(self.nx):
                n0 = node_map[(i, j)]
                n1 = node_map[(i + 1, j)]
                n2 = node_map[(i + 1, j + 1)]
                n3 = node_map[(i, j + 1)]
                x_centroid = (i + 0.5) * self.hx
                y_centroid = (j + 0.5) * self.hy

                material = self.find_material(x_centroid, y_centroid)
                self.elements.append(TriangularElement(element_id, [n0, n1, n3], material))
                element_id += 1
                self.elements.append(TriangularElement(element_id, [n1, n2, n3], material))
                element_id += 1

    def find_material(self, x, y):
        for material in self.materials:
            if material.contains(x, y):
                return material
        raise ValueError(f"No material")

    def get_boundary_nodes(self, tol=1e-9):
        boundaries = self.corners.copy()

        for node_id, node in self.nodes.items():
            if abs(node.x) < tol:
                boundaries["left"].append(node_id)
            if abs(node.x - self.width) < tol:
                boundaries["right"].append(node_id)
            if abs(node.y) < tol:
                boundaries["bottom"].append(node_id)
            if abs(node.y - self.height) < tol:
                boundaries["top"].append(node_id)

        return boundaries

    def get_boundary_edges(self):
        for i in range(self.nx):
            self.edges["bottom"].append((self.node_map[(i, 0)], self.node_map[(i + 1, 0)]))
            self.edges["top"].append((self.node_map[(i, self.ny)], self.node_map[(i + 1, self.ny)]))
        for j in range(self.ny):
            self.edges["left"].append((self.node_map[(0, j)], self.node_map[(0, j + 1)]))
            self.edges["right"].append((self.node_map[(self.nx, j)], self.node_map[(self.nx, j + 1)]))

        return self.edges