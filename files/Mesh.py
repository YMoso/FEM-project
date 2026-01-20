import numpy as np
import matplotlib.pyplot as plt
from files.Node import Node
from files.TriangularElement import TriangularElement


class Mesh:
    def __init__(self, width, height, nx, ny, materials):
        self.width = width
        self.height = height
        self.nx = nx
        self.ny = ny
        self.materials = materials
        self.x_coords = self.generate_axis_nodes(axis="x")
        self.y_coords = self.generate_axis_nodes(axis="y")
        self.nodes = {}
        self.elements = []
        self.node_map = self.create_nodes()
        self.create_elements()

    def generate_axis_nodes(self, axis):
        breakpoints = self.get_axis_breakpoints(axis)
        region_lengths, total_length = self.compute_region_lengths(breakpoints)
        element_counts = self.allocate_elements_per_region(axis, region_lengths, total_length)
        return self.assemble_axis_coordinates(breakpoints, element_counts)

    def get_axis_breakpoints(self, axis):
        breakpoints = set()
        for m in self.materials:
            if axis == "x":
                breakpoints.update([m.x_min, m.x_max])
            else:
                breakpoints.update([m.y_min, m.y_max])
        return sorted(breakpoints)

    def compute_region_lengths(self, breakpoints):
        region_length = []
        for i in range(len(breakpoints) - 1):
            region_length.append(breakpoints[i+1] - breakpoints[i])

        total_length = sum(region_length)
        return region_length, total_length


    def allocate_elements_per_region(self, axis, region_lengths, total_length):
        if axis == "x":
            total_elements = self.nx
        else:
            total_elements = self.ny
        counts = []
        for L in region_lengths:
            n = int(total_elements*L/total_length)
            counts.append(max(1, n))

        while sum(counts) < total_elements:
            for i in range(len(counts)):
                if sum(counts) < total_elements:
                    counts[i] += 1
        return counts

    def assemble_axis_coordinates(self, breakpoints, element_counts):
        coords = [breakpoints[0]]

        for i, n in enumerate(element_counts):
            start = breakpoints[i]
            end = breakpoints[i + 1]
            local = np.linspace(start, end, n + 1)[1:]
            coords.extend(local)
        return np.array(coords)

    def create_nodes(self):
        node_map = {}
        node_id = 0

        for j, y in enumerate(self.y_coords):
            for i, x in enumerate(self.x_coords):
                self.nodes[node_id] = Node(node_id, x, y)
                node_map[(i, j)] = node_id
                node_id += 1

        return node_map

    def create_elements(self):
        nx = len(self.x_coords)
        ny = len(self.y_coords)

        element_id = 0
        for j in range(ny -1):
            for i in range(nx -1):
                # its my square
                n0 = j * nx + i
                n1 = j * nx + (i+1)
                n2 = (j+1) * nx + (i +1)
                n3 = (j+1) * nx + i

                # its my two triangles in this square
                triangles = [(n0, n1, n3),(n1, n2, n3)]

                for nodes in triangles:
                    x_c = sum(self.nodes[n].x for n in nodes)/3
                    y_c = sum(self.nodes[n].y for n in nodes)/3
                    material = self.find_material(x_c, y_c)
                    self.elements.append(TriangularElement(element_id, nodes, material))
                    element_id += 1

    def find_material(self, x, y):
        for material in self.materials:
            if material.contains(x, y):
                return material
        raise ValueError(f"No material")

    def get_boundary_nodes(self, tol=1e-9):
        boundaries = {"left": [],"right": [],"top": [],"bottom": []}
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

    def boundary_edges_from_nodes(self, node_ids, side):
        if side in ("left", "right"):
            node_ids = sorted(node_ids, key=lambda i: self.nodes[i].y)
        else:
            node_ids = sorted(node_ids, key=lambda i: self.nodes[i].x)

        return list(zip(node_ids[:-1], node_ids[1:]))

    def plot_mesh_with_materials(self):
        plt.figure(figsize=(6, 6))
        colors = ["red", "green", "blue", "orange", "purple", "cyan"]
        for i, material in enumerate(self.materials):
            color = colors[i % len(colors)]
            plt.fill([material.x_min, material.x_max, material.x_max, material.x_min],
                [material.y_min, material.y_min, material.y_max, material.y_max],color=color,alpha=0.25,
                     label=material.name)

        for element in self.elements:
            ids = element.node_indices
            x = [self.nodes[i].x for i in ids] + [self.nodes[ids[0]].x]
            y = [self.nodes[i].y for i in ids] + [self.nodes[ids[0]].y]
            plt.plot(x, y, "k-", linewidth=0.4)

        x_nodes = [node.x for node in self.nodes.values()]
        y_nodes = [node.y for node in self.nodes.values()]
        plt.scatter(x_nodes, y_nodes, s=6, c="black")

        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Mesh with material regions")
        plt.gca().set_aspect("equal")
        plt.grid(True)
        plt.legend()


