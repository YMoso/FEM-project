import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


from CONFIG.config import config
from files.Material import Material
from files.Mesh import Mesh

class FEM:
    def __init__(self, mesh):
        self.mesh = mesh
        self.n_nodes = len(mesh.nodes)

        self.K = lil_matrix((self.n_nodes, self.n_nodes))
        self.F = np.zeros(self.n_nodes)
        self.solution = None

    def build_system(self):
        self.assemble_stiffness_matrix()
        self.K = self.K.tocsr()

    def assemble_stiffness_matrix(self):
        for element in self.mesh.elements:
            Ke = self.compute_element_stiffness(element)
            self.assemble_element(element, Ke)

    def assemble_element(self, element, Ke):
        for a in range(3):
            A = element.node_indices[a]
            for b in range(3):
                B = element.node_indices[b]
                self.K[A, B] += Ke[a, b]

    def compute_element_stiffness(self, element):
        nodes = self.mesh.nodes
        idx = element.node_indices

        coords = np.array([[nodes[i].x, nodes[i].y] for i in idx])
        x = coords[:, 0]
        y = coords[:, 1]

        area = 0.5 * abs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))

        B = (np.array([[y[1] - y[2], y[2] - y[0], y[0] - y[1]],[x[2] - x[1], x[0] - x[2], x[1] - x[0]]])
             /(2 * area))

        return element.material.permittivity * area * (B.T @ B)

    def apply_dirichlet_bc(self, node_ids, values):
        for node_id, value in zip(node_ids, values):
            self.K[node_id, :] = 0.0
            self.K[node_id, node_id] = 1.0
            self.F[node_id] = value

            for i in range(self.n_nodes):
                if i != node_id:
                    self.F[i] -= self.K[i, node_id] * value
                    self.K[i, node_id] = 0.0

    def solve(self):
        self.solution = spsolve(self.K, self.F)
        return self.solution


def main():
    materials = []
    for m in config["materials"]:
        materials.append(Material(m["name"],m["permittivity"],m["x_range"],m["y_range"]))

    mesh = Mesh(config["width"],config["height"],config["nx"],config["ny"],materials)

    solver = FEM(mesh)
    solver.build_system()

    boundaries = mesh.get_boundary_nodes()

    for side, bc in config["boundary_conditions"].items():
        node_ids = boundaries[side]
        if bc["type"] == "dirichlet":
            solver.apply_dirichlet_bc(node_ids,[bc["value"]] * len(node_ids))

    phi = solver.solve()

    x = np.array([node.x for node in mesh.nodes.values()])
    y = np.array([node.y for node in mesh.nodes.values()])

    plt.tricontourf(x, y, phi, levels=25, )
    plt.colorbar(label="Electric potential[V]")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("2D FEM Electric Potential")
    plt.gca().set_aspect("equal")
    plt.show()


if __name__ == "__main__":
    main()