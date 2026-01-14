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
        self.phi = None
        self.build_system()

    def build_system(self):
        self.assemble_stiffness_matrix()
        self.K = self.K.tocsr()
        self.boundaries()
        self.plots()

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
        n0, n1, n2 = element.node_indices

        x1, y1 = self.mesh.nodes[n0].x, self.mesh.nodes[n0].y
        x2, y2 = self.mesh.nodes[n1].x, self.mesh.nodes[n1].y
        x3, y3 = self.mesh.nodes[n2].x, self.mesh.nodes[n2].y

        area = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

        if area < 1e-12:
            raise ValueError(f"Element {element.id} has near-zero area")

        b1, c1 = (y2 - y3), (x3 - x2)
        b2, c2 = (y3 - y1), (x1 - x3)
        b3, c3 = (y1 - y2), (x2 - x1)

        b = np.array([b1, b2, b3], dtype=float)
        c = np.array([c1, c2, c3], dtype=float)

        Ke = (element.material.permittivity/(4.0*area)) * (np.outer(b, b) + np.outer(c, c))
        return Ke

    def apply_dirichlet_bc(self, node_ids, values):
        for node_id, value in zip(node_ids, values):
            self.K[node_id, :] = 0.0
            self.K[node_id, node_id] = 1.0
            self.F[node_id] = value

            for i in range(self.n_nodes):
                if i != node_id:
                    self.F[i] -= self.K[i, node_id] * value
                    self.K[i, node_id] = 0.0

    def apply_neumann_bc(self, edge_list, q):
        for n1, n2 in edge_list:
            x1, y1 = self.mesh.nodes[n1].x, self.mesh.nodes[n1].y
            x2, y2 = self.mesh.nodes[n2].x, self.mesh.nodes[n2].y

            L = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

            self.F[n1] += q * L / 2.0
            self.F[n2] += q * L / 2.0


    def solve(self):
        self.solution = spsolve(self.K, self.F)
        return self.solution

    def boundaries(self):
        node_boundaries = self.mesh.get_boundary_nodes()
        edge_boundaries = self.mesh.get_boundary_edges()

        for side, bc in config["boundary_conditions"].items():
            if bc["type"] == "neumann":
                self.apply_neumann_bc(edge_boundaries[side], bc["value"])

        for side, bc in config["boundary_conditions"].items():
            if bc["type"] == "dirichlet":
                self.apply_dirichlet_bc(
                    node_boundaries[side],
                    [bc["value"]] * len(node_boundaries[side]))



    def plots(self):
        self.phi = self.solve()

        x = np.array([node.x for node in self.mesh.nodes.values()])
        y = np.array([node.y for node in self.mesh.nodes.values()])

        plt.tricontourf(x, y, self.phi, levels=25)
        plt.colorbar(label="Electric potential[V]")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("2D FEM Electric Potential")
        plt.gca().set_aspect("equal")
        plt.show()



def main():
    materials = []
    for m in config["materials"]:
        materials.append(Material(m["name"],m["permittivity"],m["x_range"],m["y_range"]))

    mesh = Mesh(config["width"],config["height"],config["nx"],config["ny"],materials)

    solver = FEM(mesh)

if __name__ == "__main__":
    main()