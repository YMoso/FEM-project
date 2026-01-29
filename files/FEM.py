import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from files.Material import Material
from files.Mesh import Mesh


def initialize_fem(config):
    materials = []
    for m in config["materials"]:
        materials.append(Material(m["name"], m["permittivity"], m["x_range"], m["y_range"]))
    mesh = Mesh(config["width"], config["height"], config["nx"], config["ny"], materials)
    solver = FEM(mesh, config)


class FEM:
    def __init__(self, mesh, config):
        self.mesh = mesh
        self.config = config
        self.n_nodes = len(mesh.nodes)
        self.rows = []
        self.cols = []
        self.data = []
        self.F = np.zeros(self.n_nodes)
        self.K = None
        self.solution = None
        self.phi = None
        self.build_system()

    def build_system(self):
        self.assemble_stiffness_matrix()
        self.rows = np.array(self.rows)
        self.cols = np.array(self.cols)
        self.data = np.array(self.data)
        mask = self.rows != self.cols
        full_rows = np.concatenate([self.rows, self.cols[mask]])
        full_cols = np.concatenate([self.cols, self.rows[mask]])
        full_data = np.concatenate([self.data, self.data[mask]])
        self.K = coo_matrix((full_data, (full_rows, full_cols)),shape=(self.n_nodes, self.n_nodes)).tocsr()
        self.rows, self.cols, self.data = None, None, None
        self.boundaries()
        self.plots()
        self.plot_solution_along_line(self.config["line_values"]["x0"], self.config["line_values"]["y0"], self.config["line_values"]["x1"], self.config["line_values"]["y1"])
        self.mesh.plot_mesh_with_materials()
        self.save_nodal_results()
        plt.show()

    def assemble_stiffness_matrix(self):
        for element in self.mesh.elements:
            Ke = self.compute_element_stiffness(element)
            self.assemble_element(element, Ke)

    def assemble_element(self, element, Ke):
        for a in range(3):
            A = element.node_indices[a]
            for b in range(3):
                B = element.node_indices[b]
                if A <= B:
                    self.rows.append(A)
                    self.cols.append(B)
                    self.data.append(Ke[a, b])

    def compute_element_stiffness(self, element):
        n0, n1, n2 = element.node_indices

        x1, y1 = self.mesh.nodes[n0].x, self.mesh.nodes[n0].y
        x2, y2 = self.mesh.nodes[n1].x, self.mesh.nodes[n1].y
        x3, y3 = self.mesh.nodes[n2].x, self.mesh.nodes[n2].y

        area = 0.5 * abs((x2 - x1) * (y3 - y1)  - (x3 - x1) * (y2 - y1))

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

    def apply_neumann_bc(self, node_ids, gradient):
        for i in range(len(node_ids) - 1):
            n1 = node_ids[i]
            n2 = node_ids[i + 1]
            x1, y1 = self.mesh.nodes[n1].x, self.mesh.nodes[n1].y
            x2, y2 = self.mesh.nodes[n2].x, self.mesh.nodes[n2].y
            edge_length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            xm = 0.5 * (x1 + x2)
            ym = 0.5 * (y1 + y2)
            eps = self.mesh.find_material(xm, ym).permittivity
            q = eps * gradient
            contribution = q * edge_length/2.0
            self.F[n1] += contribution
            self.F[n2] += contribution

    def solve(self):
        self.solution = spsolve(self.K, self.F)
        return self.solution

    def boundaries(self):
        node_boundaries = self.mesh.get_boundary_nodes()
        for side, value in self.config["boundary_conditions_dirichlet"].items():
            node_ids = node_boundaries[side]
            self.apply_dirichlet_bc(node_ids, [value] * len(node_ids))

        for side, q in self.config["boundary_conditions_neumann"].items():
            self.apply_neumann_bc(node_boundaries[side],q)


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
        if self.config["line_values"]["line_visible"]:
            x0, y0, x1, y1 =  self.config["line_values"]["x0"], self.config["line_values"]["y0"], self.config["line_values"]["x1"], self.config["line_values"]["y1"]
            plt.plot([x0, x1], [y0, y1], "r--", lw=2, label="sampling line")
            plt.legend()


    def save_nodal_results(self, filename="solution_nodes.txt"):
        x = np.array([node.x for node in self.mesh.nodes.values()])
        y = np.array([node.y for node in self.mesh.nodes.values()])
        data = np.column_stack((x, y, self.solution))
        head = "x [m]\t y [m]\t phi [V]"
        np.savetxt(filename, data, header=head, comments="", fmt="%.6f")

    def plot_solution_along_line(self, x0, y0, x1, y1):
        n_points = 100
        xs = np.linspace(x0, x1, n_points)
        ys = np.linspace(y0, y1, n_points)
        values = []
        distances = []
        total_length = np.sqrt((x1-x0) ** 2 +(y1 - y0) ** 2)
        for i in range(n_points):
            x = xs[i]
            y = ys[i]
            s = total_length *i /(n_points - 1)
            phi = self.interpolate_at_point(x, y)
            values.append(phi)
            distances.append(s)

        plt.figure()
        plt.plot(distances, values, "b-")
        plt.xlabel("distance along line [m]")
        plt.ylabel("electric potential [V]")
        plt.grid(True)
        plt.title("Electric potential")

    def interpolate_at_point(self, x, y):
        for element in self.mesh.elements:
            n0, n1, n2 = element.node_indices
            x1, y1 = self.mesh.nodes[n0].x, self.mesh.nodes[n0].y
            x2, y2 = self.mesh.nodes[n1].x, self.mesh.nodes[n1].y
            x3, y3 = self.mesh.nodes[n2].x, self.mesh.nodes[n2].y
            det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
            if abs(det) < 1e-14:
                continue
            l1 = ((y2-y3)*(x-x3) + (x3 - x2) * (y -y3))/det
            l2 = ((y3 -y1) *(x -x3) +(x1-x3) *(y-y3)) / det
            l3 = 1.0 -l1 -l2
            if 0 <= l1 <= 1 and 0 <= l2 <= 1 and 0 <= l3 <= 1:
                phi = (l1 *self.solution[n0] +l2 *self.solution[n1] +l3 * self.solution[n2])
                return phi

        return np.nan
