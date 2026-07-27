# IMPORTS
import MDAnalysis as mda
import numpy as np
import networkit as nk
import networkx as nx
import pandas as pd
from itertools import combinations
from MDAnalysis.lib import distances

# CONFIGURATION
DMBA_PDB_FILES = [
    "/mnt/d/work/dmba.hso4/1dmba-h2o/dmba_hso4_1h2o_tip4e_l100ns.pdb",
    "/mnt/d/work/dmba.hso4/2dmba-h2o/dmba_hso4_2h2o_tip4e_l100ns.pdb",
    "/mnt/d/work/dmba.hso4/3dmba-h2o/dmba_hso4_3h2o_tip4e_l100ns.pdb",
    "/mnt/d/work/dmba.hso4/4dmba-h2o/dmba_hso4_4h2o_tip4e_l100ns.pdb",
    "/mnt/d/work/dmba.hso4/5dmba-h2o/dmba_hso4_5h2o_tip4e_l100ns.pdb",
    "/mnt/d/work/dmba.hso4/6dmba-h2o/dmba_hso4_6h2o_tip4e_l100ns.pdb"
]

HAM_PDB_FILES = [
    "/mnt/d/work/ha.hso4/1h2o_ha/ham_hso4_1h2o_tip4ew_l100ns.pdb",
    "/mnt/d/work/ha.hso4/2h2o_ha/ham_hso4_2h2o_tip4ew_l100ns.pdb",
    "/mnt/d/work/ha.hso4/3h2o_ha/ham_hso4_3h2o_tip4ew_l100ns.pdb",
    "/mnt/d/work/ha.hso4/4h2o_ha/ham_hso4_4h2o_tip4ew_l100ns.pdb",
    "/mnt/d/work/ha.hso4/5h2o_ha/ham_hso4_5h2o_tip4ew_l100ns.pdb",
    "/mnt/d/work/ha.hso4/6h2o_ha/ham_hso4_6h2o_tip4ew_l100ns.pdb"
]

# CUTOFFS
CUTOFFS = {

    "DMBA": {
        ("WAT", "WAT"): 3.4,
        ("WAT", "DMBA"): 5.5,
        ("WAT", "HSO"): 5.6,
        ("DMBA", "HSO"): 6.7,
        ("DMBA", "DMBA"): 5.5,
        ("HSO", "HSO"): 5.6,
    },

    "HAM": {
        ("WAT", "WAT"): 3.3,
        ("WAT", "HAM"): 4.0,
        ("WAT", "HSO"): 5.0,
        ("HAM", "HSO"): 5.2,
        ("HAM", "HAM"): 5.5,
        ("HSO", "HSO"): 5.6,
    }
}

# SOLVATION SHELL DEFINITIONS
def interaction_shell(system, t1, t2, dist):

    # DMBA - HSO
    if system == "DMBA":
        if (t1 == "DMBA" and t2 == "HSO") or (t1 == "HSO" and t2 == "DMBA"):
            if dist <= 6.7:
                return "DMBA_HSO_1"

            elif dist <= 10.5:
                return "DMBA_HSO_2"

            elif dist <= 15.8:
                return "DMBA_HSO_3"

            return None

    # HAM - HSO
    if system == "HAM":
        if (t1 == "HAM" and t2 == "HSO") or (t1 == "HSO" and t2 == "HAM"):

            if dist <= 6.3:
                return "HAM_HSO_1"

            elif dist <= 10.5:
                return "HAM_HSO_2"

            return None

    # HSO - HSO
    if t1 == "HSO" and t2 == "HSO":

        if dist <= 5.6:
            return "HSO_HSO_1"

        elif dist <= 10.5:
            return "HSO_HSO_2"

        elif dist <= 15.8:
            return "HSO_HSO_3"

        return None

    return "OTHER"

# MOLECULE GROUPS
def get_molecule_groups(u, system):

    groups = {
        "WAT": u.select_atoms("resname WAT"),
        "HSO": u.select_atoms("resname HSO")
    }

    if system == "DMBA":
        groups["DMBA"] = u.select_atoms("resname BMM")

    elif system == "HAM":
        groups["HAM"] = u.select_atoms("resname HAM")

    return groups

# NODE POSITIONS
def build_node_positions(groups):

    labels = []
    pos = []
    types = []


    for name, atoms in groups.items():
        for res in atoms.residues:
            labels.append(f"{name}_{res.resid}")

            pos.append(res.atoms.positions.mean(axis=0))

            types.append(name)


    return labels, np.array(pos), np.array(types)

# CUTOFF LOOKUP
def get_cutoff(system, t1, t2):

    cutoff = CUTOFFS[system]

    return (cutoff.get((t1, t2)) or cutoff.get((t2, t1)))

# EDGE BUILDING
def build_edges(pos, types, system, box=None):

    n = len(pos)

    if box is not None and np.all(box[:3] > 0):
        dist_matrix = distances.distance_array(pos, pos,box=box)

    else:
        dist_matrix = np.linalg.norm(
            pos[:, None, :] - pos[None, :, :],
            axis=-1
        )

    edges = []
    weighted_edges = []


    for i, j in combinations(range(n), 2):
        t1, t2 = types[i], types[j]

        # Solvation shells
        if system == "DMBA" and (
            {t1, t2} == {"DMBA", "HSO"}
            or {t1, t2} == {"HSO", "HSO"}
        ):
            cutoff = 15.8

        elif system == "HAM" and (
            {t1, t2} == {"HAM", "HSO"}
        ):
            cutoff = 10.5

        elif t1 == "HSO" and t2 == "HSO":
            cutoff = 15.8

        else:
            cutoff = get_cutoff(system, t1, t2)

        if cutoff is None or dist_matrix[i, j] > cutoff:
            continue

        shell = interaction_shell(system, t1, t2, dist_matrix[i, j])

        if shell is None:
            continue

        weight = 1 / (dist_matrix[i, j] + 1e-6)

        edges.append((i, j, shell))

        weighted_edges.append((i, j, weight))

    return edges, weighted_edges

# GRAPH
def build_graph(n_nodes, weighted_edges):
    G = nk.Graph(n_nodes, weighted=True)

    for i, j, w in weighted_edges:
        G.addEdge(i, j, w)

    return G

# CLOSENESS CENTRALITY
def compute_closeness_nx(n_nodes, edges):
    Gx = nx.Graph()
    Gx.add_nodes_from(range(n_nodes))
    Gx.add_edges_from([(i, j) for i, j, _ in edges])
    cc = nx.closeness_centrality(Gx)

    return np.array([cc[i] for i in range(n_nodes)])

# SPECIES METRICS
def species_metrics(deg_values, bc_values, ec_values, types, species):
    idx = np.where(types == species)[0]

    if len(idx) == 0:
        return 0, 0, 0

    return (
        np.mean(deg_values[idx]),
        np.mean(bc_values[idx]),
        np.mean(ec_values[idx])
    )

# SHELL COUNTS
def shell_counts(edges, system):

    if system == "DMBA":
        cation = {1: 0, 2: 0, 3: 0}

    else:
        cation = {1: 0,2: 0}

    hso = {1: 0,2: 0,3: 0}

    for _, _, shell in edges:
        if "DMBA_HSO" in shell:
            cation[int(shell[-1])] += 1

        elif "HAM_HSO" in shell:
            cation[int(shell[-1])] += 1

        elif "HSO_HSO" in shell:
            hso[int(shell[-1])] += 1

    return cation, hso

# MAIN ANALYSIS
def analyze(pdb_file, system):
    u = mda.Universe(pdb_file)
    groups = get_molecule_groups(u, system)
    results = []
    prev_edges = set()

    for ts in u.trajectory:
        labels, pos, types = build_node_positions(groups)

        if len(pos) == 0:
            continue

        box = (
            u.dimensions
            if hasattr(u, "dimensions")
            else None
        )

        edges, weighted_edges = build_edges(pos, types, system, box)

        G = build_graph(len(pos), weighted_edges)

        if G.numberOfNodes() == 0:
            continue

        # CONNECTIVITY
        cc = nk.components.ConnectedComponents(G)
        cc.run()
        percolation = (max(cc.getComponentSizes().values()) / G.numberOfNodes())

        # CENTRALITIES
        deg_values = np.array([G.degree(i) for i in G.iterNodes()])

        avg_degree = np.mean(deg_values)

        bc = nk.centrality.Betweenness(G)
        bc.run()
        bc_values = np.array(bc.scores())

        ec = nk.centrality.EigenvectorCentrality(G)
        ec.run()
        ec_values = np.array(ec.scores())

        avg_bc = np.mean(bc_values)

        avg_ec = np.mean(ec_values)

        # NETWORKX GRAPH
        Gx = nx.Graph()
        Gx.add_nodes_from(range(len(pos)))
        Gx.add_edges_from([(i, j) for i, j, _ in edges])

        # clustering
        clustering = (
            np.mean(
                list(nx.clustering(Gx).values())
            )
            if edges
            else 0
        )

        # density
        density = nx.density(Gx)

        # K-CORE
        if Gx.number_of_edges() > 0:
            core_numbers = nx.core_number(Gx)

            avg_kcore = np.mean(list(core_numbers.values()))

            wat_kcore, cat_kcore, an_kcore = [np.mean(
                    [core_numbers[i]
                        for i in range(len(pos))
                        if types[i] == t
                    ]
                )

                for t in [
                    "WAT",
                    system,
                    "HSO"
                ]
            ]

        else:
            avg_kcore = 0
            wat_kcore = 0
            cat_kcore = 0
            an_kcore = 0

        # ASSORTATIVITY
        nx.set_node_attributes(
            Gx,
            {i: t
                for i, t in enumerate(types)
            },
            "type"
        )

        assort = nx.attribute_assortativity_coefficient(Gx, "type")

        # degree variance
        deg_var = np.var(deg_values)

        # SHELL COUNTS
        cation, hso = shell_counts(edges, system)

        # EDGE PERSISTENCE
        edge_set = {(i, j)
            for i, j, _ in edges
        }

        persistence = (
            len(edge_set & prev_edges)
            /
            len(edge_set)
            if edge_set
            else 0
        )

        prev_edges = edge_set

        # SPECIES METRICS
        wat = species_metrics(deg_values, bc_values, ec_values, types, "WAT")
        cat = species_metrics(deg_values, bc_values, ec_values, types, system)

        an = species_metrics(deg_values, bc_values, ec_values, types, "HSO")

        results.append([percolation, avg_degree,
                        avg_bc, avg_cc, avg_ec,
                        clustering, density, avg_kcore,
                        persistence,
                        wat_cc, cat_cc, an_cc,
                        *wat, *cat, *an,
                        *cation.values(),
                        *hso.values(),
                        wat_kcore, cat_kcore, an_kcore,
                        assort,
                        deg_var])

        print(
            f"{system} Frame {ts.frame}: "
            f"percolation={percolation:.3f}"
        )

    return np.array(results)

# OUTPUT COLUMNS

DMBA_COLUMNS = [ "percolation","avg_degree",
                "betweenness", "closeness", "eigenvector",
                "clustering", "density", "avg_kcore",
                "persistence",
                "wat_cc", "cat_cc", "an_cc",
                "wat_deg", "wat_bc", "wat_ec",
                "cat_deg", "cat_bc", "cat_ec",
                "an_deg", "an_bc", "an_ec",
                "dmba_hso_s1", "dmba_hso_s2", "dmba_hso_s3",
                "hso_hso_s1", "hso_hso_s2", "hso_hso_s3",
                "wat_kcore", "cat_kcore", "an_kcore",
                "assortativity",
                "degree_variance"]


HAM_COLUMNS = ["percolation", "avg_degree",
               "betweenness", "closeness", "eigenvector",
               "clustering", "density", "avg_kcore",
               "persistence",
               "wat_cc", "cat_cc", "an_cc",
               "wat_deg", "wat_bc", "wat_ec",
               "cat_deg", "cat_bc", "cat_ec",
               "an_deg", "an_bc", "an_ec",
               "ham_hso_s1", "ham_hso_s2",
               "hso_hso_s1", "hso_hso_s2", "hso_hso_s3",
               "wat_kcore", "cat_kcore", "an_kcore",
               "assortativity",
               "degree_variance"]

# RUN ANALYSIS
def run_system_analysis(pdb_files, system):
    columns = (
        DMBA_COLUMNS
        if system == "DMBA"
        else HAM_COLUMNS
    )

    for pdb_file in pdb_files:
        print(f"Analyzing {system}: {pdb_file}")

        data = analyze(pdb_file, system)
        name = (pdb_file.split("/")[-1].replace(".pdb", ""))
        df = pd.DataFrame(data,columns=columns)
        output = (f"{system.lower()}_{name}_network.csv")
        df.to_csv(output, index=False)

        print(f"Saved: {output}")

# MAIN

if __name__ == "__main__":

    # DMBA-HSO4 systems
    run_system_analysis(DMBA_PDB_FILES, "DMBA")

    # HAM-HSO4 systems
    run_system_analysis(HAM_PDB_FILES, "HAM")

    print("\nAll network analyses completed")

