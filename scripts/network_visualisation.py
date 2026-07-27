# IMPORTS
import os
import MDAnalysis as mda
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from itertools import combinations
from MDAnalysis.lib import distances
import re

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

SYSTEMS = { "DMBA": DMBA_PDB_FILES,"HAM": HAM_PDB_FILES}

# OUTPUT SETTINGS
OUTPUT_DIR = "network_visualisations"


NODE_SIZE_2D = 55

NODE_SIZE_3D = 6



COLOURS = {
    "WAT":"#4C9BE8",
    "DMBA": "#FFEF00",
    "HAM": "#FFEF00",
    "HSO": "#e27602"
}

# INTERACTION CUTOFFS
CUTOFFS = {
    "DMBA": {
        ("WAT", "WAT"): 3.4,
        ("WAT", "DMBA"): 5.5,
        ("WAT", "HSO"): 5.6,
        ("DMBA", "HSO"): 6.7,
        ("DMBA", "DMBA"): 5.5,
        ("HSO", "HSO"): 5.6
    },

    "HAM": {
        ("WAT", "WAT"): 3.3,
        ("WAT", "HAM"): 4.0,
        ("WAT", "HSO"): 5.0,
        ("HAM", "HSO"): 5.2,
        ("HAM", "HAM"): 5.5,
        ("HSO", "HSO"): 5.6
    }
}

# SHELL DEFINITIONS
def interaction_shell(system, t1, t2, dist):
    # DMBA - HSO
    if system == "DMBA":
        if {t1, t2} == {"DMBA", "HSO"}:
            if dist <= 6.7:
                return "DMBA_HSO_1"

            elif dist <= 10.5:
                return "DMBA_HSO_2"

            elif dist <= 15.8:
                return "DMBA_HSO_3"

    # HAM - HSO
    if system == "HAM":
        if {t1, t2} == {"HAM", "HSO"}:
            if dist <= 6.3:
                return "HAM_HSO_1"

            elif dist <= 10.5:
                return "HAM_HSO_2"

    # HSO - HSO
    if {t1, t2} == {"HSO", "HSO"}:
        if dist <= 5.6:
            return "HSO_HSO_1"

        elif dist <= 10.5:
            return "HSO_HSO_2"

        elif dist <= 15.8:
            return "HSO_HSO_3"

    # WATER INTERACTIONS
    if {t1, t2} == {"WAT", "WAT"}:
        return "WAT_WAT"

    if "WAT" in {t1, t2}:
        return "WAT_INTERACTION"

    return None

# MOLECULE SELECTION
def get_groups(u, system):
    groups = {
        "WAT":
        u.select_atoms("resname WAT"),

        "HSO":
        u.select_atoms("resname HSO")
    }

    if system == "DMBA":
        groups["DMBA"] = u.select_atoms(
            "resname BMM"
        )

    elif system == "HAM":
        groups["HAM"] = u.select_atoms(
            "resname HAM"
        )

    return groups

# NODE GENERATION
def build_nodes(groups):
    labels = []
    species = []
    positions = []

    for mol, atoms in groups.items():
        for res in atoms.residues:
            labels.append(f"{mol}_{res.resid}")
            species.append(mol)
            positions.append(res.atoms.center_of_mass())

    return (labels, np.array(species), np.array(positions))

# BUILD NETWORK
def build_network(pdb_file, system):
    u = mda.Universe(pdb_file)
    # use final frame
    u.trajectory[-1]
    groups = get_groups(u,system)

    labels, species, positions = build_nodes(groups)

    # Distance matrix
    D = distances.distance_array(positions, positions, box=u.dimensions)
    G = nx.Graph()

    # Add nodes
    for i, label in enumerate(labels):
        G.add_node(i, species=species[i], label=label, pos=positions[i])

    # Add edges
    for i, j in combinations(range(len(labels)), 2):
        pair = (species[i], species[j])

        cutoff = (CUTOFFS[system].get(pair) or CUTOFFS[system].get(pair[::-1]))

        if cutoff is None:
            continue

        if D[i, j] <= cutoff:
            shell = interaction_shell(system, species[i], species[j], D[i, j])

            if shell is None:
                continue

            G.add_edge(i, j, distance=D[i, j], weight=1/(D[i, j] + 1e-6), shell=shell)

    print(
        "\nOriginal network:",
        G.number_of_nodes(),
        "nodes",
        G.number_of_edges(),
        "edges"
    )

    # Largest connected component
    if G.number_of_nodes() > 0:
        largest = max(nx.connected_components(G), key=len)

        G = G.subgraph(largest).copy()

    print("Largest component:",
          G.number_of_nodes(),
          "nodes",
          G.number_of_edges(),
          "edges")

    return G, positions

# NETWORK CENTRALITY
def calculate_centrality(G):
    bet = nx.betweenness_centrality(G)
    degree = nx.degree_centrality(G)
    eigen = nx.eigenvector_centrality(G, max_iter=1000)
    threshold = np.percentile(list(bet.values()), 90)

    important = [n
        for n, value in bet.items()
        if value >= threshold
]

    return (bet, degree, eigen, important)

# EDGE COLOURS
def edge_colour(shell):
    if "HSO" in shell:
        return "#e27602"

    if "WAT" in shell:
        return "#4C9BE8"

    return "0.70"

# 2D NETWORK VISUALISATION
def plot_network_2D(G, positions,system,filename,title):
    bet, degree, eigen, important = calculate_centrality(G)

    # Generate initial coordinates
    xy = np.array(
        [positions[n][:2]
            for n in G.nodes()])

    xy -= xy.mean(axis=0)

    coords = xy.copy()

    # Node separation
    min_dist = 2.5

    for _ in range(250):
        moved = False

        for i in range(len(coords)):
            for j in range(i+1, len(coords)):
                vector = (coords[j] - coords[i])
                dist = np.linalg.norm(vector)

                if dist < min_dist:
                    if dist < 1e-6:
                        vector = np.random.randn(2)
                        dist = np.linalg.norm(vector)

                    push = (min_dist - dist) / 2
                    direction = (vector / dist)
                    coords[i] -= (direction * push)
                    coords[j] += (direction * push)
                    moved = True

        if not moved:
            break

    layout = {node: coords[i]
        for i, node in enumerate(G.nodes())
    }

    fig, ax = plt.subplots(figsize=(10, 10))

    # Draw edges

    for u, v, data in G.edges(data=True):
        nx.draw_networkx_edges(G, layout,
            edgelist=[(u, v)],
            edge_color=edge_colour(data["shell"]),
            width=0.8,
            alpha=0.35,
            ax=ax)

    # Draw species nodes
    species_list = [s
        for s in ["WAT", "DMBA", "HAM", "HSO"]
        if s in set(nx.get_node_attributes(G, "species").values())
    ]

    for mol in species_list:
        nodes = [ n
            for n in G.nodes()
            if G.nodes[n]["species"] == mol
        ]

        nx.draw_networkx_nodes(G, layout, nodelist=nodes, node_size=NODE_SIZE_2D, node_color=COLOURS[mol], edgecolors="white",
                                linewidths=0.4, label=mol, ax=ax)

    # Highlight important nodes
    nx.draw_networkx_nodes(G, layout, nodelist=important, node_size=NODE_SIZE_2D+35, node_color="none",
                           edgecolors="black", linewidths=2, ax=ax)

    ax.set_title(title, fontsize=16,pad=15)
    ax.legend(frameon=False, fontsize=12)
    ax.set_aspect("equal")
    ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(filename,  dpi=800, facecolor="white", bbox_inches="tight")
    plt.close()
    print(
        "Saved:",filename)

# 3D NETWORK VISUALISATION
def plot_network_3D(G, positions, system, filename, title):
    bet, degree, eigen, important = calculate_centrality(G)

    # Edges
    edge_x = []
    edge_y = []
    edge_z = []

    edge_colours = []

    for i, j, data in G.edges(data=True):
        x1, y1, z1 = positions[i]
        x2, y2, z2 = positions[j]
        edge_x += [x1, x2, None]
        edge_y += [y1, y2,None]
        edge_z += [z1, z2, None]
        edge_colours.append(edge_colour(data["shell"]))

    edge_trace = go.Scatter3d( x=edge_x, y=edge_y, z=edge_z, mode="lines",
                              line=dict(color="black",width=1),
                              opacity=0.35, hoverinfo="none", showlegend=False)

    # Nodes
    node_traces = []
    present_species = set(nx.get_node_attributes(G, "species").values())

    for mol in ["WAT", "DMBA", "HAM", "HSO"]:
        if mol not in present_species:
            continue

        xs = []
        ys = []
        zs = []
        texts = []
        sizes = []

        for n in G.nodes():
            if G.nodes[n]["species"] != mol:
                continue

            x, y, z = positions[n]
            xs.append(x)
            ys.append(y)
            zs.append(z)
            texts.append(
                f"{G.nodes[n]['label']}<br>"
                f"Species = {mol}<br>"
                f"Betweenness = {bet[n]:.4f}<br>"
                f"Eigenvector = {eigen[n]:.4f}"
            )

            sizes.append(NODE_SIZE_3D + 30 * bet[n])

        node_traces.append(go.Scatter3d(x=xs, y=ys, z=zs,  mode="markers", name=mol,
                                        marker=dict(size=sizes, color=COLOURS[mol],
                                                    line=dict(color="black", width=0.5)),
                                                    text=texts, hoverinfo="text", showlegend=True))

    fig3d = go.Figure(data=[edge_trace] + node_traces)

    fig3d.update_layout(title=title, legend=dict(title="Species", x=1.02, y=1, bgcolor="rgba(255,255,255,0.85)",
                                                 bordercolor="black", borderwidth=1, font=dict(size=12)),
                                                scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False) , zaxis=dict(visible=False),aspectmode="data"),
                                                margin=dict(l=0, r=0, b=0, t=40))

    fig3d.write_html(filename)
    print("Saved:", filename)

# RUN SINGLE NETWORK
def analyse_network(pdb_file, system):
    name = (os.path.basename(pdb_file).replace(".pdb", ""))
    match = re.search(r"(\d+)h2o", name.lower())
    hydration = (match.group(1)
                 if match
                 else "?")
    title = (f"DMBA·HSO₄·{hydration}H₂O"
             if system == "DMBA"
             else f"HAM·HSO₄·{hydration}H₂O")
    output_folder = os.path.join(OUTPUT_DIR, system)
    os.makedirs(output_folder, exist_ok=True)
    G, positions = build_network(pdb_file, system)

    # 2D output
    plot_network_2D(G, positions, system,
                    os.path.join( output_folder,
                                  f"{name}_network_2D.png"), title)

    # 3D output
    plot_network_3D(G, positions, system,
                     os.path.join(output_folder,
                                   f"{name}_network_3D.html"), title)

# BATCH RUNNER
if __name__ == "__main__":
    for system, pdb_files in SYSTEMS.items():
        print(f"Running {system} networks")
        print("==============================")

        for pdb_file in pdb_files:
            print("\nProcessing:", pdb_file)
            analyse_network(pdb_file, system)

    print("\nAll network visualisations completed")

