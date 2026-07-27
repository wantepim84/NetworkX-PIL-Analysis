</head>

<body>

<h1>Scripts Documentation</h1>

<p>
This directory contains Python scripts used for the construction,
analysis, and visualisation of molecular interaction networks for hydrated
protic ionic liquid (PIL) systems.
</p>

<p>
The workflow uses graph-based analysis to investigate the structural
organisation of DMBA-HSO<sub>4</sub> and HA-HSO<sub>4</sub> systems at
different hydration levels.
</p>


<h2>Workflow Overview</h2>

<ol>
<li>
Construct molecular interaction networks from molecular structure files.
</li>

<li>
Calculate network topology and centrality descriptors.
</li>

<li>
Generate two-dimensional and three-dimensional network visualisations.
</li>
</ol>


<h2>Scripts</h2>

<table>

<tr>
<th>Script</th>
<th>Description</th>
</tr>

<tr>
<td><code>network_analysis.py</code></td>
<td>
Builds molecular interaction networks and calculates graph-based
properties including degree, clustering coefficient, centrality measures,
network density, k-core structure, assortativity, and percolation.
</td>
</tr>


<tr>
<td><code>network_visualisation.py</code></td>
<td>
Generates two-dimensional and interactive three-dimensional network
representations using NetworkX, Matplotlib, and Plotly.
</td>
</tr>

</table>


<h2>Network Construction</h2>

<p>
Molecules are represented as graph nodes, while intermolecular interactions
are represented as edges. Interaction criteria are defined using
distance cutoffs obtained from radial distribution function (RDF) analysis.
</p>

<p>
The scripts use:
</p>

<ul>
<li>MDAnalysis for molecular structure processing</li>
<li>NetworkX for graph construction and analysis</li>
<li>Matplotlib for 2D visualisation</li>
<li>Plotly for interactive 3D visualisation</li>
</ul>


<h2>Calculated Network Properties</h2>

<p>
The analysis scripts calculate the following graph descriptors:
</p>

<ul>
<li>Average degree</li>
<li>Degree variance</li>
<li>Network density</li>
<li>Clustering coefficient</li>
<li>Betweenness centrality</li>
<li>Closeness centrality</li>
<li>Eigenvector centrality</li>
<li>K-core values</li>
<li>Assortativity</li>
<li>Percolation behaviour</li>
<li>Species-resolved network properties</li>
</ul>


<h2>Usage</h2>

<p>
Run scripts from the repository root directory:
</p>

<pre>
python scripts/network_analysis.py

python scripts/network_visualisation.py
</pre>


<h2>Input Data</h2>

<p>
The scripts require molecular structure files generated from molecular
simulation trajectories. Interaction cutoffs used during network generation
are provided from RDF analysis.
</p>

<p>
The repository contains processed network outputs for:
</p>

<ul>
<li>DMBA-HSO<sub>4</sub> systems</li>
<li>HA-HSO<sub>4</sub> systems</li>
<li>Hydration states from 1H<sub>2</sub>O to 6H<sub>2</sub>O</li>
</ul>


<h2>Output Files</h2>

<table>

<tr>
<th>File</th>
<th>Description</th>
</tr>

<tr>
<td><code>network_metrics.csv</code></td>
<td>
Contains calculated network descriptors.
</td>
</tr>

<tr>
<td><code>network_2D.png</code></td>
<td>
Static two-dimensional network visualisation.
</td>
</tr>

<tr>
<td><code>network_3D.html</code></td>
<td>
Interactive three-dimensional network visualisation.
</td>
</tr>

</table>


<h2>Requirements</h2>

<p>
Install required Python packages using:
</p>

<pre>
pip install -r requirements.txt
</pre>


<h2>Notes</h2>

<p>
The scripts were developed for analysing molecular interaction networks in
hydrated protic ionic liquids. RDF-derived distance cutoffs define molecular
connectivity and network properties describe changes in molecular
organisation with increasing hydration.
</p>

</body>
</html>
