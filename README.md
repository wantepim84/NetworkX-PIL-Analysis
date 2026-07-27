<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>NetworkX-Hydrated-PILs</title>

<style>

body{
    font-family:Arial, Helvetica, sans-serif;
    max-width:1000px;
    margin:auto;
    padding:40px;
    background:#f7f7f7;
    color:#222;
    line-height:1.7;
}

h1{
    color:#004c78;
    border-bottom:3px solid #004c78;
    padding-bottom:10px;
}

h2{
    color:#004c78;
    margin-top:40px;
}

.card{
    background:white;
    padding:25px;
    border-radius:10px;
    box-shadow:0 2px 8px rgba(0,0,0,0.08);
    margin-bottom:30px;
}

table{
    width:100%;
    border-collapse:collapse;
    margin-top:20px;
}

th{
    background:#004c78;
    color:white;
    padding:12px;
}

td{
    padding:12px;
    border-bottom:1px solid #dddddd;
    text-align:center;
}

ul{
    margin-left:20px;
}

code{
    background:#eeeeee;
    padding:2px 5px;
    border-radius:4px;
}

.footer{
    text-align:center;
    margin-top:60px;
    color:#666;
    font-size:14px;
}

</style>

</head>

<body>

<h1>NetworkX-Hydrated-PILs</h1>

<p>
Supplementary repository accompanying the manuscript:
</p>

<p>
<b>NetworkX: A convenient tool for graph analysis of IL solution mixtures</b>
</p>

<div class="card">

<h2>Overview</h2>

<p>

This repository contains the Python scripts, graph-theoretical datasets,
and interactive molecular interaction networks supporting the published study.
Molecular interaction networks were constructed using RDF-derived interaction
cutoffs and analysed using NetworkX to investigate how hydration reorganises
the solvent environment of two protic ionic liquids.

</p>

</div>

<div class="card">

<h2>Repository Contents</h2>

<ul>

<li><b>Python analysis scripts</b> for network construction and graph-theoretical analysis.</li>

<li><b>Python visualisation scripts</b> for generating publication-quality two-dimensional and interactive three-dimensional molecular interaction networks.</li>

<li><b>Raw graph-theoretical datasets</b> generated from the NetworkX analysis.</li>

<li><b>Interactive 3D molecular interaction networks</b> for every hydration state.</li>

<li><b>2D molecular interaction networks</b>.</li>

<li><b>Supporting documentation</b>.</li>

</ul>

</div>

<div class="card">

<h2>Available Systems</h2>

<table>

<tr>

<th>System</th>

<th>Hydration States</th>

<th>Interactive 3D</th>

<th>2D Network</th>

<th>Analysis Data</th>

</tr>

<tr>

<td>DMBA</td>

<td>1–6 H₂O</td>

<td>✓</td>

<td>✓</td>

<td>✓</td>

</tr>

<tr>

<td>HA</td>

<td>1–6 H₂O</td>

<td>✓</td>

<td>✓</td>

<td>✓</td>

</tr>

</table>

</div>

<div class="card">

<h2>Representative Hydration States</h2>

<p>

The manuscript presents representative network visualisations for
<b>1</b>, <b>3</b>, and <b>6</b> water molecules per ion pair.
These systems represent the low-, intermediate-, and high-hydration
regimes discussed in the paper.

</p>

<p>

The repository contains the complete hydration series from
<b>1–6 H₂O</b> per ion pair for both protic ionic liquids.

</p>

</div>

<div class="card">

<h2>Contents of Each Hydration-State Folder</h2>

<ul>

<li><code>network_analysis.csv</code> – Graph-theoretical descriptors calculated from the production trajectory.</li>

<li><code>network_2D.png</code> – Two-dimensional molecular interaction network.</li>

<li><code>network_3D.html</code> – Interactive three-dimensional molecular interaction network.</li>

</ul>

</div>

<div class="card">

<h2>Software</h2>

<ul>

<li>Python</li>

<li>MDAnalysis</li>

<li>NetworkX</li>

<li>NumPy</li>

<li>SciPy</li>

<li>Matplotlib</li>

<li>Plotly</li>

</ul>

<p>

Install all required packages using:

</p>

<pre><code>pip install -r requirements.txt</code></pre>

</div>

<div class="card">

<h2>Citation</h2>

<p>

If you use the scripts, datasets, or visualisations provided in this repository,
please cite the accompanying publication and the archived repository DOI.

</p>

<p>

<b>DOI:</b> To be added.

</p>

</div>

<div class="footer">

NetworkX-Hydrated-PILs • Supplementary Repository

</div>

</body>

</html>
