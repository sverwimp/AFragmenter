# AFragmenter

AFragmenter is a schema-free, tunable protein domain segmentation tool for AlphaFold structures based on network analysis.


### How it works
1. **Network representation**: Each protein residue is treated as a node within a fully connected network
2. **Edge weighting**: The edges between the nodes are weighted using transformed **Predicted Aligned Error** (PAE) values from AlphaFold, reflecting relative positional confidence between residues.
    <details>
    <summary>Details on the use of PAE values</summary>

    - PAE values show differences when looking between inter- versus intra-domain residue pairs.
    - Intra-domain residue paris are expected to have lower PAE values compared to inter-domain residue pair.
    - This difference is used to distinguish well-structured regions within a protein structure from other well-structured regions and from poorly structured regions.
    - This enables us to cluster protein residue pairs of well-structured regions together.

    </details>

3. **Clustering with Leiden algorithm**: Utilizes the Leiden clustering algorithm to group residues into domains, with adjustable resolution parameters to control cluster granularity.

### Key features
- **Schema free**: AFragmenter only uses the PAE values from AlphaFold structures. No domain-segmentation scheme is learned or used for evaluation.

- **Tunable segmentation**: The 'resolution' parameter gives control over the coarseness of clustering, and thus the number of clusers / domains.
    - **Higher resolution**: Yields more, smaller clusters
    - **Lower resolution**: Yields fewer, larger clusters

<br>

<div style="display: flex; justify-content: space-around">
    <p>Resolution = 0.8</p>
    <p>Resolution = 1.1</p>
    <p>Resolution = 0.3</p>
</div>

<p style="display: flex; justify-content: center">
    <img src="images/resolution_0_8.png" width=33% />
    <img src="images/resolution_1_1.png" width=33% />
    <img src="images/resolution_0_3.png" width=33% />
</p>

<p style="display: flex; justify-content: right">protein: P15807&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</p>



## Table of contents

1. [Try it](#try-it)
2. [Installation](#installation)
3. [Overview usage](#overview-usage)
   1. [Python](#python)
   2. [Command line](#command-line)

## Try it

The recommended way to use AFragmenter is through jupyter notebooks, where visualization and fine-tuning of parameters is most easily done.
The easiest way to begin is by using our [Google colab notebook](https://colab.research.google.com/drive/1QQ3MO0kaTrJxD9EH1jghProsAe3Kmoru?usp=sharing).
- **Note**: While colab notebooks offers convenience, it can experience slower performance due to shared resources.

An alternative way to get started is by using the [webtool] (coming soon)

## Installation

### System Requirements
- **Python Version**: Ensure you have Python 3.9 or higher installed on your system.
- **Operating Systems**: The tool is compatible with Linux, macOS, and Windows.

### Installation Steps

1. **Set Up a Virtual Environment (Recommended)**:
   Creating a virtual environment helps manage dependencies effectively. Here's how to set it up:

   ```bash
   # Install virtualenv if not already installed
   pip install virtualenv
   
   # Create a new virtual environment
   virtualenv myenv
   
   # Activate the virtual environment
   # On Windows:
   myenv\Scripts\activate
   # On macOS/Linux:
   source myenv/bin/activate
   ```

   or alternatively, create and use a conda environment `conda create --name myenv pip 'python>=3.9'` & `conda activate myenv`

2. **Install AFragmenter**:
   Install the package using pip within your activated virtual environment.

   ```bash
   pip install AFragmenter
   ```

3. **Optional Dependencies**:
   - **py3Dmol**: Required for protein structure visualization.
     
     ```bash
     pip install py3Dmol
     ```

### Verification

After installation, verify that AFragmenter is correctly installed by running:

```bash
afragmenter --version
```

This command should display the installed version of AFragmenter.

## Overview usage

- [Python](#python)
- [Command line](#command-line)


### Python

Since AFragmenter is dependent on the PAE values of AlphaFold, it is a good idea to first have a look at the PAE plot.

```python
from afragmenter import AFragmenter, fetch_afdb_data

pae, structure = fetch_afdb_data('P15807')
p15807 = AFragmenter(pae) # Or bring your own files: a = AFragmenter('filename.json')
p15807.plot_pae()
```
<img src="images/P15807_pae.png" width=400 />

Here we see some collections of very low PAE values (dark green) on the PAE matrix, which could indicate different domains. However, there is still a lot of very green datapoints visible around these potential domains. So we might need to pay attention the the **PAE contrast threshold** that is used.

These PAE values are transformed into edge weights in a way that increases the contrast between high and low PAE values. The **PAE contrast theshold** that is used here can be changed. Below we can see the effect of the threshold on the weights of the graph.

<details>
<summary>Show code</summary>

```python
from afragmenter.plotting import plot_matrix

p15807 = AFragmenter(pae)
fig, ax = plt.subplots(2, 2, figsize=(10, 10))

p15807.plot_pae(ax=ax[0, 0])
plot_matrix(p15807.edge_weights_matrix, ax=ax[0, 1])

p15807 = AFragmenter(pae, threshold=3)
plot_matrix(p15807.edge_weights_matrix, ax=ax[1, 0])

p15807 = AFragmenter(pae, threshold=1)
plot_matrix(p15807.edge_weights_matrix, ax=ax[1, 1])

ax[0, 0].set_title('PAE matrix')
ax[0, 1].set_title('Edge weights matrix (threshold=5)\n[default]')
ax[1, 0].set_title('Edge weights matrix (threshold=3)')
ax[1, 1].set_title('Edge weights matrix (threshold=1)')

plt.tight_layout()
plt.show()
``` 
</details>

<img src="images/P15807_thresholds.png" width=700>

A threshold of 3 seems to give a good contrast between the higher and lower PAE values.

```python
p15807 = AFragmenter(pae, threshold=3)
p15807.cluster() # default resolution = 0.8
p15807.plot_result()
p15807.py3Dmol(structure)
```

<p style="display: flex" float="left">
    <img src="images/result_resolution_0_8.png" width=35.8% />
    <img src="images/resolution_0_8.png" width=50% />
</p>



# Continue working from here


```python
p15807.cluster(resolution=1.1)
```

<p style="display: flex" float="left">
    <img src="images/result_resolution_1_1.png" width=35.8% />
    <img src="images/resolution_1_1.png" width=50% />
</p>

```python
p15807.cluster(resolution=0.3)
```

<p style="display: flex" float="left">
    <img src="images/result_resolution_0_3.png" width=35.8% />
    <img src="images/resolution_0_3.png" width=50% />
</p>


```python
a = AFragmenter(pae, threshold=3)
a.cluster(resolution=1.1)
a.print_result()

    ┏━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━┓
    ┃ Domain ┃ Number of Residues ┃ Chopping ┃
    ┡━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━┩
    │ 1      │                137 │   10-146 │
    │ 2      │                 47 │  147-193 │
    │ 3      │                 81 │  194-274 │
    └────────┴────────────────────┴──────────┘

a.print_fasta(structure)

    >P15807_1 10-146
    QLKDKKILLIGGGEVGLTRLYKLIPTGCKLTLVSPDLHKSIIPKFGKFIQNEDQPDYRED
    AKRFINPNWDPTKNEIYEYIRSDFKDEYLDLEDENDAWYIIMTCIPDHPESARIYHLCKE
    RFGKQQLVNVADKPDLC
    >P15807_2 147-193
    DFYFGANLEIGDRLQILISTNGLSPRFGALVRDEIRNLFTQMGDLAL
    >P15807_3 194-274
    EDAVVKLGELRRGIRLLAPDDKDVKYRMDWARRCTDLFGIQHCHNIDVKRLLDLFKVMFQ
    EQNCSLQFPPRERLLSEYCSS

a.py3Dmol(structure)
```





### Command line

![help_message](images/help_message.svg)