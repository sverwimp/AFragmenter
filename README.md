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

| Resolution = 0.8 | Resolution = 1.1 | Resolution = 0.3 |
|:----------------:|:----------------:|:----------------:|
| ![Resolution 0.8](images/resolution_0_8.png) | ![Resolution 1.1](images/resolution_1_1.png) | ![Resolution 0.3](images/resolution_0_3.png) |

<p style="text-align: right">protein: P15807&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</p>



## Table of contents

1. [Try it](#try-it)
2. [Installation](#installation)
3. [Tutorial](#quick-tutorial)
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

## Quick Tutorial

In this short tutorial, we will walk through the process of using AFragmenter to segment protein domains based on AlphaFold structures. We will use the example protein P15807 (PDB: 1KYQ) to demonstrate the steps involved. This protein is classified differently by various protein domain databases, making it an interesting case for domain segmentation.

> P15807 is classified as a three-domain protein in both [CATH](https://www.cathdb.info/pdb/1kyq) and [ECOD](http://prodata.swmed.edu/ecod/af2_pdb/domain/e1kyqA1#tab-organization), a two-domain protein in [SCOPe](https://scop.berkeley.edu/pdb/code=1kyq) and [InterPro](https://www.ebi.ac.uk/interpro/protein/UniProt/P15807/), and a single-domain protein in [SCOP](https://www.ebi.ac.uk/pdbe/scop/term/8001033).

<br>

Since AFragmenter is dependent on the PAE values of AlphaFold, it is a good idea to first have a look at the PAE plot.

```python
from afragmenter import AFragmenter, fetch_afdb_data

pae, structure = fetch_afdb_data('P15807')
p15807 = AFragmenter(pae) # Or bring your own files: a = AFragmenter('filename.json')
p15807.plot_pae()
```
<img src="images/P15807_pae.png" width=400 />

Here we see some regions of very low PAE values (dark green) on the PAE matrix, which could indicate different domains. However, there are still many green (low PAE) datapoints visible around these potential domains. Therefore, it is important to consider the **PAE contrast threshold** used.

These PAE values are transformed into edge weights to increase the contrast between high and low PAE values. The **PAE contrast threshold** can be adjusted to control this contrast. Below, we can see the effect of different thresholds on the weights of the graph.



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

Next, we cluster the residues into domains using the Leiden clustering algorithm. We get a result, but the resolution parameter can be changed to explore multiple potential solutions.

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

Once a solution has been found that is satisfactory to the user, we can print the result and the FASTA file for each domain, or save them to files for further analysis.

```python
p15807 = AFragmenter(pae, threshold=3)
p15807.cluster(resolution=1.1)
p15807.print_result()

    ┏━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━┓
    ┃ Domain ┃ Number of Residues ┃ Chopping ┃
    ┡━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━┩
    │ 1      │                137 │   10-146 │
    │ 2      │                 47 │  147-193 │
    │ 3      │                 81 │  194-274 │
    └────────┴────────────────────┴──────────┘

p15807.print_fasta(structure)

    >P15807_1 10-146
    QLKDKKILLIGGGEVGLTRLYKLIPTGCKLTLVSPDLHKSIIPKFGKFIQNEDQPDYRED
    AKRFINPNWDPTKNEIYEYIRSDFKDEYLDLEDENDAWYIIMTCIPDHPESARIYHLCKE
    RFGKQQLVNVADKPDLC
    >P15807_2 147-193
    DFYFGANLEIGDRLQILISTNGLSPRFGALVRDEIRNLFTQMGDLAL
    >P15807_3 194-274
    EDAVVKLGELRRGIRLLAPDDKDVKYRMDWARRCTDLFGIQHCHNIDVKRLLDLFKVMFQ
    EQNCSLQFPPRERLLSEYCSS


# Or save it
p15807.save_result('result.csv')
p15807.save_fasta(structure, 'result.fasta')
```


## Overview usage

- [Python](#python)
- [Command line](#command-line)


### Python


### Command line

![help_message](images/help_message.svg)