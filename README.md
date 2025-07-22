<!-- markdownlint-disable MD033 -->

# AFragmenter

AFragmenter is a schema-free, tunable protein domain segmentation tool for AlphaFold structures based on network analysis.

### Key features

- **Schema free**: AFragmenter only uses the PAE values from AlphaFold structures. No domain-segmentation scheme is learned or used for evaluation.

- **Tunable segmentation**: The 'resolution' parameter gives control over the coarseness of clustering, and thus the number of clusers / domains.
  - **Higher resolution**: Yields more, smaller clusters
  - **Lower resolution**: Yields fewer, larger clusters

<br>

| Resolution = 0.8 | Resolution = 1.1 | Resolution = 0.3 |
|:----------------:|:----------------:|:----------------:|
| ![Resolution 0.8](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/resolution_0_8.png) | ![Resolution 1.1](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/resolution_1_1.png) | ![Resolution 0.3](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/resolution_0_3.png) |

<p style="text-align: right">protein: P15807&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</p>

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

## Table of contents

1. [Try it](#try-it)
2. [Installation](#installation)
3. [Tutorial](#quick-tutorial)
4. [Usage](#usage)
    - [Python](#python)
    - [Command line](#command-line)
5. [Options](#options)
    - [Threshold](#threshold)
    - [Resolution](#resolution)
    - [Objective function](#objective-function)
    - [Minimum size](#minimum-size-min_size)
    - [Merge](#merge)

## Try it

The recommended way to use AFragmenter is through jupyter notebooks, where visualization and fine-tuning of parameters is most easily done.
The easiest way to begin is by using our [Google colab notebook](https://colab.research.google.com/drive/1Fphc-SkTFsphBIL-i7Rfxg86p4znU7m4?usp=sharing).

- **Note**: While colab notebooks offers convenience, it can experience slower performance due to shared resources.

An alternative way to get started is by using the [webtool] (coming soon)

## Installation

### System Requirements

- **Python Version**: Python 3.9 or higher
- **Operating Systems**: Linux, macOS, or Windows.

### Installation

AFragmenter is available through PyPI and bioconda.

```bash
pip install afragmenter
``` 
or
```bash
conda install -c conda-forge -c bioconda afragmenter
```



#### **Optional Dependencies**:
   - **py3Dmol**: Required for protein structure visualization.

     ```bash
     pip install py3Dmol
     ```
     Or
     ```bash
     conda install -c conda-forge py3Dmol
     ```

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

<img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/P15807_pae.png" width=400 alt="P15807 PAE matrix"/>

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

<img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/P15807_thresholds.png" width=700 alt="Effect of different thresholds on edge weights">

A threshold of 3 seems to give a good contrast between the higher and lower PAE values.

Next, we cluster the residues into domains using the Leiden clustering algorithm. We get a result, but the resolution parameter can be changed to explore multiple potential solutions.

```python
p15807 = AFragmenter(pae, threshold=3)
result = p15807.cluster() # default resolution = 0.8
result.plot_result()
result.py3Dmol(structure)
```

<p style="display: flex" float="left">
    <img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/result_resolution_0_8.png" width=35.8% alt="PAE result, resolution = 0.8"/>
    <img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/resolution_0_8.png" width=50% alt="structure colored by resulting domain"/>
</p>

```python
p15807.cluster(resolution=1.1)
```

<p style="display: flex" float="left">
    <img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/result_resolution_1_1.png" width=35.8% alt="PAE result, resolution = 1.1"/>
    <img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/resolution_1_1.png" width=50% alt="structure colored by resulting domain"/>
</p>

```python
p15807.cluster(resolution=0.3)
```

<p style="display: flex" float="left">
    <img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/result_resolution_0_3.png" width=35.8% alt="PAE result, resolution = 0.3"/>
    <img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/resolution_0_3.png" width=50% alt="structure colored by resulting domain"/>
</p>

Once a solution has been found that is satisfactory to the user, we can print the result and the FASTA file for each domain, or save them to files for further analysis.

```python
p15807 = AFragmenter(pae, threshold=3)
result = p15807.cluster(resolution=1.1)
result.print_result()

    ┏━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━┓
    ┃ Domain ┃ Number of Residues ┃ Chopping ┃
    ┡━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━┩
    │ 1      │                137 │   10-146 │
    │ 2      │                 47 │  147-193 │
    │ 3      │                 81 │  194-274 │
    └────────┴────────────────────┴──────────┘

result.print_fasta(structure)

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
result.save_result('result.csv')
result.save_fasta(structure, 'result.fasta')
```

## Usage

- [Python](#python)
- [Command line](#command-line)

### Python

Docs coming soon...

### Command line

![help_message](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/help_message.svg)

## Options

- [Threshold](#threshold)
- [Resolution](#resolution)
- [Objective function](#objective-function)
- [Minimum size](#minimum-size-min_size)
- [Merge](#merge)

### <ins>Threshold</ins>

The 'contrast threshold' serves as a soft cut-off to increase the distinction between low and high PAE values. Used in calculating the edge weights of the network and will thus have a large impact on the clustering and segmentation results. It is important to consider this threshold in the context of the AlphaFold results for the protein of interest.



**Examples**:

#### [Q5VSL9](https://alphafold.ebi.ac.uk/entry/Q5VSL9)

Overall good structure with high pLDDt and low PAE scores for the majority of the protein, and lower pLDDT and high PAE scores for the disordered regions / loops, like is expected. Default threshold should be good (default PAE threshold = 5).

| AlphaFold structure | PAE plot | Edge weights |
|:-------------------:|:----------:|:------------:|
| ![Q5VSL9 AF structure](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/Q5VSL9/Q5VSL9.png) | ![Q5VSL9 PAE plot](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/Q5VSL9/Q5VSL9_pae.png) | ![Q5VSL9 edge weights](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/Q5VSL9/Q5VSL9_edge_weights.png) |

#### [P15807](https://alphafold.ebi.ac.uk/entry/P15807)

Very high pLDDT scores and low PAE scores for the AlphaFold structure indicating strong confidence, with one loop as exception. Several linkers, including the disordered N-terminal region, also show unexpectedly high pLDDT scores and low PAE scores, contrary to what would be expected for such regions. This apparent overconfidence is likely due to the inclusion of the crystal structure (1KYQ) in the AlphaFold training dataset.

Lowering the treshold can help reduce this apparent confidence, making it easier to differentiate between genuinely well-structured regions ans those that are more likely to be flexible or disordered.

| AlphaFold structure | PAE plot |
|:-------------------:|:--------:|
| ![P15807 AF structure](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/P15807.png) | <img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/P15807_pae.png" width=70% alt="P15807 PAE plot"> |
| <strong>Edge weights (default threshold = 5)</strong> | <strong>Edge weights (treshold = 3)</strong> |
| <img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/P15807_edge_weights.png" width=70% alt="edge weight using default threshold"> | <img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/P15807_edge_weights_3.png" width=67% alt="edge weights using PAE threshold = 3"> |

#### [Q9YFU8](https://alphafold.ebi.ac.uk/entry/Q9YFU8)

Q9YFU8 is a great example to remind us again that the PAE scores are not primarily intended to be used for domain segmentation,
but instead are a measure of how confident AlphaFold is in the relative position of two residues.

The PAE plot for Q9YFU8 shows two distinct parts of the protein seperated with hight PAE values, indicating uncertainty in their relative positions.
Going of off the previous examples, it would not be uncommon to assume there to be two distince domains in this protein, but this isn't necessarily the case.
Q9YFU8 has two crystal structures in the PDB: 1W5S and 1W5T. Superpositioning of these crystal structures reveals that a significant portion
of the protein overlays well, however another part shows a large deviation in orientation. AlphaFold likely learned this similarity and difference,
resulting in low PAE scores for the overlapping regions and high PAE scores between the differently oriented parts.
These structures might explain the resulting PAE scores, but this means we still need to pay attention choosing the threshold
to properly segment the remaining parts of the protein structure.

| AlphaFold structure | PAE plot | Crystal structures: 1W5S (green) and 1W5T (red) |
|:-------------------:|:--------:|:-----------------------------------------------:|
| ![Q9YFU8 AlphaFold structure](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/Q9YFU8/Q9YFU8.png) | <img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/Q9YFU8/Q9YFU8_pae.png" width=75% alt="Q9YFU8 PAE plot"> | ![1W5S (green) and 1W5T (red) structures](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/Q9YFU8/1w5s_1w5t.png) |

Lowering the treshold even if initial inspection deems it not necessary can still change the results. Without changing the threshold
we see two domains, consistent with the results from [SCOP](https://www.ebi.ac.uk/pdbe/scop/search?t=txt;q=1w5s) and [SCOPe](https://scop.berkeley.edu/pdb/code=1w5s). While lowering the threshold results in 3 domains, consistent with [ECOD](http://prodata.swmed.edu/ecod/af2_pdb/domain/e1w5sA1#tab-organization), [CATH](https://www.cathdb.info/pdb/1w5s), [Interpro](https://www.ebi.ac.uk/interpro/protein/reviewed/Q9YFU8/) and SCOP.
(SCOP can contain multiple solutions)

| Threshold = 5 | Threshold = 3 |
|:-------------:|:-------------:|
|<img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/Q9YFU8/Q9YFU8_two_domains.png" width=75% alt="Q9YFU8 two domains colored"> | <img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/Q9YFU8/Q9YFU8_three_domains.png" width=75% alt="Q9YFU8 three domains colored"> |

(Other settings kept as default values)

### <ins>Resolution</ins>

The **resolution** can be thought of as the coarseness of clustering. Increasing the resolution will result in more, smaller clusters (/domains). Decreasing the resolution will result in fewer but larger clusters.

- Default: dependent on objective_function `{"modularity": 0.7, "pm": 0.3}`

**Examples**:

#### [P15807](https://alphafold.ebi.ac.uk/entry/P15807)

| Resolution = 0.8 | Resolution = 1.1 | Resolution = 0.3 |
|:----------------:|:----------------:|:----------------:|
| ![Resolution 0.8](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/resolution_0_8.png) | ![Resolution 1.1](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/resolution_1_1.png) | ![Resolution 0.3](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P15807/resolution_0_3.png) |

#### [A0A098AQT8](https://alphafold.ebi.ac.uk/entry/A0A098AQT8)

| Resolution = 0.8 | Resolution = 1.4 |
|:----------------:|:----------------:|
|<img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/A0A098AQT8/A0A098AQT8_resolution_0_8.png" width=75% alt="Resolution 0.8"> | <img src="https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/A0A098AQT8/A0A098AQT8_resolution_1_4.png" width=75% alt="Resolution 1.4"> |

### <ins>Objective function</ins>

The objective function that is optimized during clustering, choices are either CPM (constant potts model) or Modularity. The contant potts model does not suffer from the resolution limit problem like modularity does, leading to more, smaller well-defined clusters.
This means that 'CPM' will result in more smaller, tightly connected clusters that represent specific subgroups or communities within the data. On the other hand, 'Modularity' will tend to produce fewer, larger clusters that encompass broader groups within the data.

For AFragmenter, 'CPM' translates to a more sensitive approach where we see many more smaller clusters, especially for disordered regions. 'Modularity' is less sensitive to small shifts in PAE values, and will be better at clustering residues from disordered regions together.

- Default: `modularity`
- Options: `[modularity, cpm]`

**Examples**:

#### [P04637](https://alphafold.ebi.ac.uk/entry/P04637)

![Result plots P04637](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P04637/PAE_results.png)
![Result structures P04637](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/P04637/structure_results.png)

#### [Q837X5](https://alphafold.ebi.ac.uk/entry/Q837X5)

![Result plots Q837X5](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/Q837X5/PAE_results.png)
![Result structures Q837X5](https://raw.githubusercontent.com/sverwimp/afragmenter/master/images/Q837X5/structure_results.png)

### <ins>Minimum size</ins> (`min_size`)

The `min_size` parameter specifies the minimum number of residues required for a cluster / domain to be considered valid. Clusters below this threshold are filtered out during post-processing. This helps in reducing noise and focusing on significant structural or functional regions.

- Default value: 10
- Valid range: 0 ≤ min_size ≤ Number of Residues

A higher min_size is often used for larger proteins to focus on major domains, while a lower value allows capturing smaller but potentially important regions. Adjusting this parameter can impact the granularity of your clustering results.

For examples, refer to the examples in [Objective function](#objective-function).

### <ins>Merge</ins>

The `merge` parameter, also referred to as `attempt_merge`, plays a important role in refining the clustering process. 
When enabled, it attempts to merge smaller clusters (= below min_size) with adjacent larger ones, ensuring that the resulting clusters are more meaningful and less fragmented. Resulting clusters below the min_size are first attempted to be merged with adjacent larger clusters. If merging is not possible, the small clusters are filtered out.

- Default value: `True`
