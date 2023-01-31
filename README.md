# MolCompassViewer
`MolCompassViewer` -- is a tool for the visualisation of chemical space, visual validation of QSAR/QSPR models and the analysis of chemical diversity. 
The tool is based on the molcomplib library. The core technique is the pretrained parametric t-SNE model for mapping of chemical space. The tool is implemented in Python and uses 
`plotly` `dash` ad `dash-bootstrap-components` for the web interface.
# Installation
To install run:
`pip install git+https://github.com/sergsb/mcv.git`
# Usage
At present, MCV supports binary classification analysis only. To use the tool, type the following after installation:

`mcv <path to your csv file>`

The CSV file must contain the following columns (other synonyms for each column are listed below):
* `smiles` -- SMILES representation of the molecule (mandatory), Synonyms: `smi`, `SMILES`, `Smiles`
* `labels` -- Ground truth labels (optional,reqired for visual validation), Synonyms: `label`, `target`, `y`
* `probabilities` -- Predicted probabilities (optional, required for visual validation) Synonyms: `probs`, `prob`, `p`,`preds`, `pred`
 
If no labels and probabilities are provided, the tool will only run in chemical space demonstration mode, and the visual validation mode will be disabled.



