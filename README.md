# MechFind: A Computational Framework for De Novo Prediction of Enzyme Mechanisms

## About the Project

MechFind addresses the "mechanism gap" in bioinformatics by providing a high-throughput tool to generate plausible mechanistic hypotheses. It operates in two main stages:
1.  **Parsimony Search:** Identifies the top-ten most parsimonious (fewest steps) mechanisms using a Mixed-Integer Linear Programming (MILP) formulation.
2.  **Similarity Re-ranking:** Re-ranks these ten candidates based on their mechanistic similarity to a curated database of known enzyme mechanisms from the M-CSA.

This repository provides the source code, necessary datasets, and a self-contained example to demonstrate the workflow on a sample reaction.

## Getting Started

To run the example, you will need a Python environment with the necessary scientific computing packages.

### Prerequisites

*   Python 3.8+
*   Jupyter Notebook or JupyterLab
*   Standard scientific libraries (numpy, pandas, rdkit, pulp, tqdm)

### Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/maranasgroup/MechFind.git
    cd MechFind
    ```

2.  **Install the required packages:**
    ```bash
    pip install -r requirements.txt
    ```

## How to Run the Example

The primary way to use this repository is through the provided Jupyter Notebook.

1.  **Launch Jupyter:** Open your terminal, navigate to the `MechFind` directory, and run:
    ```bash
    jupyter notebook
    ```
    or
    ```bash
    jupyter lab
    ```

2.  **Open the Notebook:** In the Jupyter interface that opens in your browser, click on **`MechFind_example.ipynb`**.

3.  **Run the Code:** Execute the cells in the notebook. The notebook is structured to:
    *   Import the core logic from the `MechFind` module.
    *   Load the required data files.
    *   Run an example prediction on a specific reaction.
    *   Display the results as a matrix of moiety changes.

## Understanding the Files

### Core Software
*   **`MechFind.py`**: The core Python module containing the main algorithm. It defines the `MechFind` function, optimization formulations (`minRules`, `OrderRules`), and helper functions for SMILES processing and similarity scoring.
*   **`MechFind_example.ipynb`**: A demonstration notebook that imports functions from `MechFind.py` and walks through a single prediction example.
*   **`Unq_Rule_Gen.py`**: A utility script used to generate the elementary rules matrix at any integer radius. It processes the raw M-CSA data files to create the CSVs found in the `Unique_Rules` folder.

### Data Files and Folders
*   **`Unique_Rules/`**: A directory containing the generated elementary rule sets. By default, it contains the radius-1 moiety-based rules (`Unique_Rules_1.csv`).
*   **`M-CSA_arrow_rules_r0.json`**: Pre-processed "arrow environment" data from the M-CSA. This file contains the fundamental chemical transformations used to calculate the similarity score for re-ranking.
*   **`MCSA_mols_rMAX_1005.json`**: Contains the reaction SMILES strings for the curated M-CSA overall reactions. Used by `Unq_Rule_Gen.py` to generate rules.
*   **`MCSA_rules_rMAX_1005.json`**: Contains the reaction SMILES strings for the specific mechanistic steps in the M-CSA. Used by `Unq_Rule_Gen.py` to generate rules.

## Understanding the Output

The example in the notebook produces two main outputs:

1.  **`solutions` (List of Mechanisms):**
    This is a list containing the top 10 predicted mechanisms, sorted by their similarity score (best first). Each mechanism is itself a list of rule IDs that, when combined, produce the overall reaction.

2.  **`Mechanism_Matrix` (DataFrame):**
    This function takes one of the predicted mechanisms (e.g., `solutions[0]`) and visualizes it as a table. The DataFrame shows the net change of each moiety for the overall reaction (`RXN`), the counts in the reactants and products, and the specific changes contributed by each elementary rule in the predicted mechanism.

## License

**Software License Agreement**

By downloading or using the Penn State MechFind software you agree to the following terms of use:

Copyright (c) 2025, The Pennsylvania State University All rights reserved.

The Software may only be used for non-profit, non-commercial purposes. Please email ottinfo@psu.edu for any queries related to licensing for commercial use. Redistribution and use in source and binary forms, with or without modification, are not permitted without prior written approval.

Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF THIS SOFTWARE.
