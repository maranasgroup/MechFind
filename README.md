# MechFind: A Computational Framework for De Novo Prediction of Enzyme Mechanisms

## About the Project

MechFind addresses the "mechanism gap" in bioinformatics by providing a high-throughput tool to generate plausible mechanistic hypotheses. It operates in two main stages:
1.  **Parsimony Search:** Identifies the top-ten most parsimonious (fewest steps) mechanisms using a Mixed-Integer Linear Programming (MILP) formulation.
2.  **Similarity Re-ranking:** Re-ranks these ten candidates based on their mechanistic similarity to a curated database of known enzyme mechanisms from the M-CSA.

This repository provides a self-contained example to demonstrate the workflow on a sample reaction.

## Getting Started

To run the example, you will need a Python environment with several scientific computing packages installed.

### Prerequisites

*   Python 3.8+
*   Jupyter Notebook or JupyterLab
*   The following Python packages:
    *   `pulp` (for MILP)
    *   `pandas`
    *   `rdkit-pypi` (for cheminformatics)

### Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/MechFind.git
    cd MechFind
    ```

2.  **Create a `requirements.txt` file** with the following content:
    ```text
    pulp
    pandas
    rdkit-pypi
    ```

3.  **Install the required packages:**
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

2.  **Open the Notebook:** In the Jupyter interface that opens in your browser, click on `MechFind_example.ipynb`.

3.  **Run the Code:** Execute the cells in the notebook from top to bottom. The notebook is structured to:
    *   Import necessary libraries.
    *   Define all helper functions.
    *   Load and pre-process the required data files.
    *   Define the main `MechFind` function.
    *   Run an example prediction and display the results.

## Understanding the Files

*   **`MechFind_example.ipynb`**: The main file containing all the code. It defines the core `MechFind` function, its helpers, loads the necessary data, and runs a demonstration.
*   **`Unique_Rules.csv`**: The elementary rules matrix. This crucial data file contains the moiety changes for each of the 4,143 unique elementary reaction rules derived from the M-CSA database.
*   **`M-CSA_arrow_rules_r0.json`**: Pre-processed "arrow environment" data from the M-CSA. This file contains the fundamental chemical transformations that are used to calculate the similarity score for re-ranking mechanisms.

## Understanding the Output

The example in the notebook (`In [39]` and `In [40]`) produces two main outputs:

1.  **`solutions` (List of Mechanisms):**
    This is a list containing the top 10 predicted mechanisms, sorted by their similarity score (best first). Each mechanism is itself a list of rule IDs that, when combined, produce the overall reaction.

2.  **`Mechanism_Matrix` (DataFrame):**
    This function takes one of the predicted mechanisms (e.g., `solutions[0]`) and visualizes it as a table. The DataFrame shows the net change of each moiety for the overall reaction (`RXN`), the counts in the reactants and products, and the specific changes contributed by each elementary rule in the predicted mechanism. This is a powerful tool for analyzing a specific prediction in detail.
