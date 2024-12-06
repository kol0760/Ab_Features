



# AbFeatures

[AbRFC](https://github.com/tbc01/AbRFC) is a model that uses a Random Forest classifier to predict the impact of protein mutations on antibody affinity. However, the authors have not fully disclosed the details of feature construction in their dataset. This project aims to build a quantitative feature set for antigen-antibody interactions, based on the previous publications by the original research group, along with my own understanding and interpretation.



## Project Requirements

### Python Environment Setup

All calculations and operations in this project are based on Python. It is recommended to use **conda** to create and manage a virtual environment for ease of use and to ensure package compatibility. You can create a new environment using the following command:

```
conda create -n project_env python=3.8
conda activate project_env
```

### Required Packages

The project relies on several scientific computing and data analysis libraries. Below is the list of required packages and their installation instructions.

#### 1. Scientific Computing Libraries

- **PyRosetta**

  - **PyRosetta** is the core of this project. It provides tools for macromolecular modeling and analysis.
  - Install **PyRosetta** using the following commands:

  ```
  pip install pyrosetta-installer
  python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
  ```

  - For more information, refer to the [PyRosetta download page](https://www.pyrosetta.org/downloads/).

- **SciPy**

- **Biopython**

#### 2. Data Analysis Libraries

- **Pandas**
- **NumPy**



