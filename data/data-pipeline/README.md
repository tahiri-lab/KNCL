
# Data pipeline for constructing datasets of biological overlapping phylogenetic trees

## 📋 Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Script Workflow](#script-workflow)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Acknowledgments](#acknowledgments)
- [License](#license)

---

## Overview

This Python script automates the process of:

- Selecting a specified number of species from a given species group.
- Creating overlapping subsets based on specified parameters.
- Retrieving phylogenetic tree data from the [VertLife](https://vertlife.org/data/) website.
- Converting the downloaded Nexus files to Newick format.
- Compiling all the data into a final dataset.

This project includes two main files that provide different ways to interact with and use the code:
1. `Biological_data_preparation.ipynb` - A Jupyter Notebook version with enhanced interactivity and additional documentation.
2. `overlap_tree_pipeline.py` - A streamlined Python script version for non-interactive use.


### `Biological_data_preparation.ipynb`
- **Description**: This Jupyter Notebook version provides an interactive way to explore and use the code.
- **Purpose**: Designed for users who want to run the code interactively and gain a deeper understanding through additional comments, suggestions, and visualizations.
- **Usage**: Open the notebook in Jupyter Notebook, JupyterLab, or Colab to start an interactive session.
- **Features**: 
  - Includes the same core functionality as `overlap_tree_pipeline.py`.
  - Provides additional comments and explanations to guide the user.
  - Includes cells for visualizing datasets, making it easier to understand the data overlapping distributions and analysis process.
  - Allows users to modify parameters and see the results in real-time.

---

## Documentation for `overlap_tree_pipeline.py`

## Prerequisites

Before running the script, ensure you have the following:

- **Python 3.6+** installed.
- Required Python packages:
  - `pandas`
  - `selenium`
  - `requests`
  - `biopython`
- **Chrome browser** installed.
- **ChromeDriver** compatible with your Chrome browser version.

---

## Installation

1. **Clone or download the script** to your local machine.

2. **Install the required Python packages**:

   ```bash
   pip install pandas selenium requests biopython
   ```

3. **Install ChromeDriver**

---

## Usage

Run the script using the command line:

```bash
python overlap_tree_pipeline.py <species_group> <n> <number_of_trees> <email>
```

- `<species_group>`: One of `amphibians`, `birds`, `mammals`, `sharks`, `squamates`.
- `<n>`: Number of unique species (integer).
- `<number_of_trees>`: Total number of trees in the resulting dataset (between 10 and 1000, divisible by 10).
- `<email>`: Your email address (used by the VertLife website to send data).

**Example**:

```bash
python overlap_tree_pipeline.py birds 135 70 your@email.com
```

---

## Script workflow

Here's an overview of what the script does:

1. **Parses command-line arguments** and validates them.
2. **Loads species data** from `all_species_lists.csv`.
3. **Selects a random sample** of `n` species from the specified group.
4. **Creates overlapping subsets** based on specified parameters.
5. **Automates data retrieval** from the VertLife website.
6. **Converts Nexus files to Newick format**.
7. **Compiles the final dataset**.
8. **Provides the necessary citations** for the data used.

---

## Examples

### Example command

```bash
python overlap_tree_pipeline.py mammals 95 100 your@email.com
```

- **Species Group**: mammals
- **Number of Unique Species (n)**: 95
- **Number of Trees**: 100
- **Email**: your@email.com

### Expected output

- `selected_species.csv`: Contains 95 randomly selected mammal species.
- `Mammals_overlapping_subsets.csv`: Contains the 10 overlapping subsets.
- Downloaded Nexus files in `mammals_nexus/`.
- Final dataset `overlapping_dataset_mammals.txt`.
- Citation information printed to the console.

Once you have the data, be sure to cite the source:

* **Amphibians**: Jetz, W., & Pyron, R. A. (2018). The interplay of past diversification and evolutionary isolation with present imperilment across the amphibian tree of life. Nature ecology & evolution, 2(5), 850-858.

* **Birds**: Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012). The global diversity of birds in space and time. Nature, 491(7424), 444-448.

* **Mammals**: Upham, N. S., Esselstyn, J. A., & Jetz, W. (2019). Inferring the mammal tree: species-level sets of phylogenies for questions in ecology, evolution, and conservation. PLoS biology, 17(12), e3000494.

* **Sharks**: Stein, R. W., Mull, C. G., Kuhn, T. S., Aschliman, N. C., Davidson, L. N., Joy, J. B., ... & Mooers, A. O. (2018). Global priorities for conserving the evolutionary history of sharks, rays and chimaeras. Nature ecology & evolution, 2(2), 288-298.

* **Squamates**: Tonini, J. F. R., Beard, K. H., Ferreira, R. B., Jetz, W., & Pyron, R. A. (2016). Fully-sampled phylogenies of squamates reveal evolutionary patterns in threat status. Biological Conservation, 204, 23-31.

* Website: https://vertlife.org/data/

---

## Troubleshooting

### 🐞 Common issues and solutions

#### 1. `all_species_lists.csv` Not Found

**Solution**: Ensure that the `all_species_lists.csv` file is present in the same directory as the script.

#### 2. ChromeDriver errors

**Solution**:

- Verify that ChromeDriver is installed and matches your Chrome browser version.
- Try `pip install chromedriver-binary`

#### 3. Selenium errors during tree requests

**Solution**:

- Check your internet connection.
- Increase the explicit wait times in the `submit_tree_request` function.
- Ensure the VertLife website is accessible.

#### 4. Download errors for subsets 9 and 10

**Solution**:

- The script includes a 60-second delay before downloading to mitigate this issue.
- If problems persist, consider increasing the delay.

#### 5. Incompatible Python version

**Solution**: Use Python 3.6 or higher.

---

## Acknowledgments

- **VertLife project**: For providing access to phylogenetic tree data.
- **Authors of the datasets**: Please refer to the citations provided for each species group.

---

## License

This script is provided under the GPL-3.0 license.

---

If you have any questions or need further assistance, feel free to reach out at Nadia.Tahiri@USherbrooke.ca.
