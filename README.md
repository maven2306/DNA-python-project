# DNA Sequence Assembler using De Bruijn Graphs

This project reconstructs a full DNA sequence from fragmented, overlapping readings. It cleans the raw data, builds a De Bruijn graph to represent the overlaps, and finds an Eulerian path through the graph to assemble the final sequence. The resulting graph is also visualized and saved as an image.

## Features

- **Data Cleaning**: Handles duplicate readings, missing data points, and invalid nucleotide entries.
- **De Bruijn Graph Construction**: Builds a `networkx` directed graph from k-mers derived from the DNA fragments.
- **Sequence Reconstruction**: Traverses the graph to find an Eulerian path, effectively reassembling the original DNA sequence.
- **Graph Visualization**: Generates and saves a `.png` image of the constructed De Bruijn graph using the `graphviz` library.

---

## Prerequisites

Before you begin, ensure you have the following installed on your system:

- **Git**: To clone the repository.
- **Anaconda or Miniconda**: To manage the Python environment and dependencies.

---

## Setup and Installation

Follow these steps to get the project running on your local machine.

### 1. Clone the Repository

Open your terminal or command prompt and clone this repository:

```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
```
*(Replace `your-username` and `your-repo-name` with your actual GitHub details)*

### 2. Create and Activate the Conda Environment

This project uses a specific set of Python packages. The included `environment.yml` file lets you install them all with a single command.

```bash
# Create the environment from the file
conda env create -f environment.yml

# Activate the new environment
conda activate dna-assembler-env
```
This command creates a new environment named `dna-assembler-env` and installs all necessary libraries (`pandas`, `networkx`, `graphviz`). You must activate this environment every time you work on the project.

---

## How to Run the Script

With the `dna-assembler-env` environment active, you can now run the main script.

### Default Usage

To run the script with the default input file (`DNA_1_5.csv`):

```bash
python project.py
```

### Specifying an Input File

You can also provide the name of a different CSV file as a command-line argument. Make sure the file is located in the same directory as the script.

For example, to run the script with `DNA_2_7.csv`:

```bash
python project.py DNA_2_7.csv
```

---

## Expected Output

After running the script (e.g., with the default input), you should see the following:

1.  **Console Output**:
    - A JSON object showing the cleaned, processed DNA fragments.
    - A step-by-step printout of the Eulerian path being traversed (e.g., `Visiting directed edge: GCA -> CAT`).

2.  **Generated Files**:
    - **`DNA_1.png`**: An image file showing the visualized De Bruijn graph.
    - **`DNA_1.text`**: A text file containing the final, reconstructed DNA sequence (or an error message if it could not be assembled).

The base name of the output files (e.g., `DNA_1`) is derived automatically from the input CSV file name.
