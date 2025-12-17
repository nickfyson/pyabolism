# Pyabolism

[![License](https://img.shields.io/badge/License-BSD-blue.svg)](LICENSE)
[![Python Version](https://img.shields.io/badge/python-2.7-blue.svg)](https://www.python.org/downloads/)

Pyabolism is a Python module for constraint-based simulation and analysis of metabolic networks. It provides tools for Flux Balance Analysis (FBA), Flux Variability Analysis (FVA), and integration of gene expression data into metabolic models.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Standard Installation](#standard-installation)
  - [Docker Installation](#docker-installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Loading Models](#loading-models)
  - [Flux Balance Analysis (FBA)](#flux-balance-analysis-fba)
  - [Flux Variability Analysis (FVA)](#flux-variability-analysis-fva)
  - [Gene Expression Integration](#gene-expression-integration)
  - [Visualization](#visualization)
- [Project Structure](#project-structure)
- [Development](#development)
  - [Setting Up Development Environment](#setting-up-development-environment)
  - [Running Tests](#running-tests)
  - [Code Style](#code-style)
- [Contributing](#contributing)
- [Examples](#examples)
- [Dependencies](#dependencies)
- [License](#license)
- [Contact](#contact)

## Features

- **Flux Balance Analysis (FBA)**: Predict steady-state metabolic flux distributions using linear programming
- **Flux Variability Analysis (FVA)**: Determine the range of possible fluxes through each reaction
- **Expression-based Flux (E-Flux)**: Integrate gene expression data to constrain flux predictions
- **Gene-Centric Flux (GC-Flux)**: Alternative method for incorporating expression data
- **SBML Support**: Load and save metabolic models in Systems Biology Markup Language (SBML) format
- **Visualization Tools**: Built-in plotting functions for flux distributions and metabolic states
- **Model Manipulation**: Tools for analyzing exchange reactions, transport reactions, and gene-protein-reaction (GPR) associations

## Installation

### Prerequisites

Pyabolism requires the following dependencies:

- **Python 2.7** (Note: Currently only compatible with Python 2.7.x)
- **Gurobi Optimizer**: Commercial optimization solver (free academic licenses available)
- **gurobipy**: Python bindings for Gurobi
- **numpy**: Numerical computing library
- **python-libsbml**: Library for reading and writing SBML files

Optional dependencies for visualization:
- **matplotlib**: For plotting flux distributions
- **seaborn**: Enhanced visualization styles

### Standard Installation

1. **Install Gurobi**:
   - Download and install Gurobi from [www.gurobi.com](http://www.gurobi.com)
   - Obtain a license (free for academic use) and activate it using `grbgetkey`

2. **Install Python dependencies**:
   ```bash
   pip install numpy python-libsbml
   ```

3. **Install Pyabolism**:
   ```bash
   git clone https://github.com/nickfyson/pyabolism.git
   cd pyabolism
   python setup.py install
   ```

### Docker Installation

A Docker container with Jupyter Notebook, Pyabolism, and Gurobi pre-installed is available. See the [jupyter-notebook README](https://github.com/nickfyson/pyabolism/tree/master/jupyter-notebook) for detailed instructions.

Quick start:
```bash
cd jupyter-notebook
cp env_make_example env_make
# Edit env_make to customize port settings
make build run
```

## Quick Start

```python
import pyabolism

# Load a metabolic model from SBML file
model = pyabolism.io.load_model('path/to/model.xml', filetype='sbml')

# Run Flux Balance Analysis
pyabolism.simulate.FBA(model, show=True)

# Access flux values
for reaction in model.reactions():
    print(f"{reaction.id}: {reaction.flux_value}")

# Visualize flux distribution
pyabolism.visualise.plot_flux_distribution(model)
```

## Usage

### Loading Models

Pyabolism supports loading models from SBML files or Python pickle format:

```python
import pyabolism

# Load from SBML
model = pyabolism.io.load_model('model.xml', filetype='sbml')

# Save to pickle for faster loading
pyabolism.io.save_model(model, 'model.pkl', filetype='pickle')

# Load from pickle
model = pyabolism.io.load_model('model.pkl', filetype='pickle')
```

### Flux Balance Analysis (FBA)

FBA predicts metabolic fluxes by optimizing an objective function (typically biomass production) subject to stoichiometric and thermodynamic constraints:

```python
import pyabolism

model = pyabolism.io.load_model('model.xml')

# Standard FBA
pyabolism.simulate.FBA(model, show=True)

# Access results
print(f"Growth rate: {model.total_objective}")
for reaction in model.reactions():
    if reaction.flux_value != 0:
        print(f"{reaction.id}: {reaction.flux_value}")
```

### Flux Variability Analysis (FVA)

FVA determines the minimum and maximum possible flux through each reaction while maintaining a specified objective value:

```python
import pyabolism

model = pyabolism.io.load_model('model.xml')

# Run FBA first to establish baseline
pyabolism.simulate.FBA(model)

# Run FVA
pyabolism.simulate.FVA(model, fraction=0.9)

# Access flux ranges
for reaction in model.reactions():
    print(f"{reaction.id}: [{reaction.flux_range[0]}, {reaction.flux_range[1]}]")
```

### Gene Expression Integration

Integrate gene expression data to constrain metabolic models:

```python
import pyabolism

model = pyabolism.io.load_model('model.xml')

# Prepare expression data as a dictionary
expression_data = {
    'gene1': 5.2,
    'gene2': 3.1,
    # ... more genes
}

# Run E-Flux
pyabolism.simulate.EFlux(model, expression_data, show=True)

# Or use GC-Flux
pyabolism.simulate.GCFlux(model, expression_data, show=True)
```

### Visualization

Pyabolism includes tools for visualizing metabolic flux distributions:

```python
import pyabolism

model = pyabolism.io.load_model('model.xml')
pyabolism.simulate.FBA(model)

# Plot flux distribution
pyabolism.visualise.plot_flux_distribution(model, title='FBA Results')

# Plot specific reactions
reactions = [r for r in model.reactions() if r.flux_value > 1.0]
pyabolism.visualise.plot_flux_distribution(model, reactions=reactions)
```

## Project Structure

```
pyabolism/
├── pyabolism/              # Main package directory
│   ├── __init__.py         # Package initialization
│   ├── model.py            # Core model classes (MetaModel, Reaction, Metabolite, etc.)
│   ├── io.py               # SBML and file I/O functions
│   ├── visualise.py        # Plotting and visualization tools
│   ├── tools.py            # Utility functions (exchange reactions, GPR parsing, etc.)
│   └── simulate/           # Simulation algorithms
│       ├── FBA.py          # Flux Balance Analysis
│       ├── FVA.py          # Flux Variability Analysis
│       ├── EFlux.py        # Expression-based flux analysis
│       ├── GCFlux.py       # Gene-centric flux analysis
│       ├── LP.py           # Linear programming utilities
│       └── simtools.py     # Simulation helper functions
├── tests/                  # Unit tests
│   ├── test_FBA.py
│   ├── test_FVA.py
│   ├── test_EFlux.py
│   ├── test_GCFlux.py
│   └── test_IO.py
├── examples/               # Jupyter notebook examples
│   ├── 01_FBA.ipynb
│   ├── 02_exploring_models.ipynb
│   ├── 03_incorporating_expression_data.ipynb
│   └── 04_Flux_Variability_Analysis.ipynb
├── jupyter-notebook/       # Docker container setup
├── setup.py                # Package installation configuration
├── requirements.txt        # Python dependencies
├── Makefile                # Build and test automation
└── README.md               # This file
```

## Development

### Setting Up Development Environment

1. **Clone the repository**:
   ```bash
   git clone https://github.com/nickfyson/pyabolism.git
   cd pyabolism
   ```

2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Install in development mode**:
   ```bash
   python setup.py develop
   ```

### Running Tests

Pyabolism uses Python's built-in unittest framework:

```bash
# Run all tests
python setup.py test

# Or use make
make test

# Run with coverage
make coverage
```

### Code Style

The project uses flake8 for linting:

```bash
# Check code style
make lint

# Or directly
flake8 pyabolism tests
```

### Available Make Commands

```bash
make help          # Show all available commands
make clean         # Remove build, test, and Python artifacts
make test          # Run tests quickly with the default Python
make coverage      # Check code coverage
make lint          # Check code style with flake8
make dist          # Build distribution packages
make install       # Install the package
```

## Contributing

We welcome contributions to Pyabolism! Here's how you can help:

### Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/pyabolism.git
   cd pyabolism
   ```
3. **Create a branch** for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

### Making Changes

1. Make your changes to the codebase
2. Add or update tests to cover your changes
3. Ensure all tests pass: `make test`
4. Check code style: `make lint`
5. Update documentation as needed

### Submitting Changes

1. **Commit your changes**:
   ```bash
   git add .
   git commit -m "Description of your changes"
   ```
2. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```
3. **Open a Pull Request** on GitHub with a clear description of your changes

### Contribution Guidelines

- Write clear, concise commit messages
- Follow existing code style and conventions
- Add tests for new functionality
- Update documentation for user-facing changes
- Keep pull requests focused on a single feature or fix
- Be respectful and constructive in discussions

### Reporting Issues

If you find a bug or have a feature request:

1. Check if the issue already exists in the issue tracker
2. If not, create a new issue with:
   - Clear description of the problem or feature
   - Steps to reproduce (for bugs)
   - Expected vs. actual behavior
   - Your environment (Python version, OS, etc.)

## Examples

The `examples/` directory contains Jupyter notebooks demonstrating various features:

1. **[01_FBA.ipynb](examples/01_FBA.ipynb)**: Basic Flux Balance Analysis
2. **[02_exploring_models.ipynb](examples/02_exploring_models.ipynb)**: Model exploration and manipulation
3. **[03_incorporating_expression_data.ipynb](examples/03_incorporating_expression_data.ipynb)**: Integrating gene expression data
4. **[04_Flux_Variability_Analysis.ipynb](examples/04_Flux_Variability_Analysis.ipynb)**: FVA and flux range analysis

To run the examples:

```bash
cd examples
jupyter notebook
```

Or use the Docker container as described in the [Docker Installation](#docker-installation) section.

## Dependencies

### Required Dependencies

- **Python 2.7**: Core language (Python 3 support is not yet available)
- **gurobipy**: Python interface to the Gurobi Optimizer
- **numpy**: Numerical computing library
- **python-libsbml**: SBML file format support

### Optional Dependencies

- **matplotlib**: For plotting and visualization
- **seaborn**: Enhanced visualization aesthetics
- **jupyter**: For running example notebooks

### Why Gurobi?

Gurobi is a high-performance commercial optimization solver that provides:
- Fast solution times for large-scale metabolic models
- Robust handling of complex linear programming problems
- Free academic licenses available at [www.gurobi.com](http://www.gurobi.com)

Alternative open-source solvers may be supported in future versions.

## License

Pyabolism is licensed under the BSD License. See the LICENSE file for details.

## Contact

- **Author**: Nick Fyson
- **Email**: nick@fyson.net
- **GitHub**: [https://github.com/nickfyson/pyabolism](https://github.com/nickfyson/pyabolism)

## Acknowledgments

Pyabolism builds upon established methods in constraint-based modeling:

- **FBA**: Orth, J. D., Thiele, I., & Palsson, B. Ø. (2010). What is flux balance analysis? *Nature Biotechnology*, 28(3), 245-248.
- **E-Flux**: Colijn, C., et al. (2009). Interpreting expression data with metabolic flux models: predicting Mycobacterium tuberculosis mycolic acid production. *PLoS Computational Biology*, 5(8), e1000489.
- **FVA**: Mahadevan, R., & Schilling, C. H. (2003). The effects of alternate optimal solutions in constraint-based genome-scale metabolic models. *Metabolic Engineering*, 5(4), 264-276.

## Additional Resources

- **SBML**: [http://sbml.org](http://sbml.org)
- **Gurobi Documentation**: [https://www.gurobi.com/documentation/](https://www.gurobi.com/documentation/)
- **Constraint-Based Reconstruction and Analysis (COBRA)**: [http://opencobra.github.io](http://opencobra.github.io)
- **BiGG Models Database**: [http://bigg.ucsd.edu](http://bigg.ucsd.edu) - Repository of genome-scale metabolic models

---

For more information, please visit the [project repository](https://github.com/nickfyson/pyabolism) or open an issue.
