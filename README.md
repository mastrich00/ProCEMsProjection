# ProCEM Projection

Implementation of metabolic network projection method based on [Marashi, David and Bockmayr's work](https://almob.biomedcentral.com/articles/10.1186/1748-7188-7-17).

## able of Contents
- [Features](#-features)
- [Prerequisites](#-prerequisites)
- [Installation](#-installation)
- [Usage](#-usage)
- [Project Structure](#-project-structure)
- [Test Cases](#-test-cases)
- [License](#-license)

## Features
- Support for both **polco** and **mplrs**
- Excel and SBML file input support

## Prerequisites
- Python 3.10
- Numerical Tools:
  - mplrs v7.2
  - polco v4.7.1
  - OpenMPI v4.1.4
- Java Runtime Environment (for polco)

## Installation
1. Clone repository:
   ```sh
   git clone git@github.com:mastrich00/ProCEMsProjection.git
   cd ProCEMsProjection
   ```

2. Install Python dependencies:
   ```sh
   pip install -r requirements.txt
   ```

3. Verify tool installations:
   ```sh
   mpirun --version
   java -jar polco.jar --version
   ```

## Usage

### General Syntax
For the A. thaliana Plastid Excel provided by Marashi:
```sh
python test_excel.py \
  --input <INPUT_FILE> \
  --tool {polco|mplrs} \
  --mplrs <MPLRS_PATH> \
  --polco <POLCO_JAR_PATH> \
  --stepsize <STEP_SIZE> \
  -n <THREAD_COUNT> \
  --outdir <OUTPUT_DIRECTORY>
```

For SBML input (`.xml` files):
```sh
python test_sbml.py \
  -f <SBML_FILE> \
  -m <MODEL_NAME> \
  --tool {polco|mplrs} \
  --mplrs <MPLRS_PATH> \
  --polco <POLCO_JAR_PATH> \
  --stepsize <STEP_SIZE> \
  -n <THREAD_COUNT>
```

## Project Structure
```
.
├── input/              # Sample input files
│   ├── 13015_2011_... # A. thaliana model
│   ├── e_coli_core.xml
│   └── mmsyn_sm00.xml
├── test_excel.py       # Excel file processor
└── test_sbml.py        # SBML file processor
```

## Test Cases

### A. thaliana Plastid Model
**Excel Input** - 120 threads  
*Using polco (stepsize 111):*
```sh
python test_excel.py \
  --input input/13015_2011_155_MOESM1_ESM.xls \
  --tool polco \
  --mplrs mplrs \
  --polco polco.jar \
  --stepsize 111 \
  -n 120 \
  --outdir testResults/measurements
```

*Using mplrs (stepsize 5):*
```sh
python test_excel.py \
  --input input/13015_2011_155_MOESM1_ESM.xls \
  --tool mplrs \
  --mplrs mplrs \
  --stepsize 5 \
  -n 120 \
  --outdir testResults/measurements
```

### E. Coli Core Model
**SBML Input** - 120 threads  
*Using mplrs (stepsize 2):*
```sh
python test_sbml.py \
  -f input/e_coli_core.xml \
  -m e_coli_core \
  --tool mplrs \
  --mplrs mplrs \
  --stepsize 2 \
  -n 120
```

*Using polco (stepsize 2):*
```sh
python test_sbml.py \
  -f input/e_coli_core.xml \
  -m e_coli_core \
  --tool polco \
  --polco polco.jar \
  --stepsize 2 \
  -n 120
```

### mmsyn_sm00 Model
**SBML Input** - 120 threads  
*Using mplrs (stepsize 2):*
```sh
python test_sbml.py \
  -f input/mmsyn_sm00.xml \
  -m mmsyn_sm00 \
  --tool mplrs \
  --mplrs mplrs \
  --stepsize 2 \
  -n 120
```