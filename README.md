# bkpir
Implementation of Keyword PIR for Private Boolean Retrieval.

## Overview

This repository provides the reference implementation of the **BKPIR** protocol proposed in our paper *"BKPIR: Keyword PIR for Private Boolean Retrieval"*. The protocol enables privacy-preserving boolean retrieval over many-to-many keyword-value databases using homomorphic encryption (based on Microsoft SEAL). It supports single-keyword and compound queries (AND, OR, NOT) while protecting query privacy.

## Features

* Privacy-preserving boolean retrieval based on homomorphic encryption
* Supports many-to-many keyword-to-value mappings
* Enables single-keyword and compound queries (AND/OR/NOT)
* Parameterizable for different database sizes and configurations

## Dependencies

- **Compiler**: GNU G++ (version >= 6.0)
- **Libraries**:
  - Microsoft SEAL (version >=4.0)
  - CMake (version >= 3.13)

## Installation

### 1. Install System Dependencies  (Recommended for Ubuntu)

```bash
# Update system packages
sudo apt update && sudo apt upgrade -y

# Install build tools and dependencies
sudo apt install build-essential cmake git libssl-dev -y
```

### 2. Install Microsoft SEAL

This project depends on [Microsoft SEAL](https://github.com/microsoft/SEAL), a homomorphic encryption library.

**To install Microsoft SEAL:**

- Clone and build the library:

   ```bash
   git clone https://github.com/microsoft/SEAL.git
   cd SEAL
   cmake -S . -B build
   cmake --build build
   sudo cmake --install build
   ```

- Ensure SEAL is installed system-wide or set `CMAKE_PREFIX_PATH` when building this project.

**Official instructions:**
[https://github.com/microsoft/SEAL#installing-microsoft-seal](https://github.com/microsoft/SEAL#installing-microsoft-seal)

### 3. Build BKPIR

```bash
git clone https://github.com/gimmehug/bkpir.git
cd bkpir
mkdir build
cd build
cmake ..
make
```

*If SEAL is installed in a custom location:*
```bash
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/seal
make
```

## Usage

### Command Line Interface

```bash
./bkpir [parameters]
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-N` | Polynomial modulus degree | 8192 |
| `-n` | Number of keywords | 128 |
| `-m` | Number of items | 8192 |
| `-kl` | Max keyword byte length | 10 |
| `-il` | Max item byte length | 2 |
| `-k` | Hamming weight for constant-weight encoding | 4 |
| `-h` | Show help message | - |

### Examples

```bash
# Basic usage with default parameters
./bkpir

# Custom parameters for larger database
./bkpir -N 8192 -n 128 -m 32768 -il 2 -k 4
```

## Interactive Menu

After starting with parameters, you'll see an interactive menu.
Select an operation by entering its number:

1. **MMKPIRL**: Basic multi-keyword PIR implementation
2. **BKPIRL Single**: Single-keyword query
3. **BKPIRL AND**: Compound query (keyword1 AND keyword2)
4. **BKPIRL OR**: Compound query (keyword1 OR keyword2)
5. **BKPIRL NOT**: Single-keyword query (NOT keyword1)

## Experimental Results
The results presented in this paper can be replicated using the implementation provided in this repository. Detailed instructions for reproducing the results are available in the /script directory.

## Citation

If you use BKPIR in your research, please cite our paper:

BKPIR: Keyword PIR for Private Boolean Retrieval
