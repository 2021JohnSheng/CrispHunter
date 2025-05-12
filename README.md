# CrispHunter

A tool for designing highly specific CRISPR-Cas12a crRNAs targeting user-defined bacterial genera.

## Overview

CrispHunter is a comprehensive bioinformatics pipeline that integrates pan-genome analysis with CRISPR-Cas12a targeting to identify highly specific guide RNAs for bacterial detection and targeting. The tool identifies conserved regions across a genus of interest and designs crRNAs that provide maximum coverage with minimal off-target effects.

## Features

- **Automated genome retrieval** from NCBI GenBank database
- **Pan-genome analysis** to identify core genes across all genomes of a genus
- **Intelligent crRNA design** that ensures 100% coverage across target genomes

## Installation

### Setup

```bash
# Clone the repository
git clone https://github.com/2021JohnSheng/CrispHunter.git
cd CrispHunter

# Create and activate the Conda environment
conda env create -f environment.yml
conda activate CrispHunter
