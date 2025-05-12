# Detailed Usage Guide

This document provides detailed instructions for using the CrispHunter tool.

## Pipeline Overview

The tool works in several distinct stages:

1. **Genome retrieval and processing**
   - Download genomes from NCBI GenBank (optional)
   - Extract and decompress genome files
   
2. **Genome annotation**
   - Annotate genomes using Prokka
   - Extract GFF files for pan-genome analysis
   
3. **Pan-genome analysis**
   - Identify core genes across all genomes using Panaroo
   - Extract core gene sequences
   
4. **crRNA design**
   - Design candidate crRNAs using CaSilico
   - Extract designed crRNAs from web interface results
   
5. **Target coverage analysis**
   - Evaluate crRNA coverage across all target genomes
   - Select crRNAs with 100% coverage
   
6. **Off-target evaluation**
   - Build database of potential off-target genomes
   - Analyze potential off-target binding sites using FlashFry
   - Filter crRNAs based on off-target analysis
   
7. **Output generation**
   - Combine target coverage and off-target information
   - Generate final list of optimal crRNAs

## Command Options Explained

### Required Arguments

- `-g, --genus`: The target bacterial genus (and species) to design crRNAs for. Example: "Burkholderia pseudomallei"

- `-offg, --offgenomes`: Path to the FASTA file containing genomes to check for off-target effects. This should include genomes of organisms that you want to ensure your crRNAs do NOT bind to (e.g., human genome, closely related bacterial species).

- `-firedriver, --firefox_driver_path`: Path to the GeckoDriver executable. This is required for web scraping the CaSilico results. Download from: https://github.com/mozilla/geckodriver/releases

- `-ffry, --flashfry_path`: Path to the FlashFry JAR file. Download from: https://github.com/mckennalab/FlashFry/releases

### Optional Arguments

- `-i, --input_path`: If you already have downloaded genome files, provide the path to the directory containing them. This will skip the download step.

- `-gff, --gff_file`: If you have already annotated your genomes with Prokka, provide the path to the directory containing the GFF files.

- `-coregene, --core_gene_folder`: If you have already performed pan-genome analysis and identified core genes, provide the path to the directory containing those results.

- `-m, --memory`: Memory allocation for FlashFry. Format: "30g" for 30 gigabytes. Increase for larger genomes.

- `-maxOff, --maximumOffTargets`: Maximum number of off-target sites to consider during FlashFry analysis. Default is 12000.

- `-PAM, --PAM_sequence`: The PAM sequence to use for crRNA design. Default is "TTTN" for Cas12a.

- `-t, --threads`: Number of CPU threads to use for parallel processing steps. Default is 8.

## Interpreting Results

The final output file "off_target_and_target_coverage.tsv" contains the following columns:

1. **gRNA ID**: Unique identifier for the guide RNA
2. **gRNA Sequence**: The sequence of the guide RNA
3. **Number of Targeted Genomes**: How many of the input genomes this gRNA can target
4. **Total Genomes Analyzed**: Total number of genomes that were analyzed
5. **Targeting Coverage (%)**: Percentage of analyzed genomes that can be targeted
6. **Potential Off-target Sites**: Number of potential off-target binding sites

The best crRNAs will have:
- 100% Targeting Coverage
- 0 Potential Off-target Sites

## Customization

### Modifying GC Content Filtering

The default GC content filter is set to 40-65%. If you need to modify this range, you can edit the R code in the `summarize_crawler_results` function.

### Changing PAM Sequence

To target different Cas enzymes, you can modify the PAM sequence using the `-PAM` parameter. Default is "TTTN" for Cas12a.

## Troubleshooting

### Common Issues

1. **Memory Errors in FlashFry**: Increase the `-m` parameter value
2. **WebDriver Issues**: Ensure GeckoDriver version is compatible with your Firefox version
3. **Pan-genome Analysis Failures**: Reduce thread count to 16 using the `-t` parameter
