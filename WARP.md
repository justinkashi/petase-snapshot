# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Project Overview

This is a protein engineering research project focused on PETase (polyethylene terephthalate hydrolase) enzymes and their variants. The project combines:
- PETase sequence dataset curation and management
- Benchmark datasets for protein stability and solubility prediction
- ESM (Evolutionary Scale Modeling) protein language model embeddings and fine-tuning
- Mutation analysis and variant generation
- Integration with external databases (UniProt, MGnify, NCBI, FireProtDB, ThermomutDB)

## Key Directories

- `petase_db/`: Core PETase dataset files including mutation tables, structural data, and reference sequences
- `data/`: External benchmark datasets (stability, solubility), FASTA files, and processed data
  - `data/benchmark/`: Consolidated benchmark datasets from multiple sources
  - `data/temp/`: Temporary processing files (gitignored)
- `notebooks/`: Jupyter notebooks for analysis and data processing
  - `main.ipynb`: Primary workflow for loading and processing benchmark datasets
  - `archive_petase.ipynb`: PETase-specific data fetching and mutation tools
  - `archive_esm.ipynb.ipynb`: ESM model embedding and inference workflows
  - `esmfinetune_tutorial.ipynb`: ESM model fine-tuning examples
- `scripts/`: Shell scripts and utilities for data fetching and processing
- `venv/`: Python virtual environment (gitignored)

## Development Commands

### Environment Setup
```bash
# Activate virtual environment
source venv/bin/activate

# The requirements.txt is currently minimal - dependencies are managed in notebooks
# Common dependencies include: pandas, biopython, requests, torch, esm, py3Dmol, biotite
```

### Running Notebooks
```bash
# Launch Jupyter from the project root
jupyter notebook

# Main analysis workflow
jupyter notebook notebooks/main.ipynb
```

### Data Fetching
```bash
# Fetch benchmark databases (from scripts/)
bash scripts/benchmarkdb_fetch.sh

# InterPro scan setup (see scripts/readme_interpro.md for detailed instructions)
bash scripts/interproscan.sh
```

### Working with FASTA Files
The project uses CD-HIT for sequence clustering at various similarity thresholds (typically 95-100% identity). Output files follow the pattern `*_cdhit*.fasta` and `*_cluster.txt`.

## Data Processing Architecture

### Multi-Source Dataset Integration
The project integrates protein data from multiple benchmark sources:

**Solubility Datasets:**
- NESG Solubility (~10k proteins, expression/solubility labels)
- Soluprot (~14k proteins, binary solubility)
- Price Solubility (~7k proteins, usability metric)
- ProtSolM (~71k proteins from HuggingFace)

**Stability Datasets:**
- FireProtDB (~53k variants with ddG, dTm, pH, mutation effects)
- ThermomutDB (~12k variants with thermodynamic measurements)
- Meltome Atlas (~1M variants with thermal stability data)
- Novozymes Kaggle (~31k enzyme variants with Tm)

**Function Datasets:**
- CAFA-5 (~142k proteins with function annotations)

### Data Loading Pipeline
The `main.ipynb` notebook defines loader functions for each dataset:
- `load_nesg()`, `load_price()`, `load_soluprot()`: Solubility loaders
- `load_fireprot()`, `load_thermomut()`, `load_meltome()`: Stability loaders
- `load_protsolm()`, `load_novozymes()`: Additional benchmark loaders

Each loader:
1. Reads CSV/JSON metadata
2. Maps sequences from FASTA files using UniProt IDs or custom identifiers
3. Standardizes column names (particularly `sequence`, `id`)
4. Returns a pandas DataFrame

### Sequence Clustering Workflow
After loading all datasets:
1. Merge DataFrames into a single collection
2. Export sequences to combined FASTA using `fasta_merger_from_dfs()`
3. Run CD-HIT for redundancy removal (typically at 100% identity)
4. Parse cluster results with `parse_cd_hit_clusters()`
5. Annotate original DataFrames with cluster membership

### PETase Mutation Tools
The `archive_petase.ipynb` contains utilities for:
- Fetching sequences from UniProt, MGnify, or NCBI given ID lists
- Applying mutations to FASTA sequences with `apply_mutations_to_fasta()`
- Handling signal peptide offsets in mutation numbering
- Generating variant FASTA files for experimental designs

## ESM Model Integration

### Model Access
The project uses ESM models via two interfaces:
1. **Local models**: ESMC (300M-6B parameters) for embeddings
2. **Forge API**: Cloud-based ESM3 inference for large-scale jobs

### Embedding Workflow
- `ESMProtein` objects encapsulate sequence, structure, secondary structure, SASA, and function annotations
- `embed_sequence()` function generates embeddings with `LogitsConfig(sequence=True, return_embeddings=True)`
- Batch processing supported via `batch_executor()` for multiple sequences
- Embeddings can be used for downstream ML tasks (stability/solubility prediction)

### Authentication
ESM Forge API requires a token (stored in notebooks - should be moved to environment variables):
```python
from esm.sdk.forge import ESM3ForgeInferenceClient
client = ESM3ForgeInferenceClient(model="esmc-6b-2024-12", token=os.getenv("ESM_TOKEN"))
```

## Data Format Standards

### FASTA Headers
- PETase master database uses `Name` field from source
- Benchmark datasets use format: `{dataframe_class}_{index}_{row_id}`
- UniProt fetches preserve original headers

### DataFrame Columns
Standardized across loaders:
- `sequence`: Amino acid sequence (string)
- `id` or `sid`: Sequence identifier
- `uniprot_id`: UniProt accession when available

Dataset-specific columns preserved (e.g., `ddG`, `dTm`, `pH`, `solubility`, `Tm`)

### Mutation Notation
Standard single-letter format: `{WT_residue}{position}{mutant_residue}`
- Example: `S121E` (Serine at position 121 to Glutamate)
- Account for signal peptide offsets when applicable

## Important Files

- `data/masterdb.tsv`: Primary PETase dataset metadata
- `data/masterdb.fasta`: Corresponding sequences
- `petase_db/fireprotdb_results.csv`: Large stability mutation dataset
- `petase_db/thermomutdb.json`: Thermodynamic mutation measurements
- `petase_db/petase_mutation_table_allseqs.csv`: Curated PETase variant table

## Notes for AI Assistants

- Python 3.12 is the target version (venv configured with 3.12.7)
- The project uses both local development (MacOS M-series with MPS) and potentially cloud compute
- Many notebooks contain hardcoded paths (e.g., `c:\\Users\\justi\\petase-1`) that may need adjustment
- API tokens for ESM Forge and HuggingFace are currently hardcoded in notebooks - should be externalized
- CD-HIT output files are large and not version controlled
- When fetching from UniProt/NCBI, implement rate limiting (e.g., 0.15s delay) to avoid 500 errors
- Some UniProt IDs in source data may be deprecated; `fetch_uniprot_fasta()` handles failures gracefully
