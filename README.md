# qPCR Primer Compatibility Analyzer

## Description

A Python GUI application that analyzes primer sequences for potential 3' end overlaps that can lead to primer-dimer formation in qPCR experiments. The tool examines all primer pairs for complementarity at their 3' ends and assigns risk levels based on overlap length and number of mismatches.

## Requirements

- Python 3.6+
- Biopython
- tkinter (included with standard Python installation)

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python overlapper_tool.py
```

## How It Works

The analyzer performs the following steps:

1. **Input Processing**: Reads primer sequences from a FASTA file
2. **3' End Analysis**: Examines the 3' ends of all primer pairs
3. **Reverse Complement Comparison**: Compares each primer's 3' end with the reverse complement of other primers' 3' ends
4. **Mismatch Calculation**: Counts base mismatches in overlapping regions
5. **Risk Assessment**: Assigns risk levels based on:
   - **HIGH**: ≥4 bases perfect match
   - **MEDIUM**: 2-3 bases perfect match OR ≥4 bases with 1 mismatch
   - **LOW**: All other overlaps

## Features

### Setup Tab
- **File Selection**: Browse and select input FASTA file
- **Output Directory**: Optional selection for saving results
- **Analysis Parameters**:
  - Overlap length range (minimum and maximum bases to check)
  - Maximum allowed mismatches
- **Preview Sequences**: View loaded primers before analysis

### Results Tab
- **Analysis Summary**: Overview of total overlaps and risk distribution
- **Detailed Results Table**: Sortable list of all detected overlaps
- **Visualization**: ASCII representation of primer alignments
- **Export**: Save results to CSV format

## Input Format

FASTA file containing primer sequences:
```
>Primer1
ATCGATCGATCGATCGATCG
>Primer2
TAGCTAGCTAGCTAGCTAGC
```

## Output

- **Interactive GUI Display**: Real-time viewing and sorting of results
- **CSV Export**: Tabular format with columns:
  - Overlap Length
  - Mismatches
  - Primer 1
  - Primer 2
  - Risk Level
- **Text Report**: Detailed analysis with visualizations (when output directory specified)

## Analysis Parameters

- **Overlap Length Range**: Defines the minimum and maximum number of bases to check for overlaps (default: 3-10)
- **Maximum Mismatches**: Number of allowed base mismatches in overlap region (default: 1)

## Technical Implementation

### Class Structure

The application is built around a single main class:
- **`PrimerCompatibilityAnalyzer`**: Main class containing all GUI and analysis logic

### Key Methods

#### Core Analysis Methods
- **`run_analysis()`**: Main analysis loop that:
  1. Iterates through overlap lengths (max to min)
  2. For each overlap length, checks all mismatch levels (0 to max)
  3. Compares all primer pairs (including self-pairs)
  4. Stores results in `self.analysis_results` list

- **`get_last_n_bases(seq_record, n)`**: Extracts the last n bases from a sequence's 3' end
- **`get_last_n_bases_rc(seq_record, n)`**: Gets the reverse complement of the last n bases
- **`count_mismatches(seq1, seq2)`**: Performs base-by-base comparison using zip()
- **`get_risk_level(overlap_length, mismatches)`**: Decision tree for risk assessment

#### GUI Methods
- **`create_widgets()`**: Initializes ttk.Notebook with two tabs
- **`create_setup_tab()`**: Builds configuration interface with ttk widgets
- **`create_results_tab()`**: Creates results display with ttk.Treeview
- **`show_overlap_detail(event)`**: Event handler for treeview selection
- **`sort_treeview(col)`**: Implements bidirectional column sorting

#### Visualization Method
- **`visualize_overlap(seq1, seq2, overlap_length)`**: 
  - Creates ASCII alignment representation
  - Calculates offset for proper alignment display
  - Generates match/mismatch indicators

### Data Structures

- **`self.sequences`**: List of Bio.SeqRecord objects from input FASTA
- **`self.analysis_results`**: List of dictionaries containing:
  ```python
  {
      'overlap_length': int,
      'mismatches': int,
      'primer1': SeqRecord,
      'primer2': SeqRecord,
      'risk_level': str,
      'visualization': str
  }
  ```
- **`self.sort_reverse`**: Dictionary tracking sort direction for each column

### Algorithm Details

1. **Overlap Detection Algorithm**:
   ```python
   for overlap_length in range(max_overlap, min_overlap - 1, -1):
       for mismatches in range(max_mismatches + 1):
           for i in range(len(sequences)):
               for j in range(i, len(sequences)):
                   # Compare 3' ends with reverse complement
   ```
   - Time complexity: O(n² × m × k) where n=sequences, m=overlap range, k=mismatch levels
   - Starts with longer overlaps (more significant) and works down

2. **Complementarity Check**:
   - Takes 3' end of primer 1: `ATCG` (last 4 bases)
   - Takes 3' end of primer 2: `GCTA` (last 4 bases)
   - Reverse complements primer 2's 3' end: `TAGC`
   - Compares primer 1's 3' end with reverse complement

3. **Risk Level Logic**:
   - Perfect matches (0 mismatches): Based purely on overlap length
   - With mismatches: Only considers MEDIUM risk if ≥4 bases with 1 mismatch

### GUI Implementation

- **Framework**: tkinter with ttk widgets for modern appearance
- **Layout**: Pack geometry manager for main sections, Grid for form elements
- **Styling**: Custom ttk.Style configurations for professional appearance
- **Event Handling**: 
  - TreeviewSelect event for detail display
  - Column header commands for sorting
  - Button commands for file operations

### File I/O

- **Input**: Uses BioPython's SeqIO.parse() for FASTA reading
- **Output**: 
  - CSV export using Python's csv module
  - Text report with formatted alignment visualizations
  - Optional timestamped output files

### Widget State Management

- Text widgets use `state='disabled'` to prevent user editing
- Progress bar uses indeterminate mode during analysis
- Status updates with color coding via ttk.Style

## Common Modifications

### Adding New Risk Categories
Modify the `get_risk_level()` method to add more granular risk assessment:
```python
def get_risk_level(self, overlap_length, mismatches):
    # Add your custom logic here
    if overlap_length >= 6 and mismatches == 0:
        return "CRITICAL"
    # ... existing logic
```

### Changing Visualization Format
The `visualize_overlap()` method can be modified to show different alignment styles:
- Add mismatch highlighting
- Include base positions
- Show melting temperature calculations

### Adding Export Formats
New export methods can be added alongside `export_results()`:
- JSON export for programmatic processing
- HTML report generation
- Integration with other bioinformatics tools

### Extending Analysis Parameters
Additional parameters can be added to the analysis:
- GC content consideration
- Temperature calculations
- Salt concentration adjustments

### GUI Customization
- Colors and fonts are defined in `setup_styles()`
- Window geometry set in `__init__()`
- Tab structure defined in `create_widgets()`

## Notes

- Self-dimers (primers forming dimers with themselves) are included in the analysis
- The tool analyzes all possible primer pair combinations including self-pairs
- Visualization shows primer alignments with match indicators (|) between complementary bases
- The treeview sorting maintains state for each column's direction