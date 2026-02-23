from Bio.Seq import Seq
from Bio import SeqIO
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
from datetime import datetime
from io import StringIO


class PrimerCompatibilityAnalyzer:
    # IUPAC ambiguity codes - maps each code to its possible bases
    IUPAC_CODES = {
        'A': {'A'},
        'C': {'C'},
        'G': {'G'},
        'T': {'T'},
        'R': {'A', 'G'},        # puRine
        'Y': {'C', 'T'},        # pYrimidine
        'S': {'G', 'C'},        # Strong
        'W': {'A', 'T'},        # Weak
        'K': {'G', 'T'},        # Keto
        'M': {'A', 'C'},        # aMino
        'B': {'C', 'G', 'T'},   # not A
        'D': {'A', 'G', 'T'},   # not C
        'H': {'A', 'C', 'T'},   # not G
        'V': {'A', 'C', 'G'},   # not T
        'N': {'A', 'C', 'G', 'T'},  # aNy
    }

    def __init__(self, root):
        self.root = root
        self.root.title("qPCR Primer Compatibility Analyzer v2.0")
        self.root.geometry("1200x800")
        self.root.configure(bg="#f0f0f0")

        # Variables
        self.input_file_location = ""
        self.sequences = []
        self.sort_reverse = {}  # Track sort direction for each column

        # Configure style
        self.setup_styles()

        # Create main interface
        self.create_widgets()

    def setup_styles(self):
        """Configure ttk styles for professional appearance"""
        style = ttk.Style()
        style.theme_use('clam')

        # Configure custom styles
        style.configure('Title.TLabel', font=('Arial', 16, 'bold'), foreground='#2c3e50')
        style.configure('Heading.TLabel', font=('Arial', 12, 'bold'), foreground='#34495e')
        style.configure('Action.TButton', font=('Arial', 10, 'bold'))
        style.configure('Success.TLabel', foreground='#27ae60')
        style.configure('Error.TLabel', foreground='#e74c3c')

    def create_widgets(self):
        """Create and arrange all GUI widgets"""
        # Main title
        title_frame = ttk.Frame(self.root)
        title_frame.pack(pady=20, padx=20, fill='x')

        ttk.Label(title_frame, text="qPCR Primer Compatibility Analyzer",
                  style='Title.TLabel').pack()
        ttk.Label(title_frame, text="Analyze 3' end overlaps to predict primer-dimer formation",
                  font=('Arial', 10)).pack(pady=(5, 0))

        # Create notebook for organized interface
        notebook = ttk.Notebook(self.root)
        notebook.pack(pady=20, padx=20, fill='both', expand=True)

        # Setup tab
        self.setup_tab = ttk.Frame(notebook)
        notebook.add(self.setup_tab, text="Setup & Configuration")
        self.create_setup_tab()

        # Results tab
        self.results_tab = ttk.Frame(notebook)
        notebook.add(self.results_tab, text="Analysis Results")
        self.create_results_tab()

    def create_setup_tab(self):
        """Create the setup and configuration tab"""
        # File selection frame
        file_frame = ttk.LabelFrame(self.setup_tab, text="File Selection", padding=15)
        file_frame.pack(pady=10, padx=20, fill='x')

        # Input file selection
        ttk.Label(file_frame, text="Input FASTA File:", style='Heading.TLabel').grid(
            row=0, column=0, sticky='w', pady=(0, 5))

        input_frame = ttk.Frame(file_frame)
        input_frame.grid(row=1, column=0, sticky='ew', pady=(0, 15))
        file_frame.columnconfigure(0, weight=1)
        input_frame.columnconfigure(1, weight=1)

        ttk.Button(input_frame, text="Browse...", command=self.select_input_file,
                   style='Action.TButton').grid(row=0, column=0, padx=(0, 10))

        self.input_file_var = tk.StringVar(value="No file selected")
        self.input_label = ttk.Label(input_frame, textvariable=self.input_file_var,
                                     foreground='#7f8c8d')
        self.input_label.grid(row=0, column=1, sticky='w')

        # Sequence input frame
        seq_input_frame = ttk.LabelFrame(self.setup_tab, text="Sequence Input (FASTA Format)", padding=15)
        seq_input_frame.pack(pady=10, padx=20, fill='both', expand=True)

        # Instructions label
        ttk.Label(seq_input_frame, text="Paste sequences in FASTA format or import from file above:",
                  style='Heading.TLabel').pack(anchor='w', pady=(0, 5))

        # Text input area with scrollbars
        seq_text_frame = ttk.Frame(seq_input_frame)
        seq_text_frame.pack(fill='both', expand=True)

        self.sequence_text = tk.Text(seq_text_frame, wrap='none', font=('Courier', 10),
                                     height=8, bg='#ffffff')

        seq_scroll_y = ttk.Scrollbar(seq_text_frame, orient='vertical', command=self.sequence_text.yview)
        seq_scroll_x = ttk.Scrollbar(seq_text_frame, orient='horizontal', command=self.sequence_text.xview)
        self.sequence_text.configure(yscrollcommand=seq_scroll_y.set, xscrollcommand=seq_scroll_x.set)

        self.sequence_text.pack(side='left', fill='both', expand=True)
        seq_scroll_y.pack(side='right', fill='y')
        seq_scroll_x.pack(side='bottom', fill='x')

        # Add placeholder text
        placeholder = """>Primer1
ATCGATCGATCGATCG
>Primer2
GCTAGCTAGCTAGCTA"""
        self.sequence_text.insert('1.0', placeholder)
        self.sequence_text.configure(foreground='#999999')

        # Bind events for placeholder behavior
        self.sequence_text.bind('<FocusIn>', self.on_sequence_focus_in)
        self.sequence_text.bind('<FocusOut>', self.on_sequence_focus_out)
        self.placeholder_active = True

        # Analysis parameters frame
        params_frame = ttk.LabelFrame(self.setup_tab, text="Analysis Parameters", padding=15)
        params_frame.pack(pady=10, padx=20, fill='x')

        # Overlap length settings
        ttk.Label(params_frame, text="Overlap Length Range:", style='Heading.TLabel').grid(
            row=0, column=0, sticky='w', pady=(0, 10))

        overlap_frame = ttk.Frame(params_frame)
        overlap_frame.grid(row=1, column=0, sticky='w', pady=(0, 15))

        ttk.Label(overlap_frame, text="Minimum:").grid(row=0, column=0, padx=(0, 5))
        self.min_overlap_var = tk.IntVar(value=3)
        ttk.Spinbox(overlap_frame, from_=2, to=15, width=5, textvariable=self.min_overlap_var).grid(
            row=0, column=1, padx=(0, 20))

        ttk.Label(overlap_frame, text="Maximum:").grid(row=0, column=2, padx=(0, 5))
        self.max_overlap_var = tk.IntVar(value=10)
        ttk.Spinbox(overlap_frame, from_=3, to=20, width=5, textvariable=self.max_overlap_var).grid(
            row=0, column=3)

        # Mismatch settings
        ttk.Label(params_frame, text="Maximum Mismatches:", style='Heading.TLabel').grid(
            row=2, column=0, sticky='w', pady=(15, 5))

        mismatch_frame = ttk.Frame(params_frame)
        mismatch_frame.grid(row=3, column=0, sticky='w', pady=(0, 15))

        self.max_mismatches_var = tk.IntVar(value=1)
        ttk.Spinbox(mismatch_frame, from_=0, to=5, width=5, textvariable=self.max_mismatches_var).grid(
            row=0, column=0, padx=(0, 10))
        ttk.Label(mismatch_frame, text="(0 = perfect matches only)").grid(row=0, column=1)

        # Ambiguity settings
        self.ambiguity_var = tk.BooleanVar(value=False)
        ambiguity_check = ttk.Checkbutton(
            params_frame,
            text="Treat ambiguous bases (R, Y, S, W, K, M, B, D, H, V, N) as matches if any variation could match",
            variable=self.ambiguity_var
        )
        ambiguity_check.grid(row=4, column=0, sticky='w', pady=(10, 0))

        # Action buttons frame
        action_frame = ttk.Frame(self.setup_tab)
        action_frame.pack(pady=20, padx=20, fill='x')

        buttons_frame = ttk.Frame(action_frame)
        buttons_frame.pack()

        ttk.Button(buttons_frame, text="Preview Sequences", command=self.preview_sequences,
                   style='Action.TButton').pack(side='left', padx=(0, 10))

        ttk.Button(buttons_frame, text="Run Analysis", command=self.run_analysis,
                   style='Action.TButton').pack(side='left', padx=(0, 10))

        ttk.Button(buttons_frame, text="Clear Results", command=self.clear_results,
                   style='Action.TButton').pack(side='left')

        # Status frame
        self.status_frame = ttk.Frame(self.setup_tab)
        self.status_frame.pack(pady=10, padx=20, fill='x')

        self.status_var = tk.StringVar(value="Ready")
        self.status_label = ttk.Label(self.status_frame, textvariable=self.status_var)
        self.status_label.pack()

        # Progress bar
        self.progress = ttk.Progressbar(self.status_frame, mode='indeterminate')
        self.progress.pack(fill='x', pady=(5, 0))

    def create_results_tab(self):
        """Create the results display tab"""
        # Results summary frame
        summary_frame = ttk.LabelFrame(self.results_tab, text="Analysis Summary", padding=10)
        summary_frame.pack(pady=10, padx=20, fill='x')

        self.summary_text = tk.Text(summary_frame, height=4, wrap='word', state='disabled',
                                    bg='#f8f9fa', font=('Arial', 10))
        self.summary_text.pack(fill='x')

        # Results display with treeview for better organization
        results_frame = ttk.LabelFrame(self.results_tab, text="Detailed Results", padding=10)
        results_frame.pack(pady=10, padx=20, fill='x')

        # Create treeview for structured results
        tree_frame = ttk.Frame(results_frame)
        tree_frame.pack(fill='x')

        # Treeview with columns - reduced height for more space for detail view
        columns = ('Overlap Length', 'Mismatches', 'Primer 1', 'Primer 2', 'Risk Level')
        self.results_tree = ttk.Treeview(tree_frame, columns=columns, show='headings', height=8)

        # Configure columns
        self.results_tree.heading('Overlap Length', text='Overlap Length')
        self.results_tree.heading('Mismatches', text='Mismatches')
        self.results_tree.heading('Primer 1', text='Primer 1')
        self.results_tree.heading('Primer 2', text='Primer 2')
        self.results_tree.heading('Risk Level', text='Risk Level')

        self.results_tree.column('Overlap Length', width=100)
        self.results_tree.column('Mismatches', width=80)
        self.results_tree.column('Primer 1', width=200)
        self.results_tree.column('Primer 2', width=200)
        self.results_tree.column('Risk Level', width=100)

        # Scrollbars for treeview
        tree_scroll_y = ttk.Scrollbar(tree_frame, orient='vertical', command=self.results_tree.yview)
        tree_scroll_x = ttk.Scrollbar(tree_frame, orient='horizontal', command=self.results_tree.xview)
        self.results_tree.configure(yscrollcommand=tree_scroll_y.set, xscrollcommand=tree_scroll_x.set)

        self.results_tree.pack(side='left', fill='both', expand=True)
        tree_scroll_y.pack(side='right', fill='y')
        tree_scroll_x.pack(side='bottom', fill='x')

        # Detailed output text area - given more space
        detail_frame = ttk.LabelFrame(self.results_tab, text="Detailed Visualization", padding=10)
        detail_frame.pack(pady=10, padx=20, fill='both', expand=True)

        detail_text_frame = ttk.Frame(detail_frame)
        detail_text_frame.pack(fill='both', expand=True)

        self.detail_text = tk.Text(detail_text_frame, wrap='none', font=('Courier', 11),
                                   state='disabled', bg='#f8f9fa', height=12)

        detail_scroll_y = ttk.Scrollbar(detail_text_frame, orient='vertical', command=self.detail_text.yview)
        detail_scroll_x = ttk.Scrollbar(detail_text_frame, orient='horizontal', command=self.detail_text.xview)
        self.detail_text.configure(yscrollcommand=detail_scroll_y.set, xscrollcommand=detail_scroll_x.set)

        self.detail_text.pack(side='left', fill='both', expand=True)
        detail_scroll_y.pack(side='right', fill='y')
        detail_scroll_x.pack(side='bottom', fill='x')

        # Bind treeview selection to show details
        self.results_tree.bind('<<TreeviewSelect>>', self.show_overlap_detail)

        # Bind column headers for sorting
        for col in columns:
            self.results_tree.heading(col, command=lambda c=col: self.sort_treeview(c))

        # Export button
        export_frame = ttk.Frame(self.results_tab)
        export_frame.pack(pady=10, padx=20)

        ttk.Button(export_frame, text="Export Results", command=self.export_results,
                   style='Action.TButton').pack()

    def on_sequence_focus_in(self, event):
        """Clear placeholder text when user focuses on the text field"""
        if self.placeholder_active:
            self.sequence_text.delete('1.0', tk.END)
            self.sequence_text.configure(foreground='#000000')
            self.placeholder_active = False

    def on_sequence_focus_out(self, event):
        """Restore placeholder text if field is empty"""
        content = self.sequence_text.get('1.0', tk.END).strip()
        if not content:
            placeholder = """>Primer1
ATCGATCGATCGATCG
>Primer2
GCTAGCTAGCTAGCTA"""
            self.sequence_text.insert('1.0', placeholder)
            self.sequence_text.configure(foreground='#999999')
            self.placeholder_active = True

    def select_input_file(self):
        """Select input FASTA file and load contents into text field"""
        filename = filedialog.askopenfilename(
            title="Select FASTA file with primer sequences",
            filetypes=[("FASTA files", "*.fasta *.fa *.fas"), ("All files", "*.*")]
        )
        if filename:
            self.input_file_location = filename
            self.input_file_var.set(os.path.basename(filename))
            self.input_label.configure(foreground='#27ae60')

            # Load file contents into the text field
            try:
                with open(filename, 'r', encoding='utf-8-sig') as f:
                    content = f.read()
                # Replace non-breaking spaces (common in files from PDFs/Word)
                content = content.replace('\xa0', ' ')

                # Clear placeholder and insert file contents
                self.sequence_text.delete('1.0', tk.END)
                self.sequence_text.insert('1.0', content)
                self.sequence_text.configure(foreground='#000000')
                self.placeholder_active = False

                # Count sequences for feedback
                sequences = list(SeqIO.parse(StringIO(content), "fasta"))
                self.update_status(f"Loaded {len(sequences)} sequences from file", "success")

            except Exception as e:
                messagebox.showerror("Error", f"Error loading file: {str(e)}")
                self.update_status("Error loading file", "error")

    def get_sequences_from_text(self):
        """Parse sequences from the text field"""
        content = self.sequence_text.get('1.0', tk.END).strip()
        # Replace non-breaking spaces (common when pasting from PDFs/Word/web)
        content = content.replace('\xa0', ' ')

        if not content or self.placeholder_active:
            return None

        try:
            sequences = list(SeqIO.parse(StringIO(content), "fasta"))
            return sequences if sequences else None
        except Exception:
            return None

    def preview_sequences(self):
        """Preview sequences from the text field"""
        sequences = self.get_sequences_from_text()

        if not sequences:
            messagebox.showerror("Error", "No valid FASTA sequences found in the text field.\n\n"
                                "Please paste sequences in FASTA format:\n"
                                ">SequenceName\n"
                                "ATCGATCG...")
            return

        try:
            self.sequences = sequences
            preview_text = f"Found {len(self.sequences)} sequences:\n\n"

            for i, seq in enumerate(self.sequences[:10]):  # Show first 10
                preview_text += f"{seq.id}: {str(seq.seq)[:50]}{'...' if len(seq.seq) > 50 else ''}\n"

            if len(self.sequences) > 10:
                preview_text += f"\n... and {len(self.sequences) - 10} more sequences"

            messagebox.showinfo("Sequence Preview", preview_text)
            self.update_status(f"Found {len(self.sequences)} sequences", "success")

        except Exception as e:
            messagebox.showerror("Error", f"Error parsing sequences: {str(e)}")
            self.update_status("Error parsing sequences", "error")

    def update_status(self, message, status_type="info"):
        """Update status message"""
        self.status_var.set(message)
        if status_type == "success":
            self.status_label.configure(style='Success.TLabel')
        elif status_type == "error":
            self.status_label.configure(style='Error.TLabel')
        else:
            self.status_label.configure(style='TLabel')
        self.root.update()

    def bases_could_match(self, base1, base2):
        """Check if two bases could match considering IUPAC ambiguity codes.

        Returns True if any possible base from base1 matches any possible base from base2.
        """
        bases1 = self.IUPAC_CODES.get(base1.upper(), {base1.upper()})
        bases2 = self.IUPAC_CODES.get(base2.upper(), {base2.upper()})
        # Check if there's any overlap between the possible bases
        return bool(bases1 & bases2)

    def count_mismatches(self, seq1, seq2, consider_ambiguity=False):
        """Count mismatches between two sequences of equal length.

        If consider_ambiguity is True, ambiguous bases that could potentially match
        are counted as matches (not mismatches).
        """
        if consider_ambiguity:
            return sum(1 for a, b in zip(seq1, seq2) if not self.bases_could_match(a, b))
        else:
            return sum(1 for a, b in zip(seq1, seq2) if a != b)

    def get_risk_level(self, overlap_length, mismatches):
        """Determine risk level based on overlap length and mismatches"""
        if mismatches == 0:
            if overlap_length >= 4:
                return "HIGH"
            elif overlap_length >= 2:
                return "MEDIUM"
            else:
                return "LOW"
        elif mismatches == 1 and overlap_length >= 4:
            return "MEDIUM"
        else:
            return "LOW"

    def get_last_n_bases(self, seq_record, n):
        """Get last n bases from 3' end"""
        return str(seq_record.seq[-n:]).upper()

    def get_last_n_bases_rc(self, seq_record, n):
        """Get last n bases from 3' end and return reverse complement"""
        trimmed = seq_record.seq[-n:].upper()
        return str(trimmed.reverse_complement())

    def visualize_overlap(self, seq1, seq2, overlap_length, consider_ambiguity=False):
        """Create ASCII visualization of overlap - original format with both full primers"""
        primer1_full = str(seq1.seq)
        primer2_full = str(seq2.seq)

        # Calculate offset for primer2 alignment
        p2_offset = len(primer1_full) - overlap_length

        # Get the overlapping regions for comparison
        primer1_3end = self.get_last_n_bases(seq1, overlap_length)
        primer2_3end_rc = self.get_last_n_bases_rc(seq2, overlap_length)

        # Create alignment visualization
        lines = []
        lines.append(primer1_full)

        # Create match/mismatch line with proper offset
        match_line = " " * p2_offset
        for i, (a, b) in enumerate(zip(primer1_3end, primer2_3end_rc)):
            if consider_ambiguity:
                # Use ambiguity-aware matching
                if self.bases_could_match(a, b):
                    match_line += "|"
                else:
                    match_line += " "  # Space for mismatch
            else:
                # Exact matching only
                if a == b:
                    match_line += "|"
                else:
                    match_line += " "  # Space for mismatch

        lines.append(match_line)

        # Create primer2 line (reversed) with proper offset
        primer2_line = " " * p2_offset + primer2_full[::-1]
        lines.append(primer2_line)
        lines.append("")

        return "\n".join(lines)

    def run_analysis(self):
        """Run the complete primer compatibility analysis"""
        # Get sequences from text field
        sequences = self.get_sequences_from_text()

        if not sequences:
            messagebox.showerror("Error", "No valid FASTA sequences found in the text field.\n\n"
                                "Please paste sequences in FASTA format or import from a file.")
            return

        # Validate parameters
        min_overlap = self.min_overlap_var.get()
        max_overlap = self.max_overlap_var.get()
        max_mismatches = self.max_mismatches_var.get()
        consider_ambiguity = self.ambiguity_var.get()

        if min_overlap >= max_overlap:
            messagebox.showerror("Error", "Minimum overlap must be less than maximum overlap.")
            return

        self.progress.start()
        self.update_status("Running analysis...", "info")

        try:
            # Use sequences from text field
            self.sequences = sequences

            # Clear previous results
            for item in self.results_tree.get_children():
                self.results_tree.delete(item)

            self.clear_text_widget(self.detail_text)
            self.clear_text_widget(self.summary_text)

            # Store results for export
            self.analysis_results = []
            total_overlaps = 0
            high_risk_count = 0

            # Analysis loop
            for overlap_length in range(max_overlap, min_overlap - 1, -1):
                for mismatches in range(max_mismatches + 1):

                    section_results = []

                    for i in range(len(self.sequences)):
                        for j in range(i, len(self.sequences)):

                            # Get 3' ends
                            primer1_3end = self.get_last_n_bases(self.sequences[i], overlap_length)
                            primer2_3end_rc = self.get_last_n_bases_rc(self.sequences[j], overlap_length)

                            # Count mismatches
                            actual_mismatches = self.count_mismatches(primer1_3end, primer2_3end_rc, consider_ambiguity)

                            if actual_mismatches == mismatches:
                                # Don't skip self-comparisons - we want to detect self-dimers!

                                risk_level = self.get_risk_level(overlap_length, mismatches)

                                # Add to treeview
                                item_id = self.results_tree.insert('', 'end', values=(
                                    overlap_length, mismatches,
                                    self.sequences[i].id, self.sequences[j].id, risk_level
                                ))

                                # Store result data
                                result_data = {
                                    'overlap_length': overlap_length,
                                    'mismatches': mismatches,
                                    'primer1': self.sequences[i],
                                    'primer2': self.sequences[j],
                                    'risk_level': risk_level,
                                    'visualization': self.visualize_overlap(
                                        self.sequences[i], self.sequences[j], overlap_length, consider_ambiguity
                                    )
                                }

                                self.analysis_results.append(result_data)
                                section_results.append(result_data)
                                total_overlaps += 1

                                if risk_level == "HIGH":
                                    high_risk_count += 1

            # Update summary
            summary = f"Analysis Complete!\n"
            summary += f"Total overlaps found: {total_overlaps}\n"
            summary += f"High-risk overlaps: {high_risk_count}\n"
            summary += f"Sequences analyzed: {len(self.sequences)}\n"
            if consider_ambiguity:
                summary += f"Ambiguity handling: ENABLED (ambiguous bases counted as matches if any variation matches)\n"
            summary += f"Risk Levels: HIGH (≥4 perfect matches), MEDIUM (2-3 perfect or ≥4 with 1 mismatch), LOW (other)"

            self.update_text_widget(self.summary_text, summary)

            self.update_status(f"Analysis complete. Found {total_overlaps} potential overlaps.", "success")

        except Exception as e:
            messagebox.showerror("Error", f"Analysis failed: {str(e)}")
            self.update_status("Analysis failed", "error")
        finally:
            self.progress.stop()

    def show_overlap_detail(self, event):
        """Show detailed visualization for selected overlap"""
        selection = self.results_tree.selection()
        if not selection:
            return

        item = selection[0]
        values = self.results_tree.item(item, 'values')

        # Find corresponding result
        overlap_length = int(values[0])
        mismatches = int(values[1])
        primer1_id = values[2]
        primer2_id = values[3]

        for result in self.analysis_results:
            if (result['overlap_length'] == overlap_length and
                    result['mismatches'] == mismatches and
                    result['primer1'].id == primer1_id and
                    result['primer2'].id == primer2_id):
                detail_text = f"Overlap Details:\n"
                detail_text += f"Primers: {primer1_id} vs {primer2_id}\n"
                detail_text += f"Overlap Length: {overlap_length} bases\n"
                detail_text += f"Mismatches: {mismatches}\n"
                detail_text += f"Risk Level: {result['risk_level']}\n\n"
                detail_text += "Visualization:\n"
                detail_text += result['visualization']

                self.update_text_widget(self.detail_text, detail_text)
                break

    def sort_treeview(self, col):
        """Sort treeview by column with alternating sort direction"""
        # Toggle sort direction for this column
        self.sort_reverse[col] = not self.sort_reverse.get(col, False)
        reverse = self.sort_reverse[col]

        # Get all items
        items = [(self.results_tree.set(k, col), k) for k in self.results_tree.get_children('')]

        # Determine sort key function based on column
        if col in ['Overlap Length', 'Mismatches']:
            # Numeric sorting
            items.sort(key=lambda x: int(x[0]), reverse=reverse)
        elif col == 'Risk Level':
            # Custom risk level sorting (HIGH > MEDIUM > LOW)
            risk_order = {'HIGH': 3, 'MEDIUM': 2, 'LOW': 1}
            items.sort(key=lambda x: risk_order.get(x[0], 0), reverse=reverse)
        else:
            # Alphabetical sorting for primer names
            items.sort(key=lambda x: x[0].lower(), reverse=reverse)

        # Rearrange items in sorted order
        for index, (val, k) in enumerate(items):
            self.results_tree.move(k, '', index)

        # Update column headers to show sort direction
        for column in ['Overlap Length', 'Mismatches', 'Primer 1', 'Primer 2', 'Risk Level']:
            if column == col:
                arrow = "↓" if reverse else "↑"
                self.results_tree.heading(column, text=f"{column} {arrow}")
            else:
                self.results_tree.heading(column, text=column)

    def clear_text_widget(self, widget):
        """Clear text widget content"""
        widget.config(state='normal')
        widget.delete(1.0, tk.END)
        widget.config(state='disabled')

    def update_text_widget(self, widget, text):
        """Update text widget content"""
        widget.config(state='normal')
        widget.delete(1.0, tk.END)
        widget.insert(1.0, text)
        widget.config(state='disabled')

    def clear_results(self):
        """Clear all results"""
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)
        self.clear_text_widget(self.detail_text)
        self.clear_text_widget(self.summary_text)
        self.analysis_results = []
        self.sort_reverse = {}  # Reset sort directions

        # Reset column headers
        for column in ['Overlap Length', 'Mismatches', 'Primer 1', 'Primer 2', 'Risk Level']:
            self.results_tree.heading(column, text=column)

        self.update_status("Results cleared", "info")

    def export_results(self):
        """Export results to CSV file"""
        if not hasattr(self, 'analysis_results') or not self.analysis_results:
            messagebox.showwarning("Warning", "No results to export.")
            return

        filename = filedialog.asksaveasfilename(
            title="Export Results",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )

        if filename:
            try:
                import csv
                with open(filename, 'w', newline='') as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerow(['Overlap Length', 'Mismatches', 'Primer 1', 'Primer 2', 'Risk Level'])

                    for result in self.analysis_results:
                        writer.writerow([
                            result['overlap_length'],
                            result['mismatches'],
                            result['primer1'].id,
                            result['primer2'].id,
                            result['risk_level']
                        ])

                messagebox.showinfo("Export Complete", f"Results exported to {filename}")

            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to export results: {str(e)}")


def main():
    root = tk.Tk()
    app = PrimerCompatibilityAnalyzer(root)
    root.mainloop()


if __name__ == "__main__":
    main()