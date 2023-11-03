import os
import gzip
import regex
import typer
from typing import Dict, Tuple
import matplotlib.pyplot as plt

app = typer.Typer()

# Define the dictionary of sequences to match with barcodes and adapters
sequence_dict = {
    'SQK-PCB111.24-FLANK1': 'ATCGCCTACCGTGA',
    'SQK-PCB111.24-FLANK2': 'TTGCCTGTCGCTCTATCTTC',
    'SQK-PCB111.24-FLANK3': 'TCTGTTGGTGCTGATATTGC',
    'SQK-PCB111.24-CRTA': 'CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG',
    'SQK-PCB111.24-RTP': 'CTTGCCTGTCGCTCTATCTTCAGAGGAG',
    'SQK-PCB111.24-RAPT': 'TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT',
    'SQK-PCB111.24-SSPII': 'TTTCTGTTGGTGCTGATATTGCTTT',
    'BP01_Forward': 'CAAGAAAGTTGTCGGTGTCTTTGTGAC',
    'BP01_Reverse': 'CAAGAAAGTTGTCGGTGTCTTTGTGTT',
    'BP02_Forward': 'CTCGATTCCGTTTGTAGTCGTCTGTAC',
    'BP02_Reverse': 'CTCGATTCCGTTTGTAGTCGTCTGTTT',
    'BP03_Forward': 'CGAGTCTTGTGTCCCAGTTACCAGGAC',
    'BP03_Reverse': 'CGAGTCTTGTGTCCCAGTTACCAGGTT',
    'BP04_Forward': 'CTTCGGATTCTATCGTGTTTCCCTAAC',
    'BP04_Reverse': 'CTTCGGATTCTATCGTGTTTCCCTATT',
}

"""
    'BP05_Forward': 'CCTTGTCCAGGGTTTGTGTAACCTTAC',
    'BP05_Reverse': 'CCTTGTCCAGGGTTTGTGTAACCTTTT',
    'BP06_Forward': 'CTTCTCGCAAAGGCAGAAAGTAGTCAC',
    'BP06_Reverse': 'CTTCTCGCAAAGGCAGAAAGTAGTCTT',
    'BP07_Forward': 'CGTGTTACCGTGGGAATGAATCCTTAC',
    'BP07_Reverse': 'CGTGTTACCGTGGGAATGAATCCTTTT',
    'BP08_Forward': 'CTTCAGGGAACAAACCAAGTTACGTAC',
    'BP08_Reverse': 'CTTCAGGGAACAAACCAAGTTACGTTT',
    'BP09_Forward': 'CAACTAGGCACAGCGAGTCTTGGTTAC',
    'BP09_Reverse': 'CAACTAGGCACAGCGAGTCTTGGTTTT',
    'BP10_Forward': 'CAAGCGTTGAAACCTTTGTCCTCTCAC',
    'BP10_Reverse': 'CAAGCGTTGAAACCTTTGTCCTCTCTT',
    'BP11_Forward': 'CGTTTCATCTATCGGAGGGAATGGAAC',
    'BP11_Reverse': 'CGTTTCATCTATCGGAGGGAATGGATT',
    'BP12_Forward': 'CCAGGTAGAAAGAAGCAGAATCGGAAC',
    'BP12_Reverse': 'CCAGGTAGAAAGAAGCAGAATCGGATT',
    'BP13_Forward': 'CAGAACGACTTCCATACTCGTGTGAAC',
    'BP13_Reverse': 'CAGAACGACTTCCATACTCGTGTGATT',
    'BP14_Forward': 'CAACGAGTCTCTTGGGACCCATAGAAC',
    'BP14_Reverse': 'CAACGAGTCTCTTGGGACCCATAGATT',
    'BP15_Forward': 'CAGGTCTACCTCGCTAACACCACTGAC',
    'BP15_Reverse': 'CAGGTCTACCTCGCTAACACCACTGTT',
    'BP16_Forward': 'CCGTCAACTGACAGTGGTTCGTACTAC',
    'BP16_Reverse': 'CCGTCAACTGACAGTGGTTCGTACTTT',
    'BP17_Forward': 'CACCCTCCAGGAAAGTACCTCTGATAC',
    'BP17_Reverse': 'CACCCTCCAGGAAAGTACCTCTGATTT',
    'BP18_Forward': 'CCCAAACCCAACAACCTAGATAGGCAC',
    'BP18_Reverse': 'CCCAACCCAACAACCTAGATAGGCTT',
    'BP19_Forward': 'CGTTCCTCGTGCAGTGTCAAGAGATAC',
    'BP19_Reverse': 'CGTTCCTCGTGCAGTGTCAAGAGATTT',
    'BP20_Forward': 'CTTGCGTCCTGTTACGAGAACTCATAC',
    'BP20_Reverse': 'CTTGCGTCCTGTTACGAGAACTCATTT',
    'BP21_Forward': 'CGAGCCTCTCATTGTCCGTTCTCTAAC',
    'BP21_Reverse': 'CGAGCCTCTCATTGTCCGTTCTCTATT',
    'BP22_Forward': 'CACCACTGCCATGTATCAAAGTACGAC',
    'BP22_Reverse': 'CACCACTGCCATGTATCAAAGTACGTT',
    'BP23_Forward': 'CCTTACTACCCAGTGAACCTCCTCGAC',
    'BP23_Reverse': 'CCTTACTACCCAGTGAACCTCCTCGTT',
    'BP24_Forward': 'CGCATAGTTCTGCATGATGGGTTAGAC',
    'BP24_Reverse': 'CGCATAGTTCTGCATGATGGGTTAGTT',
}
"""
@app.command()
def process_directory(
    directory: str = typer.Argument(..., help="Directory containing FASTQ files"),
    mismatch_threshold: float = typer.Option(0.1, help="Mismatch threshold as a decimal")
):
    report = process_files(directory, sequence_dict, mismatch_threshold)
    
    # Output the report and create plots
    for file, data in report.items():
        print(f"File: {file}")
        print(f"Total Reads: {data['total_reads']}")
        for seq_name, count in data['sequence_counts'].items():
            percent = (count / data['total_reads']) * 100 if data['total_reads'] > 0 else 0
            print(f"{seq_name}: {count} ({percent:.2f}%)")
        
        # Call the plotting function
        base_filename = os.path.splitext(os.path.basename(file))[0]
        create_and_save_plot(data['sequence_counts'], data['total_reads'], base_filename)

def find_fastq_files(path: str):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('.fastq.gz'):
                yield os.path.join(root, file)

def count_matches(
    fastq_file: str,
    sequence_dict: Dict[str, str],
    mismatch_threshold: float
) -> Tuple[int, Dict[str, int]]:
    total_reads = 0
    sequence_counts = {key: 0 for key in sequence_dict.keys()}
    with gzip.open(fastq_file, 'rt') as file:
        for line in file:
            if line.startswith('@'):
                total_reads += 1
                sequence = next(file).strip()
                for name, barcode in sequence_dict.items():
                    max_errors = max(1, round(len(barcode) * mismatch_threshold))
                    match = regex.search(f'({barcode}){{e<={max_errors}}}', sequence)
                    if match:
                        sequence_counts[name] += 1
                next(file)  # Skip the '+' line
                next(file)  # Skip the quality line
    return total_reads, sequence_counts

def process_files(
    directory: str,
    sequence_dict: Dict[str, str],
    mismatch_threshold: float
) -> Dict[str, Dict]:
    report = {}
    for fastq_file in find_fastq_files(directory):
        total_reads, sequence_counts = count_matches(fastq_file, sequence_dict, mismatch_threshold)
        report[fastq_file] = {
            'total_reads': total_reads,
            'sequence_counts': sequence_counts,
        }
    return report

def create_and_save_plot(data: Dict[str, int], total_reads: int, filename: str):
    sequences = list(data.keys())
    percents = [(count / total_reads) * 100 for count in data.values()]
    plt.figure(figsize=(10, 8))
    plt.bar(sequences, percents, color='skyblue')
    plt.title(f'{filename}')
    plt.xlabel('Component Sequences Detected')
    plt.ylabel(f'Percentage of {total_reads} Total Reads (%)')
    plt.xticks(rotation=90)
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(f"{filename}.png")
    plt.close()

if __name__ == "__main__":
    app()