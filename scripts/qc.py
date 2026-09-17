from argparse import ArgumentParser, Namespace
from collections import Counter
from csv import DictReader, DictWriter
from dataclasses import dataclass, field
from pathlib import Path
from statistics import quantiles
from typing import Final, Literal

VALID_AMINO_ACIDS: Final[frozenset[str]] = frozenset('ABCDEFGHIKLMNPQRSTVWXYZ*')
MAXIMUM_X_FRACTION: Final[float] = 0.05
MAXIMUM_FRAGMENT_PROFILE_COVERAGE: Final[float] = 0.50
OUTPUT_FIELDS: Final[tuple[str, ...]] = (
    'protein_coverage',
    'profile_coverage',
    'length_outlier',
    'x_fraction',
    'qc_flags',
    'qc_status',
)


@dataclass(frozen=True, slots=True)
class FastaRecord:
    identifier: str
    header: str
    sequence: str


@dataclass(frozen=True, slots=True)
class LengthRange:
    lower: float
    upper: float

    def contains(self, length: int) -> bool:
        """Return whether a protein length is inside the accepted range."""

        return self.lower <= length <= self.upper


@dataclass(frozen=True, slots=True)
class QcResult:
    status: Literal['exclude', 'review', 'keep']
    flags: list[str]
    protein_coverage: float
    profile_coverage: float
    x_fraction: float


@dataclass(slots=True)
class FastaReader:
    """Read FASTA records while maintaining parser state."""

    __records: dict[str, FastaRecord] = field(default_factory=dict, init=False)
    __header: str | None = field(default=None, init=False)
    __sequence_parts: list[str] = field(default_factory=list, init=False)

    def read(self, path: Path) -> dict[str, FastaRecord]:
        """Read FASTA records indexed by the first header token."""

        with path.open(encoding='utf-8') as file:
            for line in file:
                stripped_line: str = line.strip()

                if stripped_line.startswith('>'):
                    self.__store_current_record()
                    self.__header = stripped_line[1:]
                    self.__sequence_parts = []
                elif stripped_line:
                    self.__sequence_parts.append(stripped_line)

        self.__store_current_record()
        return self.__records

    def __store_current_record(self) -> None:
        if self.__header is None:
            return None
        identifier: str = self.__header.split(maxsplit=1)[0]
        if identifier in self.__records:
            raise ValueError(f'Duplicate FASTA identifier: {identifier}')
        sequence: str = ''.join(self.__sequence_parts).upper()
        self.__records[identifier] = FastaRecord(
            identifier=identifier,
            header=self.__header,
            sequence=sequence,
        )


def parse_arguments() -> Namespace:
    """Parse command-line arguments."""

    parser: ArgumentParser = ArgumentParser(
        description='Quality-control an extracted COG protein dataset.'
    )
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--sequences', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    arguments: Namespace = parser.parse_args()
    return arguments


def read_metadata(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    """Read tab-separated protein metadata."""

    with path.open(encoding='utf-8', newline='') as file:
        reader: DictReader[str] = DictReader(file, delimiter='\t')
        if reader.fieldnames is None:
            raise ValueError('Metadata file has no header.')
        field_names: list[str] = list(reader.fieldnames)
        rows: list[dict[str, str]] = [
            {key: value or '' for key, value in row.items()} for row in reader
        ]
    return rows, field_names


def coordinate_length(coordinates: str) -> int:
    """Calculate total length of a possibly segmented footprint."""

    total: int = 0
    for segment in coordinates.split('='):
        start_str, end_str = segment.split('-', maxsplit=1)
        total += int(end_str) - int(start_str) + 1
    return total


def calculate_length_range(lengths: list[int], iqr: float = 1.5) -> LengthRange:
    """Calculate Tukey `iqr`-IQR bounds."""

    if len(lengths) < 2:
        return LengthRange(lower=float('-inf'), upper=float('inf'))
    quartiles: list[float] = quantiles(lengths, n=4, method='inclusive')
    first_quartile: float = quartiles[0]
    third_quartile: float = quartiles[2]
    interquartile_range: float = third_quartile - first_quartile
    return LengthRange(
        lower=first_quartile - iqr * interquartile_range,
        upper=third_quartile + iqr * interquartile_range,
    )


def validate_sequence(sequence: str) -> list[str]:
    """Return sequence-integrity QC flags."""

    flags: list[str] = []
    if len(sequence) == 0:
        flags.append('empty_sequence')
        return flags
    invalid_residues: set[str] = set(sequence) - VALID_AMINO_ACIDS
    if len(invalid_residues) > 0:
        flags.append('invalid_residues')
    if '*' in sequence[:-1]:
        flags.append('internal_stop')
    return flags


def classify_record(
    row: dict[str, str],
    sequence: str,
    length_range: LengthRange,
) -> QcResult:
    """Calculate QC metrics and assign a QC status."""

    flags: list[str] = validate_sequence(sequence)
    protein_length: int = int(row['protein_length'])
    footprint_length: int = int(row['cog_footprint_length'])
    profile_length: int = int(row['cog_profile_length'])
    membership_class: int = int(row['membership_class'])
    profile_footprint_length: int = coordinate_length(row['profile_footprint'])
    protein_coverage: float = footprint_length / protein_length
    profile_coverage: float = profile_footprint_length / profile_length
    x_fraction: float = sequence.count('X') / len(sequence) if sequence else 0.0
    length_outlier: bool = not length_range.contains(protein_length)
    if len(sequence) != protein_length:
        flags.append('length_mismatch')
    if membership_class in {2, 3}:
        flags.append(f'membership_class_{membership_class}')
    if length_outlier:
        flags.append('length_outlier')
    if x_fraction > MAXIMUM_X_FRACTION:
        flags.append('many_unknown_residues')
    probable_fragment: bool = (
        membership_class in {2, 3}
        and profile_coverage < MAXIMUM_FRAGMENT_PROFILE_COVERAGE
        and protein_length < length_range.lower
    )
    if probable_fragment:
        flags.append('probable_fragment')
    exclusion_flags: Final[frozenset[str]] = frozenset(
        {
            'empty_sequence',
            'invalid_residues',
            'internal_stop',
            'length_mismatch',
            'probable_fragment',
        }
    )
    return QcResult(
        status=(
            'exclude' if exclusion_flags.intersection(flags) else (
                'review' if len(flags) > 0 else 'keep'
            )
        ),
        flags=flags,
        protein_coverage=protein_coverage,
        profile_coverage=profile_coverage,
        x_fraction=x_fraction,
    )


def write_fasta(
    path: Path,
    records: list[FastaRecord],
) -> None:
    """Write FASTA records."""

    with path.open('w', encoding='utf-8') as file:
        for record in records:
            file.write(f'>{record.header}\n')
            chunks: list[str] = [
                record.sequence[position : position + 80]
                for position in range(0, len(record.sequence), 80)
            ]
            for chunk in chunks:
                file.write(f'{chunk}\n')


def write_table(
    path: Path,
    rows: list[dict[str, str]],
    field_names: list[str],
) -> None:
    """Write tab-separated QC metadata."""

    with path.open('w', encoding='utf-8', newline='') as file:
        writer: DictWriter[str] = DictWriter(
            file,
            fieldnames=field_names,
            delimiter='\t',
            lineterminator='\n',
        )
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    """Write a compact QC summary."""

    statuses: Counter[str] = Counter(row['qc_status'] for row in rows)
    flags: Counter[str] = Counter(
        flag for row in rows for flag in row['qc_flags'].split(';') if flag
    )
    with path.open('w', encoding='utf-8') as file:
        file.write(f'Proteins: {len(rows)}\n\n')
        file.write('QC status:\n')
        for status, count in sorted(statuses.items()):
            file.write(f'  {status}: {count}\n')
        file.write('\nFlags:\n')
        for flag, count in sorted(flags.items()):
            file.write(f'  {flag}: {count}\n')


def process_dataset(
    metadata_path: Path,
    sequences_path: Path,
    output_path: Path,
) -> None:
    """Run quality control and write all output files."""

    metadata_result: tuple[list[dict[str, str]], list[str]] = read_metadata(
        path=metadata_path
    )
    rows: list[dict[str, str]] = metadata_result[0]
    field_names: list[str] = metadata_result[1]
    sequences: dict[str, FastaRecord] = FastaReader().read(sequences_path)
    protein_lengths: list[int] = [int(row['protein_length']) for row in rows]
    length_range: LengthRange = calculate_length_range(protein_lengths)
    output_path.mkdir(parents=True, exist_ok=True)
    grouped_sequences: dict[str, list[FastaRecord]] = {
        'keep': [],
        'review': [],
        'exclude': [],
    }
    processed_rows: list[dict[str, str]] = []
    for row in rows:
        sequence_identifier: str = row['sequence_id']
        if sequence_identifier not in sequences:
            raise ValueError(f'Protein {sequence_identifier!r} is missing from FASTA.')
        record: FastaRecord = sequences[sequence_identifier]
        qc_result: QcResult = classify_record(
            row=row,
            sequence=record.sequence,
            length_range=length_range,
        )
        processed_row: dict[str, str] = {
            **row,
            'protein_coverage': f'{qc_result.protein_coverage:.4f}',
            'profile_coverage': f'{qc_result.profile_coverage:.4f}',
            'length_outlier': str(
                not length_range.contains(int(row['protein_length']))
            ).lower(),
            'x_fraction': f'{qc_result.x_fraction:.4f}',
            'qc_flags': ';'.join(qc_result.flags),
            'qc_status': qc_result.status,
        }
        processed_rows.append(processed_row)
        grouped_sequences[qc_result.status].append(record)
    output_field_names: list[str] = [*field_names, *OUTPUT_FIELDS]
    write_table(
        path=output_path / 'qc.tsv',
        rows=processed_rows,
        field_names=output_field_names,
    )
    for status, records in grouped_sequences.items():
        selected_rows: list[dict[str, str]] = [
            row for row in processed_rows if row['qc_status'] == status
        ]
        write_table(
            path=output_path / f'{status}.tsv',
            rows=selected_rows,
            field_names=output_field_names,
        )
        write_fasta(path=output_path / f'{status}.faa', records=records)
    write_summary(path=output_path / 'summary.txt', rows=processed_rows)


def main() -> None:
    arguments: Namespace = parse_arguments()
    process_dataset(
        metadata_path=arguments.metadata,
        sequences_path=arguments.sequences,
        output_path=arguments.output,
    )
