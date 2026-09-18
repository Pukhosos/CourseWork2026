from argparse import ArgumentParser, Namespace
from collections import defaultdict
from collections.abc import Iterable, Sequence
from csv import reader
from dataclasses import dataclass
from pathlib import Path
from typing import Final


BLAST_COLUMN_COUNT: Final[int] = 13
OUTPUT_DELIMITER: Final[str] = '\t'


@dataclass(frozen=True, slots=True)
class BlastHit:
    """Represent one HSP from an NCBI BLAST hit-table export."""

    query_id: str
    subject_id: str
    identity_percentage: float
    alignment_length: int
    mismatch_count: int
    gap_opening_count: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    e_value: float
    bit_score: float
    positives_percentage: float


@dataclass(frozen=True, slots=True)
class QueryHit:
    """Represent the best HSP for one query against one chain."""

    query_id: str
    subject_id: str
    identity_percentage: float
    alignment_length: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    e_value: float
    bit_score: float
    positives_percentage: float


@dataclass(frozen=True, slots=True)
class DeduplicatedHit:
    """Represent one unique chain discovered by one or more queries."""

    subject_id: str
    query_ids: tuple[str, ...]
    best_query_id: str
    best_identity_percentage: float
    best_alignment_length: int
    best_query_start: int
    best_query_end: int
    best_subject_start: int
    best_subject_end: int
    best_e_value: float
    best_bit_score: float
    best_positives_percentage: float


def parse_blast_row(row: Sequence[str], source: Path, line_number: int) -> BlastHit:
    """Parse one row of an NCBI BLAST hit table."""

    if len(row) != BLAST_COLUMN_COUNT:
        raise ValueError(
            f'{source}:{line_number}: expected {(
                BLAST_COLUMN_COUNT
            )} columns, got {len(row)}'
        )
    return BlastHit(
        query_id=row[0],
        subject_id=row[1],
        identity_percentage=float(row[2]),
        alignment_length=int(row[3]),
        mismatch_count=int(row[4]),
        gap_opening_count=int(row[5]),
        query_start=int(row[6]),
        query_end=int(row[7]),
        subject_start=int(row[8]),
        subject_end=int(row[9]),
        e_value=float(row[10]),
        bit_score=float(row[11]),
        positives_percentage=float(row[12]),
    )


def read_blast_hits(path: Path) -> list[BlastHit]:
    """Read all HSPs from one NCBI BLAST hit-table file."""

    hits: list[BlastHit] = []
    with path.open('r', encoding='utf-8', newline='') as input_file:
        rows: Iterable[list[str]] = reader(input_file)
        for line_number, row in enumerate(rows, start=1):
            if len(row) == 0:
                continue
            hit: BlastHit = parse_blast_row(row, path, line_number)
            hits.append(hit)
    return hits


def select_best_hsp(hits: Iterable[BlastHit]) -> BlastHit:
    """Select the strongest HSP using bit score and E-value."""

    hits_tuple: tuple[BlastHit, ...] = tuple(hits)
    if len(hits_tuple) == 0:
        raise ValueError('Cannot select the best HSP from an empty collection')
    return max(
        hits_tuple,
        key=(
            lambda hit: (
                hit.bit_score,
                -hit.e_value,
                hit.alignment_length,
                hit.identity_percentage,
            )
        ),
    )


def collapse_query_hsps(hits: Iterable[BlastHit]) -> list[QueryHit]:
    """Collapse multiple HSPs for each query/chain pair."""

    grouped_hits: defaultdict[tuple[str, str], list[BlastHit]] = defaultdict(list)
    for hit in hits:
        key: tuple[str, str] = (hit.query_id, hit.subject_id)
        grouped_hits[key].append(hit)
    collapsed_hits: list[QueryHit] = []
    for grouped_hit_list in grouped_hits.values():
        best_hit: BlastHit = select_best_hsp(grouped_hit_list)
        collapsed_hit: QueryHit = QueryHit(
            query_id=best_hit.query_id,
            subject_id=best_hit.subject_id,
            identity_percentage=best_hit.identity_percentage,
            alignment_length=best_hit.alignment_length,
            query_start=best_hit.query_start,
            query_end=best_hit.query_end,
            subject_start=best_hit.subject_start,
            subject_end=best_hit.subject_end,
            e_value=best_hit.e_value,
            bit_score=best_hit.bit_score,
            positives_percentage=best_hit.positives_percentage,
        )
        collapsed_hits.append(collapsed_hit)
    return collapsed_hits


def select_best_query_hit(hits: Iterable[QueryHit]) -> QueryHit:
    """Select the strongest query hit for one chain."""

    hits_tuple: tuple[QueryHit, ...] = tuple(hits)
    if len(hits_tuple) == 0:
        raise ValueError('Cannot select the best query hit from an empty collection')
    return max(
        hits_tuple,
        key=(
            lambda hit: (
                hit.bit_score,
                -hit.e_value,
                hit.alignment_length,
                hit.identity_percentage,
            )
        ),
    )


def deduplicate_chains(hits: Iterable[QueryHit]) -> list[DeduplicatedHit]:
    """Merge occurrences of the same chain across different queries."""

    grouped_hits: defaultdict[str, list[QueryHit]] = defaultdict(list)
    for hit in hits:
        grouped_hits[hit.subject_id].append(hit)
    deduplicated_hits: list[DeduplicatedHit] = []
    for subject_id, subject_hits in grouped_hits.items():
        best_hit: QueryHit = select_best_query_hit(subject_hits)
        query_ids: tuple[str, ...] = tuple(
            sorted({hit.query_id for hit in subject_hits})
        )
        deduplicated_hit: DeduplicatedHit = DeduplicatedHit(
            subject_id=subject_id,
            query_ids=query_ids,
            best_query_id=best_hit.query_id,
            best_identity_percentage=best_hit.identity_percentage,
            best_alignment_length=best_hit.alignment_length,
            best_query_start=best_hit.query_start,
            best_query_end=best_hit.query_end,
            best_subject_start=best_hit.subject_start,
            best_subject_end=best_hit.subject_end,
            best_e_value=best_hit.e_value,
            best_bit_score=best_hit.bit_score,
            best_positives_percentage=best_hit.positives_percentage,
        )
        deduplicated_hits.append(deduplicated_hit)
    return sorted(
        deduplicated_hits,
        key=(lambda hit: (-hit.best_bit_score, hit.subject_id)),
    )


def format_number(value: float) -> str:
    """Format a floating-point BLAST statistic compactly."""

    return f'{value:g}'


def write_hits(path: Path, hits: Iterable[DeduplicatedHit]) -> None:
    """Write unique chain candidates as a TSV file."""

    header: Final[str] = OUTPUT_DELIMITER.join(
        (
            'chain',
            'queries',
            'query_count',
            'best_query',
            'identity_percentage',
            'alignment_length',
            'query_start',
            'query_end',
            'subject_start',
            'subject_end',
            'e_value',
            'bit_score',
            'positives_percentage',
        )
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as output_file:
        output_file.write(f'{header}\n')
        for hit in hits:
            values: tuple[str, ...] = (
                hit.subject_id,
                ','.join(hit.query_ids),
                str(len(hit.query_ids)),
                hit.best_query_id,
                format_number(hit.best_identity_percentage),
                str(hit.best_alignment_length),
                str(hit.best_query_start),
                str(hit.best_query_end),
                str(hit.best_subject_start),
                str(hit.best_subject_end),
                format_number(hit.best_e_value),
                format_number(hit.best_bit_score),
                format_number(hit.best_positives_percentage),
            )
            output_file.write(f'{OUTPUT_DELIMITER.join(values)}\n')


def build_argument_parser() -> ArgumentParser:
    """Build the command-line argument parser."""

    parser: ArgumentParser = ArgumentParser(
        description='Deduplicate chains from NCBI BLAST hit tables.'
    )
    parser.add_argument(
        'inputs',
        nargs='+',
        type=Path,
        help='NCBI BLAST Hit Table CSV files.',
    )
    parser.add_argument(
        '--output',
        required=True,
        type=Path,
        help='Output TSV file.',
    )
    return parser


def main() -> None:
    parser: ArgumentParser = build_argument_parser()
    arguments: Namespace = parser.parse_args()
    input_paths: list[Path] = arguments.inputs
    output_path: Path = arguments.output
    all_hits: list[BlastHit] = []
    for input_path in input_paths:
        input_hits: list[BlastHit] = read_blast_hits(input_path)
        all_hits.extend(input_hits)
    query_hits: list[QueryHit] = collapse_query_hsps(all_hits)
    deduplicated_hits: list[DeduplicatedHit] = deduplicate_chains(query_hits)
    write_hits(output_path, deduplicated_hits)
    print(
        '; '.join(
            (
                f'Read {len(all_hits)} HSPs',
                f'collapsed to {len(query_hits)} query-chain hits',
                f'found {len(deduplicated_hits)} unique chains.',
            )
        )
    )
