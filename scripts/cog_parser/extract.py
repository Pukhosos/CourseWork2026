from argparse import ArgumentParser, Namespace
from collections import defaultdict
from collections.abc import Iterable
from csv import writer
from dataclasses import dataclass
from pathlib import Path
from re import compile, Pattern  # noqa: A004
from sys import stderr
from textwrap import wrap
from typing import Final, TextIO

from .parsers import (
    CogAssignment,
    CogAssignmentFile,
    CogProtein,
    CogProteinFile,
    Footprint,
    GenomeCollection,
    GenomeRecord,
)

COG_ID_PATTERN: Final[Pattern[str]] = compile(r'^COG\d{4}$')
FASTA_LINE_WIDTH: Final[int] = 80

METADATA_FILENAME_SUFFIX: Final[str] = '_275_raw.tsv'
FASTA_FILENAME_SUFFIX: Final[str] = '_275_raw.faa'
AMBIGUOUS_FILENAME_SUFFIX: Final[str] = '_275_ambiguous.tsv'
UNMATCHED_FILENAME_SUFFIX: Final[str] = '_275_unmatched.tsv'


@dataclass(frozen=True, slots=True)
class Arguments:
    """Command-line arguments for dataset extraction."""

    cog_id: str
    genomes_path: Path
    assignments_path: Path
    sequences_path: Path
    output_directory: Path
    suppress_errors: bool


@dataclass(frozen=True, slots=True)
class SelectedAssignment:
    """A COG assignment linked to its source genome metadata."""

    assignment: CogAssignment
    genome: GenomeRecord


@dataclass(frozen=True, slots=True)
class ExtractionStatistics:
    """Summary statistics produced by the extraction pipeline."""

    representative_genomes: int
    selected_assignments: int
    selected_genomes: int
    matched_sequences: int
    ambiguous_sequences: int
    duplicate_fasta_matches: int
    unmatched_assignments: int


def require_file(path: Path, description: str) -> None:
    """Require a command-line input path to be an existing file."""

    if not path.is_file():
        raise FileNotFoundError(f'{description.capitalize()} does not exist: {path}')


def parse_arguments() -> Arguments:
    """Parse and validate command-line arguments."""

    parser: ArgumentParser = ArgumentParser(
        description=' '.join(
            (
                'Extract proteins assigned to a requested COG from the',
                'representative 275-genome dataset.',
            )
        ),
    )
    parser.add_argument(
        '--cog',
        required=True,
        help='COG identifier, for example COG0258.',
    )
    parser.add_argument(
        '--genomes',
        required=True,
        type=Path,
        help='Path to 27789-154191-1-SP.xlsx.',
    )
    parser.add_argument(
        '--assignments',
        required=True,
        type=Path,
        help='Path to cog-24.cog.csv.',
    )
    parser.add_argument(
        '--sequences',
        required=True,
        type=Path,
        help='Path to COGorg24.faa.',
    )
    parser.add_argument(
        '--output',
        required=True,
        type=Path,
        help='Output directory.',
    )
    parser.add_argument(
        '--suppress-errors',
        action='store_true',
        help='Suppress error messages.',
    )
    namespace: Namespace = parser.parse_args()
    cog_id: str = str(namespace.cog).upper()
    if COG_ID_PATTERN.fullmatch(cog_id) is None:
        parser.error(f'Invalid COG identifier {cog_id!r}; expected e.g. COG0258.')
    genomes_path: Path = Path(namespace.genomes).expanduser().resolve()
    assignments_path: Path = Path(namespace.assignments).expanduser().resolve()
    sequences_path: Path = Path(namespace.sequences).expanduser().resolve()
    output_directory: Path = Path(namespace.output).expanduser().resolve()
    require_file(genomes_path, 'genome workbook')
    require_file(assignments_path, 'COG assignment database')
    require_file(sequences_path, 'COG protein sequence database')
    return Arguments(
        cog_id=cog_id,
        genomes_path=genomes_path,
        assignments_path=assignments_path,
        sequences_path=sequences_path,
        output_directory=output_directory,
        suppress_errors=namespace.suppress_errors,
    )


def representative_genomes_by_assembly(
    genomes: GenomeCollection,
) -> dict[str, GenomeRecord]:
    """Index representative genomes by NCBI assembly accession."""

    result: dict[str, GenomeRecord] = {}
    genome: GenomeRecord
    for genome in genomes.representatives:
        assembly_accession: str = genome.assembly_accession
        if assembly_accession in result:
            raise ValueError(
                f'Duplicate representative assembly accession: {assembly_accession!r}.'
            )
        result[assembly_accession] = genome
    return result


def select_assignments(
    assignment_file: CogAssignmentFile,
    representative_genomes: dict[str, GenomeRecord],
    cog_id: str,
) -> tuple[SelectedAssignment, ...]:
    """Select requested COG assignments from representative genomes."""

    selected: list[SelectedAssignment] = []
    for assignment in assignment_file:
        if assignment.cog_id != cog_id:
            continue
        genome: GenomeRecord | None = representative_genomes.get(assignment.assembly_id)
        if genome is None:
            continue
        selected.append(SelectedAssignment(assignment=assignment, genome=genome))
    return tuple(selected)


def index_assignments_by_protein_id(
    selected_assignments: Iterable[SelectedAssignment],
) -> dict[str, tuple[SelectedAssignment, ...]]:
    """Index selected assignments without assuming protein IDs are unique."""

    temporary_index: dict[str, list[SelectedAssignment]] = defaultdict(list)
    for selected_assignment in selected_assignments:
        temporary_index[selected_assignment.assignment.protein_id].append(
            selected_assignment
        )
    return {
        protein_id: tuple(assignments)
        for protein_id, assignments in temporary_index.items()
    }


def match_protein(
    protein: CogProtein,
    candidates: tuple[SelectedAssignment, ...],
) -> SelectedAssignment | None:
    """
    Resolve a FASTA record against candidate COG assignments.

    Protein ID is used to obtain the candidate set. Protein length is then
    required to agree. If several candidates remain, the first token of the
    FASTA description is compared with the COG gene ID.
    """

    length_matches: tuple[SelectedAssignment, ...] = tuple(
        candidate
        for candidate in candidates
        if candidate.assignment.protein_length == protein.length
    )
    if len(length_matches) == 1:
        return length_matches[0]
    if len(length_matches) == 0:
        return None
    gene_id: str | None = protein_header_gene_id(protein)
    if gene_id is None:
        return None
    gene_matches: tuple[SelectedAssignment, ...] = tuple(
        candidate
        for candidate in length_matches
        if candidate.assignment.gene_id == gene_id
    )
    if len(gene_matches) == 1:
        return gene_matches[0]
    return None


def protein_header_gene_id(protein: CogProtein) -> str | None:
    """Extract the likely locus tag from the FASTA description."""

    if protein.description is None:
        return None
    description_parts: list[str] = protein.description.split(maxsplit=1)
    if len(description_parts) == 0:
        return None
    return description_parts[0]


def footprint_to_string(footprint: Footprint) -> str:
    """Serialize a COG footprint using the source database notation."""

    return '='.join(
        f'{coordinate_range.start}-{coordinate_range.end}'
        for coordinate_range in footprint.ranges
    )


def assignment_identity(
    selected_assignment: SelectedAssignment,
) -> tuple[str, str, str, str]:
    """Return an identity key for one selected assignment."""

    assignment: CogAssignment = selected_assignment.assignment
    return (
        assignment.assembly_id,
        assignment.gene_id,
        assignment.protein_id,
        assignment.cog_id,
    )


def fasta_identifier(selected_assignment: SelectedAssignment) -> str:
    """Create an unambiguous identifier for extracted FASTA records."""

    assignment: CogAssignment = selected_assignment.assignment
    return '|'.join(
        (
            assignment.assembly_id,
            assignment.gene_id,
            assignment.protein_id,
        )
    )


def write_fasta_record(
    file: TextIO,
    protein: CogProtein,
    selected_assignment: SelectedAssignment,
) -> None:
    """Write one matched protein in FASTA format."""

    genome: GenomeRecord = selected_assignment.genome
    identifier: str = fasta_identifier(selected_assignment)
    file.write(
        ' '.join(
            (
                f'>{identifier}',
                f'{selected_assignment.assignment.cog_id}',
                f'[{genome.species}]\n',
            )
        )
    )
    for sequence_line in wrap(protein.sequence, FASTA_LINE_WIDTH):
        file.write(f'{sequence_line}\n')


def write_metadata_header(metadata_writer: object) -> None:
    """Write the extracted metadata table header."""

    if not hasattr(metadata_writer, 'writerow'):
        raise TypeError('Metadata writer does not provide writerow().')
    metadata_writer.writerow(
        (
            'sequence_id',
            'protein_id',
            'gene_id',
            'assembly_id',
            'organism',
            'taxid',
            'cog_id',
            'protein_length',
            'cog_footprint',
            'cog_footprint_length',
            'membership_class',
            'membership_description',
            'bit_score',
            'e_value',
            'cog_profile_length',
            'profile_footprint',
        )
    )


def write_metadata_row(
    metadata_writer: object,
    selected_assignment: SelectedAssignment,
) -> None:
    """Write metadata for one successfully extracted sequence."""

    if not hasattr(metadata_writer, 'writerow'):
        raise TypeError('Metadata writer does not provide writerow().')
    assignment: CogAssignment = selected_assignment.assignment
    genome: GenomeRecord = selected_assignment.genome
    metadata_writer.writerow(
        (
            fasta_identifier(selected_assignment),
            assignment.protein_id,
            assignment.gene_id,
            assignment.assembly_id,
            genome.species,
            genome.taxid,
            assignment.cog_id,
            assignment.protein_length,
            footprint_to_string(assignment.cog_footprint),
            assignment.cog_footprint_length,
            int(assignment.membership_class),
            assignment.membership_class.description,
            assignment.bit_score,
            assignment.e_value,
            assignment.cog_profile_length,
            footprint_to_string(assignment.profile_footprint),
        )
    )


def candidate_description(candidates: tuple[SelectedAssignment, ...]) -> str:
    """Produce a compact diagnostic representation of candidate records."""

    return ';'.join(
        '|'.join(
            (
                candidate.assignment.assembly_id,
                candidate.assignment.gene_id,
                candidate.assignment.protein_id,
                str(candidate.assignment.protein_length),
            )
        )
        for candidate in candidates
    )


def extract_sequences(
    sequence_file: CogProteinFile,
    selected_assignments: tuple[SelectedAssignment, ...],
    metadata_path: Path,
    fasta_path: Path,
    ambiguous_path: Path,
) -> tuple[
    set[tuple[str, str, str, str]],
    int,
    int,
    int,
]:
    """Stream the FASTA database and write all unambiguously matched proteins."""

    assignments_by_protein_id: dict[
        str,
        tuple[SelectedAssignment, ...],
    ] = index_assignments_by_protein_id(selected_assignments)
    matched_assignment_ids: set[tuple[str, str, str, str]] = set()
    matched_sequences: int = 0
    ambiguous_sequences: int = 0
    duplicate_fasta_matches: int = 0
    with (
        metadata_path.open(
            mode='w',
            encoding='utf-8',
            newline='',
        ) as metadata_file,
        fasta_path.open(
            mode='w',
            encoding='ascii',
            newline='',
        ) as fasta_file,
        ambiguous_path.open(
            mode='w',
            encoding='utf-8',
            newline='',
        ) as ambiguous_file,
    ):
        metadata_writer = writer(metadata_file, delimiter='\t')
        ambiguous_writer = writer(ambiguous_file, delimiter='\t')
        write_metadata_header(metadata_writer)
        ambiguous_writer.writerow(
            (
                'protein_id',
                'description',
                'sequence_length',
                'candidate_count',
                'candidates',
                'reason',
            )
        )
        for protein in sequence_file:
            candidates: tuple[SelectedAssignment, ...] | None = (
                assignments_by_protein_id.get(protein.protein_id)
            )
            if candidates is None:
                continue
            selected_assignment: SelectedAssignment | None = match_protein(
                protein=protein,
                candidates=candidates,
            )
            if selected_assignment is None:
                ambiguous_sequences += 1
                ambiguous_writer.writerow(
                    (
                        protein.protein_id,
                        protein.description or '',
                        protein.length,
                        len(candidates),
                        candidate_description(candidates),
                        'could not resolve by protein ID, length, and gene ID',
                    )
                )
                continue
            assignment_id: tuple[str, str, str, str] = assignment_identity(
                selected_assignment
            )
            if assignment_id in matched_assignment_ids:
                duplicate_fasta_matches += 1
                ambiguous_writer.writerow(
                    (
                        protein.protein_id,
                        protein.description or '',
                        protein.length,
                        len(candidates),
                        candidate_description(candidates),
                        'assignment already matched by another FASTA record',
                    )
                )
                continue
            write_metadata_row(
                metadata_writer=metadata_writer,
                selected_assignment=selected_assignment,
            )
            write_fasta_record(
                file=fasta_file,
                protein=protein,
                selected_assignment=selected_assignment,
            )
            matched_assignment_ids.add(assignment_id)
            matched_sequences += 1
    return (
        matched_assignment_ids,
        matched_sequences,
        ambiguous_sequences,
        duplicate_fasta_matches,
    )


def write_unmatched_assignments(
    path: Path,
    selected_assignments: tuple[SelectedAssignment, ...],
    matched_assignment_ids: set[tuple[str, str, str, str]],
) -> int:
    """Write assignments for which no FASTA sequence was recovered."""

    unmatched_count: int = 0
    with path.open(mode='w', encoding='utf-8', newline=str()) as file:
        unmatched_writer = writer(file, delimiter='\t')
        unmatched_writer.writerow(
            (
                'protein_id',
                'gene_id',
                'assembly_id',
                'organism',
                'taxid',
                'protein_length',
                'membership_class',
            )
        )
        for selected_assignment in selected_assignments:
            assignment_id: tuple[str, str, str, str] = assignment_identity(
                selected_assignment
            )
            if assignment_id in matched_assignment_ids:
                continue
            assignment: CogAssignment = selected_assignment.assignment
            genome: GenomeRecord = selected_assignment.genome
            unmatched_writer.writerow(
                (
                    assignment.protein_id,
                    assignment.gene_id,
                    assignment.assembly_id,
                    genome.species,
                    genome.taxid,
                    assignment.protein_length,
                    int(assignment.membership_class),
                )
            )
            unmatched_count += 1
    return unmatched_count


def run(arguments: Arguments) -> ExtractionStatistics:
    """Execute the complete extraction pipeline."""

    arguments.output_directory.mkdir(parents=True, exist_ok=True)
    genomes: GenomeCollection = GenomeCollection.from_xlsx(arguments.genomes_path)
    representative_genomes: dict[str, GenomeRecord] = (
        representative_genomes_by_assembly(genomes)
    )
    assignment_file: CogAssignmentFile = CogAssignmentFile(
        path=arguments.assignments_path,
        suppress_errors=arguments.suppress_errors,
    )
    selected_assignments: tuple[SelectedAssignment, ...] = select_assignments(
        assignment_file=assignment_file,
        representative_genomes=representative_genomes,
        cog_id=arguments.cog_id,
    )
    selected_genome_accessions: set[str] = {
        selected_assignment.assignment.assembly_id
        for selected_assignment in selected_assignments
    }
    metadata_path: Path = (
        arguments.output_directory / f'{arguments.cog_id}{METADATA_FILENAME_SUFFIX}'
    )
    fasta_path: Path = (
        arguments.output_directory / f'{arguments.cog_id}{FASTA_FILENAME_SUFFIX}'
    )
    ambiguous_path: Path = (
        arguments.output_directory / f'{arguments.cog_id}{AMBIGUOUS_FILENAME_SUFFIX}'
    )
    unmatched_path: Path = (
        arguments.output_directory / f'{arguments.cog_id}{UNMATCHED_FILENAME_SUFFIX}'
    )
    sequence_file: CogProteinFile = CogProteinFile(path=arguments.sequences_path)
    (
        matched_assignment_ids,
        matched_sequences,
        ambiguous_sequences,
        duplicate_fasta_matches,
    ) = extract_sequences(
        sequence_file=sequence_file,
        selected_assignments=selected_assignments,
        metadata_path=metadata_path,
        fasta_path=fasta_path,
        ambiguous_path=ambiguous_path,
    )
    unmatched_assignments: int = write_unmatched_assignments(
        path=unmatched_path,
        selected_assignments=selected_assignments,
        matched_assignment_ids=matched_assignment_ids,
    )
    return ExtractionStatistics(
        representative_genomes=len(representative_genomes),
        selected_assignments=len(selected_assignments),
        selected_genomes=len(selected_genome_accessions),
        matched_sequences=matched_sequences,
        ambiguous_sequences=ambiguous_sequences,
        duplicate_fasta_matches=duplicate_fasta_matches,
        unmatched_assignments=unmatched_assignments,
    )


def main() -> None:
    try:
        arguments: Arguments = parse_arguments()
        statistics: ExtractionStatistics = run(arguments)
        print(f'Representative genomes:       {statistics.representative_genomes}')
        print(f'Selected COG assignments:     {statistics.selected_assignments}')
        print(f'Genomes containing the COG:   {statistics.selected_genomes}')
        print(f'Matched sequences:            {statistics.matched_sequences}')
        print(f'Ambiguous FASTA records:      {statistics.ambiguous_sequences}')
        print(f'Duplicate FASTA matches:      {statistics.duplicate_fasta_matches}')
        print(f'Unmatched assignments:        {statistics.unmatched_assignments}')
        print(f'Output directory:             {arguments.output_directory}')
    except (FileNotFoundError, ValueError) as error:
        print(f'error: {error}', file=stderr)
        raise SystemExit(1) from error
