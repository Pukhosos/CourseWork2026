from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from re import compile, Match, Pattern  # noqa: A004
from typing import Final, Self

FASTA_HEADER_PREFIX: Final[str] = '>'
HEADER_PATTERN: Final[Pattern[str]] = compile(
    r'^(?P<protein_id>\S+)'
    r'(?:\s+(?P<description>.*?))?'
    r'(?:\s+\[(?P<organism>[^\[\]]+)\])?$'
)


def parse_header(header: str) -> tuple[str, str | None, str | None]:
    """Parse a COGorg24.faa FASTA header."""

    m: Match[str] | None = HEADER_PATTERN.fullmatch(header)
    if m is None:
        raise ValueError(f'Invalid FASTA header: {header!r}.')
    protein_id: str = m.group('protein_id')
    description: str | None = m.group('description')
    organism: str | None = m.group('organism')
    if description is not None:
        description = description.strip() or None
    return protein_id, description, organism


@dataclass(frozen=True, slots=True)
class CogProtein:
    """One protein record from COGorg24.faa."""

    protein_id: str
    description: str | None
    organism: str | None
    sequence: str

    @property
    def length(self) -> int:
        """Return the protein sequence length."""

        return len(self.sequence)

    @classmethod
    def from_fasta_record(cls, header: str, sequence_parts: list[str]) -> Self:
        """Construct a protein from one complete FASTA record."""

        protein_id, description, organism = parse_header(header)
        sequence: str = str().join(sequence_parts)
        if len(sequence) == 0:
            raise ValueError(f'Protein {protein_id!r} has an empty sequence.')
        return cls(
            protein_id=protein_id,
            description=description,
            organism=organism,
            sequence=sequence,
        )


@dataclass(frozen=True, slots=True)
class CogProteinFile(Iterable[CogProtein]):
    """Reiterable streaming view over COGorg24.faa."""

    path: Path

    def __post_init__(self) -> None:
        if not self.path.is_file():
            raise FileNotFoundError(f'COG protein database does not exist: {self.path}')

    def __iter__(self) -> Iterator[CogProtein]:
        """Yield FASTA records one at a time."""

        with self.path.open(mode='r', encoding='ascii', newline=str()) as file:
            header: str | None = None
            sequence_parts: list[str] = []
            for line_no, line in enumerate(file, start=1):
                stripped_line: str = line.strip()
                if len(stripped_line) == 0:
                    continue
                if stripped_line.startswith(FASTA_HEADER_PREFIX):
                    if header is not None:
                        yield CogProtein.from_fasta_record(
                            header=header,
                            sequence_parts=sequence_parts,
                        )
                    header = stripped_line[1:]
                    sequence_parts = []
                    continue
                if header is None:
                    raise ValueError(f'Header before the first FASTA, line {line_no}')
                sequence_parts.append(stripped_line)
            if header is not None:
                yield CogProtein.from_fasta_record(
                    header=header,
                    sequence_parts=sequence_parts,
                )

    def with_protein_ids(
        self,
        protein_ids: set[str] | frozenset[str],
    ) -> Iterator[CogProtein]:
        """Yield proteins whose IDs belong to the requested set."""

        return (protein for protein in self if protein.protein_id in protein_ids)
