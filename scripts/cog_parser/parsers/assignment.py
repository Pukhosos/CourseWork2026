from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from enum import IntEnum
from pathlib import Path
from typing import Final, Self

EXPECTED_COLUMN_COUNT: Final[int] = 13
FOOTPRINT_SEGMENT_SEPARATOR: Final[str] = '='
COORDINATE_SEPARATOR: Final[str] = '-'


class MembershipClass(IntEnum):
    """Quality class of a protein-to-COG assignment."""

    MOST_PROTEIN_MOST_PROFILE = 0
    """Most of protein + most of COG profile."""

    MOST_PROFILE_PART_PROTEIN = 1
    """Most of COG profile + part of protein."""

    MOST_PROTEIN_PART_PROFILE = 2
    """Most of protein + part of COG profile."""

    PARTIAL_PROTEIN_PARTIAL_PROFILE = 3
    """Partial protein + partial COG profile."""

    @property
    def description(self) -> str:
        """Return a verbose description of the membership class."""

        match self:
            case MembershipClass.MOST_PROTEIN_MOST_PROFILE:
                return 'most of protein + most of COG profile'
            case MembershipClass.MOST_PROFILE_PART_PROTEIN:
                return 'most of COG profile + part of protein'
            case MembershipClass.MOST_PROTEIN_PART_PROFILE:
                return 'most of protein + part of COG profile'
            case MembershipClass.PARTIAL_PROTEIN_PARTIAL_PROFILE:
                return 'partial protein + partial COG profile'


@dataclass(frozen=True, slots=True)
class CoordinateRange:
    """Inclusive one-based coordinate range."""

    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 1:
            raise ValueError(f'Coordinate start must be positive, got {self.start}.')
        if self.end < self.start:
            raise ValueError(f'Coordinate end {self.end} precedes start {self.start}.')

    @property
    def length(self) -> int:
        """Return the inclusive length of the coordinate range."""

        return self.end - self.start + 1

    @classmethod
    def parse(cls, value: str) -> Self:
        """Parse a coordinate range *(such as '201-400')*."""

        try:
            start_text, end_text = value.split(COORDINATE_SEPARATOR, maxsplit=1)
        except ValueError as error:
            raise ValueError(f'Invalid coordinate range: {value!r}.') from error
        return cls(start=int(start_text), end=int(end_text))


@dataclass(frozen=True, slots=True)
class Footprint:
    """One or more COG footprint coordinate ranges."""

    ranges: tuple[CoordinateRange, ...]

    def __post_init__(self) -> None:
        if len(self.ranges) == 0:
            raise ValueError('A footprint must contain at least one range.')

    @property
    def length(self) -> int:
        """Return the total number of positions covered by the footprint."""

        return sum(coordinate_range.length for coordinate_range in self.ranges)

    @classmethod
    def parse(cls, value: str) -> Self:
        """Parse a COG footprint *(such as '1-100=201-300')*."""

        ranges: tuple[CoordinateRange, ...] = tuple(
            CoordinateRange.parse(segment)
            for segment in value.split(FOOTPRINT_SEGMENT_SEPARATOR)
        )
        return cls(ranges=ranges)


@dataclass(frozen=True, slots=True)
class CogAssignment:
    """One protein-to-COG assignment from cog-24.cog.csv."""

    gene_id: str
    assembly_id: str
    protein_id: str
    protein_length: int
    cog_footprint: Footprint
    cog_footprint_length: int
    cog_id: str
    membership_class: MembershipClass
    bit_score: float
    e_value: float
    cog_profile_length: int
    profile_footprint: Footprint

    def __post_init__(self) -> None:
        if self.protein_length < 1:
            raise ValueError(
                f'Protein length must be positive, got {self.protein_length}.'
            )
        if self.cog_footprint_length != self.cog_footprint.length:
            error_message: str = ' '.join(
                (
                    'Reported COG footprint length does not match its',
                    f'coordinates for protein {self.protein_id!r}:',
                    f'{self.cog_footprint_length} !=',
                    f'{self.cog_footprint.length}.',
                )
            )
            raise ValueError(error_message)
        if self.cog_profile_length < 1:
            raise ValueError(
                f'COG profile length must be positive, got {self.cog_profile_length}.'
            )

    @classmethod
    def from_csv_line(cls, line: str, separator: str) -> Self:
        """Parse one line of cog-24.cog.csv."""

        fields: list[str] = line.split(separator)
        if len(fields) != EXPECTED_COLUMN_COUNT:
            raise ValueError(
                f'expected {EXPECTED_COLUMN_COUNT} fields, got {len(fields)}.'
            )
        (
            gene_id,
            assembly_id,
            protein_id,
            protein_length_str,
            cog_footprint_str,
            cog_footprint_length_str,
            cog_id,
            reserved_cog_id,
            membership_class_str,
            bit_score_str,
            e_value_str,
            cog_profile_length_str,
            profile_footprint_str,
        ) = fields
        if reserved_cog_id != cog_id:
            raise ValueError(
                f'Reserved COG ID {reserved_cog_id!r} does not match COG ID {cog_id!r}.'
            )
        return cls(
            gene_id=gene_id,
            assembly_id=assembly_id,
            protein_id=protein_id,
            protein_length=int(protein_length_str),
            cog_footprint=Footprint.parse(cog_footprint_str),
            cog_footprint_length=int(cog_footprint_length_str),
            cog_id=cog_id,
            membership_class=MembershipClass(int(membership_class_str)),
            bit_score=float(bit_score_str),
            e_value=float(e_value_str),
            cog_profile_length=int(cog_profile_length_str),
            profile_footprint=Footprint.parse(profile_footprint_str),
        )


@dataclass(frozen=True, slots=True)
class CogAssignmentFile(Iterable[CogAssignment]):
    """Lazily iterable COG assignment database."""

    path: Path
    separator: str = ','
    suppress_errors: bool = False

    def __post_init__(self) -> None:
        if not self.path.is_file():
            raise FileNotFoundError(
                f'COG assignment database does not exist: {self.path}'
            )

    def __iter__(self) -> Iterator[CogAssignment]:
        with self.path.open(mode='r', encoding='utf-8', newline=str()) as file:
            for line_number, line in enumerate(file, start=1):
                stripped_line: str = line.strip()
                if len(stripped_line) == 0:
                    continue
                try:
                    yield CogAssignment.from_csv_line(
                        line=stripped_line,
                        separator=self.separator,
                    )
                except TypeError, ValueError:
                    if not self.suppress_errors:
                        print(f'Cannot parse {self.path} line {(
                                line_number
                            )}: {stripped_line!r}')

    def for_cog(self, cog_id: str) -> Iterator[CogAssignment]:
        """Yield only assignments belonging to the requested COG."""

        return (assignment for assignment in self if assignment.cog_id == cog_id)
