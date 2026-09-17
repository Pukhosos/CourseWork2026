from .assignment import (
    CogAssignment,
    CogAssignmentFile,
    CoordinateRange,
    Footprint,
    MembershipClass,
)
from .sequence import CogProtein, CogProteinFile
from .table import GenomeCollection, GenomeRecord

__all__: list[str] = [
    'CogAssignment',
    'CogAssignmentFile',
    'CoordinateRange',
    'Footprint',
    'MembershipClass',
    'CogProtein',
    'CogProteinFile',
    'GenomeCollection',
    'GenomeRecord',
]
