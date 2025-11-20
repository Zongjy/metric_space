"""Metric Space Data Processing System

A system for processing metric space data with support for:
- Vector data with Minkowski distance
- Protein sequences with mPAM alignment distance
- Extensible architecture for adding new data types
"""

__version__ = "0.1.0"

from metric_genhierarchy.core.metric_space import MetricObject, DistanceFunction
from metric_genhierarchy.impl.vector import VectorObject, MinkowskiDistance
from metric_genhierarchy.impl.protein import ProteinSequence, MPAMAlignmentDistance

__all__ = [
    "MetricObject",
    "DistanceFunction",
    "VectorObject",
    "MinkowskiDistance",
    "ProteinSequence",
    "MPAMAlignmentDistance",
]
