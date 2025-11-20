from abc import ABC, abstractmethod

class MetricObject(ABC):
    """metric space object abstract base class"""

    def __init__(self, raw_data, label=None):
        self.raw_data = raw_data
        self.label = label

    @abstractmethod
    def __repr__(self):
        ...

class DistanceFunction(ABC):
    """metric space distance function abstract base class"""

    @abstractmethod
    def __call__(self, obj1: MetricObject, obj2: MetricObject) -> float:
        ...