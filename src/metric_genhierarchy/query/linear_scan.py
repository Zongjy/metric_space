"""Linear scan query algorithms for metric spaces"""

import heapq
from typing import List, Tuple
from metric_genhierarchy.core.metric_space import MetricObject, DistanceFunction


class LinearScan:
    """
    Linear scan query algorithms for metric spaces.

    Supported queries:
    - Range query: Find all objects within a given distance from a query object
    - kNN query: Find k nearest neighbors
    - dkNN query: Find k nearest neighbors within a distance constraint
    """

    def __init__(self, objects: List[MetricObject], distance_fn: DistanceFunction):
        """
        Initialize linear scan with a dataset and distance function.

        Args:
            objects: List of metric objects to search
            distance_fn: Distance function to measure similarity
        """
        self.objects = objects
        self.distance_fn = distance_fn
        self.distance_computations = 0  # Track performance

    def reset_distance_counter(self):
        """Reset the distance computation counter."""
        self.distance_computations = 0

    def range_query(self, query: MetricObject, radius: float) -> List[Tuple[MetricObject, float]]:
        """
        Range query: Find all objects within radius distance from query object.

        Args:
            query: The query object
            radius: Maximum distance threshold

        Returns:
            List of (object, distance) tuples where distance <= radius
        """
        self.reset_distance_counter()
        results = []

        for obj in self.objects:
            distance = self.distance_fn(query, obj)
            self.distance_computations += 1

            if distance <= radius:
                results.append((obj, distance))

        # Sort by distance for consistent output
        results.sort(key=lambda x: x[1])
        return results

    def knn_query(self, query: MetricObject, k: int) -> List[Tuple[MetricObject, float]]:
        """
        k-Nearest Neighbor query: Find k nearest objects to query.

        Args:
            query: The query object
            k: Number of nearest neighbors to find

        Returns:
            List of (object, distance) tuples for k nearest neighbors,
            sorted by distance in ascending order
        """
        self.reset_distance_counter()

        if k <= 0:
            return []

        if k > len(self.objects):
            k = len(self.objects)

        heap = []
        counter = 0

        for obj in self.objects:
            distance = self.distance_fn(query, obj)
            self.distance_computations += 1

            if len(heap) < k:
                heapq.heappush(heap, (-distance, counter, distance, obj))
                counter += 1
            elif distance < heap[0][2]:
                heapq.heapreplace(heap, (-distance, counter, distance, obj))
                counter += 1

        # Extract results and sort by distance (ascending)
        results = [(obj, dist) for (_, _, dist, obj) in heap]
        results.sort(key=lambda x: x[1])

        return results

    def dknn_query(self, query: MetricObject, k: int, max_distance: float) -> List[Tuple[MetricObject, float]]:
        """
        Distance-constrained k-Nearest Neighbor query:
        Find up to k nearest objects that are within max_distance from query.

        Args:
            query: The query object
            k: Maximum number of nearest neighbors to find
            max_distance: Maximum distance constraint

        Returns:
            List of (object, distance) tuples for up to k nearest neighbors
            within max_distance, sorted by distance in ascending order
        """
        self.reset_distance_counter()

        if k <= 0:
            return []

        # Use max heap to maintain top k nearest neighbors within distance constraint
        # Include counter as tie-breaker to avoid comparing objects
        heap = []
        counter = 0

        for obj in self.objects:
            distance = self.distance_fn(query, obj)
            self.distance_computations += 1

            # Only consider objects within max_distance
            if distance <= max_distance:
                if len(heap) < k:
                    # Heap not full, add directly (negate distance for max heap)
                    heapq.heappush(heap, (-distance, counter, distance, obj))
                    counter += 1
                elif distance < heap[0][2]:  # heap[0][2] is the max distance in heap
                    # Found a closer object, replace the farthest one
                    heapq.heapreplace(heap, (-distance, counter, distance, obj))
                    counter += 1

        # Extract results and sort by distance (ascending)
        results = [(obj, dist) for (_, _, dist, obj) in heap]
        results.sort(key=lambda x: x[1])

        return results
