"""Pivot Table index structure for metric spaces"""

import heapq
import random
from typing import List, Tuple, Set
from metric_genhierarchy.core.metric_space import MetricObject, DistanceFunction


class PivotTable:
    """
    Pivot Table index structure for metric space queries.

    Uses pre-selected pivot points and the triangle inequality to prune
    the search space and reduce distance computations.

    For a query q, object o, and pivot p:
    |d(q,p) - d(o,p)| <= d(q,o) <= d(q,p) + d(o,p)

    We can filter out objects that violate these bounds.
    """

    def __init__(self, objects: List[MetricObject], distance_fn: DistanceFunction, num_pivots: int = 5, pivot_selection: str = 'random'):
        """
        Initialize Pivot Table index.

        Args:
            objects: List of metric objects to index
            distance_fn: Distance function to measure similarity
            num_pivots: Number of pivots to use
            pivot_selection: Pivot selection strategy ('random', 'farthest', 'incremental')
        """
        self.objects = objects
        self.distance_fn = distance_fn
        self.num_pivots = min(num_pivots, len(objects))
        self.distance_computations = 0

        # Select pivots and build the pivot table
        self.pivots = self._select_pivots(pivot_selection)
        self.pivot_table = self._build_pivot_table()

    def reset_distance_counter(self):
        """Reset the distance computation counter."""
        self.distance_computations = 0

    def _select_pivots(self, strategy: str) -> List[MetricObject]:
        """
        Select pivot points using different strategies.

        Args:
            strategy: Pivot selection strategy
                - 'random': Randomly select pivots
                - 'farthest': Select pivots that are far from each other
                - 'incremental': Incrementally select pivots maximizing minimum distance

        Returns:
            List of selected pivot objects
        """
        if self.num_pivots <= 0 or not self.objects:
            return []

        if strategy == 'random':
            return random.sample(self.objects, self.num_pivots)

        elif strategy == 'farthest':
            # Select first pivot randomly
            pivots = [random.choice(self.objects)]

            # Iteratively select the object farthest from all existing pivots
            for _ in range(self.num_pivots - 1):
                max_min_dist = -1
                farthest_obj = None

                for obj in self.objects:
                    if obj in pivots:
                        continue

                    # Find minimum distance to existing pivots
                    min_dist = min(self.distance_fn(obj, p) for p in pivots)

                    if min_dist > max_min_dist:
                        max_min_dist = min_dist
                        farthest_obj = obj

                if farthest_obj:
                    pivots.append(farthest_obj)

            return pivots

        elif strategy == 'incremental':
            # Similar to farthest, but explicitly maximize min distance to existing pivots
            pivots = [random.choice(self.objects)]

            for _ in range(self.num_pivots - 1):
                max_min_dist = -1
                next_pivot = None

                for obj in self.objects:
                    if obj in pivots:
                        continue

                    # Calculate minimum distance to existing pivots
                    min_dist = float('inf')
                    for p in pivots:
                        d = self.distance_fn(obj, p)
                        min_dist = min(min_dist, d)

                    # Select object with maximum min_dist
                    if min_dist > max_min_dist:
                        max_min_dist = min_dist
                        next_pivot = obj

                if next_pivot:
                    pivots.append(next_pivot)

            return pivots

        else:
            raise ValueError(f"Unknown pivot selection strategy: {strategy}")

    def _build_pivot_table(self) -> List[List[float]]:
        """
        Build pivot table: precompute distances from each object to each pivot.

        Returns:
            2D list where pivot_table[i][j] = distance from objects[i] to pivots[j]
        """
        table = []
        for obj in self.objects:
            distances = []
            for pivot in self.pivots:
                dist = self.distance_fn(obj, pivot)
                distances.append(dist)
            table.append(distances)

        return table

    def _can_prune(self, query_pivot_dists: List[float], obj_idx: int, threshold: float) -> bool:
        """
        Check if an object can be pruned using triangle inequality.

        For each pivot p: |d(q,p) - d(o,p)| <= d(q,o)
        If |d(q,p) - d(o,p)| > threshold, then d(q,o) > threshold, so prune.

        Args:
            query_pivot_dists: Distances from query to each pivot
            obj_idx: Index of object to check
            threshold: Distance threshold for pruning

        Returns:
            True if object can be pruned (definitely beyond threshold)
        """
        for pivot_idx, query_dist in enumerate(query_pivot_dists):
            obj_dist = self.pivot_table[obj_idx][pivot_idx]
            lower_bound = abs(query_dist - obj_dist)

            if lower_bound > threshold:
                return True

        return False

    def range_query(self, query: MetricObject, radius: float) -> List[Tuple[MetricObject, float]]:
        """
        Range query with pivot-based filtering.

        Args:
            query: The query object
            radius: Maximum distance threshold

        Returns:
            List of (object, distance) tuples where distance <= radius
        """
        self.reset_distance_counter()

        # Compute distances from query to all pivots
        query_pivot_dists = []
        for pivot in self.pivots:
            dist = self.distance_fn(query, pivot)
            self.distance_computations += 1
            query_pivot_dists.append(dist)

        results = []

        # Check each object
        for idx, obj in enumerate(self.objects):
            # Try to prune using pivot distances
            if self._can_prune(query_pivot_dists, idx, radius):
                continue

            # Cannot prune, compute actual distance
            distance = self.distance_fn(query, obj)
            self.distance_computations += 1

            if distance <= radius:
                results.append((obj, distance))

        # Sort by distance
        results.sort(key=lambda x: x[1])
        return results

    def knn_query(self, query: MetricObject, k: int) -> List[Tuple[MetricObject, float]]:
        """
        k-Nearest Neighbor query with pivot-based filtering.

        Args:
            query: The query object
            k: Number of nearest neighbors to find

        Returns:
            List of (object, distance) tuples for k nearest neighbors
        """
        self.reset_distance_counter()

        if k <= 0:
            return []

        if k > len(self.objects):
            k = len(self.objects)

        # Compute distances from query to all pivots
        query_pivot_dists = []
        for pivot in self.pivots:
            dist = self.distance_fn(query, pivot)
            self.distance_computations += 1
            query_pivot_dists.append(dist)

        # Use max heap to maintain top k nearest neighbors
        # Include counter as tie-breaker to avoid comparing objects
        heap = []
        counter = 0
        current_kth_dist = float('inf')  # Distance of k-th nearest neighbor

        # Process each object
        for idx, obj in enumerate(self.objects):
            # Try to prune using current k-th distance as threshold
            if len(heap) >= k and self._can_prune(query_pivot_dists, idx, current_kth_dist):
                continue

            # Cannot prune, compute actual distance
            distance = self.distance_fn(query, obj)
            self.distance_computations += 1

            if len(heap) < k:
                heapq.heappush(heap, (-distance, counter, distance, obj))
                counter += 1
                current_kth_dist = heap[0][2]
            elif distance < current_kth_dist:
                heapq.heapreplace(heap, (-distance, counter, distance, obj))
                counter += 1
                current_kth_dist = heap[0][2]

        # Extract and sort results
        results = [(obj, dist) for (_, _, dist, obj) in heap]
        results.sort(key=lambda x: x[1])

        return results

    def dknn_query(self, query: MetricObject, k: int, max_distance: float) -> List[Tuple[MetricObject, float]]:
        """
        Distance-constrained k-Nearest Neighbor query with pivot-based filtering.

        Args:
            query: The query object
            k: Maximum number of nearest neighbors to find
            max_distance: Maximum distance constraint

        Returns:
            List of (object, distance) tuples for up to k nearest neighbors within max_distance
        """
        self.reset_distance_counter()

        if k <= 0:
            return []

        # Compute distances from query to all pivots
        query_pivot_dists = []
        for pivot in self.pivots:
            dist = self.distance_fn(query, pivot)
            self.distance_computations += 1
            query_pivot_dists.append(dist)

        # Use max heap to maintain top k nearest neighbors
        # Include counter as tie-breaker to avoid comparing objects
        heap = []
        counter = 0

        # Process each object
        for idx, obj in enumerate(self.objects):
            # Prune using max_distance constraint
            if self._can_prune(query_pivot_dists, idx, max_distance):
                continue

            # Cannot prune, compute actual distance
            distance = self.distance_fn(query, obj)
            self.distance_computations += 1

            # Only consider objects within max_distance
            if distance <= max_distance:
                if len(heap) < k:
                    heapq.heappush(heap, (-distance, counter, distance, obj))
                    counter += 1
                elif distance < heap[0][2]:
                    heapq.heapreplace(heap, (-distance, counter, distance, obj))
                    counter += 1

        # Extract and sort results
        results = [(obj, dist) for (_, _, dist, obj) in heap]
        results.sort(key=lambda x: x[1])

        return results
