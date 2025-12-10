#!/usr/bin/env python3
"""
Pivot Table Performance Benchmark Script

Automatically tests different pivot configurations and displays results in ASCII tables.
"""

import subprocess
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def detect_data_type(data_file: str) -> str:
    """
    Detect data type from file extension.

    Returns:
        'vector' or 'protein'
    """
    ext = Path(data_file).suffix.lower()

    # Protein file extensions
    if ext in ['.fasta', '.fa', '.faa', '.aa']:
        return 'protein'

    # Vector file extensions (default)
    return 'vector'


def run_query_command(data_file: str, query_file: str, query_type: str,
                      algorithm: str, count: int, data_type: str = None, **kwargs) -> Dict[str, any]:
    """
    Run a query command and parse the results.

    Returns:
        Dict containing query results including distance computations
    """
    # Auto-detect data type if not provided
    if data_type is None:
        data_type = detect_data_type(data_file)

    cmd = [
        sys.executable, "-m", "metric_genhierarchy.main",
        data_type, "query", data_file, query_file, query_type,
        "--algorithm", algorithm,
        "--count", str(count)
    ]

    # Add query type specific parameters
    if query_type == "knn" and "k" in kwargs:
        cmd.extend(["--k", str(kwargs["k"])])
    elif query_type == "range" and "radius" in kwargs:
        cmd.extend(["--radius", str(kwargs["radius"])])
    elif query_type == "dknn":
        if "k" in kwargs:
            cmd.extend(["--k", str(kwargs["k"])])
        if "max_distance" in kwargs:
            cmd.extend(["--max-distance", str(kwargs["max_distance"])])

    # Add pivot table specific parameters
    if algorithm == "pivot_table":
        if "num_pivots" in kwargs:
            cmd.extend(["--num-pivots", str(kwargs["num_pivots"])])
        if "pivot_selection" in kwargs:
            cmd.extend(["--pivot-selection", kwargs["pivot_selection"]])

    # Add data type specific parameters
    if data_type == "vector":
        for param in ["p", "p_inf", "dim"]:
            if param in kwargs:
                if param == "p_inf" and kwargs[param]:
                    cmd.append("--p-inf")
                elif param != "p_inf":
                    cmd.extend([f"--{param}", str(kwargs[param])])
    elif data_type == "protein":
        for param in ["gap_cost", "length"]:
            if param in kwargs:
                cmd.extend([f"--{param.replace('_', '-')}", str(kwargs[param])])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        output = result.stdout

        # Parse distance computations from each query
        distance_comps = re.findall(r"Distance computations:\s+(\d+)", output)
        distance_comps = [int(dc) for dc in distance_comps]

        # Parse number of results
        results_found = re.findall(r"Results:\s+(\d+)\s+(?:objects|sequences)\s+found", output)
        results_found = [int(rf) for rf in results_found]

        return {
            "success": True,
            "distance_computations": distance_comps,
            "avg_distance_computations": sum(distance_comps) / len(distance_comps) if distance_comps else 0,
            "results_found": results_found,
            "output": output
        }
    except subprocess.TimeoutExpired:
        return {"success": False, "error": "Timeout"}
    except Exception as e:
        return {"success": False, "error": str(e)}


def print_table_header(title: str, width: int = 80):
    """Print a formatted table header."""
    print()
    print("=" * width)
    print(f"{title:^{width}}")
    print("=" * width)


def print_table(headers: List[str], rows: List[List[str]], col_widths: List[int] = None):
    """Print a formatted ASCII table."""
    if not col_widths:
        # Auto-calculate column widths
        col_widths = [len(h) for h in headers]
        for row in rows:
            for i, cell in enumerate(row):
                col_widths[i] = max(col_widths[i], len(str(cell)))

    # Add padding
    col_widths = [w + 2 for w in col_widths]

    # Print header
    header_line = "│".join(f"{h:^{col_widths[i]}}" for i, h in enumerate(headers))
    separator = "─" * len(header_line.replace("│", "─"))

    print("┌" + separator + "┐")
    print("│" + header_line + "│")
    print("├" + separator + "┤")

    # Print rows
    for row in rows:
        row_line = "│".join(f"{str(cell):^{col_widths[i]}}" for i, cell in enumerate(row))
        print("│" + row_line + "│")

    print("└" + separator + "┘")


def benchmark_pivot_counts(data_file: str, query_file: str, count: int = 5000, data_type: str = "vector"):
    """Benchmark different numbers of pivots."""
    print_table_header("Benchmark: Number of Pivots (Random Selection)", 80)

    pivot_counts = [3, 5, 10, 15, 20]
    results = []

    # Get baseline (linear scan)
    print("Running Linear Scan baseline...", flush=True)
    baseline = run_query_command(data_file, query_file, "knn", "linear_scan", count, data_type=data_type, k=10)
    baseline_comps = baseline["avg_distance_computations"] if baseline["success"] else count

    print(f"Baseline: {baseline_comps:.0f} distance computations\n")

    # Test different pivot counts
    for num_pivots in pivot_counts:
        print(f"Testing with {num_pivots} pivots...", flush=True)
        result = run_query_command(
            data_file, query_file, "knn", "pivot_table", count,
            data_type=data_type, k=10, num_pivots=num_pivots, pivot_selection="random"
        )

        if result["success"]:
            avg_comps = result["avg_distance_computations"]
            pruning_rate = (1 - avg_comps / baseline_comps) * 100 if baseline_comps > 0 else 0
            results.append([
                str(num_pivots),
                f"{avg_comps:.0f}",
                f"{pruning_rate:.2f}%",
                f"{avg_comps/baseline_comps:.4f}x"
            ])
        else:
            results.append([str(num_pivots), "ERROR", "N/A", "N/A"])

    print()
    headers = ["# Pivots", "Avg Dist Comps", "Pruning Rate", "Speedup"]
    print_table(headers, results)
    print(f"\nBaseline (Linear Scan): {baseline_comps:.0f} distance computations")


def benchmark_pivot_strategies(data_file: str, query_file: str, count: int = 5000, data_type: str = "vector"):
    """Benchmark different pivot selection strategies."""
    print_table_header("Benchmark: Pivot Selection Strategies (10 Pivots)", 80)

    strategies = ["random", "farthest", "incremental"]
    results = []

    # Get baseline
    print("Running Linear Scan baseline...", flush=True)
    baseline = run_query_command(data_file, query_file, "knn", "linear_scan", count, data_type=data_type, k=10)
    baseline_comps = baseline["avg_distance_computations"] if baseline["success"] else count

    print(f"Baseline: {baseline_comps:.0f} distance computations\n")

    # Test different strategies
    for strategy in strategies:
        print(f"Testing {strategy} strategy...", flush=True)
        result = run_query_command(
            data_file, query_file, "knn", "pivot_table", count,
            data_type=data_type, k=10, num_pivots=10, pivot_selection=strategy
        )

        if result["success"]:
            avg_comps = result["avg_distance_computations"]
            pruning_rate = (1 - avg_comps / baseline_comps) * 100 if baseline_comps > 0 else 0
            results.append([
                strategy.capitalize(),
                f"{avg_comps:.0f}",
                f"{pruning_rate:.2f}%",
                f"{avg_comps/baseline_comps:.4f}x"
            ])
        else:
            results.append([strategy.capitalize(), "ERROR", "N/A", "N/A"])

    print()
    headers = ["Strategy", "Avg Dist Comps", "Pruning Rate", "Speedup"]
    print_table(headers, results)
    print(f"\nBaseline (Linear Scan): {baseline_comps:.0f} distance computations")


def benchmark_query_types(data_file: str, query_file: str, count: int = 5000, data_type: str = "vector"):
    """Benchmark different query types."""
    print_table_header("Benchmark: Query Types (10 Pivots, Farthest)", 80)

    query_configs = [
        ("Range", "range", {"radius": 0.1}),
        ("kNN", "knn", {"k": 10}),
        ("dkNN", "dknn", {"k": 10, "max_distance": 1.0})
    ]

    results = []

    for name, qtype, params in query_configs:
        print(f"Testing {name} query...", flush=True)

        # Linear scan baseline
        baseline = run_query_command(data_file, query_file, qtype, "linear_scan", count, data_type=data_type, **params)
        baseline_comps = baseline["avg_distance_computations"] if baseline["success"] else count

        # Pivot table
        result = run_query_command(
            data_file, query_file, qtype, "pivot_table", count,
            data_type=data_type, num_pivots=10, pivot_selection="farthest", **params
        )

        if result["success"] and baseline["success"]:
            avg_comps = result["avg_distance_computations"]
            pruning_rate = (1 - avg_comps / baseline_comps) * 100 if baseline_comps > 0 else 0
            results.append([
                name,
                f"{baseline_comps:.0f}",
                f"{avg_comps:.0f}",
                f"{pruning_rate:.2f}%"
            ])
        else:
            results.append([name, "ERROR", "ERROR", "N/A"])

    print()
    headers = ["Query Type", "Linear Scan", "Pivot Table", "Pruning Rate"]
    print_table(headers, results)


def benchmark_scalability(data_file: str, query_file: str, data_type: str = "vector"):
    """Benchmark scalability with different dataset sizes."""
    print_table_header("Benchmark: Scalability (Dataset Size)", 80)

    dataset_sizes = [1000, 2000, 5000, 10000]
    results = []

    for size in dataset_sizes:
        print(f"Testing with {size} objects...", flush=True)

        # Linear scan
        baseline = run_query_command(data_file, query_file, "knn", "linear_scan", size, data_type=data_type, k=10)
        baseline_comps = baseline["avg_distance_computations"] if baseline["success"] else size

        # Pivot table
        result = run_query_command(
            data_file, query_file, "knn", "pivot_table", size,
            data_type=data_type, k=10, num_pivots=10, pivot_selection="farthest"
        )

        if result["success"] and baseline["success"]:
            avg_comps = result["avg_distance_computations"]
            pruning_rate = (1 - avg_comps / baseline_comps) * 100 if baseline_comps > 0 else 0
            results.append([
                str(size),
                f"{baseline_comps:.0f}",
                f"{avg_comps:.0f}",
                f"{pruning_rate:.2f}%"
            ])
        else:
            results.append([str(size), "ERROR", "ERROR", "N/A"])

    print()
    headers = ["Dataset Size", "Linear Scan", "Pivot Table", "Pruning Rate"]
    print_table(headers, results)


def main():
    """Main benchmark execution."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Pivot Table Performance Benchmark",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("data_file", help="Path to dataset file (vector or protein)")
    parser.add_argument("query_file", help="Path to query file (JSON)")
    parser.add_argument("--count", type=int, default=5000,
                       help="Number of objects to load from dataset (default: 5000)")
    parser.add_argument("--benchmark", choices=["all", "pivots", "strategies", "queries", "scalability"],
                       default="all", help="Which benchmark to run (default: all)")
    parser.add_argument("--data-type", choices=["vector", "protein"],
                       help="Data type (auto-detected from file extension if not specified)")

    args = parser.parse_args()

    # Auto-detect data type if not provided
    data_type = args.data_type if args.data_type else detect_data_type(args.data_file)

    print()
    print("╔" + "═" * 78 + "╗")
    print("║" + "PIVOT TABLE PERFORMANCE BENCHMARK".center(78) + "║")
    print("╚" + "═" * 78 + "╝")
    print()
    print(f"Dataset: {args.data_file}")
    print(f"Data Type: {data_type}")
    print(f"Queries: {args.query_file}")
    print(f"Count: {args.count}")
    print()

    if args.benchmark in ["all", "pivots"]:
        benchmark_pivot_counts(args.data_file, args.query_file, args.count, data_type)

    if args.benchmark in ["all", "strategies"]:
        benchmark_pivot_strategies(args.data_file, args.query_file, args.count, data_type)

    if args.benchmark in ["all", "queries"]:
        benchmark_query_types(args.data_file, args.query_file, args.count, data_type)

    if args.benchmark in ["all", "scalability"]:
        benchmark_scalability(args.data_file, args.query_file, data_type)

    print()
    print("=" * 80)
    print("Benchmark completed!")
    print("=" * 80)
    print()


if __name__ == "__main__":
    main()
