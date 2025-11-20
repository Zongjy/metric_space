"""Metric Space Data Processing System CLI"""

import argparse
import sys
from typing import List

from metric_genhierarchy.core.metric_space import MetricObject
from metric_genhierarchy.impl.vector import VectorObject, MinkowskiDistance
from metric_genhierarchy.impl.protein import ProteinSequence, MPAMAlignmentDistance


def compute_distance_matrix(objects: List[MetricObject], distance_fn):
    n = len(objects)
    if n < 2:
        print("Error: Need at least 2 objects to compute distances")
        return

    print(f"\nComputing distance matrix for {n} objects:")

    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                d = 0.0
            else:
                d = distance_fn(objects[i], objects[j])
            row.append(d)
        matrix.append(row)

    col_width = 13
    index_width = max(4, len(str(n - 1)) + 1)

    print("\nDistance Matrix:")

    header = " " * index_width + "│"
    for j in range(n):
        header += f"{j:>{col_width}}"
    print(header)
    print("─" * index_width + "┼" + "─" * (col_width * n))

    # 打印每一行
    for i in range(n):
        row_str = f"{i:>{index_width-1}} │"
        for j in range(n):
            row_str += f"{matrix[i][j]:>{col_width}.4e}"
        print(row_str)


def compute_distances(objects: List[MetricObject], distance_fn, max_pairs: int = 5):
    """计算并显示前 max_pairs 对相邻对象之间的距离"""
    n = len(objects)
    if n < 2:
        print("Error: Need at least 2 objects to compute distances")
        return

    print(f"\nComputing distances (showing up to {max_pairs} pairs):")
    pairs_shown = 0
    for i in range(n - 1):
        if pairs_shown >= max_pairs:
            break

        obj_a = objects[i]
        obj_b = objects[i + 1]

        print(f"\nPair {pairs_shown + 1}: obj[{i}] vs obj[{i+1}]")

        if isinstance(obj_a, VectorObject):
            coords_a = ", ".join(f"{x:.4e}" for x in obj_a.coords)
            coords_b = ", ".join(f"{x:.4e}" for x in obj_b.coords)
            print(f"  obj[{i}]: [{coords_a}]")
            print(f"  obj[{i+1}]: [{coords_b}]")
        elif isinstance(obj_a, ProteinSequence):
            id_a = obj_a.id or f"seq{i}"
            id_b = obj_b.id or f"seq{i+1}"
            print(f"  obj[{i}] ({id_a}):")
            print(f"    {obj_a.seq}")
            print(f"  obj[{i+1}] ({id_b}):")
            print(f"    {obj_b.seq}")

        d = distance_fn(obj_a, obj_b)

        print(f"  --> Distance = {d:.4e}")
        pairs_shown += 1


def process_vectors(args):
    """处理向量数据"""
    print(f"\n{'='*60}")
    print("Processing Vector Data")
    print(f"{'='*60}")
    print(f"File: {args.file}")
    print(f"Dimension filter: {args.dim if args.dim else 'None (all)'}")
    print(f"Count: {args.count if args.count else 'All'}")
    print(f"Minkowski p: {args.p if not args.p_inf else 'inf'}")

    # 加载向量
    vectors = VectorObject.from_file(args.file, dim=args.dim, count=args.count)

    if not vectors:
        print("Error: No vectors loaded")
        sys.exit(1)

    print(f"\nLoaded {len(vectors)} vectors")
    if vectors:
        print(f"First vector dimension: {vectors[0].dim}")

    # 计算距离
    p_value = float('inf') if args.p_inf else args.p
    dist_fn = MinkowskiDistance(p=p_value)

    if args.output_mode == 'matrix':
        compute_distance_matrix(vectors, dist_fn)
    else:  # pairs
        compute_distances(vectors, dist_fn, max_pairs=args.demo_pairs)


def process_proteins(args):
    """处理蛋白质序列数据"""
    print(f"\n{'='*60}")
    print("Processing Protein Sequence Data")
    print(f"{'='*60}")
    print(f"File: {args.file}")
    print(f"Length filter: {args.length if args.length else 'None (all)'}")
    print(f"Count: {args.count if args.count else 'All'}")
    print(f"Gap cost: {args.gap_cost}")

    # 加载蛋白质序列
    proteins = ProteinSequence.from_file(args.file, length=args.length, count=args.count)

    if not proteins:
        print("Error: No protein sequences loaded")
        sys.exit(1)

    print(f"\nLoaded {len(proteins)} protein sequences")
    if proteins:
        seq_id = getattr(proteins[0], "id", None) or "(no-id)"
        seq_len = len(proteins[0].seq)
        print(f"First sequence: id={seq_id}, length={seq_len}")

    # 计算距离
    dist_fn = MPAMAlignmentDistance(gap_cost=args.gap_cost)

    if args.output_mode == 'matrix':
        compute_distance_matrix(proteins, dist_fn)
    else:  # pairs
        compute_distances(proteins, dist_fn, max_pairs=args.demo_pairs)


def main():
    parser = argparse.ArgumentParser(
        description="Metric Space Data Processing System - Process vector or protein sequence data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest="command", help="Data type to process")
    subparsers.required = True

    # Vector 子命令
    parser_vec = subparsers.add_parser("vector", help="Process vector data")
    parser_vec.add_argument("file", help="Path to vector data file (.txt, .vec, .csv)")
    parser_vec.add_argument("--dim", type=int, help="Filter vectors by dimension (optional)")
    parser_vec.add_argument("--count", type=int, help="Maximum number of vectors to load (optional)")
    parser_vec.add_argument("--p", type=float, default=2.0, help="Minkowski p parameter (default: 2.0)")
    parser_vec.add_argument("--p-inf", action="store_true", help="Use p = infinity (L-infinity norm)")
    parser_vec.add_argument("--output-mode", choices=['pairs', 'matrix'], default='pairs',
                            help="Output mode: 'pairs' for pairwise distances, 'matrix' for full distance matrix (default: pairs)")
    parser_vec.add_argument("--demo-pairs", type=int, default=5, help="Number of pairs to demo in 'pairs' mode (default: 5)")
    parser_vec.set_defaults(func=process_vectors)

    # Protein 子命令
    parser_prot = subparsers.add_parser("protein", help="Process protein sequence data")
    parser_prot.add_argument("file", help="Path to protein sequence file (.fasta, .fa, .faa, .txt)")
    parser_prot.add_argument("--length", type=int, help="Filter sequences by length (optional)")
    parser_prot.add_argument("--count", type=int, help="Maximum number of sequences to load (optional)")
    parser_prot.add_argument("--gap-cost", type=int, default=7, help="Linear gap cost for alignment (default: 7)")
    parser_prot.add_argument("--output-mode", choices=['pairs', 'matrix'], default='pairs',
                             help="Output mode: 'pairs' for pairwise distances, 'matrix' for full distance matrix (default: pairs)")
    parser_prot.add_argument("--demo-pairs", type=int, default=3, help="Number of pairs to demo in 'pairs' mode (default: 3)")
    parser_prot.set_defaults(func=process_proteins)

    args = parser.parse_args()

    try:
        args.func(args)
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
