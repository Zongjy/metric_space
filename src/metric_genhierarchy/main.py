"""Metric Space Data Processing System CLI"""

import argparse
import sys
from typing import List

from metric_genhierarchy.core.metric_space import MetricObject
from metric_genhierarchy.impl.vector import VectorObject, MinkowskiDistance
from metric_genhierarchy.impl.protein import ProteinSequence, MPAMAlignmentDistance
from metric_genhierarchy.query import LinearScan, PivotTable
from metric_genhierarchy.utils import parse_query_file, QueryFileError


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


def query_vectors(args):
    """执行向量数据查询"""
    print(f"\n{'='*60}")
    print("Vector Query Processing")
    print(f"{'='*60}")
    print(f"Dataset: {args.data_file}")
    print(f"Query file: {args.query_file}")
    print(f"Query type: {args.query_type}")
    print(f"Algorithm: {args.algorithm}")
    print(f"Minkowski p: {args.p if not args.p_inf else 'inf'}")

    # 加载数据集
    vectors = VectorObject.from_file(args.data_file, dim=args.dim, count=args.count)
    if not vectors:
        print("Error: No vectors loaded from dataset")
        sys.exit(1)

    print(f"\nLoaded {len(vectors)} vectors from dataset")
    dataset_dim = vectors[0].dim
    print(f"Vector dimension: {dataset_dim}")

    # 加载并验证查询文件
    try:
        queries = parse_query_file(args.query_file, 'vector', expected_dim=dataset_dim)
        print(f"Loaded {len(queries)} queries from file")
    except QueryFileError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # 创建距离函数
    p_value = float('inf') if args.p_inf else args.p
    dist_fn = MinkowskiDistance(p=p_value)

    # 根据算法初始化查询引擎
    if args.algorithm == 'linear_scan':
        print("\nInitializing Linear Scan...")
        query_engine = LinearScan(vectors, dist_fn)
    elif args.algorithm == 'pivot_table':
        print(f"\nInitializing Pivot Table (pivots={args.num_pivots}, selection={args.pivot_selection})...")
        query_engine = PivotTable(vectors, dist_fn, num_pivots=args.num_pivots, pivot_selection=args.pivot_selection)
    else:
        print(f"Error: Unknown algorithm: {args.algorithm}", file=sys.stderr)
        sys.exit(1)

    # 执行查询
    print(f"\n{'='*60}")
    print(f"Executing {args.query_type.upper()} Queries")
    print(f"{'='*60}")

    for idx, query_data in enumerate(queries):
        query_coords = query_data['query_object']
        query_obj = VectorObject(query_coords)

        print(f"\nQuery #{idx+1}:")
        coords_str = ", ".join(f"{x:.8f}" for x in query_coords[:5])
        if len(query_coords) > 5:
            coords_str += ", ..."
        print(f"  Query vector: [{coords_str}]")

        # 根据查询类型执行查询
        if args.query_type == 'range':
            print(f"  Radius: {args.radius}")
            results = query_engine.range_query(query_obj, args.radius)
        elif args.query_type == 'knn':
            print(f"  k: {args.k}")
            results = query_engine.knn_query(query_obj, args.k)
        elif args.query_type == 'dknn':
            print(f"  k: {args.k}, max_distance: {args.max_distance}")
            results = query_engine.dknn_query(query_obj, args.k, args.max_distance)
        else:
            print(f"Error: Unknown query type: {args.query_type}", file=sys.stderr)
            continue

        # 显示结果
        print(f"  Results: {len(results)} objects found")
        print(f"  Distance computations: {query_engine.distance_computations}")

        if results:
            print("  Top results:")
            for i, (obj, dist) in enumerate(results[:5]):
                result_coords = ", ".join(f"{x:.8f}" for x in obj.coords[:3])
                if len(obj.coords) > 3:
                    result_coords += ", ..."
                print(f"    {i+1}. distance={dist:.8f}, vector=[{result_coords}]")
            if len(results) > 5:
                print(f"    ... and {len(results) - 5} more")


def query_proteins(args):
    """执行蛋白质序列查询"""
    print(f"\n{'='*60}")
    print("Protein Sequence Query Processing")
    print(f"{'='*60}")
    print(f"Dataset: {args.data_file}")
    print(f"Query file: {args.query_file}")
    print(f"Query type: {args.query_type}")
    print(f"Algorithm: {args.algorithm}")
    print(f"Gap cost: {args.gap_cost}")

    # 加载数据集
    proteins = ProteinSequence.from_file(args.data_file, length=args.length, count=args.count)
    if not proteins:
        print("Error: No protein sequences loaded from dataset")
        sys.exit(1)

    print(f"\nLoaded {len(proteins)} protein sequences from dataset")

    # 加载并验证查询文件
    try:
        queries = parse_query_file(args.query_file, 'protein')
        print(f"Loaded {len(queries)} queries from file")
    except QueryFileError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # 创建距离函数
    dist_fn = MPAMAlignmentDistance(gap_cost=args.gap_cost)

    # 根据算法初始化查询引擎
    if args.algorithm == 'linear_scan':
        print("\nInitializing Linear Scan...")
        query_engine = LinearScan(proteins, dist_fn)
    elif args.algorithm == 'pivot_table':
        print(f"\nInitializing Pivot Table (pivots={args.num_pivots}, selection={args.pivot_selection})...")
        query_engine = PivotTable(proteins, dist_fn, num_pivots=args.num_pivots, pivot_selection=args.pivot_selection)
    else:
        print(f"Error: Unknown algorithm: {args.algorithm}", file=sys.stderr)
        sys.exit(1)

    # 执行查询
    print(f"\n{'='*60}")
    print(f"Executing {args.query_type.upper()} Queries")
    print(f"{'='*60}")

    for idx, query_data in enumerate(queries):
        query_seq = query_data['query_object']
        query_obj = ProteinSequence(query_seq)

        print(f"\nQuery #{idx+1}:")
        seq_preview = query_seq[:50] + "..." if len(query_seq) > 50 else query_seq
        print(f"  Query sequence: {seq_preview}")
        print(f"  Length: {len(query_seq)}")

        # 根据查询类型执行查询
        if args.query_type == 'range':
            print(f"  Radius: {args.radius}")
            results = query_engine.range_query(query_obj, args.radius)
        elif args.query_type == 'knn':
            print(f"  k: {args.k}")
            results = query_engine.knn_query(query_obj, args.k)
        elif args.query_type == 'dknn':
            print(f"  k: {args.k}, max_distance: {args.max_distance}")
            results = query_engine.dknn_query(query_obj, args.k, args.max_distance)
        else:
            print(f"Error: Unknown query type: {args.query_type}", file=sys.stderr)
            continue

        # 显示结果
        print(f"  Results: {len(results)} sequences found")
        print(f"  Distance computations: {query_engine.distance_computations}")

        if results:
            print("  Top results:")
            for i, (obj, dist) in enumerate(results[:5]):
                seq_id = getattr(obj, "id", None) or f"seq_{i}"
                seq_preview = obj.seq[:30] + "..." if len(obj.seq) > 30 else obj.seq
                print(f"    {i+1}. distance={dist:.8f}, id={seq_id}, seq={seq_preview}")
            if len(results) > 5:
                print(f"    ... and {len(results) - 5} more")


def main():
    parser = argparse.ArgumentParser(
        description="Metric Space Data Processing System",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # 第一层子命令：数据类型
    type_subparsers = parser.add_subparsers(dest="data_type", help="Data type to process")
    type_subparsers.required = True

    # ============================================================
    # Vector 数据类型
    # ============================================================
    parser_vec = type_subparsers.add_parser("vector", help="Process vector data")
    vec_subparsers = parser_vec.add_subparsers(dest="operation", help="Operation to perform")
    vec_subparsers.required = True

    # Vector Distance 子命令
    parser_vec_dist = vec_subparsers.add_parser("distance", help="Compute pairwise distances")
    parser_vec_dist.add_argument("file", help="Path to vector data file (.txt, .vec, .csv)")
    parser_vec_dist.add_argument("--dim", type=int, help="Filter vectors by dimension (optional)")
    parser_vec_dist.add_argument("--count", type=int, help="Maximum number of vectors to load (optional)")
    parser_vec_dist.add_argument("--p", type=float, default=2.0, help="Minkowski p parameter (default: 2.0)")
    parser_vec_dist.add_argument("--p-inf", action="store_true", help="Use p = infinity (L-infinity norm)")
    parser_vec_dist.add_argument("--output-mode", choices=['pairs', 'matrix'], default='pairs',
                                 help="Output mode: 'pairs' for pairwise distances, 'matrix' for full distance matrix (default: pairs)")
    parser_vec_dist.add_argument("--demo-pairs", type=int, default=5, help="Number of pairs to demo in 'pairs' mode (default: 5)")
    parser_vec_dist.set_defaults(func=process_vectors)

    # Vector Query 子命令
    parser_vec_query = vec_subparsers.add_parser("query", help="Query vector dataset")
    parser_vec_query.add_argument("data_file", help="Path to vector dataset file (.txt, .vec, .csv)")
    parser_vec_query.add_argument("query_file", help="Path to query file (JSON format)")
    parser_vec_query.add_argument("query_type", choices=['range', 'knn', 'dknn'],
                                  help="Query type: 'range', 'knn', or 'dknn'")
    parser_vec_query.add_argument("--algorithm", choices=['linear_scan', 'pivot_table'], default='linear_scan',
                                  help="Query algorithm (default: linear_scan)")
    # Dataset loading options
    parser_vec_query.add_argument("--dim", type=int, help="Filter vectors by dimension (optional)")
    parser_vec_query.add_argument("--count", type=int, help="Maximum number of vectors to load (optional)")
    # Distance function options
    parser_vec_query.add_argument("--p", type=float, default=2.0, help="Minkowski p parameter (default: 2.0)")
    parser_vec_query.add_argument("--p-inf", action="store_true", help="Use p = infinity (L-infinity norm)")
    # Query parameters
    parser_vec_query.add_argument("--radius", type=float, help="Radius for range query")
    parser_vec_query.add_argument("--k", type=int, help="Number of nearest neighbors for knn/dknn query")
    parser_vec_query.add_argument("--max-distance", type=float, help="Maximum distance for dknn query")
    # Pivot table options
    parser_vec_query.add_argument("--num-pivots", type=int, default=5, help="Number of pivots for pivot_table (default: 5)")
    parser_vec_query.add_argument("--pivot-selection", choices=['random', 'farthest', 'incremental'], default='random',
                                  help="Pivot selection strategy for pivot_table (default: random)")
    parser_vec_query.set_defaults(func=query_vectors)

    # ============================================================
    # Protein 数据类型
    # ============================================================
    parser_prot = type_subparsers.add_parser("protein", help="Process protein sequence data")
    prot_subparsers = parser_prot.add_subparsers(dest="operation", help="Operation to perform")
    prot_subparsers.required = True

    # Protein Distance 子命令
    parser_prot_dist = prot_subparsers.add_parser("distance", help="Compute pairwise distances")
    parser_prot_dist.add_argument("file", help="Path to protein sequence file (.fasta, .fa, .faa, .txt)")
    parser_prot_dist.add_argument("--length", type=int, help="Filter sequences by length (optional)")
    parser_prot_dist.add_argument("--count", type=int, help="Maximum number of sequences to load (optional)")
    parser_prot_dist.add_argument("--gap-cost", type=int, default=7, help="Linear gap cost for alignment (default: 7)")
    parser_prot_dist.add_argument("--output-mode", choices=['pairs', 'matrix'], default='pairs',
                                  help="Output mode: 'pairs' for pairwise distances, 'matrix' for full distance matrix (default: pairs)")
    parser_prot_dist.add_argument("--demo-pairs", type=int, default=3, help="Number of pairs to demo in 'pairs' mode (default: 3)")
    parser_prot_dist.set_defaults(func=process_proteins)

    # Protein Query 子命令
    parser_prot_query = prot_subparsers.add_parser("query", help="Query protein sequence dataset")
    parser_prot_query.add_argument("data_file", help="Path to protein dataset file (.fasta, .fa, .faa, .txt)")
    parser_prot_query.add_argument("query_file", help="Path to query file (JSON format)")
    parser_prot_query.add_argument("query_type", choices=['range', 'knn', 'dknn'],
                                   help="Query type: 'range', 'knn', or 'dknn'")
    parser_prot_query.add_argument("--algorithm", choices=['linear_scan', 'pivot_table'], default='linear_scan',
                                   help="Query algorithm (default: linear_scan)")
    # Dataset loading options
    parser_prot_query.add_argument("--length", type=int, help="Filter sequences by length (optional)")
    parser_prot_query.add_argument("--count", type=int, help="Maximum number of sequences to load (optional)")
    # Distance function options
    parser_prot_query.add_argument("--gap-cost", type=int, default=7, help="Linear gap cost for alignment (default: 7)")
    # Query parameters
    parser_prot_query.add_argument("--radius", type=float, help="Radius for range query")
    parser_prot_query.add_argument("--k", type=int, help="Number of nearest neighbors for knn/dknn query")
    parser_prot_query.add_argument("--max-distance", type=float, help="Maximum distance for dknn query")
    # Pivot table options
    parser_prot_query.add_argument("--num-pivots", type=int, default=5, help="Number of pivots for pivot_table (default: 5)")
    parser_prot_query.add_argument("--pivot-selection", choices=['random', 'farthest', 'incremental'], default='random',
                                   help="Pivot selection strategy for pivot_table (default: random)")
    parser_prot_query.set_defaults(func=query_proteins)

    args = parser.parse_args()

    # Validate query-specific parameters
    if args.operation == 'query':
        if args.query_type == 'range' and args.radius is None:
            parser.error("--radius is required for range queries")
        if args.query_type == 'knn' and args.k is None:
            parser.error("--k is required for kNN queries")
        if args.query_type == 'dknn':
            if args.k is None:
                parser.error("--k is required for dkNN queries")
            if args.max_distance is None:
                parser.error("--max-distance is required for dkNN queries")

    try:
        args.func(args)
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
