"""Utility functions for file validation and query parsing"""

import json
from typing import List, Dict, Any
from pathlib import Path


class QueryFileError(Exception):
    """Exception raised for errors in query file parsing or validation"""
    pass


def validate_file_exists(file_path: str) -> Path:
    """
    Check if a file exists.

    Args:
        file_path: Path to the file

    Returns:
        Path object if file exists

    Raises:
        QueryFileError: If file does not exist
    """
    path = Path(file_path)
    if not path.exists():
        raise QueryFileError(f"File not found: {file_path}")
    if not path.is_file():
        raise QueryFileError(f"Path is not a file: {file_path}")
    return path


def load_query_file(file_path: str) -> List[Dict[str, Any]]:
    """
    Load and parse a query file (JSON format).

    Expected format:
    [
        {"query_object": [...]},
        {"query_object": [...]},
        ...
    ]

    Args:
        file_path: Path to the query file

    Returns:
        List of query objects

    Raises:
        QueryFileError: If file cannot be parsed or has invalid format
    """
    path = validate_file_exists(file_path)

    try:
        with open(path, 'r', encoding='utf-8') as f:
            queries = json.load(f)
    except json.JSONDecodeError as e:
        raise QueryFileError(f"Invalid JSON format in query file: {e}")
    except Exception as e:
        raise QueryFileError(f"Error reading query file: {e}")

    if not isinstance(queries, list):
        raise QueryFileError("Query file must contain a JSON array")

    if len(queries) == 0:
        raise QueryFileError("Queries list cannot be empty")

    return queries


def validate_vector_query(query: Dict[str, Any], expected_dim: int = None) -> None:
    """
    Validate a single vector query object.

    Args:
        query: Query dictionary containing 'query_object' field
        expected_dim: Expected vector dimension (optional)

    Raises:
        QueryFileError: If query is invalid
    """
    if 'query_object' not in query:
        raise QueryFileError("Each query must contain 'query_object' field")

    query_obj = query['query_object']

    if not isinstance(query_obj, list):
        raise QueryFileError(
            f"Vector query_object must be a list, got {type(query_obj).__name__}"
        )

    if len(query_obj) == 0:
        raise QueryFileError("Vector query_object cannot be empty")

    # Validate all elements are numbers
    for i, val in enumerate(query_obj):
        if not isinstance(val, (int, float)):
            raise QueryFileError(
                f"Vector query_object[{i}] must be a number, got {type(val).__name__}"
            )

    # Validate dimension if specified
    if expected_dim is not None and len(query_obj) != expected_dim:
        raise QueryFileError(
            f"Vector dimension mismatch: expected {expected_dim}, got {len(query_obj)}"
        )


def validate_protein_query(query: Dict[str, Any]) -> None:
    """
    Validate a single protein sequence query object.

    Args:
        query: Query dictionary containing 'query_object' field

    Raises:
        QueryFileError: If query is invalid
    """
    if 'query_object' not in query:
        raise QueryFileError("Each query must contain 'query_object' field")

    query_obj = query['query_object']

    if not isinstance(query_obj, str):
        raise QueryFileError(
            f"Protein query_object must be a string, got {type(query_obj).__name__}"
        )

    if len(query_obj) == 0:
        raise QueryFileError("Protein query_object (sequence) cannot be empty")

    # Validate amino acid characters
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    for i, char in enumerate(query_obj.upper()):
        if char not in valid_aa:
            raise QueryFileError(
                f"Invalid amino acid '{char}' at position {i} in query sequence"
            )


def parse_query_file(file_path: str, dataset_type: str, expected_dim: int = None) -> List[Dict[str, Any]]:
    """
    Parse and validate a query file.

    Args:
        file_path: Path to the query file
        dataset_type: Type of dataset ('vector' or 'protein')
        expected_dim: Expected vector dimension for vector datasets (optional)

    Returns:
        List of validated query objects

    Raises:
        QueryFileError: If file is invalid or queries don't match dataset format
    """
    queries = load_query_file(file_path)

    for i, query in enumerate(queries):
        try:
            if dataset_type == 'vector':
                validate_vector_query(query, expected_dim)
            elif dataset_type == 'protein':
                validate_protein_query(query)
            else:
                raise QueryFileError(f"Unknown dataset type: {dataset_type}")
        except QueryFileError as e:
            raise QueryFileError(f"Error in query #{i+1}: {e}")

    return queries
