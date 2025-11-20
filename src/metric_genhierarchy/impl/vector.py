import math
import re
from typing import Iterable, List

from metric_genhierarchy.core.metric_space import DistanceFunction, MetricObject


class VectorObject(MetricObject):
    def __init__(self, coords: Iterable[float]):
        self.coords = tuple(float(x) for x in coords)

    @property
    def dim(self) -> int:
        return len(self.coords)

    def __repr__(self) -> str:
        if self.dim <= 5:
            coords_str = ", ".join(f"{x:.4e}" for x in self.coords)
            return f"VectorObject(dim={self.dim}, [{coords_str}])"
        else:
            first_three = ", ".join(f"{x:.4e}" for x in self.coords[:3])
            return f"VectorObject(dim={self.dim}, [{first_three}, ...])"

    @staticmethod
    def _iter_vector_lines(path: str) -> Iterable[List[float]]:
        """
        兼容常见文本格式：
          - 每行空格/逗号分隔的一条向量
          - 可自动跳过含非数值的行
        """
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line: 
                    continue
                parts = re.split(r"[,\s]+", line)
                try:
                    vec = [float(x) for x in parts]
                except ValueError:
                    continue
                yield vec

    @classmethod
    def from_file(cls, file_path: str, dim: int = None, count: int = None) -> List["VectorObject"]:
        """
        从指定文件读取向量数据
        Args:
          file_path: 向量数据文件路径
          dim: 指定维度（可选，None表示接受所有维度）
          count: 指定数量（可选，None表示读取所有数据）
        """
        picked: List[VectorObject] = []
        for vec in cls._iter_vector_lines(file_path):
            if dim is None or len(vec) == dim:
                picked.append(cls(vec))
                if count is not None and len(picked) >= count:
                    return picked
        return picked

class MinkowskiDistance(DistanceFunction):
    def __init__(self, p: float = 2.0):
        assert p > 0 or math.isinf(p), "p must be >0 or inf"
        self.p = p

    def __call__(self, a: VectorObject, b: VectorObject) -> float:
        assert isinstance(a, VectorObject) and isinstance(b, VectorObject), "type mismatch"
        assert a.dim == b.dim, "dimension mismatch"
        if math.isinf(self.p):
            return max(abs(x - y) for x, y in zip(a.coords, b.coords))
        if self.p == 1:
            return sum(abs(x - y) for x, y in zip(a.coords, b.coords))
        if self.p == 2:
            s = sum((x - y) * (x - y) for x, y in zip(a.coords, b.coords))
            return math.sqrt(s)
        s = sum(abs(x - y) ** self.p for x, y in zip(a.coords, b.coords))
        return s ** (1.0 / self.p)