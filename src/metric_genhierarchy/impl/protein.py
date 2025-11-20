import os

from typing import Iterable, Optional, Tuple, List
from metric_genhierarchy.core.metric_space import DistanceFunction, MetricObject

    
_AA_ORDER = [
    "A","R","N","D","C","Q","E","G","H","I",
    "L","K","M","F","P","S","T","W","Y","V","X"  # X means other
]
_AA_SET = set(_AA_ORDER)

class ProteinSequence(MetricObject):
    def __init__(self, seq: str, id: Optional[str] = None):
        seq = seq.strip().upper()
        seq = seq.replace("B", "X").replace("Z","X").replace("U","X")
        seq = "".join([c if c in _AA_SET else "X" for c in seq])
        assert len(seq) > 0, "empty sequence"
        self.seq = seq
        self.id = id

    def __len__(self): return len(self.seq)

    def __repr__(self) -> str:
        seq_preview = self.seq[:20] + "..." if len(self.seq) > 20 else self.seq
        id_str = f"id={self.id}, " if self.id else ""
        return f"ProteinSequence({id_str}length={len(self.seq)}, seq={seq_preview})"

    @staticmethod
    def _iter_fasta(path: str) -> Iterable[Tuple[str,str]]:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            sid, buf = None, []
            for line in f:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if sid is not None and buf:
                        yield sid, "".join(buf)
                    sid, buf = line[1:].strip(), []
                else:
                    if line.strip():
                        buf.append(line.strip())
            if sid is not None and buf:
                yield sid, "".join(buf)

    @staticmethod
    def _iter_plain(path: str) -> Iterable[Tuple[str,str]]:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for i, line in enumerate(f, 1):
                s = line.strip()
                if s and not s.startswith(">"):
                    yield f"{os.path.basename(path)}:{i}", s

    @classmethod
    def from_file(cls, file_path: str, length: int = None, count: int = None) -> List["ProteinSequence"]:
        """
        从指定文件读取蛋白质序列数据
        Args:
          file_path: 蛋白质序列文件路径（支持FASTA或纯文本格式）
          length: 指定长度（可选，None表示接受所有长度）
          count: 指定数量（可选，None表示读取所有数据）
        """
        ext = os.path.splitext(file_path)[1].lower()
        exts_fa = {".fasta", ".fa", ".faa"}
        exts_txt = {".txt"}

        picked: List[ProteinSequence] = []
        try:
            if ext in exts_fa:
                it = cls._iter_fasta(file_path)
            elif ext in exts_txt:
                with open(file_path, "r", encoding="utf-8", errors="ignore") as f0:
                    head = f0.readline()
                it = cls._iter_fasta(file_path) if head.startswith(">") else cls._iter_plain(file_path)
            else:
                # 尝试自动检测格式
                with open(file_path, "r", encoding="utf-8", errors="ignore") as f0:
                    head = f0.readline()
                it = cls._iter_fasta(file_path) if head.startswith(">") else cls._iter_plain(file_path)

            for sid, s in it:
                if length is None or len(s) == length:
                    picked.append(cls(s, sid))
                    if count is not None and len(picked) >= count:
                        return picked
        except Exception as e:
            raise IOError(f"Error reading protein file {file_path}: {e}")

        return picked
    
class MPAMAlignmentDistance(DistanceFunction):
    """
    使用 UMAD 的 mPAM 替换矩阵作为替换代价；
    """
    _AA_IDX = {aa:i for i,aa in enumerate(_AA_ORDER)}

    _MPAM = [
        # A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  X
        [ 0, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 5, 4, 2, 7], # A
        [ 2, 0, 2, 2, 4, 2, 2, 2, 2, 3, 3, 2, 2, 4, 2, 2, 2, 4, 4, 3, 7], # R
        [ 2, 2, 0, 2, 4, 2, 2, 2, 2, 3, 3, 2, 2, 4, 2, 2, 2, 5, 4, 2, 7], # N
        [ 2, 2, 2, 0, 4, 2, 2, 2, 2, 3, 3, 2, 3, 4, 2, 2, 2, 6, 4, 2, 7], # D
        [ 3, 4, 4, 4, 0, 4, 4, 3, 4, 3, 4, 4, 4, 4, 3, 3, 3, 7, 3, 3, 7], # C
        [ 2, 2, 2, 2, 4, 0, 2, 2, 2, 3, 3, 2, 2, 4, 2, 2, 2, 5, 4, 3, 7], # Q
        [ 2, 2, 2, 2, 4, 2, 0, 2, 2, 3, 3, 2, 3, 4, 2, 2, 2, 6, 4, 2, 7], # E
        [ 2, 2, 2, 2, 3, 2, 2, 0, 2, 2, 3, 2, 2, 4, 2, 2, 2, 6, 4, 2, 7], # G
        [ 2, 2, 2, 2, 4, 2, 2, 2, 0, 3, 3, 2, 3, 3, 2, 2, 2, 5, 3, 3, 7], # H
        [ 2, 3, 3, 3, 3, 3, 3, 2, 3, 0, 1, 3, 2, 2, 2, 2, 2, 5, 3, 2, 7], # I
        [ 2, 3, 3, 3, 4, 3, 3, 3, 3, 1, 0, 3, 1, 2, 3, 3, 2, 4, 2, 1, 7], # L
        [ 2, 2, 2, 2, 4, 2, 2, 2, 2, 3, 3, 0, 2, 4, 2, 2, 2, 4, 4, 3, 7], # K
        [ 2, 2, 2, 3, 4, 2, 3, 2, 3, 2, 1, 2, 0, 2, 2, 2, 2, 4, 3, 2, 7], # M
        [ 3, 4, 4, 4, 4, 4, 4, 4, 3, 2, 2, 4, 2, 0, 4, 3, 3, 3, 1, 2, 7], # F
        [ 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2, 4, 0, 2, 2, 5, 4, 2, 7], # P
        [ 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2, 3, 2, 0, 2, 5, 4, 2, 7], # S
        [ 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 0, 5, 3, 2, 7], # T
        [ 5, 4, 5, 6, 7, 5, 6, 6, 5, 5, 4, 4, 4, 3, 5, 5, 5, 0, 4, 5, 7], # W
        [ 4, 4, 4, 4, 3, 4, 4, 4, 3, 3, 2, 4, 3, 1, 4, 4, 3, 4, 0, 3, 7], # Y
        [ 2, 3, 2, 2, 3, 3, 2, 2, 3, 2, 1, 3, 2, 2, 2, 2, 2, 5, 3, 0, 7], # V
        [ 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0], # X/OTHER
    ]

    def __init__(self, gap_cost: int = 7):
        self.gap = int(gap_cost)

    @classmethod
    def _sub_cost(cls, a: str, b: str) -> int:
        ia = cls._AA_IDX.get(a, cls._AA_IDX["X"])
        ib = cls._AA_IDX.get(b, cls._AA_IDX["X"])
        return cls._MPAM[ia][ib]

    def __call__(self, a: ProteinSequence, b: ProteinSequence) -> float:
        assert isinstance(a, ProteinSequence) and isinstance(b, ProteinSequence), "type mismatch"
        s, t = a.seq, b.seq
        n, m = len(s), len(t)
        dp = [ [0]*(m+1) for _ in range(n+1) ]
        for i in range(1, n+1): 
            dp[i][0] = i * self.gap
        for j in range(1, m+1): 
            dp[0][j] = j * self.gap
        for i in range(1, n+1):
            si   = s[i-1]
            di_1 = dp[i-1]
            di   = dp[i]
            dim1 = dp[i-1]
            for j in range(1, m+1):
                tj    = t[j-1]
                sub   = dim1[j-1] + self._sub_cost(si, tj)
                ins   = di[j-1]   + self.gap
                dele  = di_1[j]   + self.gap
                di[j] = min(sub, ins, dele)
        return float(dp[n][m])