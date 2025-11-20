# Metric Space Data Processing System

一个用于处理各种度量空间对象（向量、蛋白质序列等）并计算它们间距离的模块化系统。
系统遵循统一的 MetricObject / DistanceFunction 接口，支持扩展自定义对象与自定义距离。

## 项目结构

```text
metric_genhierarchy/
├── core/
│   └── metric_space.py      # 抽象基类：MetricObject / DistanceFunction
├── impl/
│   ├── vector.py            # 向量对象与 Minkowski 距离
│   └── protein.py           # 蛋白序列对象与 mPAM 比对距离
└── main.py                  # CLI 主入口
```

⸻

## 模块说明

### 1. 度量空间基类 (core/metric_space.py)

`MetricObject`：所有度量对象的基类
- `from_file(path, **kwargs)`：从文件加载对象
- `__repr__()`：打印对象信息

`DistanceFunction`：距离函数基类
- 需实现：
- `__call__(obj1, obj2)`：计算两个对象间距离


### 2. 向量模块 (impl/vector.py)

`VectorObject`：
- 多维数值向量
- 支持从 .txt, .csv加载

`MinkowskiDistance`：
- $p = 1 \to L_1$，也即曼哈顿距离；
- $p = 2 \to L_2$，也即欧几里得距离；
- $p = \infty \to L_\infty$，即切比雪夫距离；
- $p$ 为其余值，有 $d_p(\mathbf{x},\mathbf{y})=\left(\sum_{i=1}^{n}|x_i-y_i|^p\right)^{1/p}$


### 3. 蛋白序列模块 (impl/protein.py)

`ProteinSequence`：处理 FASTA 格式或纯序列文本，并自动解析 ID、序列长度等信息

`MPAMAlignmentDistance`：使用 mPAM 打分矩阵的序列比对距离，即：

$$
E(i,j)=
\begin{cases}
\displaystyle \sum_{t=1}^{i} Score(S_1(t), \text{gap}), & j = 0 \\[10pt]
\displaystyle \sum_{t=1}^{j} Score(S_2(t), \text{gap}), & i = 0 \\[12pt]
\displaystyle \min \Big(
E(i-1,j) + Score(S_1(i), \text{gap}), \\
\qquad\quad E(i,j-1) + Score(S_2(j), \text{gap}), \\
\qquad\quad E(i-1,j-1) + Score(S_1(i), S_2(j))
\Big), & i=1,\ldots,n,\ j=1,\ldots,m
\end{cases}
$$

⸻

## 使用方法（CLI）

环境部署可通过[uv](https://docs.astral.sh/uv/getting-started/installation/#standalone-installer)完成：
```shell
uv venv && source .venv/bin/activate
uv pip install -e .
```

基本命令格式：
```shell
python -m metric_genhierarchy.main <command> <file> [options]
```