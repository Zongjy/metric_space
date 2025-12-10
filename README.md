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
├── query/
│   ├── linear_scan.py       # 线性扫描查询算法
│   └── pivot_table.py       # Pivot Table 索引结构
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
\displaystyle \sum_{t=1}^{i} Score(S_1(t), \text{gap}), & j = 0 \\
\displaystyle \sum_{t=1}^{j} Score(S_2(t), \text{gap}), & i = 0 \\
\displaystyle \min \Big(
E(i-1,j) + Score(S_1(i), \text{gap}), \\
\qquad\quad E(i,j-1) + Score(S_2(j), \text{gap}), \\
\qquad\quad E(i-1,j-1) + Score(S_1(i), S_2(j))
\Big), & i=1,\ldots,n,\ j=1,\ldots,m
\end{cases}
$$

⸻

### 4. 查询模块 (query)

实现了度量空间相似性查询算法，包括线性扫描和基于 `Pivot Table` 的索引方法。

#### 线性扫描 (`LinearScan`)

直接遍历所有对象进行查询，保证精确结果：

- **范围查询** (`range_query(query, radius)`): 找出与查询对象距离小于等于 radius 的所有对象
- **kNN 查询** (`knn_query(query, k)`): 找出距离查询对象最近的 k 个对象
- **dkNN 查询** (`dknn_query(query, k, max_distance)`): 在距离约束内找出最多 k 个最近邻

#### Pivot Table 索引 (`PivotTable`)

利用预选的支撑点（pivot）和三角不等式进行剪枝，减少距离计算次数：

**三角不等式剪枝原理**：对于查询对象 $q$、数据对象 $o$ 和支撑点 $p$，有：

$$
|d(q,p) - d(o,p)| \leq d(q,o) \leq d(q,p) + d(o,p)
$$

若 $|d(q,p) - d(o,p)| > r$（查询半径），则可以直接剪枝掉对象 $o$。

**支撑点选择策略**：
- `random`: 随机选择支撑点
- `farthest`: 迭代选择与已有支撑点距离最远的对象
- `incremental`: 增量选择，最大化最小距离

**查询方法**：
- `range_query(query, radius)`: 带剪枝的范围查询
- `knn_query(query, k)`: 带剪枝的 kNN 查询
- `dknn_query(query, k, max_distance)`: 带剪枝的 dkNN 查询

**性能特点**：
- 预计算支撑点距离表，空间复杂度 $O(n \times m)$（$n$ 为数据对象数，$m$ 为支撑点数）
- 查询时可显著减少距离计算次数（通常可减少 80-95%）
- 支撑点数量和选择策略对性能有重要影响

⸻

## 测试

项目包含完整的测试套件，验证各项查询功能的正确性：

```shell
python tests/test_queries.py
```

测试覆盖：
- 线性扫描：范围查询（4个测试）、kNN 查询（4个测试）、dkNN 查询（4个测试）
- Pivot Table：正确性验证（4个测试）、性能分析（3个测试）
- 共 19 个测试用例，全面验证功能正确性和性能优化效果

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