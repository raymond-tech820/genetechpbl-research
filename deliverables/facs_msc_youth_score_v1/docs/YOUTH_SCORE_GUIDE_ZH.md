# TMS Youth Score：可复现训练与验证指南

## 1. 目标

本项目使用小鼠单细胞 RNA 测序数据建立研究级 Youth Score。分数是 0–1 之间经过校准的概率：

- `1.0`：表达状态更接近 TMS 3 月龄参考细胞。
- `0.0`：表达状态更接近 TMS 18–24 月龄参考细胞。

Transformer 同时预测月龄，作为约束表征的辅助任务；月龄输出不取代 Youth Score。

该分数不是 Risk Score，也不是因果年轻化证据、临床年龄或治疗建议。Youth Score 高并不能证明处理安全，也不能证明细胞身份得到保留。

## 2. 数据集设计

### 2.1 主模型

- 数据集：TMS FACS。
- 组织：`SCAT`。
- 细胞类型：`mesenchymal stem cell of adipose`。
- Young：3 月龄，737 个细胞、7 只小鼠。
- Old：18 和 24 月龄，833 个细胞、8 只小鼠。

SCAT 被设为主要参考，是因为其 Young/Old 细胞数和生物学重复相对平衡，并且存在公开的 adipogenic partial-reprogramming 数据。

### 2.2 次级模型

- 数据集：TMS FACS。
- 组织：`Limb_Muscle`。
- 细胞类型：`mesenchymal stem cell`。
- Young：3 月龄，264 个细胞、6 只小鼠。
- Old：18 和 24 月龄，671 个细胞、8 只小鼠。

Limb 模型与外部 MSC 实验的细胞类型更接近。由于细胞较少，而且年龄与性别、文库质量存在明显混杂，因此作为次级模型。

### 2.3 敏感性数据

TMS Droplet Limb Muscle MSC 共 13,037 个细胞，覆盖 1、3、18、21、24 和 30 月龄。主要年龄比较只有 2 只 Young 和 10 只 Old 小鼠，而且年龄与性别高度混杂，因此不用于训练，只用于检查跨技术迁移和不同月龄的分数顺序。

### 2.4 外部验证

流程自动下载两个 GSE176206 处理后文件：

- `GSE176206_adipo_sokm.h5ad.gz`：验证 SCAT 模型。
- `GSE176206_msc_sokm.h5ad.gz`：验证 Limb 模型。

外部数据绝不参与特征选择、模型训练、early stopping 或概率校准。主要检查：

1. Young control 的分数是否高于 Aged control。
2. Aged SOKM 是否相对 Aged control 向年轻方向移动。

只要 metadata 允许，结果均按真实生物学重复聚合；不会把细胞级 p 值当作生物学重复证据。

数据与软件来源：

- [TMS CELLxGENE 页面](https://cellxgene.cziscience.com/e/f16a8f4d-bc97-43c5-a2f6-bbda952e4c5c.cxg/)
- [BPCells Python `DirMatrix` 文档](https://bnprks.github.io/BPCells/python/generated/bpcells.experimental.DirMatrix.html)
- [GSE176206 GEO 页面](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176206)
- [PyTorch 2.13 发布说明](https://pytorch.org/blog/pytorch-2-13-release-blog/)

## 3. 为什么 Droplet 和 FACS 分开

Droplet 与 FACS/Smart-seq2 在捕获化学、测序深度、dropout、细胞选择、基因覆盖和 metadata 上均不同，因此不会直接合并成一个训练矩阵。两个正式模型分别使用 FACS，Droplet Limb 仅作为敏感性域。

Tokenizer 使用细胞内表达排名和分箱后的归一化表达，可以减轻但不能消除平台差异。

## 4. 数据整理

本地 BPCells 压缩包会流式解压到 `data/extracted/`，原始文件不会被修改。

每个目标细胞群都会验证：

- 矩阵维度。
- BPCells 列名与 metadata `index` 的内容及顺序完全一致。
- tissue、cell type、age、cell 和 mouse 数量符合声明。
- 文库非零且矩阵计算的检测基因数不少于 500。

整理后的训练包位于 `data/processed/youth_score/<dataset_id>/`：

- `counts.npz`：cell × gene 稀疏原始计数。
- `cells.parquet`：标签、小鼠、年龄、性别、QC 和原始位置。
- `genes.parquet`：基因和主分析排除标记。
- `donor_folds.csv`：固定的外层供体分折。
- `manifest.json`：输入输出哈希值和运行环境。

Droplet 3 月龄 Limb metadata 的 `n_genes` 缺失，因此 QC 指标统一从表达矩阵重新计算。

## 5. 特征处理与混杂控制

每个外层折都只用训练小鼠重新选择特征：

1. library-size 归一化到 10,000。
2. `log1p` 转换。
3. 稀疏方差计算。
4. 最多选择 4,096 个 highly variable genes。

主分析排除线粒体、核糖体、血红蛋白、已知 Y 连锁基因以及 `Xist/Tsix`，以减少明显的测序质量和性别捷径，但这不能去除全部性别效应。

流程不做 batch correction。TMS 中年龄、测序板、文库复杂度和性别部分混杂；如果校正变量几乎等同于年龄，既可能删除真实衰老信号，也可能制造人工信号。因此流程改为报告：

- 只使用文库量、检测基因数和性别的技术模型。
- male-only 分数方向。
- 3 月龄与 18 月龄的敏感性比较。
- 跨研究的外部验证。

若正式表达模型与技术模型的 ROC-AUC 相差不足 0.02，或任一敏感性方向反转，则结果标记为 `confound_limited`。

## 6. Transformer 架构

每个细胞会转换成一个序列，包含该折词表中表达最高的 256 个基因。

每个 token 由以下部分相加：

- 基因身份 embedding。
- 64 级表达分箱 embedding。
- 细胞内表达排名 embedding。

序列前添加可学习的 `[CLS]`。Encoder 包含：

- 2 层 pre-norm Transformer。
- hidden width 128。
- 8 个 attention heads。
- feed-forward width 512。
- GELU。
- dropout 0.2。
- 训练时随机丢弃 15% gene tokens。

两个输出头分别预测：

- Young/Old logit。
- 月龄，训练时除以 30 缩放。

损失为 donor-balanced binary cross-entropy 加 `0.25 × Huber age loss`。两个年龄类别总权重相同；类别内每只小鼠总权重相同；同一小鼠的细胞共享其权重。

## 7. 交叉验证与校准

主模型和次级模型均采用五折、按小鼠分组并按年龄分层的外层交叉验证。在一次训练中，同一只小鼠只能属于 train、validation 或 test 之一。

每个外层折中：

- 从可训练小鼠中保留一只 Young 和一只 Old 用于 validation。
- HVG、表达分箱、模型拟合、early stopping 和校准都不能使用 test 小鼠。
- 训练三个随机种子的 Transformer 并取平均。
- temperature scaling 只在 validation 小鼠上拟合。
- 二值决策阈值只根据 validation 供体均值选择，以 balanced accuracy 最大为目标；并列时优先选择最接近 0.5 的阈值。

Out-of-fold 细胞概率会按小鼠取平均。主指标是 donor-level ROC-AUC，同时报告 PR-AUC、balanced accuracy、F1、Brier score、年龄相关性和 2,000 次供体 bootstrap 区间。

## 8. 基线与正式模型选择

流程生成三个对照：

- 基于 donor-pseudobulk 的 Young/Old gene signature。
- Elastic-net logistic regression。
- 仅用于诊断的技术模型。

技术模型不能成为正式 Youth Score。在 gene signature、elastic net 和 Transformer 之间：

1. 选择 donor-level out-of-fold ROC-AUC 最高者。
2. AUC 相差小于 0.02 时，选择 Brier score 更低者。
3. Brier 仍相差小于 0.02 时，优先选择更简单模型。

即使 Transformer 未获胜，也会保留全部 checkpoint、预测和评估结果。

## 9. 安装

在项目根目录的 PowerShell 中运行：

```powershell
.\scripts\setup.ps1
```

只使用 CPU：

```powershell
.\scripts\setup.ps1 -CpuOnly
```

脚本会创建 `.venv`、安装固定依赖、以 editable 模式安装项目并执行环境检查。

## 10. 命令

检查环境：

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli doctor
```

整理全部本地数据：

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli prepare-tms --config configs/scat_primary.yaml --config configs/limb_secondary.yaml --config configs/limb_droplet_sensitivity.yaml
```

训练和评估一个模型：

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli train --config configs/scat_primary.yaml --device cuda
.\.venv\Scripts\python.exe -m youth_score.cli evaluate --config configs/scat_primary.yaml
```

运行保留全部基因的特征掩码敏感性分析：

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli train --config configs/scat_primary.yaml --feature-mode all --device cuda
```

下载并验证 GSE176206：

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli download-external
.\.venv\Scripts\python.exe -m youth_score.cli validate-external --device cuda
```

运行全部可续传阶段：

```powershell
.\scripts\run_pipeline.ps1 -Device cuda
```

为训练包或 h5ad 评分：

```powershell
.\.venv\Scripts\python.exe -m youth_score.cli score --model-dir outputs/youth_score/scat_facs --input path/to/input.h5ad --output scores.parquet --device cuda
```

## 11. 输出

每个训练数据集输出到 `outputs/youth_score/<dataset_id>/`：

- `oof_cell_scores.parquet`。
- `oof_donor_scores.csv`。
- `model_metrics.csv`。
- `selection.json`。
- `training_report.md`。
- `training_report.json`。
- `figures/`。
- `all_genes_sensitivity/`：关闭主分析排除掩码后的完整平行 OOF 结果。
- `folds/fold_*/`：tokenizer、拆分、基线、校准、Transformer checkpoint 和训练历史。

公共评分输出包括：

- `cell_id`。
- `youth_score`。
- Transformer ensemble 给出的 `predicted_age_months`。
- Youth Score 使用的 `model_id`。
- `gene_overlap`。
- `qc_status`。
- 各折的 `decision_threshold` 保存在 OOF 输出中；它只影响 balanced accuracy/F1，不改变连续 Youth Score。

## 12. 证据状态解释

- `internally_supported`：模型能区分留出 TMS 小鼠，并且未触发预定义混杂警告。
- `confound_limited`：技术预测器几乎同样强，或性别/年龄敏感性方向反转。
- 外部支持单独报告；只有 Young control 与 Aged control 的基本检查通过后，才解释 SOKM 的移动。

软件运行成功并不等于生物学验证成功。无显著结果、基线获胜或混杂警告都是有效且重要的结果。

## 13. 已完成的参考运行（2026-07-17）

本项目已经在完整目标队列上执行，而不只是完成了合成数据测试。

### 13.1 内部供体级验证

| 模型域 | 正式模型 | 供体 ROC-AUC | Brier | Transformer ROC-AUC | 纯技术模型 ROC-AUC | 状态 |
|---|---:|---:|---:|---:|---:|---|
| FACS SCAT adipose MSC | Elastic net | 0.9464 | 0.0883 | 0.8571 | 0.6964 | `internally_supported` |
| FACS Limb Muscle MSC | Gene signature | 1.0000 | 0.0467 | 0.8750 | 0.6875 | `internally_supported` |

因此，Transformer 被完整保留为比较模型和辅助月龄预测器，但没有被强制指定为正式 Youth Score。Droplet Limb 敏感性分析中，供体级 Youth Score 与实际月龄的 Spearman 相关为 -0.8564，整体年龄方向通过；但它仍然只是带有技术平台、性别和供体构成混杂的敏感性结果。

全基因敏感性运行保留了主分析排除的线粒体、核糖体、血红蛋白、Y 连锁、`Xist` 和 `Tsix` 基因。两个模型域的正式模型均未改变：SCAT elastic net 的 AUC 仍为 0.9464，Limb gene signature 的 AUC 仍为 1.0000。SCAT Transformer 的 AUC 从 0.8571 降至 0.6964，而 Limb Transformer 保持在 0.8750。因此，正式模型选择对该特征掩码敏感性是稳健的，但 SCAT Transformer 表征并不稳健。

### 13.2 GSE176206 预检与验证

| 外部模型域 | 最小基因重叠 | 与 TMS 参考表达的 Spearman | Young control - aged control | Aged SOKM - aged control | 重复情况 |
|---|---:|---:|---:|---:|---|
| Adipo / SCAT 模型 | 0.9431 | 0.5710 | +0.0498 | -0.0109 | 每个条件只有一个 pooled library；无供体置信区间 |
| MSC / Limb 模型 | 0.9421 | 0.5319 | +0.0298（95% CI +0.0135 至 +0.0453） | -0.1185（95% CI -0.1612 至 -0.0551） | 每个条件可识别 3 只小鼠 |

两个外部数据的 control 对比都通过了基本年龄方向检查。但两个 aged SOKM 对比都没有向更高的 Youth Score 移动；MSC 结果在该模型下显著向相反方向移动。因此，本次运行支持模型跨研究识别 young 与 aged control 的方向，但**不支持**用该 Youth Score 验证 SOKM 相关的年轻化。

Adipo h5ad 没有提供单只小鼠标识，所以条件级结果只能作描述性分析，并标记为 `replication_limited`。MSC h5ad 使用 `animal` 作为生物学单位；供体未知的细胞仍可获得细胞级分数，但不进入供体汇总。处理映射中，只有 `Tg+/Dox+` 被定义为 SOKM；`Tg+/Dox-`、`Tg-/Dox+` 和 `NegCtrl` 均定义为 control。

## 14. 故障排查

- BPCells 导入失败：重新运行 `scripts/setup.ps1`；本项目依赖固定的 Windows CPython 3.12 wheel。
- CUDA 不可用：运行 `doctor` 检查 NVIDIA 环境，或使用 `--device cpu` 做测试。
- 下载中断：重新运行 `download-external`；`.part` 文件会通过 HTTP Range 续传。
- 外部基因重叠低于 70%：检查 `metadata_schema.json` 及 h5ad 的 gene-symbol 列。流程会停止，不会静默生成无效分数。
- 出现 `confound_limited`：不要删除警告；在使用分数前检查技术模型和外部 control 对比。
