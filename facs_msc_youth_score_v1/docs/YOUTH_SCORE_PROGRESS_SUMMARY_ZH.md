# Youth Score 项目进展总结

**更新日期：2026-07-17**  
**研究范围：TMS MSC 年龄评分、跨技术敏感性检查及 GSE176206 外部验证**

## 1. 执行摘要

本阶段使用本地 Tabula Muris Senis（TMS）FACS 数据中的两个组织特异性 MSC 群体建立 Youth Score：SCAT adipose MSC 作为主分析，Limb Muscle MSC 作为次级分析。我们比较了 Gene-signature、Elastic Net 和 Gene-token Transformer 三种候选模型，并使用 Technical-only 逻辑回归作为技术混杂诊断对照。

所有内部性能均来自**按小鼠分组、按年龄分层的五折交叉验证**所产生的样本外预测，而不是训练集拟合成绩。SCAT 最终选择 Elastic Net（供体级 ROC-AUC 0.946），Limb Muscle 最终选择 Gene-signature（ROC-AUC 1.000）。Transformer 在两个细胞群中都学到了年龄信号，但没有优于简单模型。

将 Limb FACS 模型应用于 TMS Droplet Limb MSC 后，供体级 Youth Score 与月龄呈较强负相关（Spearman \(\rho=-0.856\)），支持年龄方向的跨技术可迁移性，但该分析仍受年龄、性别、供体数和技术平台混杂影响。

在未参与训练的 GSE176206 中，Young control 的平均 Youth Score 在 Adipo 和 MSC 中均高于 Aged control。MSC 的年龄差值为 +0.030，95% 供体 bootstrap 置信区间为 +0.0135 至 +0.0453；Adipo 在处理后元数据中每个条件只有一个混合文库，因此不能估计供体级置信区间。外部 SOKM 分析没有显示老年细胞沿当前 TMS 定义的 Youth Score 向年轻方向移动。

## 2. 数据与分析设计

### 2.1 正式训练数据

SCAT 和 Limb Muscle 不是把同一个 MSC 数据集平均切成两半，而是从 FACS 数据中独立筛选出的两个组织特异性 MSC 群体。

| 分析群体 | 技术 | 细胞数 | Young 小鼠 | Old 小鼠 | 年龄定义 | 角色 |
|---|---:|---:|---:|---:|---|---|
| SCAT adipose MSC | FACS | 1,570 | 7 | 8 | Young=3m；Old=18m/24m | 正式主分析 |
| Limb Muscle MSC | FACS | 935 | 6 | 8 | Young=3m；Old=18m/24m | 正式次级分析 |

训练标签为 Young/Old，Youth Score 定义为经过折内校准的 `P(young)`。在外部平台中，由于存在数据域偏移，分数主要用于同一外部数据内部的相对比较，不能直接解释为生物学年龄概率。

### 2.2 四种分析策略

| 方法 | 类别 | 在本项目中的作用 |
|---|---|---|
| Gene-signature | 生物信息学评分＋逻辑回归 | 每折仅用训练小鼠寻找偏年轻和偏老年的基因并建立一维评分 |
| Elastic Net | 正则化逻辑回归 | 从基因表达中学习稀疏线性权重 |
| Gene-token Transformer | 深度学习 Transformer | 学习基因身份、表达分箱和表达排序之间的非线性组合，同时辅助预测月龄 |
| Technical-only | 经典逻辑回归 | 只使用总计数、检测基因数和性别，诊断技术混杂 |

前三种是正式候选模型。Technical-only 是诊断对照，**不参与正式模型竞争**。

### 2.3 五折交叉验证的准确含义

我们使用的是**按小鼠分组、按年龄分层的五折外层交叉验证，并在每个外层训练部分保留内部验证小鼠**：

1. 以 `mouse.id` 为最小拆分单位，同一只小鼠的细胞不会出现在不同集合中。
2. 15只 SCAT 小鼠或14只 Limb 小鼠被分成5个外层折；每一折轮流作为独立测试折。
3. 其余小鼠中保留至少1只 Young 和1只 Old 小鼠，用于 early stopping、概率校准和阈值选择。
4. HVG、标准化、Gene-signature、表达分箱和概率校准只允许使用相应训练/验证小鼠。
5. 五折结束后，每只小鼠都恰好获得一次完全样本外预测；最终 AUC 由所有供体的 out-of-fold（OOF）预测合并计算。

因此，这不是一次性的“50%训练、50%验证”，也不是随机拆分细胞。

## 3. TMS FACS 内部验证

| 方法 | SCAT AUC | Limb Muscle AUC |
|---|---:|---:|
| Gene-signature | 0.786 | **1.000** |
| Elastic Net | **0.946** | 0.979 |
| Transformer | 0.857 | 0.875 |
| Technical-only | 0.696 | 0.688 |

![Internal donor-level ROC-AUC comparison](../outputs/youth_score/progress_summary/figures/figure_1_internal_auc.png)

*图1。点为合并 OOF 供体级 ROC-AUC，横线为2,000次供体 bootstrap 的95%区间，虚线为随机水平0.5，黑色描边表示最终选择的模型。矢量版本：[PDF](../outputs/youth_score/progress_summary/figures/figure_1_internal_auc.pdf)。*

### 3.1 结果解释

- SCAT 的正式模型是 Elastic Net：ROC-AUC 0.946，95%区间0.821–1.000，balanced accuracy 0.804，Brier score 0.088。
- Limb Muscle 的正式模型是 Gene-signature：ROC-AUC 1.000，balanced accuracy 1.000，Brier score 0.047。
- 三种基因表达候选模型在两个群体中的 AUC 点估计均高于 Technical-only，说明正式模型利用了超出测序深度、检测基因数和性别的信息。
- 由于供体数只有15和14，多个 AUC 区间较宽且彼此重叠；因此不应把“点估计更高”表述为所有候选模型都已被证明在统计上显著优于 Technical-only。
- Transformer 没有显示性能优势。大量细胞不能替代独立小鼠，当前独立供体数可能不足以发挥复杂神经网络的优势。
- Limb 的 AUC=1.000 表示当前14只留出小鼠被完全排序，而不是证明未来样本必然达到100%准确率。

## 4. TMS Droplet 跨技术敏感性分析

Limb FACS 的正式 Gene-signature 模型被冻结后，被直接应用于独立处理的 TMS Droplet Limb Muscle MSC；Droplet 数据没有参与正式训练、特征选择或校准。

| 月龄 | 平均 Youth Score | 小鼠数 |
|---|---:|---:|
| 1m | 0.448 | 2 |
| 3m | 0.273 | 2 |
| 18m | 0.215 | 4 |
| 21m | 0.162 | 2 |
| 24m | 0.112 | 4 |
| 30m | 0.142 | 2 |

供体级 Spearman 相关系数为：

```text
rho = -0.856
```

![Droplet age trajectory](../outputs/youth_score/progress_summary/figures/figure_2_droplet_age_trajectory.png)

*图2。小点为单只小鼠，空心点及连线为各月龄小鼠的非加权平均。矢量版本：[PDF](../outputs/youth_score/progress_summary/figures/figure_2_droplet_age_trajectory.pdf)。*

### 4.1 结果解释

Youth Score 总体随月龄下降，方向符合预期；30m 相比24m略有回升，因此不是严格单调关系。该结果支持 FACS 模型跨到 Droplet 后仍保留年龄方向，但它应被称为**TMS 内部跨技术敏感性检查**，而不是完全独立的跨数据集验证。不同月龄的性别组成和小鼠数不一致，而且每个月龄只有2–4只小鼠，因此不能把相关系数单独解释为稳定的生物年龄效应。

## 5. GSE176206 外部年龄方向验证

[GSE176206](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176206) 是独立的小鼠部分重编程单细胞 RNA-seq 研究，包含 Young/Aged、Control/SOKM 的 adipogenic cells 和 limb muscle MSC。外部数据没有参与本项目的训练、特征选择或校准。

模型对应关系为：

| TMS 模型 | GSE176206 数据 | 匹配程度 |
|---|---|---|
| SCAT adipose MSC / Elastic Net | Adipo SOKM | 脂肪组织与谱系相关，但不是完全相同的 MSC 身份；属于近似跨细胞状态验证 |
| Limb Muscle MSC / Gene-signature | MSC SOKM | 同为肢体肌肉来源 MSC，匹配更直接 |

### 5.1 Young control 与 Aged control

| 外部数据 | Young control | Aged control | Young − Aged |
|---|---:|---:|---:|
| Adipo | 0.497 | 0.447 | +0.050 |
| MSC | 0.493 | 0.463 | +0.030 |

![External age-direction validation](../outputs/youth_score/progress_summary/figures/figure_3_external_age_validation.png)

*图3。左侧显示条件分数和可用的供体点，右侧显示 Young control − Aged control 的效应值。MSC 横线为95%供体 bootstrap 区间；Adipo 空心菱形表示处理后元数据中每个条件只有一个混合文库，无法提供供体级区间。矢量版本：[PDF](../outputs/youth_score/progress_summary/figures/figure_3_external_age_validation.pdf)。*

### 5.2 置信区间的正确解释

MSC 的差值为：

```text
Young control - Aged control = +0.02985
95% donor-bootstrap CI       = [+0.01350, +0.04532]
```

`+0.0135～+0.0453` 是 **MSC 差值 +0.030 自身的不确定性区间**。点估计位于这个区间之内；不能说 Adipo 或 MSC 的差值“部分或全部大于这个区间”，也不能把 Adipo 的差值与 MSC 的置信区间直接比较。由于 MSC 区间整体高于0，数据支持 MSC Young control 的得分高于 Aged control。

Adipo 在处理后元数据中每个条件只有一个混合文库。虽然有数千个细胞，但细胞不能替代独立小鼠，因此 +0.050 只能作为方向性结果，不能估计可靠的供体级置信区间。

### 5.3 可以得出的结论

两套外部数据中的年龄方向均正确，说明模型具有初步的跨研究迁移能力。这里验证的是**相对年龄方向**，不是外部分类准确率：我们没有从这些均值计算外部 accuracy 或 ROC-AUC，而且外部分数受平台和细胞状态偏移影响，不能机械地使用0.5阈值。

MSC 是相对更强、细胞类型更匹配的外部验证；Adipo 是证据较弱的跨细胞状态检查。

## 6. GSE176206 SOKM 分析

SOKM 效应定义为：

```text
Aged SOKM - Aged control
```

正值才表示处理后沿当前 TMS 定义的 Youth Score 向年轻方向移动。

| 数据 | Aged control | Aged SOKM | SOKM − Control |
|---|---:|---:|---:|
| Adipo | 0.447 | 0.436 | −0.011 |
| MSC | 0.463 | 0.344 | −0.118 |

![External SOKM contrast](../outputs/youth_score/progress_summary/figures/figure_4_external_sokm_contrast.png)

*图4。左侧显示老年对照与老年 SOKM 条件分数；右侧显示 Aged SOKM − Aged control。MSC 横线为95%供体 bootstrap 区间；Adipo 没有供体分辨率的区间。矢量版本：[PDF](../outputs/youth_score/progress_summary/figures/figure_4_external_sokm_contrast.pdf)。*

MSC 的处理差值为 −0.118，95%供体 bootstrap 区间为 −0.161 至 −0.055，整个区间低于0。Adipo 的差值为 −0.011，但缺少可用的独立供体区间。

因此，当前结果支持的严格表述是：

> 在 GSE176206 中，SOKM 没有使老年样本沿本项目从 TMS 自然衰老数据中学习到的 Youth Score 轴向年轻方向移动；MSC 中观察到的变化方向相反。

这不等于证明 SOKM 会加速衰老或在其他指标上无效。可能原因包括自然衰老轴与重编程响应轴不同、SOKM 短期应激、细胞身份暂时受抑、测序平台与培养条件的域偏移、Adipo 与 SCAT MSC 身份不完全一致，以及供体与细胞数不平衡。原研究使用其自身定义的年龄表达程序报告了年轻化方向；本结果检验的是该处理是否与**本项目独立训练的 TMS Youth Score**一致，两者不是同一个统计问题。[Roux et al., 2022](https://doi.org/10.1016/j.cels.2022.05.002)

## 7. 当前结论分级

| 科学问题 | 当前结果 | 证据强度 |
|---|---|---|
| TMS FACS 中能否区分 Young/Old MSC？ | 支持；SCAT AUC 0.946，Limb AUC 1.000 | 较强，但供体数有限 |
| Transformer 是否优于简单模型？ | 不支持 | 清晰 |
| Limb FACS 模型跨到 Droplet 后是否保持年龄方向？ | 支持，\(\rho=-0.856\) | 中等；存在混杂 |
| 外部 MSC 中 Young control 是否高于 Aged control？ | 支持，差值+0.030，95%区间不跨0 | 初步支持；每条件3只动物 |
| 外部 Adipo 中年龄方向是否正确？ | 方向正确，差值+0.050 | 较弱；每条件一个混合文库且细胞身份近似匹配 |
| SOKM 是否提高老年样本 Youth Score？ | 不支持；MSC 中方向相反 | MSC 较明确，Adipo 不确定 |

## 8. 可声明与不可声明的内容

### 可以声明

- 我们建立并验证了两个组织特异性的 TMS MSC Youth Score。
- 简单模型在当前供体规模下优于 Transformer。
- Limb 模型在 TMS Droplet 中保持了总体年龄下降趋势。
- 外部 MSC 的年轻对照得分高于老年对照。
- 当前模型没有在 GSE176206 中检测到 SOKM 驱动的 Youth Score 增加。

### 不能声明

- 不能把 Youth Score 解释为临床年龄、治疗安全性、因果年轻化或功能恢复。
- 不能把数千个细胞当作数千个生物学重复。
- 不能从 Limb AUC=1.000 推断未来数据必然100%准确。
- 不能从负的 SOKM 差值推断 SOKM 一定有害或加速衰老。
- 不能把 Adipo 结果描述为完全匹配的 SCAT MSC 外部验证。

## 9. 文件与复现

科研图同时提供300 dpi PNG和矢量 PDF。全部图由以下英文脚本从现有结果表自动生成：

```powershell
.\.venv\Scripts\python.exe scripts\plot_youth_score_progress.py
```

- 绘图脚本：[`scripts/plot_youth_score_progress.py`](../scripts/plot_youth_score_progress.py)
- 图像目录：[`outputs/youth_score/progress_summary/figures`](../outputs/youth_score/progress_summary/figures)
- SCAT 指标：[`outputs/youth_score/scat_facs/model_metrics.csv`](../outputs/youth_score/scat_facs/model_metrics.csv)
- Limb 指标：[`outputs/youth_score/limb_facs/model_metrics.csv`](../outputs/youth_score/limb_facs/model_metrics.csv)
- Droplet 供体结果：[`outputs/youth_score/limb_droplet_sensitivity/donor_scores.csv`](../outputs/youth_score/limb_droplet_sensitivity/donor_scores.csv)
- 外部验证汇总：[`outputs/youth_score/external/GSE176206/combined_summary.json`](../outputs/youth_score/external/GSE176206/combined_summary.json)

