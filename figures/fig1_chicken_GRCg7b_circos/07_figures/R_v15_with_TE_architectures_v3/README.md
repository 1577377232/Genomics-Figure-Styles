# Figure 1 拼接 V3：五个模块与五条独立引线

本版从 V2 继续精修，保留 V15 主图、全部原始数据/窗口/变换/颜色语义及四类 TE 结构 SVG；V2 脚本和输出也保持原样。只新增独立 V3 拼接脚本、更新 TEV 模式 SVG，并生成新目录中的最终输出。

## 运行

在 `figures/fig1_chicken_GRCg7b_circos/` 下执行：

```bash
# 缺少 R 依赖时先执行；不修改 V15 的 R 库：
bash 06_scripts/install_composition_R_dependencies.sh
# 完整 SVG / PDF / 600 dpi PNG / 150 dpi preview：
bash 06_scripts/run_compose_v15_with_te_architectures_v3.sh
# 仅预览或另设输出目录：
bash 06_scripts/run_compose_v15_with_te_architectures_v3.sh --preview-only
bash 06_scripts/run_compose_v15_with_te_architectures_v3.sh --out-dir /absolute/new_output
```

运行依赖为 R、rsvg、png；沿用原拼接安装环境（安装器也保留旧版所需的 svglite）。V15 绘图脚本不执行，无需原始 VCF、基因组比对文件或重新统计。

## 主图与右侧素材

- 主图仍为 `07_figures/R_v15/chicken_GRCg7b_R_v15_1000kb_adaptive_discovery_TIP-extended_SV-all_palette-pangenome_spectrum.svg`，完整嵌入原生矢量内容及全部325个裁剪路径。
- 前四幅分别读取 `assets/TE_architectures/LINE_CR1.svg`、`LTR_ERV.svg`、`SINE_generic.svg`、`DNA_TIR.svg`，这些源文件按字节保持不变。
- 第五幅读取本轮重绘的 `assets/TE_mediated_SV/TE_mediated_SV_classification.svg`，1400×250，全部文字及形状可编辑。V2 的 >50% 专用副本不变。修改前的 >80% SVG 保存在本目录 `provenance/TE_mediated_SV_classification.before_v3.svg`，也可通过 Git 历史取得。

## 布局和引线

14×10英寸画布及左图尺寸/位置保持原版。五个模块整体向下移动0.60英寸；顶部依次为1.65、3.20、4.75、6.30、7.85英寸。所有模块宽3.95英寸，等比例缩放，主体中心间距统一1.55英寸。

LINE/LTR/SINE/DNA 的四条线从右上开口中对应第4/5/6/7轨附近引出。第五条从最内层 TIP（中心图例 TEV density，第10轨）靠右下缺口处单独引出，避开现有数据扇区和纵轴；不从其他轨道借用起点。第五条的首段沿缺口径向延伸，以安全越过主图外沿后接入 TEV 模块，其首段不强制为水平短线。

五条线均为独立三段骨架、两个二次 Bézier 圆角、小型闭合箭头，描边 #253246、0.80pt。没有共同主干、合并顶点或交叉；箭头在对应模块左端前停止，不遮住图形。圆角路径由 R 直接输出为原生 SVG `<path>`，保留可编辑节点；没有使用手工拼图或位图替代。

`composition_config()` 是位置、首段、转折、圆角、箭头等参数入口。`connector_vertices.tsv` 记录三段骨架，`connector_paths.tsv` 记录真正导出的 SVG 路径，`connector_curve_samples.tsv` 记录曲线采样点。

## TEV 示意阈值与真实数据口径（必须区分）

**用户本轮要求示意显示 >80%，但不可修改 V15 数据；因此本图右侧 >80% 仅是所要求的分类规则示意，主图 TEV 轨道仍是原来的 >50%、事件长度 >50 bp 的 extended 口径。没有按80%重新筛选或重算主图。二者不可解释为相同过滤条件。**

TEV SVG 的 SV 区间为 x=110..1000，80% 阈值位置严格为110+890×0.80=822；示例 TE 覆盖长度760/890≈85.4%，确实超过该示意阈值。珊瑚色主体与浅灰剩余序列、简短 SV bracket、>80% 刻度、分类箭头及 TEV 小模块保留；删除原有重复竖线装饰。标题与前四图使用同款58px Georgia，模块采用浅色圆角矩形。不添加底部说明句。

这是覆盖重叠分类示意，不证明 SV 的分子形成机制；形状不代表真实序列长度。`scientific_scope.tsv` 明确记录示意与数据阈值差异。若用于正文图注，需说明该差异；不得宣称当前主图为80%筛选结果。

## 输出和审计

- SVG、PDF：原生矢量，未嵌入栅格图片；SVG 保留可编辑 TE 文本。
- PNG：8400×6000，600dpi；preview：2100×1500，150dpi。
- `protected_baseline.json`：修改前 V15、V2 和四类 TE 的 SHA256 基线。
- `input_preservation.tsv`：本次运行前后文件 MD5。
- `validation_report.json`：独立解析实际 SVG M/L/Q 路径，检查曲线间距、数据及文字框避让、源图保真、PDF/PNG 格式。
- `provenance/validate_composition_v3.py`：全量审计脚本，额外依赖 Python/NumPy/SciPy/Pillow 及 Poppler；它们不属于 R 拼接必需依赖。

R 运行时先检查曲线不穿入任何数据/核型/坐标扇区、避开带1.2pt保护边距的文字框；独立 Python 审计直接解析导出的二次 Bézier 路径，以最大0.08pt步长采样，并对线间距扣除0.16pt给出保守下界。V15 和 V2 文件不被覆盖。
