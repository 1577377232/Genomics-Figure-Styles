# Figure 1：V15 + 最新四类 TE 结构 + TEV 分类示意（v2）

这是独立 R 拼接，不执行、不修改 `plot_chicken_circos_v15.R`，不改变任何 V15 统计、输入、轨道、配色或原始输出。旧拼接版也保留。

## 运行

在本 Fig1 目录下执行：

```bash
# 首次使用且未安装独立依赖时：
bash 06_scripts/install_composition_R_dependencies.sh
# 完整 SVG、PDF、600 dpi PNG 和小预览：
bash 06_scripts/run_compose_v15_with_te_architectures_v2.sh
# 也可直接用 R：
Rscript 06_scripts/compose_v15_with_te_architectures_v2.R
# 快速检查版式：
bash 06_scripts/run_compose_v15_with_te_architectures_v2.sh --preview-only
# 可使用新的输出目录：
bash 06_scripts/run_compose_v15_with_te_architectures_v2.sh --out-dir /absolute/new_output
```

依赖仍为 R/grid、svglite、rsvg、png；使用已有独立安装脚本，依赖安装到 `06_scripts/composition_R_library/`，不修改 V15 的库。已验证 R4.1.2、svglite2.1.3、rsvg2.7.0。无需重跑基因组分析。

## 输入与素材版本

主图固定为 `07_figures/R_v15/chicken_GRCg7b_R_v15_1000kb_adaptive_discovery_TIP-extended_SV-all_palette-pangenome_spectrum.svg`。PDF/PNG 作为原版保留，拼接优先读取 SVG。

本轮先获取 GitHub 最新素材，源提交为 **f926a01**。读取并保留 `TE_architectures_4up.svg`、资产 README 与原素材清单；实际拼接以下五个独立 SVG：

1. `assets/TE_architectures/LINE_CR1.svg`
2. `assets/TE_architectures/LTR_ERV.svg`
3. `assets/TE_architectures/SINE_generic.svg`
4. `assets/TE_architectures/DNA_TIR.svg`
5. `assets/TE_mediated_SV/TE_mediated_SV_classification_V15_gt50.svg`

四个结构 SVG 均为最新1400×250版本：居中衬线标题、浅色圆角模块、无底部注释小字，原始 SVG 文字、形状和样式保持。五个画布均等比缩放，统一宽3.95英寸，不拉伸。

## 四条独立右上引线

- LINE density → LINE；LTR density → LTR；SINE density → SINE；DNA density → DNA。
- 起点根据原V15的真实轨道半径和右上开口位置计算，径向位置分别位于第4、5、6、7轨。所有起点均在圆心右上方。
- 为避开原V15不可改的百分比标签，起点在各自轨道端点附近向开口内侧偏移约3–4mm；没有改动或遮盖纵轴文字。
- 每条线均为“短水平段—斜向过渡—水平箭头”三段，颜色#253246，闭合箭头。四条线的起点、折点均独立，没有共用主干，没有先汇聚再分叉。
- 连接线交叉、共享顶点、线间距离、数据区穿越及染色体/纵轴文字框均经过几何检查。详细点位见 `connector_vertices.tsv`；源轨道见 `connector_anchors.tsv`；六组线对的最短距离见 `connector_pair_QC.tsv`。
- 第五幅TEV图是分类规则示意，不是第五种TE结构类型，因此不额外增加一条连接到TEV轨道的外部引线。

## TEV 阈值副本

原仓库 `TE_mediated_SV_classification.svg` 标注 >80%。按照用户本轮明确选择，**只新增拼接用 >50% 副本**，原 >80% SVG 按字节保留。

副本只将可见阈值改为 >50%，同时将阈值线和标签从x=822移动到x=555（SV区间起点110、宽890，555=110+890×0.5），更新非可见标题/描述。其他形状、示例TE覆盖片段、颜色和尺寸不变。副本及原文件SHA256与授权范围记录在 `assets/TE_mediated_SV/TE_mediated_SV_classification_V15_gt50.provenance.json`。

该图是简化的重叠覆盖分类示意，不表示所有SV均由TE介导生成。V15原有的>50%、>50bp和extended范围、INS/DEL以及simple/low-complexity/satellite等规则继续沿用原始统计，没有使用此模式图重新分类或生成任何位点。

## 导出与保真

14×10英寸横向画布；PDF1008×720pt；PNG8400×6000、600dpi；预览2100×1500、150dpi。SVG/PDF均为矢量，无嵌入位图。

原V15 SVG整体嵌入并保留全部复杂裁剪，四个原结构图及TEV副本均保留可编辑SVG元素；仅给元素ID加前缀并添加等比外层定位。不能直接用grImport2默认重绘V15，否则会丢弃复杂扇区裁剪。箭头由grid/svglite绘制，再由rsvg导出同一SVG场景。

## 调整入口与 QC

只需修改新R脚本顶部的 `composition_config()`：结构位置、统一宽度、四条线的起点半径/角度、短水平段长度和折点。所有排版参数均为物理尺寸或图形角度，不涉及基因组坐标重算。

- `validation_report.json`：独立XML内容、原文件、交叉/遮挡、格式与DPI检查。
- `input_preservation.tsv`：原V15、原始资产及脚本在运行前后的MD5。
- `architecture_placements.tsv`：五幅SVG的位置与等比缩放。
- `connector_anchors.tsv`、`connector_vertices.tsv`、`connector_pair_QC.tsv`：独立引线证据。
- `composition_settings.R`、`source_geometry.R`、`sessionInfo.txt`：复现参数和环境。
- `connectors_only.svg`：四条独立引线及第五图标题，可单独编辑。

四种结构是代表性示意而非序列比例图：CR1保留3′重复，ERV的Env为可选，通用SINE不等同于人类Alu，DNA表示TIR型。V15着丝粒及assembly-discovery数据的所有原有限制不因本次拼接而改变。

## 独立目录复现验证

已将拼接输入和脚本复制到不含原始基因组数据的临时项目，从 `/tmp` 调用 shell 的 `--preview-only` 模式。新生成 SVG 与本目录最终 SVG 逐字节一致，详见 `provenance/portable_validation.json`。验证复用了已安装 R 依赖；新机器请先执行依赖安装脚本。Python 全量审计脚本还需 Git 源版本、Pillow 和 Poppler，仅作审计，不是 R 拼接运行依赖。
