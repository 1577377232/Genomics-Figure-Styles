# Figure 1 — V3_R3（紧凑结构栏，五条直线折线）

本版基于V3_R2，直接使用GitHub新增的五个 `composition_1_5x.svg`，保留V15原始Circos代码、数据、窗口、颜色和全部既有版本。根据本轮最后补充，结构栏再向右移一点、显著压缩纵向空隙，五个模块总高度小于Circos真实高度。

## 运行

在 `figures/fig1_chicken_GRCg7b_circos/` 下执行：

```bash
# 只在缺少R依赖时执行：
bash 06_scripts/install_composition_R_dependencies.sh
# 完整导出：
bash 06_scripts/run_compose_v15_with_te_architectures_v3_r3.sh
# 快速预览或另设新输出目录：
bash 06_scripts/run_compose_v15_with_te_architectures_v3_r3.sh --preview-only
bash 06_scripts/run_compose_v15_with_te_architectures_v3_r3.sh --out-dir /absolute/new_output
```

完整、独立R代码位于 `06_scripts/compose_v15_with_te_architectures_v3_r3.R`，不执行V15或旧拼接代码。运行依赖R/rsvg/png；无需原始VCF或全基因组分析。Python独立审计额外需要NumPy/Pillow及Poppler，不属于拼接运行依赖。

## 冻结素材与1.5倍放大

源素材提交为 `06e8112`。实际读取：

1. `assets/TE_architectures/LINE_CR1_composition_1_5x.svg`
2. `assets/TE_architectures/LTR_ERV_composition_1_5x.svg`
3. `assets/TE_architectures/SINE_generic_composition_1_5x.svg`
4. `assets/TE_architectures/DNA_TIR_composition_1_5x.svg`
5. `assets/TE_mediated_SV/TEV_classification_composition_1_5x.svg`

五个SVG显示宽度统一为7.40625in，是R2的4.9375in的1.5倍。保持原始1400×250画布比例，每个高约1.32254in；内部文字、框和边线同步放大1.5倍。素材标题已经从58改为39源单位，最终约14.855pt，相比R2的14.728pt仅大0.86%，因此视觉上基本保持原字号。没有二次修改这些SVG。

TEV保留标题、SV括线、TE珊瑚条、浅灰非TE区及 >80% 阈值；无内部箭头、无独立结果框。五个原生SVG内容直接嵌入，仅给ID/引用加前缀以避免冲突，结构文字仍可编辑。

## 最终结构栏位置（以最后补充为准）

从实际扇区角度、轨道外半径、文字的有向边界框计算Circos右沿，增加0.5pt描边余量，得到右沿8.24012in。结构栏左侧8.94012in，离右沿0.70in；本轮初始预览的0.45in已按最新要求加大。

结构栏内部的SVG画布间距从初始0.52746in压至0.08in；模块尺寸保持不变。总高6.93272in，对比Circos实际高7.64824in，为其90.64%，并按两者实际边界的垂直中心对齐。顶端依次约1.48506、2.88760、4.29015、5.69269、7.09524in，均由脚本计算，不再使用初始试验的0.80/2.65/4.50/6.35/8.20in。

页面宽度动态计算为8.9401169+7.40625+0.20=16.5463669in，高10in。V15仍固定左上位置0.15/0.15in和显示高度9.70in、缩放0.97；没有缩小左图来容纳结构。

几何细节见 `layout_geometry.tsv`、`architecture_placements.tsv`。

## 真实锚点、留白和三段折线

源轨道表与染色体映射均来自原 `07_figures/R_v15/`。五个真实锚点在Chr1/0bp的上方24°端面，内部轨道索引4/5/6/7/10分别是LINE/LTR/SINE/DNA/TIP。内部索引10对应最内侧TEV density，不是Chr7或其他外轨。TEV真锚点半径0.333，和R2一致，已由原SVG珊瑚填充多边形验证。

**按照用户已明确确认的方向，留白沿轨道端面法向进入右侧缺口0.14in。** 对本版上方24°端面，页面坐标中的单位法向为(sin24°,cos24°)，

```
true_anchor = circle_center + R * track_radius * (cos24°, -sin24°)
connector_start = true_anchor + 0.14in * (sin24°, cos24°)
```

真锚点不移动，正式绘线从偏移后的空白区开始。`connector_anchors.tsv` 同时记录真锚点、实际绘线起点、法向、偏移长度和源轨道。

每条线严格为三段直线，命令仅 `M L L L`；不含C/Q，不调用任何平滑函数。颜色#253246、线宽0.85pt、round linejoin仅影响描边接头。每条线独立，不共用主干、顶点或线段。

第二折点横坐标按LINE/LTR/SINE/DNA/TEV逐步增加，分别为Circos右沿加0.42/0.49/0.56/0.63/0.70in。第一折点及前三段的方向按源轨道和标签避让设置。SINE/DNA/TEV首段带有向下分量，以避开原图百分比标签并保持不交叉；最内TEV需要先经缺口向外走，不能直接斜穿外轨。这些段均为直线。

末段水平，小闭合箭头尖端距对应结构**可见描边外缘**0.09in，已扣除原生SVG边框的半线宽，而非只按矩形坐标估计。实际路径、四个顶点和线对最短距离分别在 `connector_paths.tsv`、`connector_vertices.tsv`、`connector_pair_QC.tsv` 中。

## DEBUG 与正式输出

正式 `.svg` / `.pdf` 不含DEBUG标记。DEBUG全图及起点放大图中，空心圈是真实端点，带编号实心圈是0.14in留白之后的绘线起点，编号5对应最内侧TEV；仅DEBUG显示两点间的虚线。

- 正式PNG：9928×6000，600dpi。
- Preview与DEBUG全图：2482×1500，150dpi。
- DEBUG起点放大：2130×1380，300dpi。
- SVG为16.54637×10in，保留全部325个V15裁剪路径，未嵌入栅格。
- PDF为矢量。当前rsvg会把PDF页尺寸取整数点，脚本在右侧只补0.6616pt白边，保证内容按原点坐标1:1导出，未缩放Circos；见 `PDF_geometry.tsv`。

## QC与数据口径

`protected_baseline.json` 校验V15、V2、旧V3、R2及所有已有素材；`input_preservation.tsv` 校验运行前后输入。独立Python直接解析实际SVG的M/L顶点，按线段解析几何计算线对距离；同时验证两个端部留白、模块间距、标题补偿、数据扇区/文字框避让及不穿其他SVG画布。`validation_report.json` 和 `provenance/validate_composition_v3_r3.py` 提供完整证据。

主图数据完全保留原来的TEV >50%、事件 >50bp extended口径；右侧分类示意沿用用户指定的 >80%。没有重新筛选或统计主图。此模式图是覆盖分类示意，不是分子形成机制证据。
