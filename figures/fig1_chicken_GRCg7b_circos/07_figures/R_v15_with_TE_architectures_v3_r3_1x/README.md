# Figure 1 — original-size modules, compact column (V3_R3_1x)

本版按最新要求恢复原版SVG的1倍显示大小，再缩小五部分纵向间距。保留V15与R3/R2等全部旧版；不执行V15绘图，也不重新计算数据。

## 运行

在 `figures/fig1_chicken_GRCg7b_circos/` 下执行：

```bash
# 缺依赖时才需要：
bash 06_scripts/install_composition_R_dependencies.sh
bash 06_scripts/run_compose_v15_with_te_architectures_v3_r3_1x.sh
# 只生成预览：
bash 06_scripts/run_compose_v15_with_te_architectures_v3_r3_1x.sh --preview-only
```

完整代码：`06_scripts/compose_v15_with_te_architectures_v3_r3_1x.R`。输入为已冻结的R3 SVG、R3几何表与下列5个原版SVG（不使用标题补偿版）：

- `assets/TE_architectures/LINE_CR1.svg`
- `assets/TE_architectures/LTR_ERV.svg`
- `assets/TE_architectures/SINE_generic.svg`
- `assets/TE_architectures/DNA_TIR.svg`
- `assets/TE_mediated_SV/TE_mediated_SV_classification_v3_r2.svg`

最后一个是已简化的原尺寸TEV图，保留SV括线、TE条、>80%阈值，无内部箭头或结果框。所有源SVG均未修改。

## 本轮调整

五个模块统一宽4.9375in（R3的7.40625in / 1.5），高0.881696in，保持1400×250原始比例；恢复未补偿的原版标题。SVG画布之间的空隙为0.08in，中心间距由1.402545in缩至0.961696in，减少约31.4%。整列高4.728482in，以实际Circos上下边界中点4.951418in居中；Circos自身高7.648243in。

五个结构的可见条带左端横坐标不动。Circos、页面宽高、颜色、源轨道、真实锚点、0.14in端面法向留白、引线前两个顶点与全部横坐标不变；只同步调整末两个顶点的纵坐标及箭头纵坐标。所有引线仍为三段M/L折线、0.85pt深灰蓝，箭头与结构可见边框保持0.09in间隙。TEV源仍为上方Chr1端面的最内TIP/TEV轨道。

SVG、PDF、600dpi PNG、Preview、DEBUG及起点放大图均已生成并查看。正式PNG为9928×6000。PDF仅补不足1pt右侧白边以适应本机rsvg尺寸取整，内容保持1:1。全部内容仍为可编辑矢量SVG；主图没有栅格化。

## QC

`QC_summary.json` 为独立检查结果；`connector_pair_QC.tsv` 记录10对引线距离，0交叉、0共用线段，最短间隔9.992263pt。独立校验原始V15/R3/源SVG等95个冻结文件不变，并逐字比较除5个结构块和对应引线末端几何外的所有SVG内容。`provenance/validate_original_size.py` 需Python/Pillow，可在完整仓库中重跑。

主图TEV数据仍为>50%、事件>50bp的既有extended统计；示意图保留用户指定>80%，两者均未重定义。原始基因组/VCF无需用于本次拼接。

便携包包含运行所需的冻结R3输入、5个原版SVG、R/shell及最终输出；完整95文件历史校验须在完整仓库运行。R依赖库不打入包。
