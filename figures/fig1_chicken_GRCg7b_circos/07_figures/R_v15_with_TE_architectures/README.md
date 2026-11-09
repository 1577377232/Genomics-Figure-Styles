# Figure 1：V15 Circos + 四类 TE 结构（独立 R 拼接版）

本目录只做矢量拼接，不执行 V15 主图脚本，不重新读取或统计 SNP/SV/TIP。原 V15 代码、SVG/PDF/PNG、统计数据及中心图例均保留。

## 运行

在 `chicken_7b_TE_circos/` 下执行：

```bash
# 首次使用新环境时安装独立依赖；不会修改 V15 的 R_library：
bash 06_scripts/install_composition_R_dependencies.sh
# 完整 SVG / PDF / 600 dpi PNG / 小预览：
bash 06_scripts/run_compose_v15_with_te_architectures.sh
# 也可以直接运行 R：
Rscript 06_scripts/compose_v15_with_te_architectures.R
# 只生成 SVG 与小预览，便于调整排版：
bash 06_scripts/run_compose_v15_with_te_architectures.sh --preview-only
# 可另选输出目录，禁止写入原 R_v15：
bash 06_scripts/run_compose_v15_with_te_architectures.sh --out-dir /absolute/new_output
```

新环境需要 R、C/C++ 编译器和 librsvg>=2.52、libpng、fontconfig、freetype 系统开发库。安装脚本固定了本机验证版本：R4.1.2；rsvg2.7.0、svglite2.1.3、png0.1-9，以及cpp11/systemfonts依赖，均安装到 `06_scripts/composition_R_library/`。

## 输入

1. `07_figures/R_v15/chicken_GRCg7b_R_v15_1000kb_adaptive_discovery_TIP-extended_SV-all_palette-pangenome_spectrum.svg`。
2. 同名的 `center_legend.tsv`、`figure_geometry.tsv`、`chromosome_display_mapping.tsv`：只用于定位圆心、缩放和右下开口锚点，不改变图形数据。
3. `assets/TE_architectures/LINE_CR1.svg`、`LTR_ERV.svg`、`SINE_generic.svg`、`DNA_TIR.svg`。来自 GitHub `1577377232/Genomics-Figure-Styles` 的 e9072dd 提交；路径及逐文件SHA256见该目录 `source_manifest.json`。原始结构文件按字节保留。

## 实现与排版

- 左侧完整嵌入 V15 的原 SVG 内容，保持等比例缩放；原有全部裁剪路径、七轨、着丝粒及 SNV density / SV density / TEV density 图例不变。
- 右侧依次 LINE、LTR、SINE、DNA，使用仓库原始结构图及配色。结构 SVG 的文字保留为 `<text>` 元素，可以继续编辑。V15 本身的字体轮廓属性沿用原 SVG。
- 黑色短总引导线从右下开口、邻近 TE 轨道的空白处引出，再扇出四条折线箭头；仅四个箭头，没有第五类 TE。连接线是类型说明，不对应特定位点或序列长度。
- 箭头、折点和箭头头部通过 grid/svglite 产生独立矢量图层。SVG 元素 ID 加入各图层前缀以防重名；没有栅格嵌入。
- 没有直接用 grImport2 重绘 V15：该包的默认 `ext="none"` 不保留复杂裁剪，`clipbbox` 只使用边界矩形。原主图依赖真实扇区裁剪，完整 SVG 嵌入可以同时保留裁剪和结构文字。参见 [grImport2 手册](https://cran.r-project.org/web/packages/grImport2/grImport2.pdf)；格式导出使用 [rsvg](https://docs.ropensci.org/rsvg/reference/rsvg.html)。
- 全部矢量/栅格输出使用同一 SVG 场景，物理画布14×10英寸。PDF1008×720pt；PNG8400×6000、600dpi；预览2100×1500、150dpi。示例中4200×3000@600dpi对应7×5英寸，本实现已统一尺寸。

## 最常修改的入口

只修改 `compose_v15_with_te_architectures.R` 顶部的 `composition_config()`：

- `main`：主图位置和高度，宽度按源 viewBox 自动推导，禁止独立拉伸。
- `architectures`：结构图的 `left_in`、`top_in`、`width_in`。
- `connectors`：总引导线、四条折线的折点、箭头长度和线宽。

排版参数的单位为英寸、原点在左上。生物坐标不参与排版修改。代码检查连接线不得穿越主数据环；独立 QC 进一步检查编号/物理坐标/纵轴文字框。

## 科学含义

这些是代表性结构、非序列比例图：CR1的3′重复不替换成人类L1的通用Poly(A)；ERV的Env可选；SINE示意为PolIII内部启动子例子而非Alu双单体；DNA为TIR型，不将Helitron/RC混入。详见资产README。原V15的TEV仍为包含INS/DEL和其他重复的extended面板；近似着丝粒的限制不因拼接而改变。

## QC 与交付

- `validation_report.json`：独立内容与格式检查。
- `input_preservation.tsv`：原V15、实际结构输入与脚本的MD5；重绘前后不变。
- `circos_placement.tsv`、`architecture_placements.tsv`、`connector_vertices.tsv`：精确排版与折点。
- `composition_settings.R`、`source_geometry.R`、`sessionInfo.txt`：可复现参数与环境。
- `connectors_only.svg`：可单独编辑的箭头层。
- `output_manifest.tsv`：绘图生成时的输出清单，不包含自身和运行日志；后续独立QC另有记录。

检查包括原V15所有SVG子元素逐项相同、四个结构图所有子元素相同、裁剪保留、ID引用有效、无连线遮盖数据/注释、PDF/SVG无嵌入位图、PDF页尺寸与PNG真实DPI。所有原文件保留。
