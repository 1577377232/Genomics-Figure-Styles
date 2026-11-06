# V11：蓝紫色 Circos 与 Z/W SNP 补算

基于 V10 精修，右侧开口 64°、Chr1 右上开始逆时针、1–39/Z/W 顺序、1.1° 染色体间缝隙与分段框结构保持。四类参考 TE 与 SNP/SV/TIP 多态性仍是不同含义的轨道。保留旧版和原始文件；仅在 V11 窗口表中补充真实 Z/W adaptive-discovery SNP。

## 直接运行

在项目根目录 `chicken_7b_TE_circos/`：

```bash
bash 06_scripts/run_circos_R_v11.sh

# 检查 100 kb / 500 kb / 1 Mb 的输入，不出图
bash 06_scripts/run_circos_R_v11.sh --all-windows --validate-only

# 修改美化后另存，建议每次使用新目录
bash 06_scripts/run_circos_R_v11.sh --out-dir 07_figures/R_v11_modified

# 可选：绘制其他窗口大小
bash 06_scripts/run_circos_R_v11.sh --window 500000 --out-dir 07_figures/R_v11_500kb
```

绘图入口为 `06_scripts/plot_chicken_circos_v11.R`，完整源码可直接交给 ChatGPT 修改。需要 R、circlize 及其依赖；新环境可先运行 `bash 06_scripts/install_circos_R_dependencies.sh`。本机验证环境 R 4.1.2、circlize 0.4.18。绘图不依赖 Python、BAM、原始 FASTA 或重新进行 SNP 提取。

## 输入与修改入口

- `04_window_statistics/v11/windows_1000kb.tsv`：默认真实统计输入；另有 100 kb、500 kb。新增字段 `SNP_v11_ZW_callable_bp` 为 Z/W 窗口内至少一个非参考组装可调用的并集长度。其他染色体该辅助字段为 NA，不改变其原有 SNP 数据。
- `01_reference/main.chrom.sizes`：GRCg7b 41 条主要染色体真实长度。
- `01_reference/v2/GRCg7b.centromeres.approximate.tsv`：沿用已核验的 17 个转移位置；本轮仅改变画法。
- `default_config()`：几何尺寸、字体、TE 分位数与平滑参数。
- `get_palette_colors()`：主方案颜色映射。LTR `#69A7AF` 是仓库原青绿 `#7BB6AC` 与青蓝 `#5899B2` 的等权 RGB 混合；DNA `#8F8FAE` 来自仓库 botanical_genome 的灰紫色。其余采用 pangenome_spectrum；边框 `#253246`、文字 `#273746`。
- `coordinate_ticks()`：1–9/Z 主刻度 50 Mb，10–39/W 主刻度 10 Mb；标签只有数字，统一单位在副标题中说明。小染色体数字加粗。为避免重叠，允许省略个别坐标文字，刻度和数据始终保留。
- `draw_ideogram_v11()`：17 个标记全部为 `#E75A46`，y=0–1，填满 ideogram 带宽。真实转移 evidence 仍在输入和标记表；展示角宽不代表着丝粒真实长度。
- `draw_grid_center()`：浅蓝紫圆底、短线色标、Outer/Inner 双向箭头。

## 独立 TE 纵轴与平滑

| 轨道 | 分位数 | 1 Mb 显示上限 | 原始最大覆盖率 | 原始超限窗口 |
| --- | ---: | ---: | ---: | ---: |
| LINE | q99 | 32% | 71.148463% | 7 |
| LTR | q95 | 8.2% | 69.758300% | 53 |
| SINE | q99 | 0.23% | 0.322200% | 10 |
| DNA | q99 | 15% | 24.495012% | 11 |

各类原始 coverage fraction 分位数向上取两位有效数字；每条轨道分别标注百分比轴。同样高度不代表同样的覆盖率。LTR q99=41.965393%，q95=8.116210%；极少数高覆盖窗口压低了旧版多数区域的可见变化，因此仅将 LTR 显示上限改成 q95。原始表与原始最大值全部保留；53/1063 个原始 LTR 窗口超限，47 个平滑节点超限。100 kb、500 kb 按各自窗口表重新确定上限。

沿用 V10 的展示平滑：每条染色体内进行长度加权的三窗口三角平均，再使用不超出局部节点范围的三次 Hermite 曲线。少于三个窗口时不作滚动平均。统计仍以原始窗口为准，不跨染色体连接。`TE_profile_nodes.tsv`、`TE_profile_curves.tsv` 和 `TE_smoothing_QC.pdf` 同时保留原始值及显示曲线；截顶发生在最终绘图阶段。

SNP 仍为线性、上限 50,000/Mb（图例显示 SNP/kb）；SV/TIP 仍为原来的 log1p、原 q99 上限，位点总数保持 123,901 / 66,804。TIP extended 包含 INS/DEL 与简单重复、低复杂度、卫星序列，不应描述成纯 TE-INS 集合。SNP 总数由 19,572,468 增至 20,552,376，仅来自新补入的 Z/W。

## Z/W 补算来源与方法

此前 Z/W 为 NA 是因为原 `extract_assembly_snps.cpp` 仅循环 1–39。已有 51 份组装到 GRCg7b 的带索引 BAM 中有 Z/W 比对，本轮没有重新比对。

补算脚本 `06_scripts/recover_sexchr_snps_v11.py` 复用原 C++，只将染色体范围改为 Z/W，过滤逻辑不变：竞争比对 MAPQ≥5、长度≥10 kb；可调用比对 MAPQ≥20、长度≥50 kb、一致性≥95%；参考覆盖必须恰好一个竞争比对；末端裁去 100 bp；indel 两侧屏蔽 5 bp；仅 A/C/G/T SNP。与原程序相同，排除 unmapped、secondary、QC-fail、duplicate；保留 supplementary 参与上述规则。

在 51 份结果间按参考位置合并，排除跨组装出现多种 ALT 的位点，输出 PASS 双等位 sites-only VCF。密度为独立位点数/实际窗口长度×1 Mb；窗口存在非参考组装 callable bp 才能输出计数，否则为 NA。GRCg7b 提供 REF 等位；51 是组装条目数，不等于 51 个独立、性别均衡的个体。**此补充是探索性 assembly discovery，没有联合基因型、性别分层 call rate 或群体等位基因频率。** W 同源/误组装序列仍可能影响结果，不能把有 W 比对的条目直接判为雌性。

| 染色体 | 长度 | 有 callable bp 的组装 | callable 并集 bp | 保留 SNP | 排除多等位位置 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Z | 86,044,486 | 51 | 84,005,525 | 955,908 | 4,925 |
| W | 9,109,940 | 40 | 8,696,547 | 24,000 | 54 |

1 Mb 的 Z 87/W 10 个窗口均有观测。100 kb 仍保留 Z 11/W 2 个缺失窗口，500 kb 保留 Z 1 个缺失窗口；没有为了显示而填零。常染色体、TE/SV/TIP 及其他 SNP 面板保持旧表值不变。只有默认 adaptive_discovery 面板纳入该补充；population95/discovery50 仍保持原始定义与缺失。

可追溯结果在 `03_variants/v11/sexchr_snps/`：sites VCF 及索引、callable union BED 及索引、逐样本提取结果、BAM 清单、逐染色体和逐窗口 QC、`method.json`。仅重新绘图不需要重算；如需在原工作区重现提取：

```bash
/home/danding/softwares/svim-asm-1.0.3-20261103/bin/python 06_scripts/recover_sexchr_snps_v11.py
/home/danding/softwares/svim-asm-1.0.3-20261103/bin/python 06_scripts/validate_sexchr_snps_v11.py
```

提取需要现有样本清单/BAM/FASTA、pysam、g++/HTSlib，路径集中在 Python 文件顶部；无需重新比对。交接包包含合并位点及 QC，未复制巨大的原始 BAM。

## 输出及核查

默认输出于 `07_figures/R_v11/`：PNG 600 dpi（5880×6000）、160 dpi 预览、PDF/SVG、参数、输入哈希、独立尺度表、颜色映射、坐标及标签计划、328 个 sector×track 面板记录。PDF/SVG 无嵌入位图，PDF 字体已嵌入。保护 V10 目录与原始表，不允许 V11 渲染写入旧 R_v3–R_v10。

`validation_report.json` 检查三种窗口数据一致性、保留原值、布局映射、标签碰撞、独立纵轴、平滑边界、着丝粒、位图 DPI 与实际 SVG/PDF 栅格化结果。`independent_CIGAR_validation.json` 使用 BAM CIGAR 和 query 碱基独立核对均匀抽取的真实 SNP，没有重复使用 cs 解析。

GitHub 风格来源为 `1577377232/Genomics-Figure-Styles`，本轮读取 commit `13f1fbb7a8d67da1c335934ce8ca839a5c6be6de`。原始风格文件快照保存在 `provenance/style/`。这里适配的是 Circos 分类轨道，不沿用 Manhattan 图的染色体配色语义；用户明确要求的浅色圆底、密度渐变、分位数显示截顶和微染色体展开优先于仓库针对 Manhattan 的默认条款。

建议图注注明：四类 TE 独立 coverage 纵轴；LTR q95、其他 q99 显示截顶；三窗口平滑仅用于展示；后三层为独立位点密度；SNP 为 assembly discovery、含本轮 Z/W 补充；微染色体展示角度放大但染色体内部坐标为真实 bp；着丝粒为转移坐标的示意标记。
