# 鸡 GRCg7b TE、SNP、SV 与重复序列多态性 Circos：数据准备与预览

本轮完成真实数据盘点、三种窗口统计、QC、可编辑矢量预览。所有新增文件位于本目录。未修改原始 FASTA、VCF、BAM 或注释，未启动全基因组比对、变异调用或 RepeatMasker 重跑。源文件大小/mtime 前后核对见 `05_QC/original_input_integrity.json`；SHA256 见 `00_data_inventory/input_fingerprints.json`。

## 本次主图使用的数据

| 外→内 | 数据与统计单位 | 主染色体范围的真实结果 |
|---|---|---|
| 1 LINE | 同类注释并集覆盖 bp / 窗口实际长度 | 见 `02_TE_annotation/class_summary.tsv` |
| 2 LTR | 同上 | 同上 |
| 3 SINE | 同上 | 同上 |
| 4 DNA | 同上；MITE 归 DNA | 同上 |
| 5 SNP | 独立 SNP 位点/Mb；未做 LD pruning | **9,741,807**，仅 1–39；Z/W 为 NA |
| 6 SV | PASS 独立非 BND SV 位点/Mb | Jasmine **123,901**；SURVIVOR 独立对照 **207,470** |
| 7 TE / repeat polymorphism | INS+DEL，事件>50 bp、重复命中覆盖严格>50% | 扩展 **66,804**；严格 TE **28,295** |

推荐首先检查 `07_figures/chicken_GRCg7b_1000kb_extended_SV-all.png`。对应 SVG 的文字与图形可编辑；PDF 字体嵌入，PNG 为600 dpi。另提供100 kb、500 kb预览，以及1 Mb严格TE版。图中的结构示意按任务指定的通用模板绘制，不代表鸡每个家族都具备所示结构（如双单体 SINE 或 Env）。

## 参考坐标及排除范围

- 实际参考：`/data/CWD_1/00_Ref_TERef/7b_genome.fa`；序列SHA256见指纹表，与原 SV/SNP 工作流记录的参考核对一致。41条主要染色体：1–39、Z、W，总长 **1,041,122,857 bp**。完整FASTA有214条序列，总长 **1,053,332,251 bp**。
- 名称保持FASTA中的 `1`、`2`、…、`39`、`Z`、`W`；显示标签可加Chr。`01_reference/chromosomes.tsv`记录全部序列；MT 16,784 bp 与172条其他contig不纳入主图。后者总长 **12,192,610 bp**。
- 主染色体顺序固定，所有窗口从0开始、BED为0-based half-open，终点裁到真实染色体长度。VCF POS 仍保留1-based字段。三种窗口完全使用同一参考；不同分辨率的总覆盖bp和位点数通过守恒检查。
- `05_QC/reference_dictionary_checks.json`核对各VCF声明的contig长度，SNP/原始组装BAM的参考溯源文件保存在 `00_data_inventory/source_methods/`。

## TE注释与拷贝口径

使用 `/data/CWD_1/Pan_TE_chicken/02_Genome_annotation/GRCg7b/RepeatMasker_out/GRCg7b.fa.out`，RepeatMasker4.1.5 / RMBlast2.14.1，原命令含 `-poly -xsmall -lib ... -engine rmblast`。该注释用的是 `7b_genome.fa.mod.EDTA.TElib.fa`，不能仅因项目名把它称为之后构建的Pan-TE库。无效目录 `RepeatMasker_out.invalid_GRCg6a_copy_20260919` 已排除；根目录旧27样本汇总不作为本次输入。

分类规则：LINE/*→LINE、LTR/*→LTR、SINE/*→SINE、DNA/*及MITE→DNA；RC/Helitron、Unknown/unknown、Simple_repeat、Low_complexity、Satellite及其他非TE类别分别保留。原始分类到统一分类的完整映射见 `classification_mapping.tsv`。

覆盖采用同一类别区间并集；跨类别重叠各自保留，因此四类覆盖不可简单相加解释为全TE无重复覆盖。拷贝数采用 `(chrom, RM_ID, repeat_name, raw_class)` 唯一组合，将一个RepeatMasker ID的碎片视作操作性拷贝，以该拷贝跨度中点计入唯一窗口。它不等价于经完整结构验证的生物学完整拷贝。输出还保留原始片段数、每类union BED、copy BED和统一注释BED。

## SNP的筛选、完整性和限制

来源 `/data/CWD_1/Pan_TE_chicken/02_Genome_annotation/15_Homologous_SNP_phylogeny/06_merged/assembly52.union.snps.vcf.gz`。母表 **19,588,797** 条位点记录。52列是51个其他组装+参考GRCg7b；一些组装是同一个体的不同单倍型，不能写为52个独立个体。原始SNP由已存在的asm5 BAM/cs标签提取，过滤包括MAPQ、比对长度、identity、竞争比对、indel附近及末端屏蔽，详见复制的原始方法README。

本次重新从母VCF筛选：`FILTER=PASS`、ACGT SNP、严格双等位、`NS/52>=0.95`、`min(AC,AN-AC)/AN>=0.01`，不做LD pruning。独立统计与现有 `call95.ld02/qc.bim` 的全部有序位置SHA256一致，逐位点REF与GRCg7b核对；重复位置数与REF不一致数均为0。`dataset.vcf.gz`是LD剪枝后系统发育输入，未用于密度。

现有调用限定常染色体。Z/W的SNP count/density为 **NA**，图中灰色；不能将其当0。补齐需要对已有BAM增加性染色体callable区域和SNP提取，并处理性别/倍性；本轮未启动。这里的PASS来自组装比对过滤，QUAL缺失，不能描述为有测序深度/GQ证明的重测序SNP。单倍体GT也不意味着所有组装均已完整相位。

缺失率直接从未剪枝 `qc.bed`解码，见 `05_QC/SNP_unpruned_sample_missingness.tsv`。位点call rate>=95%不保证每个样本缺失率<=5%，也不消除组装完整度差异对密度的影响。

进一步复用了51个已有callable BED，加参考的ACGT区域，逐碱基统计至少50/52个组装可判定的区域（与95%位点call rate一致）。**全部9,741,807个保留SNP均落在该区域内**。常染色体中完全没有达到该阈值可判定区域的是 **16, 29, 30, 31, 32, 33, 35, 36, 37, 38, 39**；这些染色体以及窗口内共同可判定bp=0的情况，主SNP count/density为NA、图示灰色，不能说没有变异。`SNP_retained_count`另保留原始筛选通过数（这些窗口为0），联合callable BED、每窗口callable bp/比例和平均样本可判定比例均已输出。仍按物理Mb显示密度，未用极小callable分母放大数值；partial-callable窗口应结合比例列解释。

## SV面板选择及去重

新Jasmine母VCF168,236条，SURVIVOR母VCF256,639条，均为同一批51个组装、GRCg7b参考，但合并规则不同。本轮分别检查、分别统计；主图选新Jasmine，因为本次全量INS/DEL注释与它逐ID对应。没有将两套面板或旧TE筛选面板合并。

原合并命令及源报告已存入 `source_methods/`。Jasmine命令带 `--use_edit_dist --keep_var_ids --output_genotypes`；原审计指出默认 `min_seq_id=0` 未实际施加序列相似度阈值。SURVIVOR旧版本存在跨类型合并，主图选择不表示已证明Jasmine每个簇都更真实。

密度过滤：严格PASS、主染色体、51–100,000 bp；BND单列排除。INS用VCF断点基座 `POS-1`；DEL及其他区间事件去除左锚定基座后，以半开区间 `[POS,END)` 的整数中点 `(POS+END-1)//2` 计数。跨窗口事件只计一次。相同chrom/POS/type/END/REF/ALT的精确重复只保留第一条；相邻或同位置的不同等位SV不被强行合并。这里的“独立位点”为合并目录操作性位点，不能声称所有邻近多等位事件都已归一为同一个生物学位点。

Jasmine主范围精确重复移除 **5,659** 条，别名与原SUPP_VEC见 `05_QC/Jasmine.exact_duplicate_aliases.tsv`。BND端点、排除数、SV类型和FILTER分布见 `05_QC/sv_panels.json`。源VCF与所有主范围记录均保留原文件/派生表，未改变原始合并记录。

原51列GT是组装发现调用，未发现的记录通常为`./.`，不能当作已证实参考型；单列ref/1|0的图构建文件是目录编码，不能用来计算群体AF或缺失率。派生TIP VCF保留真实51列发现调用，GT缺失率表仅反映被选代表记录的发现矩阵，不是群体重测序分型质量。重复位点的发现支持未重新合并，不推断频率。

## TE / repeat多态性的重新注释

复用新Jasmine全部 **154,228** 个INS/DEL查询序列的原始RepeatMasker命中，没有仅从旧80%候选中下调阈值。原GraffiTE脚本提交 `f37674bd5fcd9e92e08a470911813a1c5807d084` 与参数记录见检测配置；本轮以独立脚本替换筛选逻辑，无需修改原GraffiTE main或重跑比对。

1. 核对每批VCF记录、查询FASTA数与序列长度；全部原REF与GRCg7b一致。
2. INS以ALT、DEL以REF的原查询命中为基础，移除左端1 bp VCF锚定碱基。以实际事件长度作为分母。
3. 使用**原始RepeatMasker真实命中区间并集**，重叠不重复计算，OneCode跨片段间隙不算作匹配碱基。这与旧 `total_match_span` 的部分跨度统计有差别，80%回算数量不应期待完全相等。
4. 严格TE包含LINE/LTR/SINE/DNA/MITE/RC/Helitron；扩展集合加Simple_repeat、Low_complexity、Satellite。Unknown不自动当作TE或允许的其他重复，另保留all_RM覆盖列。要求**事件长度>50 bp且允许类别并集覆盖>50%**，不额外要求匹配碱基绝对长度>50；如需后一条件，可用 `expanded_match_bp_gt50` 再筛。原调用上限100 kb仍存在。

全部contig的严格TE为 **35,163**，扩展为 **87,349**；主染色体/PASS/精确去重后为上表数值。严格主集合 INS **13,602**、DEL **14,693**；扩展INS **36,904**、DEL **29,900**。DEL表示相对参考的TE/重复缺失，不能称为全部都是非参考新插入。

每条记录完整覆盖比例、主类别、事件长度和筛选标记见 `INS_DEL.repeat_reclassified.tsv.gz`；每类统计见 `TIP_class_counts.tsv`。短片段可能来自非完整TE，采用宽松50%阈值的目录是候选多态性集。重复序列注释只能支持序列组成，不能单独证明事件的转座发生机制。

SV总密度与TIP本来存在包含关系，不应解释为互斥。另提供按同一50%口径未达到TE阈值的已检测INS/DEL **94,873**，未达到扩展重复阈值 **56,364**，以及尚未按此方法检测的DUP/INV **733**。列名 `SV_nonTE_assessed` 只表示“已检测且TE覆盖<=50%”，不证明它们完全无TE来源；不足50%的TE片段仍可能存在。未检测的SV不能归为非TE。

## 窗口、单位与尺度

`04_window_statistics/windows_100kb.tsv`、`windows_500kb.tsv`、`windows_1000kb.tsv`具有相同列定义、同一染色体顺序与各分辨率窗口边界；每类有covered_bp、coverage_fraction、copy_count、copies_per_Mb。后三层有count和per_Mb，并附SURVIVOR、严格/扩展TIP、INS/DEL子集和SV补集。

`effective_window_bp=End-Start`是**该窗口实际物理长度**，末端短窗不按整1 Mb除；同时记录nonACGT_bp和ACGT_bp供检查。主密度单位为位点/实际Mb，主覆盖率除以实际物理窗口长度。没有将N碱基排除后的密度称为群体callable校正密度；真正的群体可判定区域随样本变化，不能只靠参考非N长度替代。N比例高的小染色体需结合QC解读。

1 Mb适合主图概览；100 kb保留更细节但噪声与空窗更多，500 kb为折中。三个分辨率的空窗比例、median/p95/max见 `05_QC/window_size_comparison.tsv`。小于1 Mb的染色体仍按真实长度绘制，通常只有一个1 Mb窗口，因此另提供更细窗口版本，未人为放大染色体弧长。

实际SINE空窗比例从100 kb的68.0%降为1 Mb的10.6%，严格TIP从48.5%降为1.41%，因此1 Mb更适合七轨主图。`SNP_per_callable_Mb`提供按95%共同可判定长度归一的补充密度；当callable长度很小时它可能不稳定，不与物理Mb密度混为同一单位。QC诊断图 `05_QC/SNP_callability_and_resolution.*`展示共同可判定比例与窗口分辨率的影响。

前四层采用相同的线性面积轨道，显示各自参考覆盖率，按类别单独定上限；后三层为独立 `log1p(位点/Mb)` 连续色阶。绘图上限取各轨道99分位并截断，**原始统计表不截断**，每轨道实际上限及被截断窗口数在对应 `.settings.json` 中。不同轨道颜色深浅与高度不可跨量纲直接比较。

## 着丝粒：公开数据和现有比对的检查

检查了UCSC GRCg7b GenArk (GCF_016699485.2) 官方轨道API，结果已完整保存：cytoBandIdeo **214**条、acen **0**条；这些只是全长gneg染色体方向示意，不能用作着丝粒。AGP gap轨道并不提供可直接确认的着丝粒注释。本轮没有找到可直接使用、覆盖41条染色体且可核验的公开GRCg7b着丝粒坐标；不等于宣称网上绝不存在。相关GGswu文献与GRCg7b局部CENP-A/Capture-C资料列在文献说明中。

您给出的35条GGswu区间已原样保留，缺15、27、33、37、38、39。输入0和1起点混用，来源坐标制及版本尚待确认。现有BAM对应GCA_024206055.2；Chr25实际长度4,095,597，所给终点4,096,356越界759 bp，未裁剪。部分小染色体编号之间也不对应，例如GGswu Chr29的主映射可能落到7b Chr33。

`transfer_centromeres.py`复用已有BAM的CIGAR，先从原FASTA重建breakgaps段编号（后缀是片段编号而不是原始偏移），只投射MAPQ>=30的非secondary M/=/X块，正确处理反链、硬剪切及supplementary。没有用包围跨度填补未比对区域。无query长度不一致记录。

得到 **6**条高覆盖序列同源性候选（W、Z、1、2、3、5），标准为源覆盖>=80%、主目标/方向占比>=95%、重复映射比<=1.05、目标跨度/源长度0.25–4。详细块、覆盖率及歧义见 `01_reference/centromere_transfer_qc.tsv`。**来源版本未确认，候选不等价于已证明的GRCg7b活性着丝粒，默认图中不绘制。**下一步应核实输入来源版本，再用一致版本区间转移，结合两侧唯一序列共线性与CENP-A验证；低覆盖区域不能自动补全。

## 文献和可复用代码核查

人类论文官方代码声明指向 [Zenodo 19220934](https://doi.org/10.5281/zenodo.19220934)，本次API解析到记录19220935。7个官方分析脚本和元数据已保存于 `00_data_inventory/literature/`；未找到其中的Circos绘图脚本、参数文件或结构示意源码，也未找到可确认的官方GitHub绘图仓库，因此没有声称复用了原绘图代码。Fig.6参考图已保存供本项目核查。其展示的是**癌症中表达TE**，本图四层是**参考组装全部已注释TE的覆盖**，生物学含义不同。[原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13666008/)

稻论文代码包 [Zenodo 16625002](https://doi.org/10.5281/zenodo.16625002) 的小型rice_pangenome.zip已保存，包含流程shell脚本；没有下载不需要的数GB组装。本图使用Matplotlib重新实现开放环、四种同色结构示意和引导线，保留脚本及矢量文字。[原文](https://www.nature.com/articles/s41467-025-62887-4)

## 目录和复现

- `00_data_inventory/`：2734项文件清单、644项VCF头检查、输入指纹、面板核查、样本/方法快照、官方文献代码。
- `01_reference/`：完整染色体长度和主图顺序、非ACGT窗口统计、GGswu原区间及探索性转移。
- `02_TE_annotation/`：统一分类、同类union、操作性copy BED与汇总。
- `03_variants/`：SNP位置表、两套SV逐位点表、全部INS/DEL重筛表及51列TIP VCF。
- `04_window_statistics/`：三种分辨率的真实数据宽表与窗口BED。
- `05_QC/`：参考字典、缺失率、各染色体总量、守恒/分区检查、来源文件完整性、阈值和着丝粒状态。
- `06_scripts/`：本轮完整统计、转移、面板核查、绘图和报告脚本。
- `07_figures/`：预览稿及可编辑SVG、PDF、600dpi PNG。该轮图可评审，不宣称已补齐性染色体SNP或已验证全部着丝粒。

现有Python环境：`/home/danding/softwares/svim-asm-1.0.3-20261103/bin/python`，依赖版本在 `software_versions.json` 和 `requirements.txt`。无需pandas、R或Circos软件。

```bash
# 在项目工作目录执行：完整重建，只写本独立目录
bash chicken_7b_TE_circos/06_scripts/run_all.sh

# 仅重画1 Mb扩展重复集合主图
/home/danding/softwares/svim-asm-1.0.3-20261103/bin/python \
  chicken_7b_TE_circos/06_scripts/plot_circos.py --window 1000000 --tip extended --dpi 600

# 严格TE版本，或以已检测非TE阈值集合替换SV层进行比较
/home/danding/softwares/svim-asm-1.0.3-20261103/bin/python \
  chicken_7b_TE_circos/06_scripts/plot_circos.py --window 1000000 --tip strict --sv nonTE
```

本轮不以模拟数据补足任何轨道。剩余实质性缺口是Z/W的SNP、GGswu着丝粒来源版本及低覆盖区的可靠坐标证据；另需注意宽松重复候选、组装发现GT和100 kb上限这些实际限制。
