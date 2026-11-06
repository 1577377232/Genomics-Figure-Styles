# 52 个基因组 TE 注释核对（2026-11-04）

后续选择（2026-11-04）：用户确认沿用 RIR、TB 旧注释，并选择以旧 WD 注释代表 Wuding；当前完整的 GCA_030914275.2 保留。已按该选择整理为 [14_TE_annotation_52](../14_TE_annotation_52/README.md)，包括 52 个样本的结果链接、来源清单和重新生成的统计。本文件下方保留首次核对的差异记录；这些差异记录不代表当前结果选择。

范围：以 /data/CWD_1/Pan_TE_chicken/01_Pan_TE_Lib/Construct_Pan_TE/Assemblies_51/raw_assemblies_51.tsv 中的 51 个原始组装为准，另加 GRCg7b 参考。逐样本路径见 [status_52.tsv](status_52.tsv)，序列核验原始记录见 [sequence_checks.txt](sequence_checks.txt)。

## 结论

以前确实进行了整基因组 RepeatMasker 注释；后续补跑的 31 个任务全部成功结束，最后一个 Ross_2 的完成时间为 2026-09-22 18:42:04。

按当前指定的 FASTA 核对，不能把所有 52 个输入都当作版本一致的已完成结果：

| 类别 | 数量 | 判断 |
|---|---:|---|
| 后续重注释，输入链接与当前清单一致 | 31 | 四项核心结果齐全，批处理状态成功；含 GRCg7b |
| 旧注释，当前输入的序列数、总长度及首段序列匹配 | 16 | 找到可对应的完整旧结果 |
| 旧注释完整覆盖当前输入，但含额外序列 | 2 | GRCg6a、Cornish_1；复用时需过滤额外序列及核对 contig 名称 |
| 旧注释中部分实际序列与当前输入不同 | 2 | Rhode_Island_Red_1、Tibetan2 |
| 旧结果对应另一组装 | 1 | Wuding |

因此，有 49 个输入找到了可对应的完整结果（其中 2 个含额外序列）；其余 3 个输入不能直接沿用现有结果。这里的“可对应”使用上述检查深度，不表示其余 47 个均已做全基因组逐碱基 MD5。

## 需要处理的 3 个输入

| 样本 | 当前 FASTA | 旧结果 | 问题 |
|---|---|---|---|
| Wuding | GCA_030914275.2，127 条，1,025,968,684 bp | WD.fasta，287 条，1,121,566,479 bp | 明确是不同组装；本目录未找到当前组装的匹配注释，需按当前输入重注释 |
| Rhode_Island_Red_1 | 1,339 条，1,023,494,120 bp | RIR.fasta，1,339 条，1,023,494,118 bp | Chr2 与 Scaffold1 对应序列的 MD5 不同；各净少 1 bp，不能归结为简单的单碱基删除 |
| Tibetan2 | 1,342 条，1,031,193,845 bp | TB.fasta，1,342 条，1,031,193,843 bp | Chr2 与 ChrW 对应序列的 MD5 不同；各净少 1 bp，不能归结为简单的单碱基删除 |

后两个样本的差异已通过实际序列确认，不是 .tbl 的长度统计误差。按当前 FASTA 统一分析前，需要重新注释，或验证差异并完成可靠的坐标转换；本次未重跑或修改原注释。

## 完整覆盖但含额外序列的 2 份旧结果

- GRCg6a：当前 FASTA 的 463 条序列全部在旧 .masked 中找到大小写归一后的逐序列 MD5 匹配；旧结果另含 MT，16,775 bp。
- Cornish_1：当前 FASTA 的 1,001 条序列全部匹配；旧结果另含 Scaffold420（7,575 bp）及 Scaffold803（20,074 bp）。

这些样本的当前输入序列没有漏注释。旧结果的 contig 标识可能经过重命名，后续与原始组装坐标联用时需恢复对应关系。

## GRCg7b

有效结果位于 ../GRCg7b/RepeatMasker_out/GRCg7b.fa.*，2026-09-20 00:21:06 完成。输入为 /data/CWD_1/00_Ref_TERef/7b_genome.fa，与 Assemblies_51 的参考是同一文件；214 条序列、1,053,332,251 bp，匹配参考 .fai。

历史上错误复制的 GRCg6a 结果已保存在 ../GRCg7b/RepeatMasker_out.invalid_GRCg6a_copy_20260919，本次未将其计入有效 GRCg7b 注释。

## 方法及旧文件状态

- 软件：RepeatMasker 4.1.5、RMBlast 2.14.1+。
- 使用的自定义库：/data/CWD_1/00_Ref_TERef/7b_genome.fa.mod.EDTA.TElib.fa。库身份来自运行脚本及日志，没有依据当前项目目录名将其视为后来构建的 Pan-TE 库。
- 完整输出检查：同一前缀的 .tbl、.out、.out.gff、.masked 四项均存在且非空。
- 31 个补跑样本：检查当前输入链接、当前样本清单与成功批处理状态一致。
- 21 份沿用旧结果：实际读取当前原始 FASTA，统计序列数及总长度；对 20 个非 GRCg6a 的旧结果比较首条序列的前 8,000 bp。Wuding 首段也不匹配。
- 对 GRCg6a、Cornish_1、Rhode_Island_Red_1、Tibetan2 另做完整逐序列 MD5 比较，忽略大小写、contig 名及顺序；差异证据见附表。

2026-09-16 的旧 audit 报告早于补跑任务，不能用于判断当前完成情况。根目录 RepeatMasker_tbl_summary.csv、RepeatMasker_out_class_family_wide.csv 以及 05_TE_annotation_plots_v2 的主分类统计仍只含 27 个旧结果，尚未更新为当前 52 输入的统一汇总。

本次仅新增核对表和核对说明，保留全部原始输入、旧报告、注释结果和汇总文件。
