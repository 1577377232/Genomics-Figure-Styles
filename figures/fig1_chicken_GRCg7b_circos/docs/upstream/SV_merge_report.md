**51 个鸡个体 SV 合并比较：SURVIVOR 与指定 Jasmine 命令**

两种方法使用同一批新调用的 51 个 VCF。合并前保留 `50 < abs(SVLEN) <= 100000`，也就是 51–100,000 bp；BND 没有统一事件长度，全部保留并单列统计。原始 821,869 条调用中，排除 53 条 >100 kb 记录，进入合并共 821,816 条。保留所有 FILTER 值，没有 TE 筛选。

比较的是同一批新输入下两种合并规则的差异。旧历史合并文件的 151,677 条不能直接作为本次对照，因为发现下限、Xinghua2 来源及 BND 调用范围已改变。

**实际命令与默认规则**

```bash
# 在本目录运行；完整路径、时间、退出码见 logs/*.run.json
SURVIVOR merge all_sample_vcf_list.txt 0.1 0 1 0 0 50 raw/merged_survivor_previous.vcf

jasmine file_list=all_sample_vcf_list.txt \
  --use_edit_dist --keep_var_ids --output_genotypes \
  threads=16 out_dir=tmp_jasmine_paper out_file=raw/merged_jasmine_paper.vcf
```

SURVIVOR 1.0.7，提交 `07404b74d3fe42b20f1362f4d8625f284b9d536c`；Jasmine 1.1.5，提交 `a4f6988eeb191c1f8cf2e8d73869b6a5ffc14665`。SURVIVOR 沿用旧规则，仅把最后的长度门槛 100 改为 50。Jasmine 使用所给算法参数，没有额外添加序列阈值、类型转换或样本内合并选项。

| 规则 | SURVIVOR 旧规则 | 指定 Jasmine 命令 |
|---|---|---|
| 距离 | 同染色体按当前输入断点跨度的 10%，要求两个断点均通过；INS 构造比较终点 | 默认在 POS–长度二维空间用欧氏距离；阈值约为 `max(100, round(0.5×SV长度))`，本命令没有额外固定上限 |
| 跨染色体 BND | 比例模式下使用 1000 bp | 无 SVLEN 时长度取 0，默认阈值为 100 bp，比较两个断点坐标 |
| 类型与方向 | 启用内部类型检查，但实测存在绕过；忽略方向；DUP 亚型归为 DUP，同染色体 BND 会归为 INV | 默认区分类型与方向；DUP 亚型分开；BND 作为 TRA 类内部聚类 |
| 同一样本多个调用 | 允许并入同一簇 | 默认禁止同一簇中出现同一样本的多个调用 |
| 插入序列比较 | 不比较序列 | `min_seq_id` 默认是 0，因此仅写 `--use_edit_dist` 不会实际施加序列相似度约束 |

上述距离、序列比较和样本约束依据本次运行版本的官方源码：[SURVIVOR 断点比较](https://github.com/fritzsedlazeck/SURVIVOR/blob/07404b74d3fe42b20f1362f4d8625f284b9d536c/src/merge_vcf/IntervallTree.cpp)、[Jasmine 参数与默认值](https://github.com/mkirsche/Jasmine/blob/a4f6988eeb191c1f8cf2e8d73869b6a5ffc14665/src/Settings.java)、[Jasmine 距离阈值构造](https://github.com/mkirsche/Jasmine/blob/a4f6988eeb191c1f8cf2e8d73869b6a5ffc14665/src/VariantInput.java)、[Jasmine 序列比较](https://github.com/mkirsche/Jasmine/blob/a4f6988eeb191c1f8cf2e8d73869b6a5ffc14665/src/Variant.java)。

后续补充核查发现：当前版本启用正值 `min_seq_id` 时，编辑距离返回值与相似度判据方向不一致；同位置候选还存在绕过序列检查的路径，均已用当前 JAR 的最小 VCF 测试复现。本次全量合并未启用正值序列门槛，不应把现有比较解释为序列筛选的收益，也不应直接添加 `--use_edit_dist min_seq_id=0.8/0.9` 并假定改进。见[序列门槛核查与复现](diagnostics/jasmine_sequence_threshold_check/核查说明.md)。

**合并后的数量**

下表按保留的代表原始调用类型统计；所有簇成员保存在 IDLIST 和簇对应表中。BND 数量是 VCF 断点记录簇数，不能直接解释为独立易位事件数。

| 类型 | 输入调用 | SURVIVOR | Jasmine | Jasmine 相对变化 |
|---|---:|---:|---:|---:|
| INS | 417,351 | 125,903 | 76,203 | -39.47% |
| DEL | 387,339 | 116,543 | 78,025 | -33.05% |
| INV | 1,107 | 465 | 378 | -18.71% |
| DUP:TANDEM | 2,148 | 1,452 | 1,197 | -17.56% |
| DUP:INT | 9 | 6 | 9 | +50.00% |
| BND | 13,862 | 12,270 | 12,424 | +1.26% |
| 合计 | 821,816 | **256,639** | **168,236** | **-34.45%** |

Jasmine 相比 SURVIVOR 减少 88,403 个记录簇。总数变化只能描述合并程度，不能证明准确率更高。两者默认距离模型、样本约束和 BND 类型处理均不同，因此这不是单独评估“是否启用编辑距离”的实验。

**同一批原始调用的分簇一致性**

- SURVIVOR 覆盖 821,816 条输入，缺失 0 条；Jasmine 覆盖 821,816 条，缺失 0 条。
- 成员完全相同的簇有 107,966 个，分别占 SURVIVOR 的 42.07%、Jasmine 的 64.18%；其中多调用簇 35,377 个。
- 359,886 条输入调用位于完全一致的簇中，占全部输入的 43.79%。
- 12,975 个 SURVIVOR 簇被分散到多个 Jasmine 簇；44,192 个 Jasmine 簇汇入了多个 SURVIVOR 簇。两类变化可以同时发生。
- SURVIVOR 判为同簇的原始调用对中，61.89% 在 Jasmine 中仍同簇；Jasmine 同簇调用对中，74.50% 在 SURVIVOR 中也同簇。调整兰德指数 ARI=0.6761。这些是方法一致性指标，没有使用真值集。

**样本支持数**

| 支持样本数 | SURVIVOR | Jasmine |
|---|---:|---:|
| 1 | 161,695 | 79,066 |
| 2–5 | 67,073 | 48,130 |
| 6–10 | 13,979 | 19,949 |
| 11–25 | 9,725 | 15,421 |
| 26–50 | 4,133 | 5,627 |
| 51 | 34 | 43 |

支持数表示多少份组装调用进入这个簇。合并输出的 `1/1` 来自原始 haploid 组装调用；没有调用的样本保留 `./.`，不能把它当作真实参考型。这一步没有做 vg 或 reads 联合分型。

**类型和输出元数据核查**

SURVIVOR 有 10,876 个簇包含同一个体的多条调用；Jasmine 为 0。混合原始 SV 类型的簇分别为 1,160 和 0；类型组合见 comparison.json。

其中 SURVIVOR 有 1,085 个 INS/DEL 混合簇，不能用“已经设置 type_save=1”推断其严格按类型分开。从真实数据抽取 GuangxiYellow 的 DEL 与 Venda 的 INS：两者 POS 均为 1:389731、长度均为 52 bp；使用同样的 `0.1 0 1 0 0 50` 参数，仍被合为一个 SUPP=2 的簇。复现输入与结果保存在 [diagnostics/type_mix_reproducer/](diagnostics/type_mix_reproducer/)。结合源码可解释为：overlap() 在类型判据不通过后的排序分支仍可能返回 0，insert() 又把 0 当成应合并；这不是两个调用在生物学上等价的证据。见 [SURVIVOR overlap 与 insert 实现](https://github.com/fritzsedlazeck/SURVIVOR/blob/07404b74d3fe42b20f1362f4d8625f284b9d536c/src/merge_vcf/IntervallTree.cpp)。因此 SURVIVOR 结果适合作为旧方法比较基线；对本次要求严格区分全部 SV 类型的后续分析，不建议直接沿用该版本的混合类型簇。

Jasmine 原始输出存在 378 条显式序列 INV 被默认 `fixImprecision()` 写成 `SVLEN=0`、`END=POS`。SURVIVOR 原始输出会改变部分 BND 和 DUP 的类型标签，并用簇平均长度搭配代表等位基因。若直接对这些原始输出按 SVLEN 再过滤，会产生类型误判或误删。对应实现见 [Jasmine VcfEntry](https://github.com/mkirsche/Jasmine/blob/a4f6988eeb191c1f8cf2e8d73869b6a5ffc14665/src/VcfEntry.java) 与 [SURVIVOR VCF 解析](https://github.com/fritzsedlazeck/SURVIVOR/blob/07404b74d3fe42b20f1362f4d8625f284b9d536c/src/vcfs/Merge_VCF.cpp)。

因此保留 `raw/` 下两份原始算法输出，同时发布两份整理后的目录 VCF：保留合并器选出的原始 ID，并恢复该 ID 在输入中的 CHROM、POS、REF、ALT、FILTER、INFO；增加 SUPP、SUPP_VEC、全部成员 IDLIST、MERGE_CLUSTER 和 CLUSTER_TYPES。按原输入样本名输出 51 列 GT；CN 一致时保留，同一样本同簇 CN 冲突时写缺失。整理过程不改变任何簇成员，所有比较均依据原始调用的分簇关系。代表 ID 回退次数为 SURVIVOR 0、Jasmine 0。

SURVIVOR 1.0.7 在样本内合并时会覆盖已有 FORMAT/ID，无法仅从普通输出完整恢复所有成员。为准确比较，额外构建了只旁路记录成员 ID 的版本；没有修改聚类、代表选择或 VCF 输出代码。该版本生成的 VCF 与未修改 1.0.7 的 VCF 经 `cmp` 检查逐字节相同后，才使用旁路成员表。修改补丁和一致性证明均保留。

输入中有 [107 条 INS/DEL](diagnostics/input_indel_length_mismatches.tsv) 的原生 SVLEN 与序列长度差不完全一致，但这些记录按两种长度算法都在所需区间内。本次为保持调用与比较的可追溯性保留原表示；进入 vg 建图准备时还应结合等位基因检查这些记录。

**已完成的验证与文件**

两份整理后 VCF 均已排序、BGZF 压缩并建立 tabix 索引。独立 bcftools 解析核对 51 个样本、全部类型、记录数、唯一 ID、支持向量和长度范围，警告为 0。非 BND 长度越界记录为 0。完整结果见 [validation.json](validation.json)。

- [Jasmine VCF](merged_jasmine_paper.vcf.gz) / [未压缩 VCF](merged_jasmine_paper.vcf)
- [SURVIVOR VCF](merged_survivor_previous.vcf.gz) / [未压缩 VCF](merged_survivor_previous.vcf)
- [类型数量](type_comparison.tsv)、[逐支持数统计](support_comparison.tsv)、[完整指标](comparison.json)
- [比较图 PNG](merge_comparison.png) / [PDF](merge_comparison.pdf)
- [簇之间的对应关系](cluster_crosswalk.tsv.gz)、[每条输入调用对应的两个簇](input_to_clusters.tsv.gz)、[差异例子](difference_examples.tsv)
- [SURVIVOR 完整成员](survivor.clusters.tsv.gz)、[Jasmine 完整成员](jasmine.clusters.tsv.gz)
- [原始 SURVIVOR 输出](raw/merged_survivor_previous.vcf)、[原始 Jasmine 输出](raw/merged_jasmine_paper.vcf)
- [输入过滤统计](input_filter_summary.json)、[软件指纹](software_versions.json)、[成员旁路补丁](scripts/survivor_provenance_only.patch)、[VCF 字节一致性证明](logs/survivor_trace_verified.txt)

要复现本次结果，可依次运行 `scripts/prepare_inputs.py`、`scripts/run_merges.py`、`scripts/build_survivor_trace.py`、`scripts/compare_merges.py`、`scripts/publish_and_validate.py`、`scripts/make_report.py`。已有原始合并输出和旁路成员表时，从 compare_merges.py 开始即可。源文件和 51 份原始 VCF 均未改写。
