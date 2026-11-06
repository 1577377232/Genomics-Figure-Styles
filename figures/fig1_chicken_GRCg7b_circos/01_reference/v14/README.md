# V14 着丝粒转移审计

供体 GCA_024206055.2；输入为用户提供的35个区域。只读取已有 BAM，不重新比对。

`GRCg7b.centromeres.approximate.tsv` 是主图的34条 GRCg7b BED 0-based half-open 候选。包括6条既有较高覆盖转移、10条既有部分转移、13条低覆盖同源匹配簇、5条纯侧翼限定的终端候选；不等于 CENP-A 验证位置。

`centromere_status_all35.tsv` 逐条记录供体输入去向，`main_chromosome_centromere_status.tsv` 记录全部41条参考主染色体，`centromere_changes_from_V13.tsv` 列16条保持、17条新增、1条Chr31修正。未定位供体36的候选单列，不复制到7b36。7b15/27/34/35/36/38/39 无支持标记。

`whole_chromosome_correspondence.tsv`/`primary_correspondence.tsv` 是染色体对应证据；`centromere_flanking_M_blocks.tsv` 是侧翼CIGAR M块，`eligible_alignment_spans.tsv` 保存可用全比对跨度；`anchor_method.json` 记录过滤与输入。供体29→7b33、34→7b29、35→7b37，小染色体不能按编号直接搬运。

读取已有准备表复现候选：`python3 06_scripts/select_centromeres_v14.py`（依赖pysam）。原始 BAM 锚点重建脚本为 `06_scripts/centromere_anchors_v14.py`，原始输入绝对路径见该脚本与JSON。完整科学限制与图形标记规则见 `07_figures/R_v14/README.md` 及 `QC_report.md`。
