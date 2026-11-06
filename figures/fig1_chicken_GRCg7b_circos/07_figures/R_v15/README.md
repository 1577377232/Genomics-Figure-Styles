# V15 最终文字调整与 Fig1 上传

基于 V14，只修改标签：删除开口处 SNP/MB (linear)、SV/MB (log1p)、TIP/MB (log1p) 标题，保留数值刻度与短渐变条；中心三个名称为 **SNV density、SV density、TEV density**。

窗口、原始数据、统计定义、TE独立尺度、平滑、着丝粒、配色和48°开口均不变。中心 TEV density 仍指原 extended TIP/repeat面板，包括INS/DEL及简单重复、低复杂度和卫星重复，改名不代表数据范围变化。密度数值仍为每1,000,000参考碱基，TE为覆盖率百分比。

运行：

```bash
bash 06_scripts/run_circos_R_v15.sh
```

完整脚本 `06_scripts/plot_chicken_circos_v15.R`；PNG、PDF、SVG在本目录。旧版均保留。

GitHub目标目录：`https://github.com/1577377232/Genomics-Figure-Styles/tree/main/figures/fig1_chicken_GRCg7b_circos`。

该目录包含最终绘图代码、shell入口、实际窗口表、小型BED/位点目录及着丝粒依据；不上传个体基因型矩阵或大型原始FASTA/BAM/VCF，原路径及大小见 MATERIALS.tsv。可直接使用包内表重绘，无需原始大文件。
