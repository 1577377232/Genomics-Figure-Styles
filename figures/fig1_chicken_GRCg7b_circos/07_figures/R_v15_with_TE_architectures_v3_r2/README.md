# Figure 1 — V3 修订版（v3_r2）

基于已发布 V3，按本轮要求制作独立修订版。不覆盖 V2，也不覆盖先前 V3 的脚本、素材或输出。V15 原始脚本、真实数据、窗口、轨道、颜色及左侧主图比例保持不变。

## 运行

在仓库的 `figures/fig1_chicken_GRCg7b_circos/` 目录执行：

```bash
# 新机器缺少依赖时，沿用原安装器：
bash 06_scripts/install_composition_R_dependencies.sh
# 完整矢量、高清 PNG 和 DEBUG：
bash 06_scripts/run_compose_v15_with_te_architectures_v3_r2.sh
# 仅检查预览：
bash 06_scripts/run_compose_v15_with_te_architectures_v3_r2.sh --preview-only
# 输出到新的目录：
bash 06_scripts/run_compose_v15_with_te_architectures_v3_r2.sh --out-dir /absolute/new_output
```

新脚本为完整、独立 R 脚本，不调用 V15 原绘图流程，也不执行旧拼接脚本。运行依赖 R/rsvg/png；无需 VCF、原始基因组或重算统计。全量 Python 审计额外需要 NumPy/SciPy/Pillow 和 Poppler，不是 R 拼接的运行依赖。

## 五个模块放大

五个 SVG 的显示宽度均从3.95in增加为4.9375in，精确为1.25倍；标题、内部文字、模块、边线一起等比放大。保留1400×250原始画布比例。标题均为原生58px Georgia，在本版最终物理尺寸上约14.73pt；常见35px内部标签约8.89pt。

画布从14×10in加宽为17×10in；左侧主图继续使用原来的左上位置0.15in/0.15in和9.70in高度，等比缩放值0.97保持。没有为了右侧模块缩小左图。

右侧模块顺序为LINE/LTR/SINE/DNA/TEV，统一左侧11.75in，顶部依次1.50/2.95/4.40/5.85/7.30in；主体中心间距1.45in，整体中心约4.989in，与圈图中心4.9515in相近。

前四个结构文件继续直接读取 `assets/TE_architectures/` 的原SVG，按字节不变。本次 TEV 使用独立副本 `assets/TE_mediated_SV/TE_mediated_SV_classification_v3_r2.svg`；旧版原SVG和 >50% 副本都保留。

TEV 副本删除内部 `classification-arrow`、`tev-result` 与结果框文字，只保留标题TEV、SV括线和标签、TE覆盖条、浅灰非TE区及 >80% 阈值。剩余横条居中展开至x=110..1290，阈值严格为110+1180×0.8=1054，示例TE长1008/1180≈85.4%。这属于示意排版调整，不代表真实序列长度变化。

## 五个真实起点与平滑路径

本版不是把旧TEV起点简单向上挪动，而是从原V15的 `track_geometry.tsv`、`sector_panels.tsv`、`chromosome_display_mapping.tsv`、`figure_geometry.tsv` 及中心坐标重算：

| 外部模块 | V15实际轨道 | 内部轨道索引 | 起点参考 |
|---|---|---:|---|
| LINE | LINE | 4 | Chr1 x=0，右上端点24° |
| LTR | LTR | 5 | Chr1 x=0，右上端点24° |
| SINE | SINE | 6 | Chr1 x=0，右上端点24° |
| DNA | DNA | 7 | Chr1 x=0，右上端点24° |
| TEV | TIP（中心标签TEV density） | 10 | Chr1 x=0，右上端点24° |

V15 的前两个内部轨道是名称/坐标，因此内部索引10对应生物学列表中第8层、最内层TEV；索引7是DNA轨道，不是Chr7。

TEV内/外半径为0.3185/0.3765；选取带内25%位置得到r=0.333。用圈心(cx,cy)、圈图半径R及θ=24°计算 SVG 页面坐标：x=cx+R×0.333×cosθ，y=cy−R×0.333×sinθ。最终起点为页面(5.9058679455in,4.4322391293in)，y从上计；它位于圈心右上方。

实际SVG检查不止依赖表格：在同一半径向Chr1内部偏0.05°的验证点，确实落入原始V15珊瑚色填充多边形；真实起点与该多边形边界距离小于0.005个原SVG点，记录于 `validation_report.json`。

五条路径均从真实扇区边界离开。短三次Bézier起始段绕过不可变的百分比/密度轴标签，独立水平引段到达圈图外侧后以平滑三次Bézier扇出，再水平接小闭合箭头。TEV起始段沿右侧缺口内侧自然转出，避免从上方直接斜穿外层轨道和标签；不用共同主干，也不合并顶点。线色#253246、宽0.80pt，箭头长0.075in。

`connector_anchors.tsv` 保存轨道、Chr1/0bp、真实24°端点、r和物理坐标；`connector_paths.tsv` 保存真实SVG M/C/L/C/L路径；`connector_curve_samples.tsv` 保存曲线采样，修改入口为 R 脚本顶部 `composition_config()`。

## DEBUG 与正式文件

- `.svg`、`.pdf`：正式矢量，V15全部325个裁剪路径保留，无嵌入栅格图。
- `.png`：10200×6000、600dpi。
- `.preview.png`：2550×1500、150dpi。
- `.DEBUG.svg` / `.DEBUG.preview.png`：五个起点依次标1..5，顶部键明确轨道索引与坐标。
- `.DEBUG.anchor_zoom.svg` / `.DEBUG.anchor_zoom.png`：起点局部放大，5号标记对应珊瑚TEV最内轨；PNG2130×1380。

DEBUG层独立附加，正式输出不含起点标记。起点圆心逐一核对为实际路径起点，不能用DEBUG标记掩饰错位。

## 审计与数据口径

独立Python直接解析真实SVG三次曲线路径，以不超过0.08pt步长采样，检查源轨道、标签框、数据扇区和所有线对。最小间距扣除0.16pt采样界限后仍大于3.5pt。边界接触只允许真实起点，轨道内部穿越仍判失败。

`protected_baseline.json` 保存原V15、V2、旧V3和全部旧素材的SHA256；`input_preservation.tsv` 保存运行前后MD5。`validation_report.json` 和 `provenance/validate_composition_v3_r2.py` 给出完整审计，`QC_report.md` 为简版。

**TEV模式图沿用用户指定的 >80% 示意阈值；V15真实TEV density仍沿用 >50%、事件 >50bp 的原extended统计。未重新过滤或计算主图数据，二者不是同一统计口径。** 示例不是形成机制证据，也不是按真实序列比例绘制。
