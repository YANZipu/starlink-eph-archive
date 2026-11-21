# Starlink 公开星历档案（协方差真实性审计数据基础设施）

自动抓取 SpaceX 公开星历（[api.starlink.com/public-files/ephemerides](https://www.starlink.com/orbit-files)，
Modified ITC 格式，含 6×6 协方差），在云端即时处理成紧凑派生数据后归档。
**服务器只保留当前 vintage**（已实测旧文件 404），所以时间序列从本仓库第一次运行起算——
这是协方差真实性审计（TraCSS 时代星历质量对账）的唯一可行数据路径。

## 抓什么 / 存什么

| 内容 | 位置 | 体量估算 |
|---|---|---|
| 每个 vintage 的 MANIFEST（全量文件清单） | `manifests/*.txt.gz` | ~50 KB × 4/天 → ~6 MB/月 |
| 全星座逐卫星 cross-vintage 对比统计（协方差真实性核心表） | `pairs/*.csv.gz` | ~150 KB × 4/天 → ~18 MB/月 |
| 抽样 ~200 星的详细状态（30 分钟采样 r,v,P6） | `states_subset/YYYYMM/*.npz` | ~2 MB × 4/天 → ~250 MB/月 |
| 滚动状态 refs.npz（全星稀疏参考历元，供下次对比） | GitHub Release `latest`（滚动覆盖，不进 git） | ~5 MB 常量 |
| 月度冻结包 | GitHub Release `vYYYYMM` | ~250 MB/月 |

原始 2 MB/星文件**不存**（14 GB/vintage 不可行）；抽样星的逐 30 分钟状态足以做全部
协方差真实性分析，全星座逐星统计在 `pairs/` 里 100% 保留。

## `pairs/` 字段

每行 = 一颗卫星在相邻两个 vintage 间的对比：新 vintage 的状态（Hermite 插值到旧 vintage
的参考历元）视为"实现"，旧 vintage 的预测 + 协方差视为"预测"。
`med_dr_m`（中位位置偏差）、`med_chi2`（3 自由度马氏距离，~3 为协方差真实，>>3 为过自信）、
`b_0_6h…b_36_48h`（按预测龄期分桶）、`flag`（|dr|>1 km，机动样事件）。

## 运行节奏

- GitHub Actions 每天 4 次（01:23/07:23/13:23/19:23 UTC，vintage 约 8 小时一轮换），
  公开仓库 Actions 分钟数无限，你的电脑可以永久关机。
- 每次 ~10–20 分钟（下载 ~14 GB 流式处理即弃）。
- 月度工作流把上月 `states_subset` 打成 Release 冻结包。

## 首次部署

1. 浏览器建一个**空 Public 仓库**（不要勾选 README/license）。
2. 本目录下：
   ```
   git remote add origin https://github.com/<你的用户名>/<仓库名>.git
   git push -u origin main
   ```
3. Actions → fetch-vintage → Run workflow 手动跑一次，看日志末尾的 JSON 报告
   （`sats_total`、`med_dr_nonmaneuver_m`、体积）。
4. 之后全自动；`runs.log` 每次运行追加一行摘要，可随时检查健康状态。

## 注意

- 文件头标注 `UVW` 框架、`ephemeris_source: blend`；跨源对比（vs TLE）时需先确认框架语义。
- 仓库 60 天无提交时 GitHub 会自动停用定时器——本仓库每天自动提交，自续期。
- 数据为 SpaceX 公开发布，公开再分发无碍；引用时注明来源 starlink.com/orbit-files。
