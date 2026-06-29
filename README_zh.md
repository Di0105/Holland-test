# Holland-test: Holland 参数化台风风场 → AtmoSurge DL 风暴潮 hindcast（香港）

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.6+](https://img.shields.io/badge/python-3.6%2B-blue.svg)](https://www.python.org/)

**用 Holland (1980) 参数化热带气旋风场重建方法，为 AtmoSurge 深度学习风暴潮模型生成香港 6 个潮位站的历史台风风压驱动场。**

---

## 项目概述

本仓库实现 **STORM-lite Holland 参数风场模型**，加入科氏力、移动不对称性和 Willoughby-Rahn RMW 估算，对 **8 个影响香港的历史热带气旋**（2001–2023）生成逐站点的风速、风向和局地气压序列。

输出格式兼容 **AtmoSurge** 深度学习风暴潮预测模型。

### 研究范围

- **8 个 TC**: Utor (2001), Hagupit (2008), Vicente (2012), Usagi (2013), Hato (2017), Mangkhut (2018), Kompasu (2021), Saola (2023)
- **6 个香港潮位站**: 启德 (SE), 长洲 (CCH), 流浮山 (LFS), 西贡 (SKG), 大美督 (PLC), 横澜岛 (WGL)
- **三种方法对比**: Holland-driven AtmoSurge → observed-meteorology-driven AtmoSurge → SLOSH 历史 report error
- **核心发现**: Raw Holland forcing 在 $r \gg RMW$ 时系统性低估站点风速；校准工具已包含

---

## 目录结构

```
Holland-test/
├── README.md                          ← English README
├── README_zh.md                       ← 中文说明（本文件）
├── .gitignore
├── scripts/
│   ├── holland_atmosurge_pipeline_calculation.py   ← 主流程脚本
│   └── calibrate_holland_wind.py                   ← 风速校准工具
├── outputs/
│   ├── README.md                      ← 输出文件说明
│   └── wind_unit_audit.csv            ← MaxWind 单位核查
├── docs/
│   └── formula_explanation.md         ← 完整物理公式
└── data/
    └── README.md                      ← 输入数据格式说明
```

---

## 快速开始

### 环境要求

- Python 3.6+
- 零外部依赖（纯标准库）

### 运行主流程

```bash
python scripts/holland_atmosurge_pipeline_calculation.py \
    --track /path/to/AtmoSurge_test1_track.csv \
    --judy-table /path/to/JUDYTABLE3.2.csv \
    --tc-nontc /path/to/TC_NonTC_Comparison_Table.csv \
    --output-dir ./outputs
```

### 打印公式说明

```bash
python scripts/holland_atmosurge_pipeline_calculation.py --explain-formula
```

### 风速校准

```bash
python scripts/calibrate_holland_wind.py
```

---

## 关键设计原则

### 1. Rmax ≠ R34（严格区分）

| 字段 | 含义 | Holland 中用法 |
|------|------|---------------|
| `Rmax_average_km` | 最大风速半径 (RMW) | **Holland RMW 输入** |
| `R34_*` | 烈风圈半径 | **仅用于 QC，绝不作为 RMW** |

### 2. MaxWind 单位 = kt（海里/小时），不是 km/h

track CSV 的 `MaxWind` 单位为 **nautical miles per hour (knots)**。转换: `Vmax_ms = MaxWind × 0.514444`。已通过单位核查（`outputs/wind_unit_audit.csv`）。

### 3. AtmoSurge pressure 输入 = 站点局地气压 P(r)，不是中心气压 Pc

$$P(r) = P_c + \Delta P \cdot \exp(-(RMW/r)^B)$$

### 4. AtmoSurge wind 输入 = 站点风，不是 TC 中心 Vmax

站点风 = 对称 Holland 梯度风 + 移动背景流，降至 10 m 高度。

---

## RMW 缺失处理策略

1. **观测**: track CSV 中的 `Rmax_average_km`
2. **插值**: 同一 TC 事件内线性插值（禁止跨 TC 插值）
3. **公式**: Willoughby & Rahn (2004) — $R_{max} = 51.6 \cdot \exp(-0.0223 \cdot V_{max} + 0.0281 \cdot |lat|)$

RMW 来源记录在 `holland_inputs.csv` 的 `rmw_source` 列中。

---

## 已知局限与改进路线

| 问题 | 状态 | 解决方案 |
|------|------|----------|
| $r \gg RMW$ 时风速低估 | 🔴 进行中 | 放开 B clamp, R34 约束, 距离校正, ERA5 混合 |
| 固定 Penv = 1013.25 hPa | 🟡 计划中 | 逐事件 ERA5 外圈 MSLP |
| B clamp [0.8, 2.5] 太窄 | 🔴 进行中 | 允许 B ∈ [0.3, 3.0] + R34 验证 |
| 无 AtmoSurge 时序适配器 | 🟡 计划中 | 10 min 到 1 hr forcing 序列 |
| 单 Holland profile 局限 | 🟡 计划中 | 双 B 或 GAHM 方案 |

---

## 输出文件

| 文件 | 行×列 | 说明 |
|------|-------|------|
| `holland_inputs.csv` | 358 × 35 | 逐时步 B, RMW, Pc, Vmax, dP, Coriolis f, 移动速度 |
| `holland_station_forcing.csv` | 2148 × 22 | 站点级风速/风向、气压 P(r)、距 TC 中心距离 |
| `event_pattern_coverage.csv` | 8 × 12 | 每 TC 的 RMW 覆盖率和 WR 公式使用次数 |
| `largest_surge_reference_comparison.csv` | — | 实测最大 storm surge vs 最大 Holland 风 |

---

## 相关仓库

| 仓库 | 用途 |
|------|------|
| `Di0105/AtmoSurge_8TC_consolidated_package` | 最终合并 CSV + SLOSH + 图件 |
| `Di0105/CiteAgent-Copilot` | VS Code Copilot AI 文献引用助手 |

---

## 参考文献

- Jelesnianski, C.P. et al. (1992). *SLOSH: Sea, Lake, and Overland Surges from Hurricanes*. NOAA TR NWS 48.
- Holland, G.J. (1980). An analytic model of the wind and pressure profiles in hurricanes. *MWR*, 108(8).
- Holland, G.J. et al. (2010). A revised model for radial profiles of hurricane winds. *MWR*, 138(12).
- Willoughby, H.E. & Rahn, M.E. (2004). Parametric representation of the primary hurricane vortex. *MWR*, 132(12).
- Shashank, V.G. et al. (2022). 五种参数风场模型对比. *Applied Ocean Research*.
- Liu, F. & Sasaki, J. (2019). ERA5 + 参数台风混合方法. *Scientific Reports*, 9.

---

## 许可证

MIT License — 见 LICENSE 文件。

## 作者

Judy Zhu (Di0105) — 香港 TC 风暴潮 hindcast 项目, 2025–2026.
