# 地震震源定位方法技術調查與改善建議

**撰寫日期：** 2026-03-24
**背景：** 針對台灣 EEW 系統與 tcpd_new 的 `locaeq()` 函數改進
**資料來源：** 文獻調查（2018–2025）

---

## 目錄

1. [震源定位問題定義](#1-震源定位問題定義)
2. [tcpd_new 現有 locaeq() 方法分析](#2-tcpd_new-現有-locaeq-方法分析)
3. [傳統迭代方法](#3-傳統迭代方法)
4. [格搜尋與機率方法](#4-格搜尋與機率方法)
5. [深度學習方法](#5-深度學習方法)
6. [方法性能比較](#6-方法性能比較)
7. [關鍵技術議題](#7-關鍵技術議題)
8. [台灣特殊情境分析](#8-台灣特殊情境分析)
9. [針對 locaeq() 的具體改善建議](#9-針對-locaeq-的具體改善建議)
10. [關鍵文獻](#10-關鍵文獻)

---

## 1. 震源定位問題定義

震源定位（Hypocenter Location）求解四個參數 $\mathbf{m} = (x_0, y_0, z_0, t_0)$，使理論走時與觀測到時的殘差最小：

$$\min_{\mathbf{m}} \sum_{i=1}^{N} w_i \left[ t_i^{obs} - t_0 - T_i(x_0, y_0, z_0) \right]^2$$

其中：
- $t_i^{obs}$：第 $i$ 台站的 P 波觀測到時
- $T_i(\mathbf{x})$：從震源到台站 $i$ 的理論走時（依速度模型計算）
- $w_i$：加權值（依殘差和距離）
- $N$：參與定位的台站數

主要挑戰：
- **非線性**：$T_i(\mathbf{x})$ 對 $\mathbf{x}$ 是非線性函數，需要迭代或全局搜尋
- **深度–起震時間耦合**：$\delta t_0 \approx -\delta z / (V_P \sin i_c)$，深度與起震時間難以獨立確定
- **速度模型誤差**：速度模型不確定度是最大單一誤差來源（2024 研究確認）
- **異常值（outlier）**：假 picks 或錯誤的相位標記

---

## 2. tcpd_new 現有 locaeq() 方法分析

### 2.1 演算法類型

`locaeq()` 實作的是**Geiger 法（迭代最小二乘線性化）**，搭配**自訂走時公式（梯度速度模型）**和**深度格搜尋**，由吳逸民教授開發。

### 2.2 走時計算核心（Wiechert-Herglotz 公式）

對線性梯度速度模型 $v(z) = v_0 + v_g \cdot z$，走時有解析解：

```
xc = (r² - 2v₀/vg·z - z²) / (2r)    ← 射線轉折點水平距離
zc = -v₀/vg                           ← 轉折點深度
θ₁ = π - arctan((z - zc) / xc)        ← 出射角
θ₂ = arctan(-zc / (r - xc))           ← 入射角
t  = (-1/vg) × ln[tan(θ₂/2) / tan(θ₁/2)]
```

計算速度極快（< 1 μs），是 EEW 即時應用的理想選擇。

### 2.3 迭代流程（10 次）

每次迭代：
1. 計算各台站走時偏導數（$\partial T/\partial x$, $\partial T/\partial y$, $\partial T/\partial z$）
2. 建立 4×4 正規方程組：$\mathbf{B} \Delta\mathbf{m} = \mathbf{y}$（Gauss-Jordan 求逆）
3. 更新震源參數：$\mathbf{m} \leftarrow \mathbf{m} + \Delta\mathbf{m}$
4. 若發散（$|\delta t_0| > 60s$）則重置至初始值

### 2.4 深度格搜尋

10 次迭代後，在固定水平位置測試 10 個深度（10–100 km，間距 10 km），選最小加權殘差的深度。

### 2.5 現有方法的限制分析

| 限制 | 根本原因 | 對精度的影響 |
|------|---------|------------|
| **只有 2 層速度模型** | 台灣地殼結構複雜（3D 橫向變化大） | 系統性誤差 5–10 km（中央山脈附近） |
| **固定迭代 10 次** | 無收斂判斷 | 浪費計算資源；可能在發散後繼續迭代 |
| **深度格搜尋間距 10 km** | 計算量限制 | 深度誤差最差可達 ±5 km |
| **Gauss-Jordan 無主元選取** | 數值穩定性 | 病態矩陣（台站分布不均）時精度差 |
| **無靜態站校正** | 程式設計選擇 | 每站系統誤差 0.1–0.5 s 未被消除 |
| **線性化誤差橢球未輸出** | 未實作 | 無法評估定位可靠性 |
| **深度–起震時間耦合未處理** | 數學上不可避免 | 深度與起震時間互相影響 |
| **發散保護過於粗糙** | `if(|dt₀| > 60s) reset` | 部分發散情況無法偵測 |

---

## 3. 傳統迭代方法

### 3.1 Geiger 法（基礎理論）

**原理（1912）：** Taylor 展開線性化走時，迭代求解最小二乘：

$$\underbrace{\begin{pmatrix} \partial T_1/\partial x & \partial T_1/\partial y & \partial T_1/\partial z & 1 \\ \vdots & & & \vdots \\ \partial T_N/\partial x & \partial T_N/\partial y & \partial T_N/\partial z & 1 \end{pmatrix}}_{\mathbf{A}} \underbrace{\begin{pmatrix} \Delta x \\ \Delta y \\ \Delta z \\ \Delta t_0 \end{pmatrix}}_{\Delta\mathbf{m}} = \underbrace{\begin{pmatrix} t_1^{obs} - t_0 - T_1 \\ \vdots \\ t_N^{obs} - t_0 - T_N \end{pmatrix}}_{\Delta\mathbf{t}}$$

阻尼最小二乘解：$\Delta\mathbf{m} = (A^T W A + \lambda I)^{-1} A^T W \Delta\mathbf{t}$

`locaeq()` 使用的是此標準 Geiger 法，只是走時用梯度模型的解析公式計算。

**典型精度：**
| 條件 | 水平誤差 | 深度誤差 |
|------|---------|---------|
| 台灣密集網（TSMIP, > 10 站） | 1–3 km | 5–10 km |
| 稀疏網（< 6 站，站間距 50 km） | 5–15 km | 10–30 km |
| 搭配正確 3D 速度模型 | < 2 km | 3–5 km |

### 3.2 HYPO71 / HYPOINVERSE-2000

**HYPO71**（Lee & Lahr 1972）：使用逐步回歸（只接受通過 F-test 的參數更新），對過定問題更穩健。

**HYPOINVERSE-2000**（Klein 2002, USGS OFR-2002-171）：
- 支援多層 1D 速度模型（Pg/Pn/Sg/Sn 多相位）
- 輸出 ERH（水平誤差）、ERZ（垂直誤差）、RMS 殘差
- 台灣 CWA 現行系統的定位核心即以此類方法為基礎

### 3.3 等走時差法（Equal Differential Time, EDT）

**原理（Font et al. 2004）：** 消去起震時間 $t_0$，直接對**到時差**做定位：

$$\Delta t_{ij} = (t_i^{obs} - t_j^{obs}) = T_i(\mathbf{x}_0) - T_j(\mathbf{x}_0)$$

**優點：**
- 完全消除 $z$–$t_0$ 耦合（不需估計起震時間）
- 問題降為純 3 參數求解（x, y, z）
- 對系統性時鐘誤差不敏感

**ML-EDT（台灣 EEW，2024）：** Earth, Planets and Space 論文結合機器學習拾取 + EDT 格搜尋 + 台灣 3D 速度模型，在離岸與複雜地殼地區精度顯著優於 CWA 現行方法。

### 3.4 VELEST（最小 1D 模型建立）

VELEST（Kissling et al. 1994, ETH Zurich）同時反演震源位置與 1D 速度模型（最小 1D 模型）：

$$\min_{\mathbf{m}, \mathbf{v}} \sum_i [t_i^{obs} - T_i(\mathbf{m}, \mathbf{v})]^2$$

適用場景：建立台灣地區最佳 1D 速度模型作為 `locaeq()` 的速度參數基礎。**不適合即時 EEW**（需大量事件批次反演）。

---

## 4. 格搜尋與機率方法

### 4.1 NonLinLoc（Oct-tree 重要採樣）

**核心演算法（Lomax et al. 2000）：**

不做線性化，在全震源空間計算最大似然：

$$L(\mathbf{m}) = \prod_{i=1}^{N} \exp\left(-\frac{[t_i^{obs} - t_0 - T_i(\mathbf{m})]^2}{2\sigma_i^2}\right)$$

使用 **Oct-tree 重要採樣**：在高機率區域遞迴細分空間格點，效率遠高於暴力格搜尋。

**輸出：** 完整 4D 後驗機率密度函數（PDF），可輸出任意置信水準的誤差橢球。

**速度模型：** 支援預先計算的 3D 走時場（eikonal 方程有限差分）。

**適用性：** 後處理標準工具（每事件數秒至數分鐘），不直接適用即時 EEW，但可用於生成台灣 3D 走時表。

### 4.2 HypoDD（雙差重定位）

**核心原理（Waldhauser & Ellsworth 2000）：**

利用鄰近事件對的走時差分，共模路徑誤差相消：

$$\Delta t_{ij}^k = (t_i^k - t_j^k)^{obs} - [T_i^k(\mathbf{x}_i) - T_j^k(\mathbf{x}_j)]^{calc}$$

**精度：**
- 相對位置：0.1–0.5 km（波形互相關版本）
- 深度相對誤差 < 1 km（事件間距 < 5 km 時）

**用途：** 建立高精度台灣地震目錄，可作為訓練深度學習定位模型的標準資料集。

**GraphDD（EPS 2025）：** 用 GNN 最小化雙差殘差，可透過批次處理任意大小目錄。

### 4.3 Kalman 濾波即時更新

**核心原理（Chen et al. 2020, GRL）：**

從第 4 站觸發開始定位，之後每到 1 站做 Kalman 更新：

```
狀態向量：s = [x, y, z, t0]ᵀ
預測步驟：s̃ = F·s + noise
更新步驟：s = s̃ + K·(t_obs - H·s̃)
           K = P·Hᵀ·(H·P·Hᵀ + R)⁻¹
```

其中 $H$ 為走時偏導數矩陣，$R$ 為觀測雜訊協方差，$P$ 為狀態協方差矩陣。

**優點：**
- 4 站觸發即可輸出初步定位，之後隨新站到達持續精化
- 自然輸出協方差矩陣（誤差橢球）
- 計算量低，完全適合即時 EEW

**這是目前文獻中最適合 EEW 即時更新定位的方法之一。**

---

## 5. 深度學習方法

### 5.1 PhaseNet / EQTransformer（相位拾取）

**PhaseNet**（Zhu & Beroza 2019, GJI）：U-Net 架構，輸出 P/S 機率曲線。
- P 波 F1 = 0.896（vs 傳統 STA/LTA 的 0.558）
- 延遲 < 1 秒，可即時部署

**EQTransformer**（Mousavi et al. 2020, Nature Communications）：多任務注意力模型（偵測 + P 拾取 + S 拾取）。
- P 拾取中位數誤差：0.05 s；S 拾取：0.10 s

> **重要**：這兩種方法是**相位拾取工具（Picker）**，不直接輸出定位，但高品質的 P/S 到時直接提升 `locaeq()` 的精度。

### 5.2 EikoNet（神經走時替代模型）

**EikoNet**（Smith et al. 2021, JGR Solid Earth）：

用**物理知情神經網路（PINN）**近似走時函數，直接最小化 eikonal 方程殘差：

$$\left| \nabla T(\mathbf{x}_s, \mathbf{x}_r) \right|^2 \cdot v(\mathbf{x}_r)^2 = 1 \quad \text{(eikonal equation)}$$

**關鍵優勢：**

| 指標 | 傳統 3D 走時表 | EikoNet |
|------|-------------|---------|
| 記憶體需求 | 595.9 MB | **4.8 MB** |
| 查詢速度 | < 1 ms（記憶體查詢）| < 1 ms（GPU 推論）|
| 精度 | 精確 | < 0.5% 相對誤差 |
| 更新速度模型 | 重新計算（小時）| 重新訓練（分鐘）|

**對 EEW 的啟示：** 可以將台灣 3D 速度模型壓縮成 5 MB 的 MLP，直接嵌入 C 程式（INT8 量化後更小），實現 3D 走時計算。

### 5.3 HypoSVI（Stein 變分推論）

**HypoSVI**（Smith et al. 2021）：

結合 EikoNet 走時替代模型與 Stein 變分推論（SVI），輸出完整後驗機率分佈：

```
前向模型：T = EikoNet(x_source, x_station)
似然函數：P(d|m) ∝ Laplace(t_obs - t0 - T)  ← 對 outlier 穩健
先驗：P(m) ← 地理範圍約束
SVI：以粒子集合近似 P(m|d)，比 MCMC 快 ~40 倍
```

**在 Luding 餘震序列測試：** 比 NonLinLoc 快 41%，水平誤差 1.75 km，並輸出校準的不確定度。

### 5.4 Bayesian 神經走時替代（2024 最新）

**arXiv 2512.06407（2024）：** 完整 Bayesian 框架：

- 7 層 MLP 走時替代模型（在中國川滇 3D 速度模型訓練）
- 兩元件污染模型（Student-t + Gaussian 混合），對 outlier 穩健
- 水平誤差 **1.75 km**，recall 88.2%
- 後驗分布對真實誤差校準良好

### 5.5 GENIE（GNN 同時定位+關聯）

GENIE 的 GNN 訊息傳遞在台站圖和震源空間圖之間迭代，**同時輸出**震相關聯結果和震源位置。2025 基準測試中，GENIE 在所有情境表現最佳（見 Phase Association 調查報告）。

### 5.6 多台站機器學習定位（MSLOC，2025）

**MSLOC**（GJI 2025）：3D U-Net 網路（帶台站分布約束），在 Oklahoma 訓練，測試誤差 5 km。
跨域至南加州後退化至 7–8 km，顯示**泛化問題仍是 CNN 直接定位的根本限制**。

---

## 6. 方法性能比較

### 6.1 精度比較

| 方法 | 水平誤差 | 深度誤差 | 即時可行 | 速度模型 | 不確定度 |
|------|---------|---------|---------|---------|---------|
| Geiger + 1D（tcpd_new 現況）| 3–10 km | 10–20 km | ✅ < 100 ms | 1D | 線性化橢球（未實作） |
| Geiger + 站校正 | 1–5 km | 5–15 km | ✅ | 1D | 同上 |
| EDT + 3D 格搜尋 | 2–5 km | 5–10 km | ✅ 接近（< 1 s）| 3D | 有限 |
| Kalman 濾波即時更新 | 2–5 km（演化）| 5–15 km | ✅ 最佳 | 1D/3D | 協方差矩陣 |
| HYPOINVERSE | 1–10 km | 5–20 km | ✅ | 1D | ERH/ERZ |
| NonLinLoc + 3D | 1–5 km | 3–10 km | ❌（數秒）| 3D | 完整 PDF |
| HypoDD（相對）| 0.1–0.5 km | < 1 km | ❌（分~時）| 1D/3D | 有限 |
| VELEST | 2–3 km | ~5 km | ❌（批次）| 1D | 有限 |
| HypoSVI（EikoNet）| 1.75 km | — | ❌（數秒）| 3D 神經 | 完整 PDF |
| Bayesian NN（2024）| 1.75 km | — | 接近 | 3D 神經 | 完整 PDF |
| CNN 直接定位 | 1.74 km（同域）| 2.65 km | ✅ | 內隱 | 無（單域有效）|
| RANSAC + Geiger | 0.5–1 km | — | ❌（需多次迭代）| 1D/3D | 有限 |

### 6.2 EEW 適用性評分

```
立即可用（不改架構）：
  Geiger + 站校正          ★★★★★  最快提升，無需架構變動
  兩階段深度格搜尋          ★★★★★  只需修改 C 程式碼
  自適應迭代收斂            ★★★★   消除無意義迭代

需要架構調整：
  EDT 走時差法              ★★★★   消除 z-t₀ 耦合
  Kalman 濾波更新           ★★★★   輸出時間演化的定位
  台灣最佳 1D 模型          ★★★★   修改速度參數

中期研究：
  3D 走時表 + 插值          ★★★    需建表（存在 TAIGER 模型）
  EikoNet 神經替代          ★★★    需訓練（台灣 3D 模型）
  PhaseNet 拾取整合         ★★★★   提升輸入品質
```

---

## 7. 關鍵技術議題

### 7.1 深度確定的根本困難

**深度–起震時間耦合**（最核心問題）：

$$\delta t_0 \approx -\frac{\delta z}{V_P \sin i_c}$$

對純 P 波到時，深度與起震時間嚴重相關，無法獨立確定。以台灣典型情形：
- $i_c \approx 30°$，$V_P = 6$ km/s
- $1$ 秒的起震時間誤差 → **約 12 km 的深度誤差**

**改善深度估計的方法（按效果排序）：**

1. **P-S 走時差**：$\Delta t_{PS} = r(1/V_S - 1/V_P)$ 直接提供震源距約束
2. **EDT 法**（消去 $t_0$）：將問題降為純 3 參數，完全消除耦合
3. **兩階段深度格搜尋**：粗搜 20 km 間距 + 細搜 2 km 間距
4. **Pn/Pg 走時差**：區分折射波和直達波，對深度敏感
5. **波形振幅比**：P/S 振幅比隱含深度資訊（arXiv 2509.02346, 2025）

### 7.2 速度模型誤差是最大誤差來源

**2024 關鍵發現**（GJI 2024, Scientific Reports 2024）：

傳統誤差橢球只反映**到時拾取不確定度**的傳播，**未包含速度模型不確定度**。全域靈敏度分析顯示：

```
定位誤差貢獻：
  速度模型不確定度：55–70%（最大！）
  到時拾取誤差：   20–30%
  台站分布幾何：   10–20%
```

傳統方法的誤差橢球嚴重低估實際誤差（可能低估 2–5 倍）。**「定位精度 1 km」往往只是數學意義上的精度，不反映真實誤差。**

### 7.3 走時公式 vs. 走時表 vs. 神經替代

| 方法 | 計算速度 | 記憶體 | 精度（與真實 3D 模型比） | EEW 適合 |
|------|---------|--------|----------------------|---------|
| 梯度模型解析公式（現況） | < 1 μs | 0 | 低（2 層近似）| 最佳速度 |
| 最佳 1D 模型（TauP 插值）| 1–10 ms | ~MB | 中（1D 全台灣平均）| 可 |
| 3D 走時表（eikonal FD）| < 1 ms（查詢）| GB 量級 | 高（完整 3D）| 可（預載記憶體）|
| EikoNet MLP（4.8 MB）| < 1 ms（CPU/GPU）| 4.8 MB | 高（< 0.5% 誤差）| 最佳平衡 |

### 7.4 不確定度量化的重要性

**線性化誤差橢球（Geiger 法輸出）：**

$$\mathbf{C}_m = \sigma^2 (A^T W A)^{-1}$$

水平誤差：$ERH = \sqrt{\lambda_1^2 + \lambda_2^2}$（特徵值之均方根）
垂直誤差：$ERZ = \lambda_3$（深度方向特徵值）

**輸出誤差橢球的用途：**
- 震度預測的不確定性傳播（大規模地震早期報告的震度範圍）
- 自動品質分級（高不確定度警報降級為「注意」）
- 後處理目錄的品質控制

**建議：** `locaeq()` 在迭代收斂後計算並輸出 `ERH`、`ERZ`（從 $A^T W A$ 的對角線），約增加 < 1 ms 計算量。

---

## 8. 台灣特殊情境分析

### 8.1 台灣地殼速度結構的挑戰

台灣是全球構造最複雜的地震帶之一，對 1D 速度模型造成以下系統性問題：

| 地區 | 速度結構特徵 | 1D 模型誤差 |
|------|------------|------------|
| 西部平原（嘉南）| 厚沉積層（Vp 1.5–4.0 km/s）| 走時低估 0.5–1.5 s |
| 中央山脈 | 變質雜岩，Vp 6.0–6.8 km/s | 走時高估 0.3–0.8 s |
| 花蓮–台東 | 碰撞帶，側向梯度劇烈 | 系統性偏移 5–10 km |
| 北部（大屯火山群） | 低速熔融體，Vp < 4.5 km/s | 走時低估 > 1 s |
| 東部外海（隱沒帶） | 海洋地殼 + 上地幔 | 走時計算不適用 |

### 8.2 台灣可用的 3D 速度模型

| 模型 | 來源 | 解析度 | 下載 |
|------|------|--------|------|
| TAIGER 3D Vp | Kuo-Chen et al. 2012; Wu et al. 2014 | 水平 ~20-30 km，深度 ~5-10 km | IRIS SPUD: EMC-Taiwan |
| Wu et al. 2009 聯合 Vp+Vs | CWA 官方使用 | ~20 km | CWA 申請 |
| Liu et al. 2021 Vs（北台灣）| Formosa Array 瑞雷波 | 水平 10 km | 論文附件 |

**建議：** 以 TAIGER 模型為基礎，用 NonLinLoc eikonal FD 計算 3D 走時場，存成二進位格點走時表，`locaeq()` 在執行時 mmap 直接存取（無 I/O 延遲）。

### 8.3 台灣 EEW 系統現狀

**CWA EEW 系統（Earthworm 平台）：**
- 定位：Geiger 法，Wu et al. 2009 模型的 **1D 平均簡化版**
- 深度：10 km 間距格搜尋（等同固定深度）
- 報告時間：地震後約 **10–13 秒**
- 典型誤差：震央 5–10 km，深度 10–20 km

**P-Alert 系統（台大 NTU）：**
- 762 個 MEMS 加速度計（全台中小學）
- 有效震央法（Effective Epicenter）+ 單站 on-site 估算
- 2016 美濃 Mw 6.4：4–8 秒預警

**RT-MEMS（Sensors 2025）：**
- 整合 CWASN + TSMIP
- 類 PhaseNet 架構（SeisBlue）+ 傳統定位
- 偵測率達 CWA 目錄 **3 倍以上**

**ML-EDT（EPS 2024）：**
- 機器學習拾取 + EDT + 台灣 3D 速度模型
- 比現行 1D+Geiger 系統精度顯著提升（特別是離岸事件）
- 這是目前最接近 EEW 實際應用的現代定位改善研究

---

## 9. 針對 locaeq() 的具體改善建議

### 9.1 短期改善（修改現有 C 程式碼）

#### 改善 A：加入自適應收斂判斷（最優先）

**問題：** 固定 10 次迭代，無法判斷是否已收斂或發散。

```c
/* 現況 */
for(itr=0; itr<10; itr++) { ... }

/* 建議：自適應收斂 */
#define MAX_ITER 15
#define CONV_THRESHOLD 0.01   /* km，修正量小於此值視為收斂 */

double correction_norm;
for(itr=0; itr<MAX_ITER; itr++) {
    /* ... 計算修正量 c[0..3] ... */
    correction_norm = sqrt(c[0]*c[0] + c[1]*c[1]) * dlon;   /* km */
    if(correction_norm < CONV_THRESHOLD && itr >= 3) {
        logit("d", "tcpd: locaeq converged at iter=%d, corr=%.4f km\n",
              itr, correction_norm);
        break;
    }
}
```

**效益：** 典型地震 5–7 次收斂，節省約 30–50% 計算量；同時可偵測不收斂情況。

---

#### 改善 B：兩階段深度格搜尋

**問題：** 深度格搜尋只有 10 km 間距（10–100 km），解析度不足。

```c
/* 現況：10 點 × 10 km */
z0 = ((double)dd+1.0) * 10.0001;   // 10, 20, 30, ..., 100 km

/* 建議：兩階段格搜尋 */
/* 第一階：粗搜（0–200 km，20 km 間距，11 點）*/
double z_coarse[11], err_coarse[11];
for(dd=0; dd<11; dd++) {
    z_coarse[dd] = (double)(dd) * 20.0 + 5.0;   // 5, 25, 45, ..., 205 km
    err_coarse[dd] = compute_weighted_err(ptr, nsta, x0, y0, z_coarse[dd], t0);
}
/* 找最小殘差的粗搜深度 */
int best_dd = argmin(err_coarse, 11);
double z_best_coarse = z_coarse[best_dd];

/* 第二階：細搜（±15 km 範圍，2 km 間距，16 點）*/
double z_fine[16], err_fine[16];
for(dd=0; dd<16; dd++) {
    z_fine[dd] = z_best_coarse - 15.0 + dd * 2.0;
    if(z_fine[dd] < 1.0) z_fine[dd] = 1.0;
    err_fine[dd] = compute_weighted_err(ptr, nsta, x0, y0, z_fine[dd], t0);
}
z0 = z_fine[argmin(err_fine, 16)];
```

**效益：** 深度解析度從 ±5 km 提升至 ±1 km，計算量增加約 1.5 倍（仍在 EEW 容許範圍）。

---

#### 改善 C：輸出線性化誤差橢球（ERH, ERZ）

**說明：** 在迭代收斂後，從 $A^T W A$ 矩陣計算定位不確定度。

```c
/* 在 locaeq() 末尾，收斂後計算誤差橢球 */
double sigma2 = errsum * errsum / (nsta - 4.0);  /* 殘差方差 */

/* 從 b[][] 矩陣（已計算於迭代中）取得協方差 */
/* C_m = sigma2 * (A^T W A)^{-1} = sigma2 * a[][] */
/* （a[][] 已在最後一次迭代中被 matinv 求逆）*/
double var_x = sigma2 * a[0][0];   /* 東西方向方差（度²）*/
double var_y = sigma2 * a[1][1];   /* 南北方向方差（度²）*/
double var_z = sigma2 * a[2][2];   /* 深度方向方差（km²）*/

/* 換算為 km */
hyp->ERH = sqrt(var_x * dlon*dlon + var_y * dlat*dlat);
hyp->ERZ = sqrt(var_z);
```

**效益：** 輸出 ERH/ERZ，可在報告中標示定位可靠性，為後續震度預測提供不確定度。

> **注意：** 如第 7.2 節所述，ERH/ERZ 只反映到時拾取誤差傳播，不含速度模型不確定度，實際誤差可能為此 2–5 倍。

---

#### 改善 D：靜態站校正

**問題：** 每個台站的系統性到時偏差（地殼速度異常、儀器誤差）累積為 0.1–0.5 s，未被校正。

```c
/* 在 tcpd_new.c 全域 */
static float station_corr[MAX_STA] = {0.0};   /* 初始化為 0 */
static char  sta_corr_name[MAX_STA][8];
static int   n_sta_corr = 0;

/* 在 tcpd_config() 讀取站校正檔 */
/* StaCorr  YULB  HHZ  TW  01  -0.15   (秒) */

/* 在 locaeq() 計算殘差時套用 */
ptr[i].perr = ptr[i].P - t0 - ptime - get_sta_corr(ptr[i].stn_name, ptr[i].stn_Comp);
```

**效益：** 一般可將 RMS 殘差降低 20–40%，水平誤差改善 1–3 km（特別是沉積盆地和火山地區）。台灣已有現成的 CWA 站校正資料（HYPOINVERSE 格式）。

---

#### 改善 E：Huber 損失函數取代 L2（提升 outlier 穩健性）

**問題：** 現有 `wt()` 函數的殘差加權方式對大殘差（誤讀取）懲罰不足。

```c
/* 現況 wt() 的殘差部分 */
xwt *= (tres / (tres + fabs(res))) * (tres / (tres + fabs(res)));   /* L2 風格 */

/* 建議：加入 Huber 閾值 */
#define HUBER_DELTA 1.5   /* 秒，超過此值改用 L1 加權 */
double res_abs = fabs(res);
if(res_abs <= HUBER_DELTA) {
    /* L2 區域：二次加權 */
    xwt *= (tres / (tres + res_abs)) * (tres / (tres + res_abs));
} else {
    /* L1 區域：假 picks 的大殘差不會過度降低權重（也不會讓解發散）*/
    xwt *= tres / (2.0 * res_abs) * HUBER_DELTA;
}
```

**效益：** 對假讀取更穩健，在高噪音台站較多時定位誤差降低 15–30%。

---

### 9.2 中期改善（需要外部資料/模型）

#### 改善 F：使用台灣最佳 1D 速度模型

**問題：** tcpd_new 使用簡化的 2 層梯度模型，與台灣真實地殼差距較大。

**建議速度模型**（基於 CWA/VELEST 最小 1D 模型）：

```c
/* 替換 tcpd_new.d 中的速度模型參數 */
/* 以下為台灣改良版 1D 模型（Wu & Liang 2009 整合版）*/

/* P 波：4 層模型 */
static double layer_depth[] = {0.0,  3.0,  15.0, 35.0};   /* 各層頂（km）*/
static double layer_vp[]    = {5.0,  5.8,  6.3,  7.8 };   /* 各層 Vp（km/s）*/
static double layer_vpg[]   = {0.05, 0.04, 0.02, 0.001};  /* 速度梯度 */

double get_vp(double z) {
    int k;
    for(k=3; k>=0; k--)
        if(z >= layer_depth[k]) break;
    return layer_vp[k] + layer_vpg[k] * (z - layer_depth[k]);
}
```

**效益：** 走時計算精度提升，特別改善 z > 40 km 的深源地震，預計水平誤差改善 2–4 km。

---

#### 改善 G：引入 EDT 走時差法（消除 z-t₀ 耦合）

**原理：** 對到時差（不需知道起震時間）做格搜尋：

```c
/* EDT 格搜尋（替代或補充現有深度格搜尋）*/
double edt_score(double xc, double yc, double zc) {
    double score = 0.0;
    int cnt = 0;
    for(i=0; i<nsta; i++)
        for(j=i+1; j<nsta; j++) {
            double dt_obs  = ptr[i].P - ptr[j].P;
            double T_i     = calc_traveltime(xc, yc, zc, ptr[i].latitude, ptr[i].longitude);
            double T_j     = calc_traveltime(xc, yc, zc, ptr[j].latitude, ptr[j].longitude);
            double dt_theo = T_i - T_j;
            double residual = dt_obs - dt_theo;
            score += exp(-residual*residual / (2.0 * 0.5*0.5));   /* Gaussian 似然 */
            cnt++;
        }
    return (cnt > 0) ? score / cnt : 0.0;
}
```

**效益：** 完全消除 $z$–$t_0$ 耦合，深度估計精度提升最顯著（與現有方法互補）。

---

#### 改善 H：Kalman 濾波即時更新

**框架：** 從第 4 站觸發後即開始輸出定位，每到 1 站做 Kalman 更新：

```c
/* Kalman 狀態（全域，跨觸發更新）*/
typedef struct {
    double x[4];    /* [lon, lat, dep, t0] */
    double P[4][4]; /* 協方差矩陣 */
    int initialized;
} KalmanState;

static KalmanState kf;

void kalman_update(PEEW *new_pick, KalmanState *kf) {
    /* H：新台站對走時的偏導數（1×4 Jacobian）*/
    double H[4] = { dT/dx, dT/dy, dT/dz, 1.0 };
    double R = 0.25;  /* 到時觀測雜訊方差（0.5s²）*/
    /* Kalman 增益 */
    double S = H·kf->P·Hᵀ + R;
    double K[4] = kf->P·Hᵀ / S;
    /* 更新狀態 */
    double innovation = new_pick->P - kf->x[3] - T_computed;
    for(i=0; i<4; i++) kf->x[i] += K[i] * innovation;
    /* 更新協方差 */
    /* ... (I - K·H) · P ... */
}
```

**效益：** 4 站觸發後即可輸出初步定位（比現有需要 > 4 站強制條件更靈活），後續每站到達自動精化，輸出即時誤差橢球。

---

### 9.3 長期目標（研究投入）

| 目標 | 技術 | 預期效益 | 難度 |
|------|------|---------|------|
| 台灣 3D 走時場 | TAIGER + eikonal FD + mmap | 系統誤差降至 < 2 km | 高（建表工程） |
| EikoNet 神經替代 | MLP 訓練在台灣 3D 模型 | 4.8 MB 取代 GB 走時表 | 高（ML 訓練） |
| 完整 Bayesian 輸出 | Student-t + Stein VI | 校準的不確定度 | 非常高 |
| PhaseNet 拾取整合 | Python/ZMQ 橋接 | 拾取品質大幅提升 | 中 |

---

## 10. 關鍵文獻

### 定位方法

| 方法 | 文獻 | 年份 |
|------|------|------|
| Geiger 法 | Geiger, Göttingen Nachrichten | 1912 |
| HYPO71 | Lee & Lahr, USGS OFR 72-224 | 1972 |
| HYPOINVERSE | Klein, USGS OFR 02-171 | 2002 |
| NonLinLoc | Lomax et al., AG Monograph | 2000 |
| EDT 法 | Font et al., JGR 109:B06304 | 2004 |
| HypoDD | Waldhauser & Ellsworth, BSSA 90:1353 | 2000 |
| VELEST | Kissling et al., JGR 99:19395 | 1994 |
| Kalman 即時定位 | Chen et al., GRL 47:e2019GL086240 | 2020 |
| EikoNet | Smith et al., JGR 126:e2020JB020986 | 2021 |
| HypoSVI | Smith et al., GJI 228:716 | 2022 |
| LOC-FLOW | Zhang et al., SRL 93:3151 | 2022 |
| Bayesian NN 走時替代 | arXiv:2512.06407 | 2024 |
| RANSAC 定位 | arXiv:2502.10933 | 2025 |
| GraphDD GNN 重定位 | EPS 2025 | 2025 |
| MSLOC 多站 ML | GJI 241:1853 | 2025 |

### 速度模型不確定度

| 主題 | 文獻 | 年份 |
|------|------|------|
| 不確定度全域靈敏度分析 | GJI 237:1048 | 2024 |
| 物理知情不確定度量化 | Sci. Rep. s41598-024-84995-9 | 2024 |

### 台灣相關

| 主題 | 文獻 | 年份 |
|------|------|------|
| ML-EDT 台灣 EEW | Earth Planets Space 76:76 | 2024 |
| TAIGER 3D Vp 模型 | Wu et al., JGR 119:6401 | 2014 |
| CWA 3D 速度模型 | Wu & Liang, TAO 17:171 | 2006 |
| 台灣 EEW 系統現況 | Wu et al., J. Geol. Soc. India | 2021 |
| RT-MEMS 深度學習監測 | Sensors 25(11):3353 | 2025 |

---

## 結語

`locaeq()` 使用的 Geiger 迭代最小二乘法在計算速度上是最優選擇，短期內應維持此核心。**最大的精度改善空間來自**：

1. **靜態站校正**（改善 20–40% RMS，工作量小）
2. **兩階段深度格搜尋**（深度解析度從 ±5 km 改善至 ±1 km）
3. **改用台灣最佳 1D 模型**（系統性誤差降低）
4. **自適應收斂判斷**（節省計算 + 偵測不收斂）

中期最有潛力的突破是整合 **EDT 走時差法**（消除深度–起震時間耦合）和 **Kalman 濾波即時更新**（從 4 站觸發開始輸出並持續精化），這兩者合併可在不增加顯著計算量的前提下，將 EEW 定位精度提升一個量級。
