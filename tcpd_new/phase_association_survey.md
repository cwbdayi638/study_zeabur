# 地震震相關聯方法技術調查與改善建議

**撰寫日期：** 2026-03-24
**背景：** 針對台灣 TSMIP/CWASN 地震早期預警系統（EEW）之應用
**對象系統：** tcpd_new（Earthworm EEW 定位與規模估算模組）

---

## 目錄

1. [什麼是震相關聯](#1-什麼是震相關聯)
2. [tcpd_new 目前的做法](#2-tcpd_new-目前的做法)
3. [傳統方法](#3-傳統方法)
4. [現代統計/機率方法](#4-現代統計機率方法)
5. [機器學習方法](#5-機器學習方法)
6. [2025 基準測試比較](#6-2025-基準測試比較)
7. [EEW 特殊挑戰](#7-eew-特殊挑戰)
8. [針對 tcpd_new 的改善建議](#8-針對-tcpd_new-的改善建議)
9. [台灣特殊情境建議](#9-台灣特殊情境建議)
10. [關鍵文獻](#10-關鍵文獻)

---

## 1. 什麼是震相關聯

**震相關聯（Seismic Phase Association）** 是地震監測流程中的核心步驟：

```
地震波形
    │
    ▼ 震相拾取（Picker：PhaseNet, EQTransformer...）
P 波到時序列（來自 N 個台站的 picks）
    │
    ▼ 震相關聯（Associator：本文主題）
歸屬到同一地震事件的 picks 集合
    │
    ▼ 定位與規模估算（Locator：tcpd_new 的 locaeq + Magnitude）
震央座標、深度、規模
```

關聯問題的困難在於：
- 同一時間可能有多個地震同時發生（主震 + 餘震、遠地 + 近地）
- 噪音誤觸發（風、海浪、工業噪音）產生的假 picks
- 台站分布不均（台灣東部海域無台站）
- EEW 要求必須在 **1–3 秒內**完成關聯

---

## 2. tcpd_new 目前的做法

### 2.1 演算法描述

tcpd_new 採用的是一種**單事件、單次均值距離過濾（Single-Event Centroid Filtering）**方法：

```
輸入：ptr[] 緩衝區中所有有效 P 波到時

步驟一：計算全部有效站的空間重心與時間重心
    avg_lat   = Σlat_i / N
    avg_lon   = Σlon_i / N
    avg_Parr  = ΣP_i   / N

步驟二：逐站過濾
    對每個站 i：
        d_i  = delaz(站_i 位置, 重心位置)
        dt_i = |P_i - avg_Parr|
        若 d_i > Trig_dis_win(100km) OR dt_i > Trig_tm_win(15s)
            → 剔除（不加入 vsn_ntri[]）
        否則 → 保留

步驟三：觸發條件
    若 vsn_trigger > 4 且 vsn_trigger > max_sta
        → 呼叫 processTrigger()

步驟四：定位內異常值排除（processTrigger 內）
    while(averr > 0.8 or avwei < 0.3):
        找最大殘差站 → 標記 flag=0
        若剩餘站 < 5 → 放棄
        重新 locaeq()
```

### 2.2 現有方法的限制

| 限制 | 具體問題 | 影響 |
|------|---------|------|
| **單次過濾** | 重心受離群站影響後，鄰近正常站可能被誤剔除 | 低召回率 |
| **單一事件假設** | 若同時有兩個地震，所有 picks 混在一起計算重心，兩個事件都無法正確關聯 | EEW 失效 |
| **固定閾值** | `Trig_dis_win=100km`、`Trig_tm_win=15s` 不依地震規模調整 | M>6 大震 picks 時間差可能超過 15s |
| **無統計推斷** | 無機率模型，無法量化關聯的信心水準 | 無法評估警報可靠性 |
| **均值計算含離群站** | 待過濾的離群站也參與重心計算（Masking Effect） | 過濾效果不穩定 |
| **二維時空過濾** | 只考慮水平位置與到時，不考慮深度 | 深震與淺震 picks 可能混淆 |

---

## 3. 傳統方法

### 3.1 格點搜尋法（Grid-Search Association）

**原理：**
在三維空間格點上，對每個格點計算理論 P 波走時，統計與實際 picks 吻合的台站數，找最大吻合格點作為震源假說。

```
for each candidate (x, y, z, t):
    score = 0
    for each pick P_i:
        theo_T = travel_time(x,y,z,station_i)
        if |P_i - t - theo_T| < tolerance:
            score += 1
    if score > threshold:
        → 記錄為候選事件
```

**優點：** 直觀、不需訓練資料、能處理稀疏網路
**缺點：** 計算量 O(G × N)，格點數 G 爆炸性增長；對假 picks 敏感
**EEW 適用性：** 粗格點時可用，但精度低

---

### 3.2 REAL（Rapid Earthquake Association and Location）

**核心演算法（Zhang et al., 2019, SRL）：**
兩階段關聯：
1. **讀取計數（Pick Counting）**：格點上統計落入走時窗內的 P+S picks 數
2. **走時殘差篩選**：多個格點通過閾值時，選走時殘差總和最小者

```
Phase 1: 讀取計數
  for each (x, y, z):
    n_P = count(|P_observed - t0 - T_P(x,y,z)| < dt_P)
    n_S = count(|S_observed - t0 - T_S(x,y,z)| < dt_S)
    if n_P >= N_P_min AND (n_P + n_S) >= N_total_min:
      → 候選事件

Phase 2: 殘差最小化
  best = argmin Σ|P_i - t0 - T_P(xi,yi,zi)| over candidates
```

**優點：** 計算效率比純格點搜尋高；已在全球多個地震序列驗證；開源（GitHub: Dal-mzhang/REAL）
**缺點：** 2025 基準測試中速度慢於 PyOcto 約 10 倍；高讀取率下召回率下降
**EEW 適用性：** 有限；計算效能不足以應付高速串流

---

### 3.3 HYPOINVERSE-2000

**本質：** 震源定位工具，**非**自動關聯工具
在 tcpd_new 流程中對應的角色是 `locaeq()`，而非前置的關聯步驟
通常需搭配 STA/LTA 觸發或人工篩選後才能使用
**EEW 適用性：** 不直接適用

---

## 4. 現代統計/機率方法

### 4.1 GaMMA（Gaussian Mixture Model Association）

**核心演算法（Zhu et al., 2022, JGR）：**
將震相關聯視為**無監督叢集問題**，使用**期望最大化（EM）演算法**求解。

每個地震對應一個高斯混合模型分量，每個 pick 被軟分配給各分量：

```
模型：
  每個事件 k 的 picks 在 (時間, 振幅) 空間服從高斯分佈：
    P(pick_i | event_k) ∝ exp(-½[(t_i - μ_t_k)²/σ_t² + (A_i - μ_A_k)²/σ_A²])

  其中：
    μ_t_k = t0_k + T(x_k, y_k, z_k, station_i)   ← 理論走時
    μ_A_k = 振幅衰減關係估計值                      ← 含規模資訊

EM 演算法：
  E-step：計算每個 pick 屬於各事件的後驗機率
  M-step：用後驗機率加權更新震源參數
  重複直到收斂
```

**優點：**
- 同時輸出震源位置、發震時間與**規模估計**（對 EEW 有附加價值）
- 不需要訓練資料
- 在 QuakeFlow 雲端流程中已有串流處理模式
- 可處理時空密集的地震序列（已用於 2019 Ridgecrest 序列）

**缺點：**
- 高噪音、高密度情境下準確度下降
- EM 演算法對初始化敏感，需預設最大事件數
- 2025 基準測試表現遜於 GENIE 和 PyOcto

**EEW 適用性：** 串流模式可達準即時；延遲約 1–5 秒

---

### 4.2 BAYESLOC（LLNL，Myers et al., 2007）

**原理：** 貝葉斯層級模型（Bayesian Hierarchical Model），以 MCMC 採樣求解多事件聯合定位
**用途：** 核爆監測、事後目錄精化
**EEW 適用性：** 完全不適用（計算量龐大，MCMC 需分鐘至小時）

---

## 5. 機器學習方法

### 5.1 PhaseNet（Zhu & Beroza, 2019, GJI）

**性質：** 震相**拾取**工具（Picker），非關聯工具，但是現代關聯流程的前置步驟

**架構：** 修改版 U-Net，輸入三分量波形，輸出 P/S 機率時間序列

```
輸入：[E, N, Z] 波形（3 × T 陣列）
    ↓  U-Net 編碼器（下採樣 + 特徵提取）
    ↓  U-Net 解碼器（上採樣 + 跳接連結）
輸出：[P 機率, S 機率, 噪音機率]（3 × T 陣列）
```

**性能：**
P 波 F1 = 0.896，S 波 F1 = 0.801（vs 傳統 AR picker：0.558 / 0.165）

**EEW 適用性：** 延遲 < 1 秒；已在多個即時系統部署
**與 tcpd_new 的關係：** tcpd_new 的輸入來自 pick_eew 模組，若以 PhaseNet 取代傳統 STA/LTA 拾取，可直接提升輸入品質

---

### 5.2 PhaseLink（Ross et al., 2019, JGR）

**架構：** LSTM 循環神經網路，輸入多台站 picks 序列，輸出各 picks 對的同源機率

```
輸入（滑動時窗內）：
  [ta, ta_type, lat_a, lon_a, tb, tb_type, lat_b, lon_b, ...]

LSTM 處理序列資訊
↓
輸出：P(picks i, j 來自同一事件)
```

**訓練資料：** 以一維速度模型生成數千萬組合成讀取序列（無需真實目錄）

**缺點：** 2025 基準測試在高噪音與隱沒帶情境下性能急劇下降，最惡劣情境接近 F1=0
**EEW 適用性：** 架構支援即時，但高噪音環境可靠性不足

---

### 5.3 GENIE（Graph Earthquake Neural Interpretation Engine）

**核心演算法（McBrearty & Beroza, 2023, BSSA）：**
使用**圖神經網路（GNN）**，同時進行震相關聯與震源定位

```
架構：兩個 K-NN 圖

台站圖（Station Graph）：
  節點 = 台站，邊 = K 個最近鄰台站
  節點特徵：[台站位置, 到時, P/S 機率, 振幅]

震源空間圖（Source Graph）：
  節點 = 候選震源格點
  邊 = K 個最近鄰格點

訊息傳遞（Message Passing）：
  台站圖 → 震源圖：將 picks 資訊投影到震源空間
  震源圖 → 台站圖：以震源假說反向驗證 picks
  迭代 L 次

輸出：
  - 每個格點是否為事件的機率
  - 每個 pick 歸屬於哪個事件的機率
```

**優點：**
- 同時輸出關聯與定位，端到端處理
- 對台網幾何變化具有強健性（GNN 自動學習幾何關係）
- 可擴展至全球規模（北加州 100 天連續處理：偵測 USGS 目錄 4 倍事件數）
- **2025 基準測試：幾乎所有情境近乎完美，最惡劣 F1 > 0.8**

**缺點：**
- 需要訓練資料（合成 + 真實資料混合）
- 遷移到新地區需要重新訓練（或微調）
- 部署需要 GPU

**EEW 適用性：** 架構支援準即時；需要 GPU 加速

---

### 5.4 PyOcto（Münchmeyer, 2024, Seismica）

**核心演算法：** 四維時空分割（4D Space-Time Partitioning）

```
核心思想：OcTree 遞迴空間分割

1. 建立走時查找表（travel time lookup table）
2. 以 OcTree 將 (x, y, z, t) 空間分割為區塊
3. 對每個時空區塊，快速計算：
   n_picks_consistent = count(picks 落入理論走時窗內)
4. n > threshold 的區塊進一步細分（遞迴）
5. 葉節點即為候選事件

複雜度：O(N × log(G))，G 為格點數（比 O(N×G) 快得多）
```

**開源：** GitHub: yetinam/pyocto（PyPI 可直接安裝）

**優點：**
- **速度比 REAL 快 ≥ 10 倍**；2014 年 Iquique 序列加速因子 ≥ 70
- 不需要訓練資料，新地區立即可用
- 穩定的高精準率（precision）
- **2025 基準測試：與 GENIE 並列最佳**
- 小規模資料集（EEW 典型情境）計算速度最快

**缺點：**
- 依賴一維速度模型假設
- 純讀取值輸入，不使用波形振幅（無內建規模估算）

**EEW 適用性：** ★★★★★ 目前最適合 EEW 即時應用的關聯方法

---

### 5.5 EQNet / QuakeFlow（Zhu et al., 2022）

**EQNet：** 端對端架構，直接從波形輸出事件偵測，不需獨立關聯步驟：
```
波形 → ResNet-18 特徵提取 → Shift-and-Stack → 事件偵測
```

**QuakeFlow：** 雲端流程，Docker + Kubernetes 部署：
```
PhaseNet（picks）→ GaMMA（關聯+定位）→ 後處理
```

已在 2022 年池上地震（Mw 6.9）序列建立高品質目錄。

---

### 5.6 HARPA（High-Rate Phase Association）

**核心：** 深度生成模型 + 走時神經場（travel time neural fields）
使用最優傳輸（Optimal Transport）度量比較到時分佈
特別適合超高事件率情境（如密集餘震序列）
**狀態：** 2023–2024 年預印本，實際部署案例尚少

---

## 6. 2025 基準測試比較

基於 Münchmeyer et al. (2025)，*Seismica*（arXiv:2501.03621）

### 測試情境

| 情境 | 描述 |
|------|------|
| 標準 | 正常台網密度、低噪音率 |
| 高噪音 | 30–50% 假讀取率 |
| 高事件率 | 密集餘震序列 |
| 隱沒帶 | 深震 + 複雜速度結構 |

### 性能比較

| 方法 | 類型 | 標準 F1 | 高噪音 F1 | 高事件率 F1 | 速度 | 訓練需求 | EEW 適用 |
|------|------|---------|----------|------------|------|---------|---------|
| **PyOcto** | 時空分割 | ~1.00 | >0.85 | >0.85 | ★★★★★ | 不需要 | ★★★★★ |
| **GENIE** | GNN | ~1.00 | >0.85 | >0.85 | ★★★★（GPU）| 合成+真實 | ★★★★ |
| **GaMMA** | EM/GMM | 良好 | 中等 | 困難 | ★★★ | 不需要 | ★★★ |
| **REAL** | 計數+殘差 | 良好 | 尚可 | 召回率下降 | ★★（慢10×）| 不需要 | ★★ |
| **PhaseLink** | LSTM | 尚可 | 差 | 極差 | ★★★★（小規模）| 合成 | ★★ |
| **tcpd_new 現況** | 均值距離 | 中等 | 差 | 完全失效 | ★★★★★ | 不需要 | ★★★ |

> **tcpd_new 現況** 是作者對本系統的評估，非基準測試直接測量。

---

## 7. EEW 特殊挑戰

### 7.1 延遲需求

```
地震發生
  │  telemetry delay（1–3 s）
  ▼
picks 到達
  │  phase association（目標：< 2 s）
  ▼
location + magnitude（< 1 s）
  │
  ▼
警報發出（目標：地震後 10–17 s 內）
```

台灣 CWA 目前系統：地震後約 **10–13 秒**發出首報
**瓶頸**：telemetry 延遲通常大於演算法本身

### 7.2 台灣特殊挑戰

| 挑戰 | 說明 |
|------|------|
| **東部海域無台站** | 宜蘭外海、台東外海為主要震源區，幾乎無台站 |
| **複雜速度結構** | 西部沉積盆地 vs 中央山脈（速度差異 > 30%），1D 模型系統誤差大 |
| **主震＋餘震序列** | 2024 花蓮 Mw 7.4 後數小時內多個 M>5 餘震，傳統單事件假設完全失效 |
| **海浪/強風噪音** | 沿海台站（台東、花蓮、蘇澳）假讀取率高 |
| **台網幾何不對稱** | 東部台站密度遠低於西部，方位角覆蓋不均 |

### 7.3 多個同時地震

| 方法 | 多事件處理能力 |
|------|-------------|
| tcpd_new（現況） | 完全無法（單事件假設） |
| REAL | 差（格點可能混淆） |
| GaMMA | 好（EM 設計支援多叢集） |
| GENIE | 最好（GNN 同時輸出多震源） |
| PyOcto | 好（OcTree 可找多個候選） |

---

## 8. 針對 tcpd_new 的改善建議

### 8.1 短期改善（可直接修改現有程式碼）

#### 改善 A：迭代式群聚過濾（解決 Masking Effect）

**問題：** 當前做法用所有站（含離群站）計算重心，離群站會拉偏重心，導致附近正常站被誤剔除。

**改善方案：**

```c
/* 現況：單次過濾 */
avg_lon = sum_lon / num_tt;
avg_lat = sum_lat / num_tt;
// 過濾一次，done

/* 建議：迭代式過濾直到收斂 */
int changed = 1;
while(changed) {
    changed = 0;
    // 重新計算重心（只用目前標記為有效的站）
    sum_lat = 0.0; sum_lon = 0.0; num_valid = 0;
    for(i=0; i<number; i++)
        if(ptr[i].flag > 0) {
            sum_lat += ptr[i].latitude;
            sum_lon += ptr[i].longitude;
            num_valid++;
        }
    if(num_valid <= 0) break;
    avg_lat = sum_lat / num_valid;
    avg_lon = sum_lon / num_valid;

    // 過濾：若有新的站被剔除，設 changed=1
    for(i=0; i<number; i++) {
        if(ptr[i].flag > 0) {
            dis = delaz(ptr[i].latitude, ptr[i].longitude, avg_lat, avg_lon);
            dt  = fabs(ptr[i].P - avg_ptime);
            if(dis > Trig_dis_win || dt > Trig_tm_win) {
                if(!is_exempt_station(ptr[i].stn_name)) {
                    ptr[i].flag = 0;
                    changed = 1;
                }
            }
        }
    }
}
/* 最多迭代 5 次防止死循環 */
```

**預期效果：** 提升過濾穩健性，降低正常站誤剔除率約 20–40%（估算）

---

#### 改善 B：動態觸發時間窗（依地震規模調整）

**問題：** 固定 `Trig_tm_win = 15s` 對大規模地震（M≥6）可能過緊，因大地震的 P 波可能跨越更大範圍。

**改善方案：**

```c
/* 在計算群聚過濾前，依當前觸發站的時間跨度動態調整窗口 */
double P_span = P_max - P_min;   // 最晚與最早 P 到時差
double dynamic_trig_tm = Trig_tm_win;

// 若觀測到的 P 時間跨度已超過設定窗口的 80%
// 可能是大地震，放寬時間窗
if(P_span > Trig_tm_win * 0.8 && vsn_trigger >= 5) {
    dynamic_trig_tm = P_span * 1.5;   // 寬鬆 50%
    logit("o", "tcpd: Dynamic Trig_tm_win=%.1fs (P_span=%.1fs)\n",
          dynamic_trig_tm, P_span);
}
```

---

#### 改善 C：加入 P 波振幅一致性檢查

**問題：** 當前只用到時和位置做關聯，未使用振幅資訊。

**改善方案：** 在群聚過濾階段加入振幅一致性檢查：

```c
/* 計算所有有效站的中位數 Pa（P 波加速度振幅）*/
double median_Pa = compute_median_Pa(ptr, number);
double Pa_ratio;

for(i=0; i<number; i++) {
    if(ptr[i].flag > 0 && ptr[i].Pa > 0) {
        // 理論振幅衰減：與震央距成反比
        // 若觀測振幅偏離理論值超過 2 個數量級 → 疑為異常站
        dis = delaz(ptr[i].latitude, ptr[i].longitude, avg_lat, avg_lon);
        double expected_ratio = dis / avg_dis;  // 簡化距離修正
        Pa_ratio = ptr[i].Pa / (median_Pa / expected_ratio);
        if(Pa_ratio > 100.0 || Pa_ratio < 0.01) {
            logit("o", "tcpd: Amplitude outlier: %s Pa=%.6f\n",
                  ptr[i].stn_name, ptr[i].Pa);
            ptr[i].flag = 0;   // 疑似錯誤讀取
        }
    }
}
```

---

### 8.2 中期改善（需要架構調整）

#### 改善 D：整合 PyOcto 作為前置關聯器

**方案：** 在 Earthworm 環境中，新增一個 Python 橋接模組，使用 PyOcto 作為關聯前置步驟：

```
PICK_RING（raw picks）
    │
    ▼ 新增：pyocto_bridge 模組（Python）
    ├── 收集 2 秒內所有 picks
    ├── 呼叫 PyOcto 關聯
    └── 輸出：每個關聯事件的 picks 子集

PICK_RING_ASSOC（已關聯 picks）
    │
    ▼ tcpd_new（現有定位+規模模組）
    └── 只需做定位，不需再做群聚過濾
```

**PyOcto 配置（針對台灣）：**

```python
from pyocto import OctoAssociator
import pyocto

velocity_model = pyocto.VelocityModel1D(
    path="taiwan_1d_model.csv",   # 台灣 1D 速度模型
    phase="P",
    tolerance=2.0,                # 走時殘差容忍度（秒）
)

associator = OctoAssociator(
    velocity_model=velocity_model,
    lat=center_lat,
    lon=center_lon,
    depth=(-2, 300),             # 深度範圍（km）
    time=60.0,                   # 時間窗口（秒）
    n_picks=5,                   # 最少觸發站數
    n_p_picks=3,                 # 最少 P 波站數
    pick_match_tolerance=1.5,    # 關聯走時容忍度（秒）
)

# 即時處理
for batch_picks in stream_picks(PICK_RING):
    events = associator.associate(batch_picks)
    for event in events:
        send_to_PICK_RING_ASSOC(event.picks)
```

**優點：** PyOcto 速度最快，不需訓練，可立即部署
**工作量：** 中等（需實作 Python/C 橋接層）

---

#### 改善 E：多事件並行處理

**問題：** 現有架構假設同一時間只有一個地震，若同時有兩個事件（如主震後餘震），系統只能追蹤其中之一。

**改善方案：** 維護多個並行的 `ptr[]` 緩衝區，每個對應一個潛在事件：

```c
#define MAX_EVENTS 5

PEEW ptr[MAX_EVENTS][number];   // 最多同時追蹤 5 個事件
int  event_active[MAX_EVENTS];  // 此事件是否活躍
int  num_active = 0;

// 新 picks 到達時，嘗試分配到最匹配的事件
// 若無匹配事件 → 開啟新事件
int assign_to_event(PEEW *new_pick) {
    double best_score = -1;
    int    best_event = -1;
    for(e=0; e<MAX_EVENTS; e++) {
        if(!event_active[e]) continue;
        score = compute_association_score(new_pick, ptr[e]);
        if(score > best_score) {
            best_score = score;
            best_event = e;
        }
    }
    if(best_score < ASSOCIATION_THRESHOLD)
        best_event = open_new_event();  // 開新事件
    return best_event;
}
```

**工作量：** 大（需要重構核心資料結構）

---

### 8.3 長期改善（需要研究投入）

#### 改善 F：以 PhaseNet 取代傳統 P 波拾取

**現況：** tcpd_new 的輸入來自傳統 STA/LTA 或 pick_eew 模組
**建議：** 在 Earthworm 中整合 PhaseNet 拾取：

```
波形串流（Wave Ring）
    │
    ▼ PhaseNet（Python + TensorFlow/PyTorch）
    ├── 延遲 < 1 秒
    ├── P 波 F1 = 0.896（vs 傳統 0.558）
    └── 輸出機率值（可用作 pick weight）

PICK_RING（高品質 picks）
    │
    ▼ tcpd_new
```

**參考：** 台灣 RT-MEMS 系統已整合類 PhaseNet 架構（SeisBlue），在台灣地震序列偵測率達 CWA 目錄 3 倍以上（Sensors, 2025）

---

#### 改善 G：以 3D 速度模型改善走時計算

**現況：** tcpd_new 使用 2 層梯度 1D 速度模型
**建議：** 引入 Lin et al. (2014) 台灣 3D Vp 模型

```c
/* 修改 locaeq() 內的走時計算 */

// 現況：解析公式（適用 1D 梯度模型）
ptime = (-1./vg) * log(fabs(tan(ang2/2.) / tan(ang1/2.)));

// 建議：預先計算 3D 走時查找表（Travel Time Table）
// 使用 SeisArray 或 NonLinLoc 生成格點走時表
// 執行時做三線性插值
ptime = interpolate_3D_TT(x0, y0, z0, ptr[i].longitude, ptr[i].latitude);
```

**預期效果：** 西部平原定位誤差可從 ~5 km 降至 ~2 km（估算，基於 Lin et al. 研究）

---

## 9. 台灣特殊情境建議

### 9.1 東部海域地震

**問題：** 台灣東部外海為主要震源區，但幾乎無海底台站，導致：
- 方位角覆蓋嚴重不足（Gap > 270°）
- 定位偏差（向有台站的西側偏移）
- EEW 延遲增加

**建議：**
1. **整合海底光纖 DAS 資料**：台灣東部海纜（如 TAIGER 海底地震儀）可提供 DAS 感測，dEPIC 框架（Scientific Reports, 2025）已示範以 GPU 加速的 DAS+陸上台網整合 EEW
2. **On-site 方法備援**：當 Gap > 200° 時，自動切換到單台站 EEW 估算（P-alert 系統）
3. **方位角缺口懲罰**：在觸發邏輯中加入 Gap 閾值：

```c
// 建議：Gap > 250° 時降低警報信心等級
if(hyp.gap > 250.0) {
    logit("o", "tcpd: Large azimuthal gap=%.0f, location uncertainty HIGH\n", hyp.gap);
    // 輸出降級警報或附加不確定性標記
}
```

### 9.2 高噪音沿海台站

**建議：** 對靠近海岸的台站（台東、花蓮、蘇澳）設定更嚴格的 pick weight 閾值，或加入動態黑名單機制（連續假觸發超過閾值次數時暫時排除）：

```c
/* 台站動態信任分數 */
double sta_trust[MAXSTA];  // 初始化為 1.0

// 每次關聯後更新信任分數
// 若某台站的殘差持續 > 2σ → 降低信任分數
// 信任分數 < 0.3 時不納入關聯
```

### 9.3 即時建議優先級

| 改善項目 | 實施難度 | 預期效益 | 優先級 |
|---------|---------|---------|--------|
| A：迭代群聚過濾 | 低（修改現有 C 程式） | 中（穩健性提升） | ★★★★★ |
| B：動態觸發時間窗 | 低 | 中（大震改善） | ★★★★ |
| C：振幅一致性檢查 | 低-中 | 中（假讀取過濾） | ★★★★ |
| D：PyOcto 橋接 | 中（Python 橋接） | 高（根本性改善） | ★★★★ |
| E：多事件並行 | 高（架構重構） | 高（餘震序列） | ★★★ |
| F：PhaseNet 整合 | 中（需 ML 環境） | 高（拾取品質） | ★★★★ |
| G：3D 速度模型 | 高（需研究資源） | 高（系統性誤差） | ★★★ |
| DAS 整合 | 非常高（硬體建設） | 非常高（東部覆蓋） | ★★★★ |

---

## 10. 關鍵文獻

### 方法論文

| 方法 | 文獻 | 年份 | 連結 |
|------|------|------|------|
| REAL | Zhang et al., SRL 90(6):2276 | 2019 | GitHub: Dal-mzhang/REAL |
| PhaseNet | Zhu & Beroza, GJI 216:261 | 2019 | GitHub: AI4EPS/seisbench |
| PhaseLink | Ross et al., JGR 124:3082 | 2019 | GitHub: interseismic/PhaseLink |
| GaMMA | Zhu et al., JGR 127:e2021JB023249 | 2022 | GitHub: AI4EPS/GaMMA |
| QuakeFlow | Zhu et al., GJI 232(1):684 | 2022 | GitHub: AI4EPS/QuakeFlow |
| EQNet | Zhu et al., JGR 127:e2021JB023283 | 2022 | arXiv:2109.09911 |
| GENIE | McBrearty & Beroza, BSSA 113(2):524 | 2023 | GitHub: imcbrearty/GENIE |
| HARPA | arXiv:2307.07572 | 2023 | arXiv |
| PyOcto | Münchmeyer, Seismica | 2024 | GitHub: yetinam/pyocto |
| **基準測試** | arXiv:2501.03621, Seismica | **2025** | arXiv |

### 台灣相關

| 主題 | 文獻 | 年份 |
|------|------|------|
| 台灣 EEW 系統現況 | Wu et al., J. Geol. Soc. India | 2021 |
| CWA 地震資料集 | Tang et al., SRL | 2024 |
| RT-MEMS 深度學習監測 | Sensors 25(11):3353 | 2025 |
| DAS 整合 EEW（dEPIC）| Scientific Reports | 2025 |
| ML+EDT 台灣 EEW | Earth Planets Space | 2024 |
| 2025 ML 6.4 大埔即時監測 | The Seismic Record 5(3):320 | 2025 |

### 速度模型

| 主題 | 文獻 |
|------|------|
| 台灣 3D Vp 模型 | Lin et al., JGR 119(8):6401 (2014) |
| 台灣地殼速度結構 | Kuo-Chen et al., GRL 39:L02313 (2012) |

---

## 結語

tcpd_new 目前採用的群聚過濾方法屬於**第一代關聯方法**，在正常情況下運作良好，但在**多個同時事件**、**高噪音**、**大規模地震**等情境下有顯著限制。

基於 2025 年最新基準測試，**PyOcto** 是目前最適合 EEW 即時應用的關聯工具，兼具速度（最快）、精準率（最高之一）和無訓練需求等優點。短期內最具實效的改善是對現有群聚過濾邏輯進行**迭代化改良（改善 A）**，並在中期整合 **PyOcto 作為前置關聯器（改善 D）**，這兩步即可大幅提升台灣 EEW 系統對複雜地震情境的應對能力。
