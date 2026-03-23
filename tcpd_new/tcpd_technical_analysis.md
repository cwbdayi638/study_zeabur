# tcpd_new.c — 深度技術分析與演算法改善方向

**程式版本：** v2（優化版）
**分析日期：** 2026-03-23
**原始碼行數：** 2986 行
**作者：** Prof. Yih-Min Wu（吳逸民教授）定位核心；EEW 系統整合

---

## 目錄

1. [系統定位與架構](#1-系統定位與架構)
2. [主程式流程（main）](#2-主程式流程-main)
3. [P 波緩衝區管理](#3-p-波緩衝區管理)
4. [群聚過濾演算法](#4-群聚過濾演算法)
5. [震源定位演算法（locaeq）](#5-震源定位演算法-locaeq)
6. [規模估算演算法（Magnitude）](#6-規模估算演算法-magnitude)
7. [報告產出（Report_seq）](#7-報告產出-report_seq)
8. [輔助數學函式](#8-輔助數學函式)
9. [速度模型](#9-速度模型)
10. [衰減關係](#10-衰減關係)
11. [未來可改善的演算法方向](#11-未來可改善的演算法方向)

---

## 1. 系統定位與架構

### 1.1 在 Earthworm 框架中的角色

```
外部測站
   │  地震波形（即時串流）
   ▼
pick_eew / pick_ew          ← P 波拾取模組
   │  TYPE_EEW 訊息 → PICK_RING
   ▼
tcpd_new（本程式）          ← 定位與規模估算
   │  TYPE_EEW 訊息 → EEW_RING
   ▼
EEW 分發模組                ← 警報發布
```

`tcpd_new` 是 EEW 管線中最核心的計算模組。它從共享記憶體環（PICK_RING）消費 P 波到時資料，執行震源定位與規模估算，再將結果寫回共享記憶體環（EEW_RING）並存成 `.rep` 報告檔。

### 1.2 程式架構圖

```
main()
├── 初始化
│   ├── 讀取 num_eew_status（事件流水號）
│   ├── 設定預設參數（硬編碼初始值）
│   ├── tcpd_config()  ← 讀取 .d 設定檔（覆蓋預設值）
│   └── tcpd_lookup()  ← 查詢 Earthworm 訊息類型 ID
│
└── while(1) 主循環
    ├── 心跳機制（每 15s）← tcpd_status()
    ├── 讀取訊息 ← tport_getmsg()
    │   ├── 格式驗證（14 欄位）
    │   ├── P 波品質過濾
    │   ├── 緩衝區管理（ptr[1000]）
    │   ├── 去重與 Pd 合併
    │   ├── 群聚過濾 → vsn_ntri[]
    │   └── 觸發判斷 → processTrigger()
    │       ├── locaeq()       ← 迭代最小二乘定位
    │       ├── Magnitude()    ← Mpd + Mtc 規模估算
    │       └── Report_seq()   ← 寫 .rep + 發送 EEW_RING
    └── 60s 超時重置（等待下次地震）
```

### 1.3 記憶體資料流

| 結構 | 大小 | 生命週期 | 用途 |
|------|------|---------|------|
| `ptr[1000]` | 全域 | 跨事件持續 | 全部已知 P 波到時（含過期尚未清除者） |
| `vsn_ntri[1000]` | 全域 | 每次訊息處理重建 | 群聚過濾後的有效觸發站 |
| `new_ntri[1000]` | processTrigger 本地 | 每次定位週期 | 去異常後的定位用站列 |
| `PEEW.Pd[15]` | 每站 15 個 double | 隨 upd_sec 累積 | P 波位移振幅時間序列 |

---

## 2. 主程式流程（main）

### 2.1 初始化邏輯

程式啟動時先讀取 `num_eew_status` 檔案取得事件流水號，使得重啟後能延續上次計數：

```c
fp = fopen("num_eew_status", "r");
// 若檔案不存在或內容損壞 → 重設為 1
// 有效範圍：1 – 10000，超過則循環
num_eew = atoi(tmp);
```

參數優先順序：**硬編碼預設值** → **tcpd_config() 設定檔覆蓋**。這意味著若設定檔缺少某項參數，程式仍能以預設值運行（但預設值與設定檔值可能相差甚遠，需特別注意）。

### 2.2 心跳機制

Earthworm 的 `startstop` 守護程序透過心跳監控所有模組存活：

```c
if (time(&timeNow) - timeLastBeat >= HeartBeatInterval) {
    timeLastBeat = timeNow;
    tcpd_status(TypeHeartBeat, 0, "");
}
```

若 `tcpd_new` 超過 `HeartBeatInterval`（15 秒）未發送心跳，`startstop` 會自動重啟該模組。這是 Earthworm 系統的核心可靠性機制。

### 2.3 觸發計數器邏輯

程式維護三個關鍵計數器來控制報告節奏：

| 變數 | 作用 |
|------|------|
| `vsn_trigger` | 本次處理週期的有效觸發站數 |
| `max_sta` | 上次產出報告時的觸發站數 |
| `count_max` | 連續沒有新觸發的計數 |

決策邏輯：
```
vsn_trigger > max_sta  →  新觸發 → processTrigger()，count_max = 0
vsn_trigger == max_sta →  無變化 → 跳過
vsn_trigger < max_sta  →  站數減少（過期）→ count_max++，跳過
count_max > Report_Limit → 事件結束，重置全部狀態
```

`Report_Limit = 1`（硬編碼），意味著觸發站數一旦停止增加，只允許再嘗試 1 次後就認定事件結束。

---

## 3. P 波緩衝區管理

### 3.1 TYPE_EEW 訊息格式

每筆訊息為 14 個空白分隔欄位：

```
YULB  HHZ  TW  01  121.2971  23.3924  0.006514  0.000187  0.000074  4.297026  1328085502.958  0  2  1
[sta] [chn][net][loc][ lon ]  [ lat ]  [  Pa  ]  [  Pv  ]  [  Pd  ]  [ Tc  ]  [ P_epoch   ] [w][inst][upd]
```

欄位 13（`upd_sec`）是關鍵：每秒更新一次，指示 `Pd[upd_sec]` 為最新位移振幅。此設計允許規模估算隨時間累積更多 P 波能量資訊。

### 3.2 Pd 陣列累積機制

```c
ptr[j].Pd[upd_sec] = atof(out_ss[8]);   // 寫入當秒的 Pd 值
ptr[j].upd_sec     = upd_sec;           // 更新最新秒數
```

去重時，若同一 SCNL 有兩筆記錄（`fabs(P[i]-P[k]) < 1e-6`），保留較新的 `upd_sec`，並將舊記錄的較早 Pd 值複製進來：

```c
if(ptr[i].upd_sec > ptr[k].upd_sec) {
    for(jj=2; jj<=ptr[k].upd_sec; jj++)
        ptr[i].Pd[jj] = ptr[k].Pd[jj];  // 合併歷史 Pd
    ptr[k].flag = 0;                     // 刪除舊記錄
}
```

這確保了每個測站的 `Pd[2..upd_sec]` 是完整的時間序列，用於後續選用最佳 Pd 秒數。

### 3.3 有效期管理

```c
if(fabs(now_time - ptr[i].P) > Active_parr_win)  // 80 秒
    ptr[i].flag = 0;   // 標記為無效（下次可被覆寫）
```

`Active_parr_win = 80s` 的設計考量：臺灣島長約 400 km，P 波以 ~6 km/s 傳播，全島傳遞時間約 67 秒。80 秒的窗口確保所有測站的 P 波到時都在有效期內。

---

## 4. 群聚過濾演算法

### 4.1 設計原理

群聚過濾的目的是排除「偶發雜訊」或「遠距地震」造成的誤觸發站。核心思想：真實地震的觸發站在**空間上集中**、在**時間上接近**。

### 4.2 過濾邏輯

```
步驟一：計算所有有效站的均值
    avg_lat = Σlat / num_tt
    avg_lon = Σlon / num_tt
    avg_P   = ΣP   / num_tt

步驟二：逐站判斷
    dis  = delaz(站位置, 均值位置)     # 球面距離
    dt   = |P_到時 - avg_P|

    if (dis > Trig_dis_win OR dt > Trig_tm_win):
        → 剔除此站（不納入 vsn_ntri[]）
    else:
        → 保留（納入 vsn_ntri[]）

步驟三：例外清單（直接保留，不做過濾）
    EOS*, YOJ, JMJ, PCY, LAY, TWH,
    PNG, PHU, PTM, PTT, VCH, VWU, WLC
```

**注意**：均值是所有有效站（含可能要被剔除的站）計算的，並非迭代式。此設計在偏遠雜訊站存在時會偏移均值，可能影響過濾效果。

### 4.3 觸發閾值

| 參數 | 值 | 物理意義 |
|------|-----|---------|
| `Trig_dis_win` | 100 km | 台灣島寬約 150 km，設 100 km 確保核心站群密集 |
| `Trig_tm_win` | 15 s | P 波以 6 km/s 傳播 100 km 需 ~17 s，略偏保守 |

---

## 5. 震源定位演算法（locaeq）

### 5.1 演算法概述

`locaeq()` 使用**吳逸民教授開發的迭代最小二乘法（Iterative Least-Squares）**，在梯度速度模型下求解震源位置 (x₀, y₀, z₀, t₀)。

### 5.2 速度模型

採用**線性梯度（Linear Gradient）P 波速度模型**：

```
v(z) = v₀ + vg × z
```

| 深度範圍 | v₀ (km/s) | vg (km/s/km) |
|---------|-----------|--------------|
| z < 40 km（淺層） | 5.10298 | 0.06659 |
| z ≥ 40 km（深層） | 7.80479 | 0.00457 |

在此速度模型下，**走時解析解**為（Wiechert-Herglotz 公式）：

```
t(r) = (-1/vg) × ln[ tan(θ₂/2) / tan(θ₁/2) ]
```

其中：
- `xc = (r² - 2v₀/vg × z - z²) / (2r)`  ← 轉折點水平距離
- `zc = -v₀/vg`                            ← 射線轉折點深度
- `θ₁ = π - arctan((z-zc)/xc)`            ← 出射角
- `θ₂ = arctan(-zc/(r-xc))`              ← 入射角
- `r` = 震央距（km）

### 5.3 初始震源估計

```c
// 選用早於平均 P 到時 0.5s 的測站做初始位置
for(i=0; i<nsta; i++)
    if(ptr[i].P < (ave_parr - 0.5))
        { 累加 lon, lat, P }

t0 = 平均P到時 - 2.0   // 假設 P 波從震源到最近站約 2 秒
x0 = 平均經度
y0 = 平均緯度
z0 = 10.0 km           // 初始深度假設 10 km
```

### 5.4 迭代線性化（10 次）

每次迭代：

**1. 計算偏導數**（走時對位置的梯度）：

```
∂T/∂x = -sin(θ₁)/(v₀+vg×z) × Δx/r
∂T/∂y = -sin(θ₁)/(v₀+vg×z) × Δy/r
∂T/∂z = -cos(θ₁)/(v₀+vg×z)
∂T/∂t = 1.0
```

**2. 殘差加權最小二乘**：

```
建立方程組：B × Δ = y
  B[j][k] = Σᵢ (∂T/∂xⱼ)(∂T/∂xₖ) × wᵢ
  y[j]    = Σᵢ (∂T/∂xⱼ) × residᵢ × wᵢ
  residᵢ  = P_observed[i] - t0 - T_computed[i]

Δ = B⁻¹ × y    (4×4 矩陣求逆)

x0 += Δ[0]/dlon    (換算回經度)
y0 += Δ[1]/dlat    (換算回緯度)
z0 += Δ[2]
t0 += Δ[3]
```

**3. 加權函式 wt()**：

```c
H = √(depth² + dist² + ε)   // 震源距
if (H > 100 km):
    xwt *= (600-100) / (9×H + 600 - 10×100)   // 遠站降權
xwt *= (1/(1+|res|))²                           // 殘差大者降權
```

距離越遠、殘差越大的測站權重越低，提升定位穩健性。

**4. 發散保護**：
```c
if (|t0 - t_ini| > 60s OR errsum > 6.0)
    → 重置回初始值（防止發散）
```

### 5.5 深度格搜（Grid Search）

10 次迭代後，在固定水平位置下，測試 10 個深度候選：

```c
z_candidates = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100] km
```

對每個候選深度計算**加權走時殘差總和**，選最小者為最終深度。

### 5.6 P-S 時差計算

定位完成後，以最終震源計算每站的理論 P-S 時差：

```c
ptr[i].P_S_time = stime - ptime   // 單位：秒
```

此值用於後續規模估算時，判斷 Pd 時間序列是否已受 S 波污染。

### 5.7 非同址站計數

剔除站間距 < 0.5 km 的「同址站」（同一地點的不同儀器），計算真正獨立方位的測站數 `G_sta_num`。觸發報告門檻為 `G_sta_num >= 5`。

---

## 6. 規模估算演算法（Magnitude）

### 6.1 雙重規模方法

程式同時估算兩種規模並取平均：

| 方法 | 指標 | 原理 | 優點 |
|------|------|------|------|
| **Mpd** | P 波最大位移 (Pd) | 峰值振幅衰減關係 | 到站後 1-3 秒即可估算 |
| **Mtc** | 卓越週期 (Tc) | 頻率內容對應震級 | 對大規模地震較穩定 |

整合規模：
```c
if (Mpd > 0 && Mtc > 0):
    Mall = (Mpd + Mtc) / 2.0
elif Mpd > 0:
    Mall = Mpd
```

### 6.2 Pd 選用邏輯（防 S 波污染）

```c
// 計算理論 P-S 時差
if (P_S_time - upd_sec < 0):
    // S 波已到達：使用 P_S_time 秒前的 Pd
    ind = floor(P_S_time)
else:
    // S 波未到：使用最新 Pd
    ind = upd_sec

Pd_used = Pd[ind]
```

此設計確保 Pd 不受 S 波污染，是 EEW 系統中對規模估算的重要保護機制。

### 6.3 儀器類型衰減關係

依 `inst` 欄位選用對應衰減關係：

**Mpd（P 波位移規模）**：
```
HH（超寬頻，inst=2）: M = 5.000 + 1.102×log₁₀(Pd) + 1.737×log₁₀(R)
HL（強震儀，inst=1）: M = 5.067 + 1.281×log₁₀(Pd) + 1.760×log₁₀(R)
HS（短週期，inst=3）: M = 4.811 + 1.089×log₁₀(Pd) + 1.738×log₁₀(R)
```

其中 `R = √(dep² + epi_dis²)`（震源距，km），`Pd` 單位為 cm。

**Mtc（卓越週期規模）**：
```c
Mtc = cal_Tc(Tc)   // Tc 為測量的卓越週期
```

### 6.4 Z-score 異常值濾除

在計算加權平均之前，先對各站規模估計值做標準化濾除：

```c
cal_avg_std_mag(Mpd_array, dd, &avg, &std);

for(i=0; i<dd; i++) {
    z = fabs(Mpd[i] - avg) / std
    if (z < 1.0 && Mpd[i] > 0):
        → 保留（|z| < 1σ）
}

// 加權平均（v2 修正：sum_wei > 0 防除以零）
if(sum_wei > 0.0)
    Mpd_final = Σ (wei[i]/sum_wei) × Mpd[i]
```

Z-score 門檻為 1.0σ（較嚴格），可能在站數少時過度排除。

### 6.5 P 波振幅調整因子（Padj）

理論 vs 觀測 P 波振幅比值，用於評估規模估算的一致性：

```c
Theo_pa_ratio[i] = Pa_observed / pa_theoretical(Mpd, R, inst)
Padj = weighted_mean(Theo_pa_ratio)   // Z-score 過濾後加權平均
```

---

## 7. 報告產出（Report_seq）

### 7.1 報告控制邏輯

```c
rep_fin++              // 報告計數
ccount = rep_fin       // 本事件已產出報告數

if(ccount > Term_num)  // 超過 50 次 → 停止
    return;

// 第一份報告時記錄時間（用於 60s 超時重置）
if(rep_fin == 0) first_rp_time = t_now;
```

### 7.2 proc_time 計算

```c
pro_time = t_now - hyp.time0   // 目前時間 - 發震時刻
if(pro_time > 28800) pro_time -= 28800   // 8小時偏移修正（UTC+8 錯誤處理？）
```

> **注意**：`28800 = 8×3600`（UTC+8 秒數）。此修正意味若 `t_now - hyp.time0 > 8小時`，程式認為時間計算有誤並強制減去時區差。此邏輯存在隱患（見第 11 節）。

### 7.3 .rep 檔案格式

報告檔名：`YYYYMMDDHHMMSS_n<N>.rep`

```
Reporting time   YYYY/MM/DD HH:MM:SS.ss  averr=X.X Q=XX Gap=XXX n=X n_c=X

year  month  day  hour  min  sec     lat       lon      dep   Mall  Mpd   Mtc  proc_time
XXXX  X      XX   XX    XX   XX.XX  XX.XXXX  XXX.XXXX  XX.XX  X.XX  X.XX  X.XX  XX.XX

Sta   C  N  L   lat    lon     pa       pv       pd      tc    Mtc  MPd  Perr  Dis  H_Wei  Parr  [各站詳細資料]
```

---

## 8. 輔助數學函式

### 8.1 delaz() — 球面距離計算

採用橢球體近似公式，非簡單球面三角：

```c
avlat = (elat + slat) / 2.0
a = 1.840708 + avlat×(0.0015269 + avlat×(-0.00034 + avlat×1.02337e-6))
b = 1.843404 + avlat×(-6.93799e-5 + avlat×(8.79993e-6 + avlat×-6.47527e-8))

dx = a × Δlon × 60    (經度差換算 km)
dy = b × Δlat × 60    (緯度差換算 km)
distance = √(dx² + dy²)
```

係數 `a`, `b` 是以台灣緯度範圍（22°–26°N）擬合的多項式，精度優於標準球面公式。

### 8.2 matinv() — 4×4 矩陣求逆

使用 **Gauss-Jordan 消去法**（不帶主元選取）：

```c
for(k=0; k<n; k++) {
    pivot = 1.0 / a[0][0]
    // 消去第 0 行、第 0 列
    for(i=0; i<nm1; i++) {
        yy = -v[i] × pivot
        a[i][n-1] = yy
        for(j=0; j<nm1; j++) a[i][j] = a[i+1][j+1] + v[j]×yy
    }
    a[n-1][n-1] = -pivot
}
// 對稱化
```

> **注意**：不帶主元選取的 Gauss-Jordan 在矩陣接近奇異時數值不穩定。程式透過加 `0.00001` 正則化項（`a[i][i] += 0.00001`）緩解此問題，但不能完全避免。

### 8.3 wt() — 殘差加權函式

```
H = √(dep² + dist²)               ← 震源距
若 H > 100 km：
    xwt × = 500 / (9H - 400)      ← 遠站降權（線性）
xwt × = (1 / (1 + |res|))²        ← 殘差降權（Huber-like）
```

### 8.4 cal_Tc() 與 cal_pgv()

```c
cal_Tc(Tc):   // Tc → Mtc（線性關係，具體係數需確認）
cal_pgv(pd):  // Pd → PGV 換算（經驗公式）
```

---

## 9. 速度模型

### 9.1 P 波梯度速度模型

| 層次 | 深度 | v₀ | vg |
|------|------|----|----|
| 淺層 | z < 40 km | 5.10298 km/s | 0.06659 km/s/km |
| 深層 | z ≥ 40 km | 7.80479 km/s | 0.00457 km/s/km |

40 km 深度的 P 波速度：`v(40) = 5.103 + 0.06659×40 = 7.767 km/s`（與深層初速 7.805 km/s 近似連續）

### 9.2 S 波梯度速度模型

| 層次 | 深度 | v₀ | vg |
|------|------|----|----|
| 淺層 | z < 50 km | 2.9105 km/s | 0.0365 km/s/km |
| 深層 | z ≥ 50 km | 4.5374 km/s | 0.0023 km/s/km |

### 9.3 Vp/Vs 比值

淺層 Vp/Vs（z=20 km）：`(5.103 + 0.06659×20) / (2.9105 + 0.0365×20) = 6.435 / 3.640 ≈ 1.77`

此值接近台灣地殼典型值（1.73–1.80）。

---

## 10. 衰減關係

衰減關係形式：**M = a + b×log₁₀(A) + c×log₁₀(R)**

| 儀器 | 方法 | a | b | c |
|------|------|---|---|---|
| HH（超寬頻） | Mpd | 5.000 | 1.102 | 1.737 |
| HL（強震儀） | Mpd | 5.067 | 1.281 | 1.760 |
| HS（短週期） | Mpd | 4.811 | 1.089 | 1.738 |

其中 `R` 為震源距（km），`A = Pd`（P 波位移，cm）。
衰減係數 `c ≈ 1.74` 介於體波幾何擴散（1.0）與面波（0.5）之間，反映了近場體波加上非彈性衰減。

---

## 11. 未來可改善的演算法方向

### 11.1 定位演算法改進

#### 現況問題
- 10 次迭代為固定次數，無收斂判斷（可能在已收斂後繼續無意義的迭代）
- 深度格搜只有 10 個點（10–100 km，間距 10 km），解析度不足
- 初始震源選取偏向較早到時站，在觸發站分布不均時初始位置可能偏離
- 矩陣求逆使用無主元 Gauss-Jordan，在病態矩陣時數值不穩

#### 改善方向

**a. 自適應迭代收斂判斷**
```c
// 現況：固定 10 次
for(itr=0; itr<10; itr++) { ... }

// 建議：加入收斂條件
double correction_norm;
for(itr=0; itr<MAX_ITER; itr++) {
    ...計算修正量...
    correction_norm = fabs(c[0]) + fabs(c[1]) + fabs(c[2]) + fabs(c[3]);
    if(correction_norm < CONVERGENCE_THRESHOLD) break;
}
```

**b. 深度格搜細化**
```c
// 現況：10–100 km，10 km 間距
// 建議：兩階段格搜
// 第一階：0–200 km，20 km 間距（粗搜）
// 第二階：最佳粗搜深度 ± 10 km，1 km 間距（細搜）
```

**c. 以 SVD 取代 Gauss-Jordan**
```c
// 使用奇異值分解（Singular Value Decomposition）
// 優點：對病態矩陣數值穩定，可識別欠定方程
// 可使用 LAPACK 的 dgelss 或 dgels 函式
```

**d. 加入臺灣 3D 速度模型**
```
現況：2 層梯度模型（各向同性）
建議：引入 3D Vp 模型（例如 Lin et al. 2014 台灣地殼速度模型）
    - 可減少系統性走時誤差
    - 特別改善西部平原（沉積層）與中央山脈交界地區的定位精度
```

---

### 11.2 規模估算改進

#### 現況問題
- Z-score 門檻 1.0σ 較嚴格，站數少時可能排除過多
- Mpd 與 Mtc 簡單算術平均，未考慮各自的不確定性
- 衰減關係為固定係數，不考慮路徑效應（如台灣東西向的速度差異）
- 大規模地震（M > 7）的 Pd 飽和問題未處理

#### 改善方向

**a. 自適應 Z-score 門檻**
```c
// 現況：固定 z < 1.0
// 建議：依站數自適應
double z_threshold;
if(dd <= 4)       z_threshold = 1.5;   // 站數少 → 寬鬆
else if(dd <= 8)  z_threshold = 1.2;
else              z_threshold = 1.0;   // 站數多 → 嚴格
```

**b. 考慮 Pd 時間窗口大小的不確定性加權**
```
越晚到的 Pd（upd_sec 越大）通常越穩定 → 可賦予較高權重
P_S_time 越小的測站（S 波較快到）→ Pd 可用時間越短 → 降低權重
```

**c. 引入 τ_c（卓越週期）飽和修正**
```
大規模地震（M > 7）的 τ_c 會飽和在 ~3 秒（Kanamori 2005）
建議：對 τ_c > 2.5s 的估計值加入飽和修正或切換到其他指標
```

**d. 即時規模更新（Rolling Magnitude）**
```
現況：每次觸發才重新計算
建議：隨著新 upd_sec 到達，只更新有變化的站，避免重複計算
    - 可將規模計算分解為「累積統計量更新」形式
    - 降低計算延遲
```

---

### 11.3 群聚過濾改進

#### 現況問題
- 均值計算包含待過濾站，造成均值偏移（Masking Effect）
- 只做一輪過濾，偏遠雜訊站可能拉偏均值後影響附近正常站的判定
- 時間窗（15s）和距離窗（100km）為固定值，不依地震規模調整

#### 改善方向

**a. 迭代式群聚過濾（DBSCAN-like）**
```
第一輪：使用全部站計算均值，剔除離群站
第二輪：用剩餘站重新計算均值，再次剔除離群站
重複直到收斂（沒有站被剔除）
```

**b. 基於最近觸發時間的動態觸發窗**
```c
// 現況：固定 Trig_tm_win = 15s
// 建議：依預估震央距調整
// 若預估震央距 ~50 km，P 波傳播 50 km 需 ~8s
// 動態窗口 = 預估傳播時間 × 1.5（容錯係數）
```

---

### 11.4 架構層面改進

#### a. 震源參數不確定性估計

```
現況：只輸出單一最佳解（xla0, xlo0, depth0）
建議：輸出信賴區間（可用共變異數矩陣對角線元素估計）
    σ_lat = √(A⁻¹[1][1] × errsum)
    σ_lon = √(A⁻¹[0][0] × errsum)
    σ_dep = √(A⁻¹[2][2] × errsum)
```

#### b. 連續定位（Associator）

```
現況：每次訊息到達都重新做群聚過濾
建議：引入事件關聯器（Associator），維護多個並行事件假說
    - 同時有兩個地震發生時（如主震+餘震）可分別追蹤
    - 使用 REAL（Rapid Earthquake Association and Location）等方法
```

#### c. 機器學習增強

```
定位階段：
  - 以 PhaseNet 或 EQTransformer 替代傳統 P 波拾取，提升 Pd 品質
  - 使用 CNN 直接從波形估算規模（不依賴單一特徵）

群聚過濾：
  - 以 DBSCAN 或 HDBSCAN 替代均值距離過濾，
    對不規則形狀的觸發群（如斜斷層）更魯棒

震源定位：
  - GNN（Graph Neural Network）以測站網絡為圖，
    直接輸出震源座標，端到端訓練
```

#### d. proc_time 計算修正

```c
// 現況（有潛在問題）：
if(pro_time > 28800) pro_time -= 28800;

// 問題：若實際處理時間剛好超過 8 小時（不可能但防禦性）
//       或 hyp.time0 計算錯誤，此修正會產生負的 proc_time

// 建議：直接判斷合理範圍
if(pro_time < 0 || pro_time > 300) {
    logit("e", "tcpd: Abnormal proc_time=%.1f, check hyp.time0\n", pro_time);
    pro_time = -1.0;  // 標示無效
}
```

---

## 總結

| 模組 | 現有技術 | 成熟度 | 改善優先級 |
|------|---------|--------|-----------|
| 速度模型 | 2 層梯度 | 中等 | 高（引入 3D 模型） |
| 迭代定位 | Gauss-Jordan 最小二乘 | 高 | 中（SVD + 自適應收斂） |
| 深度估算 | 10 點格搜 | 低 | 高（兩階段細化） |
| 群聚過濾 | 單次均值距離 | 低 | 高（迭代式 DBSCAN） |
| 規模估算 | Mpd + Mtc 平均 | 中等 | 中（加權不確定性） |
| Z-score 過濾 | 固定 1.0σ | 中等 | 低（自適應門檻） |
| 異常值排除 | 最大殘差逐一移除 | 高 | 低（成熟穩健） |
| 整體架構 | 單事件順序處理 | 中等 | 高（多事件並行） |

`tcpd_new` 作為即時 EEW 系統，在計算速度與演算法精度之間已做了良好的平衡。主要瓶頸在於**速度模型的簡化**（單一 1D 梯度模型）和**深度估算解析度不足**，這兩點對定位精度影響最大，也是最值得優先投入改進的方向。
