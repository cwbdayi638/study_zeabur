# tcpd_new — Earthquake Early Warning Location & Magnitude Module

**版本：** v2（優化版）
**最後編譯：** 2026-03-23
**容器：** `earthworm_eew` (`cwadayi/earthworm_ubuntu22.04_eew:v1`)
**部署路徑：** `/opt/earthworm/earthworm_8.0/bin/tcpd`

---

## 概述

`tcpd_new` 是 Earthworm EEW（地震早期預警）系統的**核心定位與規模估算模組**。
從 `PICK_RING` 讀取 P 波到時，進行震源定位與規模估算，輸出 EEW `.rep` 報告並寫入 `EEW_RING`。

本目錄為改良版，相較原始 `tcpd/` 模組進行了以下修正與優化：
- **Bug 修正**：讀取格式錯誤時意外終止程式的錯誤（`return 0` → `continue`）
- **廢棄 API 替換**：`ftime()` → `gettimeofday()`（微秒精度）
- **全面 logit 遷移**：24 處 `printf` → `logit()`，訊息正確寫入 EW 日誌
- **未使用變數清除**：Magnitude() 與 Report_seq() 中 11 個未使用變數

詳細說明見 [`tcpd_optimization_report.md`](tcpd_optimization_report.md)。

---

## 目錄結構

```
tcpd_new/
├── tcpd_new.c                  # 主程式原始碼（3008 行，優化版）
├── tcpd_new.d                  # 設定檔
├── makefile.unix               # Linux 編譯腳本
├── locate.h                    # 資料結構：PEEW, HYP, MAG
├── dayi_time.h                 # 時間轉換工具（epoch ↔ 日曆）
├── time_ew.h                   # Earthworm 時間函式介面
├── num_eew_status              # EEW 事件編號持久化計數器（1–10000）
├── tcpd_new                    # 已編譯 Linux 二進位（231 KB）
├── README.md                   # 本文件
├── tcpd_new_summary.md         # 第一版改進摘要（v1）
├── tcpd_optimization_report.md # 優化報告（v2，42 項修改）
├── tcpd_flow_analysis.md       # 程式流程中文詳細解析
└── optimize_tcpd.py            # 自動化優化腳本
```

---

## 處理流程

```
PICK_RING (TYPE_EEW)
        │
        ▼
1. P 波品質過濾
   ├── 僅處理垂直向（??Z）
   ├── 丟棄 pick weight ≥ Ignore_weight_P (2) 的資料
   └── 特殊站（YOJ, EOS2-4）門檻寬鬆（weight < 3）
        │
        ▼
2. 管理 P 波緩衝區（最多 1000 筆）
   ├── 清除過期到時（age > Active_parr_win = 80s）
   └── 去重（同 SCNL：保留最新 upd_sec，合併 Pd 歷史）
        │
        ▼
3. 群聚過濾（篩出有效觸發站）
   ├── 計算所有有效站的平均位置與平均 P 到時
   ├── 剔除：距平均位置 > Trig_dis_win (100 km)
   │      或   P 到時差   > Trig_tm_win (15 s)
   └── 離島/偏遠站（EOS, YOJ, JMJ, PCY, LAY...）免篩選
        │
        ▼
4. 觸發條件：有效觸發站數 > 4 且比上次報告增加
        │
        ▼
5. processTrigger() — 定位 + 規模
   ├── locaeq()：迭代最小二乘定位（10 次迭代 + 深度格搜）
   └── Magnitude()：Mpd + Mtc → Mall
        │
        ▼
6. Report_seq() — 輸出報告
   ├── 寫入 .rep 檔案（$EW_PARAMS/）
   └── 發送 TYPE_EEW 訊息 → EEW_RING
```

---

## 設定檔（`tcpd_new.d`）

| 參數 | 值 | 說明 |
|------|----|------|
| `MyModuleId` | `MOD_TCPD` | Earthworm 模組 ID |
| `RingName` | `PICK_RING` | 輸入 ring（P 波到時） |
| `RingName_out` | `EEW_RING` | 輸出 ring（EEW 訊息） |
| `HeartBeatInterval` | 15 s | 心跳間隔 |
| `LogFile` | 0 | 日誌（0=關閉，1=開啟） |
| `MagMin` / `MagMax` | 0.5 / 10.0 | 規模接受範圍 |
| `Ignore_weight_P` | **2** | P 波品質碼門檻（≥ 此值丟棄） |
| `Ignore_weight_S` | **2** | S 波品質碼門檻 |
| `Trig_tm_win` | **15.0 s** | P 到時差觸發窗 |
| `Trig_dis_win` | **100.0 km** | 距離觸發窗 |
| `Active_parr_win` | **80.0 s** | P 到時有效期 |
| `Term_num` | **50** | 每次地震最大更新報數 |
| `Show_Report` | 1 | 報告輸出（0=關閉，1=開啟） |

### P 波速度模型

| 層次 | 參數 | 值 |
|------|------|----|
| 淺層（z < 40 km） | `Boundary_P` | 40.0 km |
| | `SwP_V` | 5.10298 km/s |
| | `SwP_VG` | 0.06659 km/s/km |
| 深層（z ≥ 40 km） | `DpP_V` | 7.80479 km/s |
| | `DpP_VG` | 0.00457 km/s/km |

### S 波速度模型

| 層次 | 參數 | 值 |
|------|------|----|
| 淺層（z < 50 km） | `Boundary_S` | 50.0 km |
| | `SwS_V` | 2.9105 km/s |
| | `SwS_VG` | 0.0365 km/s/km |
| 深層（z ≥ 50 km） | `DpS_V` | 4.5374 km/s |
| | `DpS_VG` | 0.0023 km/s/km |

---

## 與原始 `tcpd/` 模組比較

| 項目 | 原始 `tcpd/`（生產環境） | `tcpd_new/`（本版） |
|------|--------------------------|---------------------|
| `Trig_tm_win` | 40.0 s | **15.0 s**（更緊） |
| `Trig_dis_win` | 220.0 km | **100.0 km**（更緊） |
| `Active_parr_win` | 45.0 s | **80.0 s**（更長） |
| `Ignore_weight_P` | 3 | **2**（更嚴格） |
| `Term_num` | 15 | **50**（更多更新報） |
| `Intensity_thr` | 存在（啟動時產生警告） | **接受但忽略**（向後相容） |
| `Pv_thr` | 存在（啟動時產生警告） | **接受但忽略** |
| `Max_Epi_dis` | 存在（啟動時產生警告） | **接受但忽略** |
| `return 0` bug | ✗ 存在 | **✅ 已修正**（改為 `continue`） |
| 時間精度 | 毫秒（`ftime`） | **微秒**（`gettimeofday`） |
| 編譯警告 | 4 個 | **0 個** |
| `printf` 呼叫 | ~30 個 | **0 個**（全改為 `logit`） |

---

## 編譯

### 環境需求

```bash
# 在容器內執行
sudo docker exec -it earthworm_eew bash
source /opt/earthworm/earthworm_8.0/environment/ew_linux.bash
```

需設定：`$EW_HOME`、`$EW_VERSION`、`$PLATFORM`、`$GLOBALFLAGS`

### 編譯指令

```bash
cd /opt/earthworm/EEW_src/tcpd_new
make -f makefile.unix        # 編譯
make -f makefile.unix clean  # 清除中間檔
```

**編譯結果：** `tcpd_new`（231 KB，零警告零錯誤）

---

## 部署

```bash
# 1. 備份舊版二進位
cp /opt/earthworm/earthworm_8.0/bin/tcpd \
   /opt/earthworm/earthworm_8.0/bin/tcpd.bak_$(date +%Y%m%d)

# 2. 部署新版二進位
cp /opt/earthworm/EEW_src/tcpd_new/tcpd_new \
   /opt/earthworm/earthworm_8.0/bin/tcpd

# 3. （可選）部署新版設定檔
cp /opt/earthworm/EEW_src/tcpd_new/tcpd_new.d \
   /opt/earthworm/run_working/params/tcpd.d

# 4. 重啟 Earthworm
source /opt/earthworm/earthworm_8.0/environment/ew_linux.bash
cd $EW_PARAMS
stopmodule MOD_STARTSTOP
sleep 5
startstop &
```

---

## EEW 報告格式

報告檔名：`YYYYMMDDHHMMSS_n<N>.rep`

```
Reporting time   2026/03/23 02:38:47.XX  averr=0.5 Q=-6 Gap=271 n=8 n_c=6

year  month  day  hour min  sec      lat       lon       dep   Mall  Mpd   Mtc  proc_time
2026  3      23   02   38   35.XX  24.1856  120.9307   30.00  6.19  6.19  0.00   12.34

Sta     C  N  L    lat       lon      pa       pv       pd       tc    Mtc  MPd  Perr  Dis  H_Wei  Parr  ...
```

| 欄位 | 說明 |
|------|------|
| `lat` / `lon` | 震央（WGS84，十進位度） |
| `dep` | 震源深度（km） |
| `Mall` | 整合規模（Mpd 與 Mtc 平均） |
| `Mpd` | Pd 法規模（P 波位移） |
| `Mtc` | Tc 法規模（主要週期） |
| `averr` | 平均走時殘差（s） |
| `Gap` | 方位角缺口（度） |
| `n` / `n_c` | 使用站數 / 非同址站數 |
| `proc_time` | 第一個 P 波到時至報告產出的秒數 |

---

## 主要資料結構（`locate.h`）

### `PEEW` — 單站 P 波資料
```c
char   stn_name[8], stn_Comp[8], stn_Net[8], stn_Loc[8];
double latitude, longitude;
double P;           // P 波到時（epoch s）
double Pa;          // P 波加速度振幅
double Pv;          // P 波速度振幅
double Pd[20];      // P 波位移振幅（各秒更新值）
double Tc;          // 主要週期（s）
double P_S_time;    // 理論 P-S 波時差（s）
int    weight;      // P 波品質碼
int    upd_sec;     // 已累積 Pd 更新秒數
int    flag;        // 0=空, 1=有效, 2=不用於規模估算
```

### `HYP` — 震源解
```c
double xla0, xlo0;  // 震央緯經度
double depth0;       // 震源深度（km）
double time0;        // 發震時刻（epoch s）
double averr;        // 平均走時殘差（s）
double avwei;        // 平均加權值
double gap;          // 方位角缺口（度）
int    Q;            // 定位品質指標
```

### `MAG` — 規模估算
```c
double xMpd;     // Pd 法規模
double mtc;      // Tc 法規模
double ALL_Mag;  // 整合規模
double Padj;     // P 波振幅調整因子
```

---

## 注意事項

- `num_eew_status` 儲存事件流水號（1–10000 循環），跨次啟動保存。
- 模組每 15 秒發送心跳給 `startstop`；若逾時未收到，`startstop` 會嘗試重啟模組。
- 離島/偏遠站（EOS, YOJ, JMJ, PCY, LAY, TWH, PNG, PHU, PTM, PTT, VCH, VWU, WLC）免距離/時間群聚過濾，但仍參與定位與規模估算。
- 規模衰減關係依儀器類型選用：HH（超寬頻）、HL（強震儀）、HS（短週期）。
- 所有時間為 UTC（台灣時間 = UTC+8）。
