# tcpd_new — Earthworm EEW 定位與規模估算模組

**版本：** v2（優化版）
**最後編譯：** 2026-03-23
**容器：** `earthworm_eew` (`cwadayi/earthworm_ubuntu22.04_eew:v1`)
**原始碼路徑：** `/opt/earthworm/EEW_src/tcpd_new/`
**部署路徑：** `/opt/earthworm/earthworm_8.0/bin/tcpd`

---

## 概述

`tcpd_new` 是 Earthworm EEW（地震早期預警）系統的**核心定位與規模估算模組**。從 `PICK_RING` 讀取 P 波到時（TYPE_EEW 格式），進行群聚過濾、震源定位（`locaeq()`）、規模估算（`Magnitude()`），最後輸出 `.rep` 報告並寫入 `EEW_RING`。

本目錄為 v2 改良版，相較原始 `tcpd/` 模組進行了 42 項修正與優化：
- **Bug 修正**：讀取格式錯誤時意外終止程式（`return 0` → `continue`）
- **廢棄 API 替換**：`ftime()` + `struct timeb` → `gettimeofday()` + `struct timeval`（微秒精度）
- **全面 logit 遷移**：24 處 `printf` → `logit()`，訊息正確寫入 EW 日誌
- **未使用變數清除**：`Magnitude()` 中 `j, l, vv, s_num`；`Report_seq()` 中 `p, cc, tmp, tp`

---

## 目錄結構

```
tcpd_new/
├── tcpd_new.c                  # 主程式原始碼（2986 行，優化版）
├── tcpd_new.d                  # 模組設定檔
├── makefile.unix               # Linux 編譯腳本
├── locate.h                    # 資料結構定義：PEEW, HYP, MAG 及函式原型
├── dayi_time.h                 # 時間工具函式原型（epoch ↔ 日曆時間）
├── time_ew.h                   # Earthworm 時間函式介面（官方標頭檔）
├── num_eew_status              # EEW 事件流水號持久化檔（1–10000 循環）
├── tcpd_new                    # 已編譯 Linux 二進位（236 KB，零警告）
├── README.md                   # 本文件
├── tcpd_new_summary.md         # 第一版（v1）改進摘要
├── tcpd_optimization_report.md # v2 優化報告（42 項修改完整清單）
├── tcpd_flow_analysis.md       # 程式流程中文詳細解析（8 章）
└── optimize_tcpd.py            # 自動化優化腳本（Python，逐行匹配法）
```

---

## 各檔案詳細說明

### `tcpd_new.c` — 主程式（2986 行）

#### 全域資料結構

| 變數 | 類型 | 說明 |
|------|------|------|
| `ptr[1000]` | `PEEW[]` | 全域 P 波緩衝區，最多容納 1000 筆測站資料 |
| `vsn_ntri[1000]` | `PEEW[]` | 群聚過濾後的有效觸發站暫存 |
| `vsn_trigger` | `int` | 當前有效觸發站數 |
| `max_sta` | `int` | 上次報告時的觸發站數（比較用） |
| `count_max` | `int` | 連續無新觸發計數器（達 `Report_Limit` 後重置） |
| `num_eew` | `int` | 當前事件流水號（1–10000） |
| `rep_fin` | `int` | 已產出報告旗標（>0 代表有報告產出，60 秒後重置） |
| `first_rp_time` | `double` | 第一個 P 波到時（epoch s，計算 proc_time 基準） |

#### 函式清單

| 函式 | 行號 | 說明 |
|------|------|------|
| `main()` | 151 | 主循環：讀取訊息 → 過濾 → 群聚 → 觸發判斷 |
| `tcpd_config()` | 747 | 解析 `.d` 設定檔（使用 kom.c） |
| `tcpd_lookup()` | 1032 | 查詢 Earthworm 訊息類型（TYPE_EEW 等） |
| `tcpd_status()` | 1090 | 發送心跳 / 錯誤訊息至 transport ring |
| `locaeq()` | 1136 | 迭代最小二乘震源定位（10 次 + 深度格搜） |
| `hyp_cmp()` | 1590 | qsort 比較函式（震源解排序） |
| `wt()` | 1598 | 計算走時殘差加權值 |
| `matinv()` | 1618 | 4×4 矩陣求逆（定位線性化用） |
| `delaz()` | 1647 | 計算兩點球面距離（km） |
| `processTrigger()` | 1666 | 呼叫 locaeq() + Magnitude() + Report_seq() |
| `Magnitude()` | 1848 | 估算 Mpd、Mtc、Mall，Z-score 濾除異常值 |
| `Report_seq()` | 2160 | 寫入 `.rep` 報告檔並發送 TYPE_EEW 至 EEW_RING |
| `cal_pgv()` | 2428 | Pd → PGV 換算 |
| `cal_Tc()` | 2440 | Tc 值換算（主要週期 → 規模指標） |
| `time_transe()` | 2451 | 日曆時間 → epoch 秒（台灣時區 UTC+8） |
| `disp_time()` / `disp_time1()` | 2519/2476 | 格式化時間字串輸出 |
| `split_c()` | 2576 | 以空白分割字串（解析 TYPE_EEW 訊息） |
| `Mpd1_HH/HL/HS()` | 2645+ | Pd 法規模衰減關係（超寬頻/強震/短週期） |
| `Mpv1_HH/HL/HS()` | 2655+ | Pv 法規模衰減關係 |
| `cal_avg_std()` | 2711 | 計算平均值與標準差（含最大最小值） |
| `cal_avg_std_mag()` | 2748 | MAG_DATA 陣列平均/標準差 |
| `cal_z()` | 2828 | 計算 Z-score |
| `sort_array()` | 2858 | double 陣列排序 |
| `sort_array_mag()` | 2876 | MAG_DATA 陣列排序 |
| `sort_array_P_S()` | 2894 | PEEW 陣列依 P 到時排序 |
| `ReportEEW()` | 2912 | 將 EEW 訊息寫入 EEW_RING |
| `ReportEEW_record()` | 2932 | 將完整紀錄訊息寫入 EEW_RING |
| `pa_HS/HL/HH()` | 2954+ | 理論 P 波振幅衰減模型（各儀器類型） |

#### main() 主循環邏輯

```
啟動
├── 讀取 num_eew_status（事件流水號）
├── 初始化 ptr[1000]（flag=0），設定預設參數
├── tcpd_config() 讀取設定檔
├── tcpd_lookup() 查詢訊息類型
├── 連接 PICK_RING（輸入）與 EEW_RING（輸出）
└── while(1) 主循環
    ├── 每 15s 發送心跳
    ├── tport_getmsg() 從 PICK_RING 取訊息
    │
    │  ── 訊息處理（TYPE_EEW）──
    ├── split_c() 分割 14 欄位
    │   └── 格式錯誤 → logit("e",...) + continue（Bug 修正點）
    ├── 品質過濾：
    │   ├── 僅保留垂直向（out_ss[1][2]=='Z'）
    │   ├── 一般站：weight < Ignore_weight_P (2)
    │   └── 特殊站（YOJ/EOS2-4）：weight < 3
    ├── 存入 ptr[] 緩衝區（空位插入，upd_sec 更新 Pd）
    │
    ├── gettimeofday() 取得目前時間
    ├── 清除過期資料（age > Active_parr_win = 80s）
    ├── 去重（同 SCNL 保留最新 upd_sec，合併 Pd 歷史）
    │
    ├── 群聚過濾：
    │   ├── 計算所有有效站的平均位置與平均 P 到時
    │   ├── 剔除：dis > Trig_dis_win(100km) 或 dt > Trig_tm_win(15s)
    │   └── 13 個離島/偏遠站免篩選
    │
    ├── vsn_trigger > 4 時觸發：
    │   ├── count_max > Report_Limit → 重置（新事件）
    │   ├── vsn_trigger == max_sta → 無新觸發，跳過
    │   ├── vsn_trigger > max_sta → 新觸發 → processTrigger()
    │   └── vsn_trigger < max_sta → count_max++，跳過
    │
    └── 60 秒後若 rep_fin>0 → 重置狀態等待下次地震
```

#### TYPE_EEW 訊息格式（14 欄位）

```
YULB HHZ TW 01 121.2971 23.3924 0.006514 0.000187 0.000074 4.297026 1328085502.958 0 2 1
[0]  [1] [2][3]  [4]       [5]     [6]      [7]      [8]       [9]         [10]     [11][12][13]
 sta  chn net loc  lon       lat     Pa       Pv       Pd        Tc        P_epoch  wt  inst upd_sec
```

#### 關鍵 Bug 修正（行 ~419）

```c
/* 原始版本 — BUG：從 while(1) 主循環意外退出 */
printf("Read Buffer Error !");
return 0;

/* 修正版本 — 正確行為：記錄錯誤並繼續循環 */
logit("e", "tcpd: Read Buffer Error - expected 14 fields, got %d\n", num_split);
continue;
```

#### 時間 API 替換（3 處）

```c
/* 原始：毫秒精度，POSIX 廢棄 API */
struct timeb tp;
ftime(&tp);
now_time = (double)tp.time + (double)tp.millitm/1000;

/* 修正：微秒精度，POSIX 現行 API */
struct timeval tv;
gettimeofday(&tv, NULL);
now_time = (double)tv.tv_sec + (double)tv.tv_usec/1000000.0;
```

---

### `locate.h` — 資料結構定義

#### `PEEW` — 單站 P 波資料（完整欄位）

```c
typedef struct {
    int    flag;          // 0=空, 1=有效, 2=有效但不用於規模估算
    int    serial;        // 序號
    char   stn_name[8];   // 測站代碼
    char   stn_Loc[8];    // 位置碼（LC）
    char   stn_Comp[8];   // 分量（HHZ/HLZ/HSZ）
    char   stn_Net[8];    // 網絡碼
    double latitude;      // 測站緯度（WGS84）
    double longitude;     // 測站經度（WGS84）
    double altitude;      // 測站高度（m）
    double P;             // P 波到時（epoch 秒）
    double Pa;            // P 波加速度振幅（gal）
    double Pv;            // P 波速度振幅（kine）
    double Pd[15];        // P 波位移振幅序列（各秒更新，upd_sec 索引）
    double Tc;            // 主要週期（s）
    double dura;          // 持時
    double report_time;   // 報告時間
    int    npoints;       // 波形點數
    double perr;          // P 到時殘差（s，locaeq 計算）
    double wei;           // 定位加權值（locaeq 計算）
    int    weight;        // P 波品質碼（0=最佳，數字越大越差）
    int    inst;          // 儀器類型（1=HL強震, 2=HH超寬頻, 3=HS短週期）
    int    upd_sec;       // 已接收的 Pd 更新秒數
    int    usd_sec;       // 實際用於規模估算的 Pd 秒數（避免 S 波污染）
    double P_S_time;      // 理論 P-S 波時差（s，由震央距推算）
    int    pin;           // SCNL 唯一識別碼
} PEEW;
```

> **Pd 選用邏輯**：若 `P_S_time - upd_sec < 0`（S 波已到），使用 `Pd[floor(P_S_time)]`；否則使用 `Pd[upd_sec]`，確保不受 S 波污染。

#### `HYP` — 震源解

```c
typedef struct {
    double xla0;    // 震央緯度
    double xlo0;    // 震央經度
    double depth0;  // 震源深度（km）
    double time0;   // 發震時刻（epoch s）
    int    Q;       // 定位品質（-6 = 最高品質）
    double averr;   // 平均走時殘差（s）
    double gap;     // 方位角缺口（度，0–360）
    double avwei;   // 加權平均值
} HYP;
```

#### `MAG` — 規模估算

```c
typedef struct {
    double xMpd;    // Pd 法規模（P 波位移，依儀器類型選用衰減關係）
    double mtc;     // Tc 法規模（主要週期法）
    double ALL_Mag; // 整合規模（Mpd 與 Mtc 的加權平均）
    double Padj;    // P 波振幅調整因子（理論振幅比）
} MAG;
```

#### `MAG_DATA` — 規模資料點（Z-score 濾除用）

```c
typedef struct {
    double mag;  // 個別測站規模估計值
    double wei;  // 對應加權值
} MAG_DATA;
```

#### `stat` — 統計累積

```c
typedef struct {
    double avg;      // 平均值
    double std;      // 標準差
    double new_sum;  // 加權總和
    int    new_num;  // 有效點數
} stat;
```

#### 函式原型（locate.h 宣告）

```c
void   locaeq(PEEW *ptr, int nsta, HYP *hyp);          // 震源定位
void   locaeq_grid(PEEW *ptr, int nsta, HYP *hyp);     // 深度格搜
int    hyp_cmp(const void *x1, const void *x2);        // qsort 比較
void   matinv(double a[4][4], int n);                  // 矩陣求逆
double wt(double dist, double res, int iter, double depth, double weight); // 加權
double delaz(double elat, double elon, double slat, double slon); // 球面距離

void   processTrigger(int vsn_trigger, PEEW *vsn_ntri);
void   Magnitude(PEEW *vsn_ntri, int qtrigger, HYP hyp, MAG *mag, double *pPgv, int *ntc);
void   Report_seq(PEEW *vsn_ntri, int qtrigger, HYP hyp, MAG mag, int ntc);

double Mpd1_HH(double pd, double dis);  // 超寬頻 Pd 規模
double Mpd1_HL(double pd, double dis);  // 強震儀 Pd 規模
double Mpd1_HS(double pd, double dis);  // 短週期 Pd 規模
```

---

### `tcpd_new.d` — 設定檔

#### 基本參數

| 參數 | 值 | 說明 |
|------|----|------|
| `MyModuleId` | `MOD_TCPD` | Earthworm 模組 ID |
| `RingName` | `PICK_RING` | 輸入 ring（P 波到時訊息） |
| `RingName_out` | `EEW_RING` | 輸出 ring（EEW 報告訊息） |
| `HeartBeatInterval` | `15` s | 心跳間隔（超時則被 startstop 重啟） |
| `LogFile` | `0` | 日誌開關（0=關閉，1=開啟） |
| `MagMin` | `0.5` | 最小有效規模 |
| `MagMax` | `10.0` | 最大有效規模 |
| `Ignore_weight_P` | `2` | P 波品質碼門檻（≥ 此值丟棄；原版為 3） |
| `Ignore_weight_S` | `2` | S 波品質碼門檻 |
| `Term_num` | `50` | 每次地震最大更新報數（原版為 15） |
| `Show_Report` | `1` | 報告輸出（0=關閉，1=開啟） |

#### 觸發窗參數

| 參數 | 值 | 說明 | 原版對比 |
|------|----|------|---------|
| `Trig_tm_win` | `15.0 s` | P 到時差觸發窗 | 原版 40.0 s（更嚴格） |
| `Trig_dis_win` | `100.0 km` | 距離觸發窗 | 原版 220.0 km（更嚴格） |
| `Active_parr_win` | `80.0 s` | P 到時有效期 | 原版 45.0 s（更長） |

#### P 波速度模型（兩層梯度）

| 層次 | 參數 | 值 |
|------|------|----|
| 淺層（z < 40 km） | `Boundary_P` | 40.0 km |
| | `SwP_V` | 5.10298 km/s（層頂速度） |
| | `SwP_VG` | 0.06659 km/s/km（速度梯度） |
| 深層（z ≥ 40 km） | `DpP_V` | 7.80479 km/s |
| | `DpP_VG` | 0.00457 km/s/km |

#### S 波速度模型

| 層次 | 參數 | 值 |
|------|------|----|
| 淺層（z < 50 km） | `Boundary_S` | 50.0 km |
| | `SwS_V` | 2.9105 km/s |
| | `SwS_VG` | 0.0365 km/s/km |
| 深層（z ≥ 50 km） | `DpS_V` | 4.5374 km/s |
| | `DpS_VG` | 0.0023 km/s/km |

---

### `makefile.unix` — 編譯腳本

**編譯器：** `$(CC)`（GCC），由 Earthworm 環境變數 `$PLATFORM` 決定平台旗標

**關鍵變數：**

| 變數 | 值/說明 |
|------|---------|
| `APP` | `tcpd_new`（輸出二進位名稱） |
| `CFLAGS` | `-D_REENTRANT $(GLOBALFLAGS) -Wall -g` |
| `EW_LIBS` | `swap.o trheadconv.o -lew_mt`（Earthworm 函式庫） |
| `BINARIES` | `tcpd_new.o kom.o getutil.o sleep_ew.o time_ew.o transport.o lockfile_ew.o lockfile.o` |

**編譯指令：**
```bash
source /opt/earthworm/earthworm_8.0/environment/ew_linux.bash
cd /opt/earthworm/EEW_src/tcpd_new
make -f makefile.unix        # 編譯 → tcpd_new（236 KB，零警告）
make -f makefile.unix clean  # 清除 .o 中間檔
```

---

### `dayi_time.h` — 時間工具函式

此標頭檔宣告四個時間轉換函式（實作位於 `tcpd_new.c` 末段，2451–2574 行）：

```c
// 日曆時間 → epoch 秒（基準：1970/01/01 08:00:00 UTC+8）
void time_transe(int yr, int mo, int dy, int hr, int mn, int se, char *s, int *tt);

// 格式化時間顯示（詳細格式）
void disp_time(int yr, int mo, int dy, int hr, int mn, int se, char *ss, int type);

// 格式化時間顯示（簡潔格式）
void disp_time1(int yr, int mo, int dy, int hr, int mn, int se, char *ss, int type);

// 取得目前時間字串（可加 plus_min 分鐘偏移）
void time_now(char *ss, int plus_min, int type);
```

---

### `time_ew.h` — Earthworm 時間介面

官方 Earthworm 標頭檔（RCS 版本管理，v1.4 2009/06/15），宣告執行緒安全的時間函式：

```c
time_t     timegm_ew(struct tm *);               // UTC 轉 time_t（跨平台）
struct tm *localtime_ew(const time_t *, struct tm *);  // 本地時間（執行緒安全）
struct tm *gmtime_ew(const time_t *, struct tm *);     // UTC 時間（執行緒安全）
double     hrtime_ew(double *);                  // 高精度時間
char      *datestr23(double, char *, int);       // epoch → 23 字元時間字串
char      *datestr23_local(double, char *, int); // epoch → 本地時間字串
char      *utc_ctime_ew(time_t *);               // UTC 格式時間字串
```

> **注意：** 此標頭檔對 Windows（`_WINNT`）定義了 `gettimeofday()` 包裝函式及 `struct timezone`，在 Linux 環境下這些定義不啟用。`tcpd_new.c` 使用 `#include <sys/time.h>` 原生提供的 `gettimeofday()`。

---

### `num_eew_status` — 事件流水號

純文字檔，記錄當前 EEW 事件流水號（整數，範圍 1–10000）。

- 程式啟動時讀取；若檔案不存在或損壞，重置為 1
- 每次地震處理完成後自動遞增並寫回
- 超過 10000 時循環回 1
- 目前值：`1`（上次啟動後重置或首次使用）

---

### `tcpd_new`（二進位）

- **格式：** ELF 64-bit LSB executable（Linux x86_64）
- **大小：** 236,256 bytes（236 KB）
- **編譯時間：** 2026-03-23 07:23（UTC）
- **編譯旗標：** `-Wall -g -D_REENTRANT`，零警告零錯誤
- **執行方式：** `tcpd_new tcpd_new.d`（需先設定 Earthworm 環境變數）

---

## 完整處理流程

```
PICK_RING (TYPE_EEW，14 欄位)
        │
        ▼
1. 訊息格式驗證
   ├── split_c() 分割 → 必須恰好 14 欄位
   └── 格式錯誤 → logit("e",...) + continue（Bug 修正：原為 return 0）

        │
        ▼
2. P 波品質過濾
   ├── 僅處理垂直向（stn_Comp[2]=='Z'）
   ├── 一般站：丟棄 weight ≥ 2
   └── 特殊站（YOJ, EOS2/3/4）：丟棄 weight ≥ 3

        │
        ▼
3. 緩衝區管理（ptr[1000]）
   ├── 空位插入新測站（flag=0 → flag=1）
   ├── gettimeofday() 取得目前時間
   ├── 清除過期（|now - P| > 80s，flag → 0）
   └── 去重（同 SCNL 且 P 相同：保留較新 upd_sec，合併 Pd 歷史）

        │
        ▼
4. 群聚過濾（建立 vsn_ntri[]）
   ├── 計算所有有效站的平均經緯度、平均 P 到時
   ├── 剔除條件（一般站）：
   │   ├── delaz(站, 均值) > 100 km
   │   └── |P - 平均P| > 15 s
   └── 免篩選站：EOS*/YOJ/JMJ/PCY/LAY/TWH/PNG/PHU/PTM/PTT/VCH/VWU/WLC

        │
        ▼
5. 觸發判斷
   ├── vsn_trigger ≤ 4 → 等待更多測站
   ├── count_max > Report_Limit → 重置（本次地震結束）
   ├── vsn_trigger == max_sta → 無新觸發，跳過
   ├── vsn_trigger < max_sta → count_max++，跳過
   └── vsn_trigger > max_sta → 新觸發！→ processTrigger()

        │
        ▼
6. processTrigger() — 定位與規模
   ├── locaeq()：迭代最小二乘定位
   │   ├── 初始震源：以最早 P 到時測站為中心
   │   ├── 10 次迭代：速度模型 + 矩陣求逆 + 殘差加權
   │   ├── locaeq_grid()：深度格搜（0/10/20/30/50/70/100/150 km）
   │   └── 輸出 HYP（lat, lon, depth, time0, averr, gap, Q）
   └── Magnitude()：規模估算
       ├── 排除：BH 網、非垂直向、Tc<0、P_S_time<2s、upd_sec<2s
       ├── 選用 Pd 秒數：min(upd_sec, floor(P_S_time))
       ├── Mpd：依儀器類型（HH/HL/HS）套用衰減關係
       ├── Mtc：cal_Tc(Tc) 換算
       ├── Z-score 濾除異常值（|z| > 1.5）
       └── 輸出 MAG（xMpd, mtc, ALL_Mag）

        │
        ▼
7. Report_seq() — 輸出報告
   ├── 寫入 .rep 檔（$EW_PARAMS/YYYYMMDDHHMMSS_n<N>.rep）
   └── ReportEEW() → tport_putmsg() → EEW_RING（TYPE_EEW）
```

---

## EEW 報告格式

**報告檔名：** `YYYYMMDDHHMMSS_n<N>.rep`

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
| `Mall` | 整合規模（Mpd 與 Mtc 加權平均） |
| `Mpd` | Pd 法規模（P 波位移） |
| `Mtc` | Tc 法規模（主要週期） |
| `averr` | 平均走時殘差（s） |
| `Gap` | 方位角缺口（度，越小代表方位覆蓋越好） |
| `n` / `n_c` | 定位用站數 / 非同址站數 |
| `proc_time` | 第一個 P 到時至產出報告的秒數 |
| `Q` | 品質指標（-6 = 最高） |

---

## 與原始 `tcpd/` 模組比較

| 項目 | 原始 `tcpd/` | `tcpd_new/`（本版） |
|------|-------------|---------------------|
| `Trig_tm_win` | 40.0 s | **15.0 s**（更嚴格，減少誤報） |
| `Trig_dis_win` | 220.0 km | **100.0 km**（更嚴格，避免遠站誤觸發） |
| `Active_parr_win` | 45.0 s | **80.0 s**（更長，允許更多測站加入） |
| `Ignore_weight_P` | 3 | **2**（更嚴格，排除品質較差拾取） |
| `Term_num` | 15 | **50**（更多更新報，追蹤慢演化事件） |
| 讀取格式錯誤處理 | `return 0`（程式意外終止） | **`continue`**（跳過並繼續） |
| 時間精度 | 毫秒（`ftime`，已廢棄） | **微秒**（`gettimeofday`） |
| 編譯警告 | 4 個 | **0 個** |
| `printf` 呼叫 | ~30 個 | **0 個**（全改為 `logit`） |
| 未使用變數 | `j, l, vv, s_num, p, cc, tmp, tp` | **已全部清除** |
| `Intensity_thr` 警告 | 存在 | **接受但忽略**（向後相容） |
| `Pv_thr` 警告 | 存在 | **接受但忽略** |
| `Max_Epi_dis` 警告 | 存在 | **接受但忽略** |

---

## 編譯與部署

### 環境準備

```bash
# 進入容器
sudo docker exec -it earthworm_eew bash

# 設定 Earthworm 環境變數
source /opt/earthworm/earthworm_8.0/environment/ew_linux.bash
# 設定變數：$EW_HOME, $EW_VERSION, $PLATFORM, $GLOBALFLAGS
```

### 編譯

```bash
cd /opt/earthworm/EEW_src/tcpd_new
make -f makefile.unix        # 輸出：tcpd_new（236 KB，零警告）
make -f makefile.unix clean  # 清除 *.o 中間檔
```

### 部署

```bash
# 1. 備份現有二進位
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

## 注意事項

- `num_eew_status` 儲存事件流水號（1–10000 循環），跨次啟動保存，修改需謹慎
- 模組每 15 秒發送心跳；若超時未收到，`startstop` 會嘗試自動重啟
- 13 個離島/偏遠站（EOS\*, YOJ, JMJ, PCY, LAY, TWH, PNG, PHU, PTM, PTT, VCH, VWU, WLC）免距離/時間群聚過濾，但仍參與定位與規模估算
- 規模衰減關係依儀器類型選用：HH（超寬頻）= inst 2、HL（強震儀）= inst 1、HS（短週期）= inst 3
- BH 網測站自動標記為 flag=2（不用於規模估算）
- 所有時間戳記為 UTC（台灣標準時間 = UTC+8）
- `Report_Limit = 1`（程式碼硬編碼）：每次事件觸發後，若觸發站數不再增加，連續 1 次後重置
