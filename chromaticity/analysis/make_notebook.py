# make_notebook.py - generates analysis.ipynb (then: jupyter nbconvert --execute --inplace analysis.ipynb)
import nbformat as nbf, pathlib
nb = nbf.v4.new_notebook(); c = []
md = lambda s: c.append(nbf.v4.new_markdown_cell(s)); code = lambda s: c.append(nbf.v4.new_code_cell(s))

md(r"""# Хроматичность и Twiss-функции magnetic_2_md: OptiM, BELL, COSY, MAD-X
Решётка `C:\GIT_REPOS\SCT_2\lattice\magnetic_2_md.opt`, дейтрон m = 1875.61294257 МэВ, W = 270 МэВ,
γ0 = 1.143952941, β0 = 0.485634729, Bρ = 3.47568504 Т·м.

**Общее определение:** ξ = dQ/dδ, δ = dp/p, натуральная (секступоли = 0), первый порядок; eta = ΔT/T/δ = alpha_c − 1/γ².
Пересчёты: COSY δ_K = ΔK/K = (1+1/γ)δ ⇒ ξ = (1+1/γ)·dQ/dδ_K = 1.874162·raw; MAD-X TWISS: PT = β0·δ ⇒ ξ = β0·DQ = 0.485635·DQ.

Структура: (1) каждый запуск отдельно — команда, сырой вывод, что берём, пересчёт; (2) сводная таблица;
(3) Twiss-функции всех программ; (4) почему ξy различается.""")

code(r'''import re, json, pathlib, subprocess, sys
import numpy as np, pandas as pd, matplotlib.pyplot as plt
R = pathlib.Path(r"C:\GIT_REPOS\SCT_2\chromaticity")
G0 = 1 + 270/1875.61294257; B0 = np.sqrt(1 - 1/G0**2); F = 1 + 1/G0
pd.set_option("display.precision", 6); pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 70)
def show(path, pat=None, first=None):
    """print raw lines of an output file (optionally only lines matching regex pat)"""
    p = R/path; print(f"----- {p} -----")
    lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
    if pat: lines = [l for l in lines if re.search(pat, l)]
    for l in lines[:first]: print(l)
def kv(path):
    d = {}
    for l in (R/path).read_text(encoding="utf-8-sig").splitlines():
        m = re.match(r"\s*(\w+)\s*=\s*([-+0-9.Ee]+)\s*$", l)
        if m: d[m[1]] = float(m[2])
    return d
def tfs(path):
    lines = (R/path).read_text().splitlines()
    hdr = next(l for l in lines if l.startswith("*")).split()[1:]
    rows = [l.split() for l in lines if not l.startswith(("@", "*", "$"))]
    return pd.DataFrame(rows, columns=hdr)
ROWS = []   # (program, method, Qx, Qy, alpha_c, eta, xi_x, xi_y)
def add(*r):
    ROWS.append(r); print("-> row:", r)
print(f"gamma0={G0:.9f} beta0={B0:.9f} 1+1/gamma0={F:.9f}")
from IPython.display import Code, display
def src(path, lang, a=None, b=None):
    """show source file (lines a..b, 1-based) with syntax highlighting"""
    p = pathlib.Path(path); p = p if p.is_absolute() else R/p
    L = p.read_text(encoding="utf-8", errors="replace").splitlines()
    print(f"{p}  (lines {a or 1}-{b or len(L)})")
    display(Code("\n".join(L[(a or 1)-1:b]), language=lang))''')

# ---------------- OptiM
md(r"""## 1. OptiM (эталон, GUI)
**Запуск:** `OptiM\run_optim.bat` открывает `magnetic_2_md.opt` в OptiMX (только GUI, без пакетного режима), Twiss — вручную.
**Вывод:** окно OptiM, переписано пользователем в `OptiM/output/optim_reference.txt`.
**Берём:** Tunes, Momentum compaction, Phase Slip Factor, Chromaticities. OptiM уже даёт ξ = dQ/d(dp/p) — пересчёт не нужен.""")
code(r'''show("OptiM/output/optim_reference.txt")
t = (R/"OptiM/output/optim_reference.txt").read_text(encoding="utf-8-sig")
g = lambda p: float(re.search(p, t)[1])
add("OptiM", "GUI Twiss (linear formulas, hard-edge bends)", g(r"Qx =\s*([-\d.]+)"), g(r"Qy =\s*([-\d.]+)"),
    g(r"compaction :\s*([-\d.]+)"), g(r"Slip Factor :\s*([-\d.]+)"), g(r"Chromaticities :\s*([-\d.]+)"), g(r"\(hor\)\s*([-\d.]+)"))''')

# ---------------- BELL
md(r"""## 2. BELL (Julia)
**Запуск:** `BELL\run_bell.bat` = `julia --project=C:\GIT_REPOS\BELL chrom_bell.jl` (BELL v0.10.0 как есть).
Читает .opt через `convert_optim` → `magnetic_2_md.bell`. Вывод: `BELL/output/bell_results.txt`, `bell_twiss.csv`.
δ в BELL = dp/p (канонический), поэтому ξ берётся без пересчёта. Модель диполя: только тело, без краевого члена.""")
code(r'''show("BELL/output/bell_results.txt")
bell = kv("BELL/output/bell_results.txt")''')
md(r"""**Код BELL.** Скрипт запуска `chrom_bell.jl` и функции ядра BELL, которые реально считают хроматичность
(`C:\GIT_REPOS\BELL\src\higher_order.jl`: `closed_orbit`, `_frac_tune`, `tunes_at_delta`, `path_expansion`, `chromaticity_orders`; только чтение).""")
code(r'''src("BELL/chrom_bell.jl", "julia")
src(r"C:\GIT_REPOS\BELL\src\higher_order.jl", "julia", 34, 106)''')

md(r"""### 2a. BELL, одно-оборотная матрица (линейно)
Строки `linear_Qx/Qy` (целая часть из набега фазы), `linear_alpha_c`, `linear_eta`. Хроматичность этим методом не считается.""")
code(r'''add("BELL", "one-turn matrix (linear)", bell["linear_Qx"], bell["linear_Qy"], bell["linear_alpha_c"], bell["linear_eta"], np.nan, np.nan)''')
md(r"""### 2b. BELL, отображение Тейлора 3-го порядка (`ring_taylor_map`)
`taylor_Qx_frac` — дробный тюн (+2 целая часть), `taylor_alpha0`, `taylor_eta0` — из `path_expansion`,
`taylor_xi_x/y` — `chromaticity_orders`: центральная разность ±1e-4 тюнов якобиана на замкнутой орбите при δ.""")
code(r'''add("BELL", "Taylor map o3, Q(δ) central diff", 2+bell["taylor_Qx_frac"], 2+bell["taylor_Qy_frac"], bell["taylor_alpha0"], bell["taylor_eta0"], bell["taylor_xi_x"], bell["taylor_xi_y"])''')
md(r"""### 2c. BELL, Ли-фактор (`observables`: `ring_lie_factor` → `factor_to_map`)
Те же читалки, другое построение одно-оборотного отображения. Строки `lie_*`.""")
code(r'''add("BELL", "Lie factor map, Q(δ) central diff", 2+bell["lie_Qx"], 2+bell["lie_Qy"], bell["lie_alpha0"], bell["lie_eta0"], bell["lie_xi_x"], bell["lie_xi_y"])''')

# ---------------- COSY
md(r"""## 3. COSY INFINITY 10.2
**Запуск:** `COSY\run_cosy.bat`: `lattice_md.fox` (решётка, один раз → `lattice_md.bin`), `ov321.fox`, `ov330.fox`, `twiss.fox`,
(`rf_test.fox` по аргументу), затем `python ov330_map_analysis.py`. Элементы: `DL`, `MS` (секторный, 22.5°, без краевых углов),
`MQ`, `MH` (0), `FR 0` (жёсткий край без fringe). Переменная энергии COSY: δ_K = ΔK/K.""")
code(r'''src("COSY/lattice_md.fox", "pascal")''')
code(r'''src("COSY/run_cosy.bat", "bat")''')

md(r"""### 3a. COSY OV 3 2 1, TP
`ov321.fox`: `OV 3 2 1`, `PART 270*PARA(1)` (энергия — DA-параметр), `TP MU`. Берём `CONS(MU)` = дробный тюн и
`MU|(0,0,0,0,1)` = dQ/dδ_K (сырое), ξ = (1+1/γ)·raw.""")
code(r'''src("COSY/ov321.fox", "pascal"); show("COSY/output/cosy_ov321.txt")
c1 = kv("COSY/output/cosy_ov321.txt")
print(f"xi_x = {c1['factor_1p1g']:.9f} * {c1['dQx_ddK_raw']:.9f} = {c1['xi_x_p']:.9f};  xi_y = {c1['factor_1p1g']:.9f} * {c1['dQy_ddK_raw']:.9f} = {c1['xi_y_p']:.9f}")
add("COSY", "OV 3 2 1, TP, ×(1+1/γ)", 2+c1["Qx_frac"], 2+c1["Qy_frac"], np.nan, np.nan, c1["xi_x_p"], c1["xi_y_p"])''')
md(r"""### 3b. COSY OV 3 3 0: дамп 6D-отображения и коэффициенты
`ov330.fox`: `OV 3 3 0`, `PART 270`, `PM` → `map_ov330.txt`, коэффициенты MAP(1..6) полной точности → `map_ov330_full.txt`.
`ov330_map_analysis.py` строит из них: линейный блок, хроматические члены 2-го порядка, (i) замкнутую орбиту при δ_K и тюны якобиана,
(ii) 4D-отображение с δ_K как параметром, собственные числа 6×6, eta и alpha_c из (l|δ_K) на замкнутой орбите.
TP здесь не работает: продольный блок [[1, 17.836],[0,1]] — жорданова клетка (собственные числа 1,1), 6D нормальной формы нет.""")
code(r'''src("COSY/ov330.fox", "pascal"); src("COSY/ov330_map_analysis.py", "python"); show("COSY/output/map_ov330.txt", first=24)''')
code(r'''show("COSY/output/ov330_map_analysis.txt")
oa = json.load(open(R/"COSY/output/ov330_map_analysis.json"))
add("COSY", "OV 3 3 0 map → closed orbit(δK) + 2×2 Jacobian tunes, FD, ×(1+1/γ)", 2+oa["Qx_frac"], 2+oa["Qy_frac"], oa["alpha_c"], oa["eta"], oa["closed_orbit_FD"][0]*F, oa["closed_orbit_FD"][1]*F)
add("COSY", "OV 3 3 0 map → 4D map, δK as parameter, dTr/dδK, ×(1+1/γ)", 2+oa["Qx_frac"], 2+oa["Qy_frac"], oa["alpha_c"], oa["eta"], oa["param4D_x"]*F, oa["param4D_y"]*F)''')
md(r"""### 3c. COSY OV 3 3 0 + тонкий ВЧ (`rf_test.fox`)
h = 1, фаза 0 (устойчива; ниже критической). TP работает: Q1 = Qx, Q2 = Qs ∝ √V, Q3 = Qy. Линейной хроматичности в TP нет:
δ — координата синхротронного движения, а не параметр. Для таблицы не используется.""")
code(r'''src("COSY/rf_test.fox", "pascal"); show("COSY/output/rf_tp_ov330.txt")''')

# ---------------- MAD-X
md(r"""## 4. MAD-X 5.09.01
**Запуск:** `MADX\run_madx.bat`: `chrom_plain.madx` и `chrom_ptc.madx` (madx.exe из `C:\GIT_REPOS\NICA\MADX\2.Programms\Twiss`).
Решётка `magnetic_2_md.seq`: SBEND angle=22.5°, E1=E2=0; K1 = G/Bρ. В MAD-X TWISS DQ = dQ/dPT (PT = β0δ) ⇒ ξ = β0·DQ; DELTAP — это dp/p.""")
code(r'''src("MADX/magnetic_2_md.seq", "text"); src("MADX/chrom_plain.madx", "text"); show("MADX/output/madx_plain_results.txt")
mp = kv("MADX/output/madx_plain_results.txt")''')
md(r"""### 4a. TWISS (без CHROM): Q1, Q2, ALFA из SUMM; DQ1, DQ2 сырые (по PT) ⇒ ×β0""")
code(r'''print(f"xi = beta0*DQ: {B0:.9f}*{mp['DQ1_pt']:.6f} = {B0*mp['DQ1_pt']:.6f};  {B0:.9f}*{mp['DQ2_pt']:.6f} = {B0*mp['DQ2_pt']:.6f}")
add("MAD-X", "TWISS, β0·DQ(PT)", mp["Q1"], mp["Q2"], mp["alfa"], mp["eta"], mp["DQ1_dp"], mp["DQ2_dp"])''')
md(r"""### 4b. TWISS, CHROM (DQ конечной разностью внутри MAD-X) ⇒ ×β0""")
code(r'''add("MAD-X", "TWISS CHROM, β0·DQ(PT)", mp["Q1"], mp["Q2"], mp["alfa"], mp["eta"], mp["DQ1_chrom_dp"], mp["DQ2_chrom_dp"])''')
md(r"""### 4c. TWISS при DELTAP = ±1e-4 (dp/p): ξ = (Q(+) − Q(−))/2e-4 — без пересчёта""")
code(r'''add("MAD-X", "TWISS DELTAP=±1e-4 (dp/p) FD", mp["Q1"], mp["Q2"], mp["alfa"], mp["eta"], mp["DQ1_fd_dp"], mp["DQ2_fd_dp"])''')
md(r"""### 4d. PTC_TWISS (метод `Compaction.madx`: model=2 method=6 nst=20 exact, time=false ⇒ δ = dp/p), icase=56 no=5
Берём из `ptc_twiss_summary.tfs`: Q1, Q2 (дробные), ALPHA_C, ETA_C, DQ1, DQ2.""")
code(r'''src("MADX/chrom_ptc.madx", "text")
ps = tfs("MADX/output/ptc_twiss_summary.tfs").iloc[0]
print(ps[["LENGTH","ALPHA_C","ALPHA_C_P","ALPHA_C_P2","ETA_C","GAMMA_TR","Q1","Q2","DQ1","DQ2"]].to_string())
ps = ps.apply(pd.to_numeric, errors="ignore")
add("MAD-X PTC", "ptc_twiss icase=56 no=5 exact", 2+float(ps["Q1"]), 2+float(ps["Q2"]), float(ps["ALPHA_C"]), float(ps["ETA_C"]), float(ps["DQ1"]), float(ps["DQ2"]))''')
md(r"""### 4e. PTC_NORMAL icase=5 no=5: DQ порядков 1–3 (n-е производные по dp/p), ANHX/ANHY""")
code(r'''show("MADX/output/ptc_normal.tfs", pat=r'^ "')
pn = {(x[0].strip('"'), int(x[1])): float(x[-1]) for x in tfs("MADX/output/ptc_normal.tfs").values.tolist()}
add("MAD-X PTC", "ptc_normal icase=5 no=5", 2+pn[("Q1",0)], 2+pn[("Q2",0)], np.nan, np.nan, pn[("DQ1",1)], pn[("DQ2",1)])''')

# ---------------- Python exact
md(r"""## 5. Python: точный гамильтониан, RK4 (`analysis/exact_hamiltonian_check.py`)
Секторный диполь с фиксированным полем, квадруполи, дрейфы; замкнутая орбита Ньютоном при δ = ±1e-4, тюны якобиана.
Вариант 1 — только тело; вариант 2 — плюс краевые удары жёсткого края Δp_y = ∓h·x'·y (вход/выход).""")
code(r'''src("analysis/exact_hamiltonian_check.py", "python")
show("analysis/exact_hamiltonian_check.out")
ex = (R/"analysis/exact_hamiltonian_check.out").read_text().splitlines()
nums = lambda s: [float(v) for v in re.findall(r"-?\d+\.\d+", s)]
q0 = nums(ex[0]); xb = nums(ex[1]); xe = nums(ex[2])[-2:]
add("Python", "exact H, RK4, body only (no edge kick), FD", 2+q0[0], 2+q0[1], np.nan, np.nan, *xb)
add("Python", "exact H, RK4, body + hard-edge kick Δp_y=∓h x' y, FD", 2+q0[0], 2+q0[1], np.nan, np.nan, *xe)''')

# ---------------- summary
md(r"""## 5b. Метод пользователя (`COSY\jobs\delta_eq.fox`, ноутбук `COSY\notebooks\delta_eq.ipynb`) против нашего

**Как получается «ξ процедура TP» = −1.8408 / −1.8477.** `delta_eq.fox` (только чтение): `OV 3 2 1` → `SET_FOR_DEUTERONS_CHROM GAMMA`
(`header.fox`: масса 2.01410 а.е.м. — масса *атома* дейтерия, γ = 1.143914, `RP EKIN*PARA(1)`) → `LATTICE` (`structures/magnetic_2.fox`,
конвертер OptiM→COSY) → `TP MU_TP` → ξ = (MU_TP|5-я переменная)·(1+1/γ). Пересчёт (1+1/γ) верный (5-я переменная в OV 3 2 1 — параметр δK).
FR не задаётся (по умолчанию = `FR 0`, проверено: тот же результат). ВЧ в этом шаге нет (RF = 0; полость добавляется позже только для
трекинга и DANF — отсюда Qs).

**Причина:** в `magnetic_2.fox` квадруполи вызываются `QUAD L_QF A QF`, а `elements.fox` объявляет `PROCEDURE QUAD L TILT K1`.
Апертура A = 0.05 попадает в TILT ⇒ **все квадруполи повёрнуты на 0.05 рад (2.86°)**. Отображение связано (M13 = 0.052, M14 = 1.004,
det 2×2 = 0.991). TP даёт собственные тюны связанного движения 0.44628 / 0.47668, а хроматичности двух нормальных мод почти равны
(−1.84 / −1.85) — смесь ξx и ξy при близких тюнах. «Линейные Qx, Qy» (0.45992 / 0.45741) взяты из следов 2×2 блоков связанного
отображения, поэтому тоже не совпадают с Julia; «ξ из MAP(5) на эллипсе» (−1.4726 / −2.2303) близки к правильным, но тоже сдвинуты связью.

**Проверка** (`COSY/user_lattice_check/`, его bin-файлы только читаются): его расчёт на его решётке −1.840801 / −1.847679 (воспроизведено);
та же решётка с `QUAD L 0 K1` (копия `magnetic_2_notilt.fox`) → Q = 0.461720 / 0.461337, ξ = **−1.461536 / −2.225485** — ровно наш результат.
`COSY/ov321_user_method.fox` — его процедура на нашей `lattice_md` → ξ = −1.460889 / −2.224733; отличие от −1.461536 / −2.225485 только из-за
массы атома (1876.12 вместо 1875.61 МэВ ⇒ другой импульс при 270 МэВ).

**Наш `ov321.fox`:** `OV 3 2 1; FR 0; PART 270*PARA(1)` (масса 1875.61294257/AMUMEV) → `RING` → `TP MU` → ξ = (1+1/γ)·(MU|(0,0,0,0,1)) — тот же метод,
но решётка без поворота квадруполей и точная масса дейтрона.""")
code(r'''src(r"C:\GIT_REPOS\SCT_2\COSY\jobs\delta_eq.fox", "pascal", 187, 196)
src(r"C:\GIT_REPOS\SCT_2\COSY\src\header.fox", "pascal", 53, 67)
src(r"C:\GIT_REPOS\SCT_2\COSY\src\elements.fox", "pascal", 61, 76)
print([l.strip() for l in open(r"C:\GIT_REPOS\SCT_2\COSY\structures\magnetic_2\magnetic_2.fox") if "QUAD L_QF " in l][:1])
src("COSY/user_lattice_check/ov321_user_lattice.fox", "pascal")
show("COSY/user_lattice_check/user_lattice.out", pat="RES"); show("COSY/user_lattice_check/user_lattice_notilt.out", pat="RES")
src("COSY/ov321_user_method.fox", "pascal"); show("COSY/output/cosy_ov321_user_method.txt")
um = kv("COSY/output/cosy_ov321_user_method.txt")
add("COSY", "user method (mass 2.01410 amu, EKIN*PARA(1), TP) on lattice_md", 2+um["Qx_frac"], 2+um["Qy_frac"], np.nan, np.nan, um["xi_x_p"], um["xi_y_p"])
pd.DataFrame([("user table: TP on his magnetic_2 (quads tilted 0.05 rad)", -1.8408, -1.84768),
  ("his magnetic_2, tilt removed, same TP", -1.461536, -2.225485),
  ("user method on our lattice_md", um["xi_x_p"], um["xi_y_p"]), ("our ov321.fox", -1.461536, -2.225485),
  ("user table: MAP(5) on ellipse (his lattice)", -1.47261, -2.23029), ("user table: Julia/BELL", -1.4624, -2.2267)], columns=["case", "xi_x", "xi_y"])''')

md(r"""## 6. Сводная таблица (ξ по dp/p)""")
code(r'''df = pd.DataFrame(ROWS, columns=["program", "model / method", "Qx", "Qy", "alpha_c", "eta", "xi_x", "xi_y"])
df''')
md(r"""**Итог:** Qx, Qy, alpha_c, eta, ξx совпадают везде. ξy: −1.8171 у моделей с жёстким краем (OptiM, MAD-X, PTC, Python + край)
и −2.2255 у моделей «только тело» (BELL, COSY FR 0, Python без края). Разница +0.408 — краевой член (раздел 8).""")

# ---------------- Twiss
md(r"""## 7. Twiss-функции: то же ли это кольцо?
Источники (D всюду приведена к dp/p):
* **BELL** `BELL/output/bell_twiss.csv` — `optics_along`, шаг 2 см, D на dp/p.
* **COSY** `COSY/output/cosy_twiss.txt` — `twiss.fox`: периодические β, α, D из одно-оборотной матрицы, затем после каждого элемента
  β(s) = M11²β0 − 2M11M12α0 + M12²γ0, D(s) = M11D0 + M12D0' + M16; D на δ_K ⇒ ×(1+1/γ).
* **MAD-X TWISS** `MADX/output/twiss_plain.tfs` — концы элементов; DX по PT ⇒ ×β0.
* **PTC_TWISS** `MADX/output/ptc_twiss.tfs` — BETA11, BETA22, DISP1 (time=false ⇒ по dp/p).
* **OptiM** — таблицы Twiss вдоль кольца нет: OptiM только GUI, экспорта в репозитории нет. Есть лишь начальные значения из .opt
  (βx0 = 19.3675 м, βy0 = 4.05447 м, Dx0 = 0 — это начальное условие из файла, а не периодическое решение; периодическое D(0) = 0.0111 м).
  Чтобы добавить OptiM: в GUI вывести Twiss-таблицу (окно Twiss → сохранить в файл) в `OptiM/output/optim_twiss.txt`.""")
code(r'''bt = pd.read_csv(R/"BELL/output/bell_twiss.csv")
ct = pd.read_csv(R/"COSY/output/cosy_twiss.txt", sep=r"\s+", comment="#", header=None, names=["name","s","betx","bety","dx"]); ct["dx"] *= F
mt = tfs("MADX/output/twiss_plain.tfs"); mt[["S","BETX","BETY","DX"]] = mt[["S","BETX","BETY","DX"]].astype(float); mt["DXp"] = mt["DX"]*B0
pt = tfs("MADX/output/ptc_twiss.tfs"); pt[["S","BETA11","BETA22","DISP1"]] = pt[["S","BETA11","BETA22","DISP1"]].astype(float)
print("raw first rows:"); print(bt.head(3).to_string()); print(ct.head(3).to_string()); print(mt[["NAME","S","BETX","BETY","DX"]].head(3).to_string()); print(pt[["NAME","S","BETA11","BETA22","DISP1"]].head(3).to_string())
print(f"MAD-X DX(0)={mt.DX[0]:.6e} (per PT) * beta0 = {mt.DXp[0]:.6e};  PTC DISP1(0)={pt.DISP1[0]:.6e};  BELL D(0)={bt.Dx_m[0]:.6e};  COSY D(0)={ct.dx[0]:.6e}")''')
md(r"""### Рис. 1. β_x, β_y, D_x вдоль кольца (все программы)
Линии — BELL (мелкий шаг), маркеры — COSY / MAD-X / PTC в концах элементов; × — начальные значения OptiM из .opt.
Оси в метрах. Если кольцо одно и то же, все маркеры лежат на линиях BELL.""")
code(r'''fig, ax = plt.subplots(3, 1, figsize=(13, 10), sharex=True)
for a, (col_b, col_c, col_m, col_p, lab) in zip(ax, [("betx_m","betx","BETX","BETA11","β_x [m]"), ("bety_m","bety","BETY","BETA22","β_y [m]"), ("Dx_m","dx","DXp","DISP1","D_x [m] (per dp/p)")]):
    a.plot(bt.s_m, bt[col_b], "k-", lw=1, label="BELL")
    a.plot(ct.s, ct[col_c], "o", ms=5, mfc="none", label="COSY")
    a.plot(mt.S, mt[col_m], "s", ms=3, label="MAD-X TWISS")
    a.plot(pt.S, pt[col_p], "+", ms=6, label="PTC_TWISS")
    a.set_ylabel(lab); a.grid(alpha=.3)
ax[0].plot([0], [19.3675], "rx", ms=10, label="OptiM initial (.opt)"); ax[1].plot([0], [4.05447], "rx", ms=10); ax[2].plot([0], [0.0], "rx", ms=10)
ax[0].legend(ncol=5); ax[2].set_xlabel("s [m]"); ax[0].set_title("magnetic_2_md: Twiss functions, all programs"); plt.tight_layout(); plt.show()''')
md(r"""### Рис. 2. Разности относительно MAD-X TWISS (в точках MAD-X, концы элементов)
BELL интерполирован линейно на s MAD-X (точки MAD-X совпадают с границами элементов BELL); COSY и PTC также интерполированы по s (их точки — те же концы элементов, поэтому интерполяция фактически не нужна, кроме совпадающих s у маркеров).""")
code(r'''s = mt.S.values
def at(sx, y): return np.interp(s, sx, y)
cmp = {
 "BELL": (at(bt.s_m, bt.betx_m), at(bt.s_m, bt.bety_m), at(bt.s_m, bt.Dx_m)),
 "COSY": (at(ct.s, ct.betx), at(ct.s, ct.bety), at(ct.s, ct.dx)),
 "PTC":  (at(pt.S, pt.BETA11), at(pt.S, pt.BETA22), at(pt.S, pt.DISP1)),
}
ref = (mt.BETX.values, mt.BETY.values, mt.DXp.values)
print("lengths: MAD-X", len(s), "COSY", len(ct), "PTC", len(pt))
fig, ax = plt.subplots(3, 1, figsize=(13, 8), sharex=True); stat = []
for k, (bx, by, dx) in cmp.items():
    d = (bx-ref[0], by-ref[1], dx-ref[2])
    for a, dd, lab in zip(ax, d, ["Δβ_x [m]", "Δβ_y [m]", "ΔD_x [m]"]): a.plot(s, dd, ".-", ms=3, label=k); a.set_ylabel(lab); a.grid(alpha=.3)
    stat.append((k, *[np.max(np.abs(x)) for x in d], np.max(np.abs(d[0]/ref[0])), np.max(np.abs(d[1]/ref[1]))))
ax[0].legend(); ax[2].set_xlabel("s [m]"); ax[0].set_title("difference vs MAD-X TWISS"); plt.tight_layout(); plt.show()
pd.DataFrame(stat, columns=["program", "max|Δβx| [m]", "max|Δβy| [m]", "max|ΔDx| [m]", "max rel Δβx", "max rel Δβy"])''')
md(r"""**Что следует:** если максимальные разности на уровне ~1e-6 м и меньше (округление вывода/шаг), то BELL, COSY, MAD-X и PTC
считают одно и то же кольцо с одной и той же линейной оптикой; расхождение ξy не связано с геометрией или фокусировкой,
а только с моделью края диполя. OptiM: совпадение βx0 = 19.3675 м, βy0 = 4.05447 м с периодическими значениями (19.36752, 4.05447) подтверждает ту же
оптику в начале кольца; полной таблицы OptiM нет (см. выше).""")

md(r"""## 8. Почему ξy: OptiM/MAD-X = −1.817, а BELL/COSY = −2.2255

### Гамильтониан секторного магнита (поле B0 однородно, h = 1/ρ0, δ = dp/p, канонические x, p_x, y, p_y)
$$H = -(1+hx)\sqrt{(1+\delta)^2-p_x^2-p_y^2} + hx + \tfrac12 h^2x^2 .$$
Разложение до 3-го порядка по (x, p, δ):
$$H \approx -(1+\delta)(1+hx) + \frac{(1+hx)(p_x^2+p_y^2)}{2(1+\delta)} + hx + \tfrac12h^2x^2 .$$
* $-(1+\delta)hx + hx = -hx\delta$ — источник дисперсии.
* $\frac{p^2}{2(1+\delta)} \approx \frac{p^2}{2} - \delta\frac{p^2}{2}$ — «обычная» хроматичность (как в дрейфах и квадруполях).
  Так как $\oint\gamma\,ds = \oint K\beta\,ds$ для периодических функций, это даёт член **A** $= -\frac1{4\pi}\oint\gamma\,ds = -\frac1{4\pi}\oint K\beta\,ds$.
  Член $\tfrac12h^2x^2$ — от векторного потенциала, от δ не зависит (в канонических переменных).
* $\frac{hx\,p^2}{2}$ на дисперсионной орбите $x = D\delta + x_\beta$, $p_x = D'\delta + p_\beta$ даёт
  $\delta\big[\tfrac12 hD\,p_{x\beta}^2 + hD'\,x_\beta p_{x\beta} + \tfrac12 hD\,p_{y\beta}^2\big]$.
  Усреднение по фазе: $\langle p^2\rangle = \gamma J$, $\langle xp\rangle = -\alpha J$, $\Delta Q = \frac1{2\pi}\oint\partial\langle H_1\rangle/\partial J$:
  **C** $= +\frac1{4\pi}\int hD\gamma\,ds$ (x и y), **E** $= -\frac1{2\pi}\int hD'\alpha_x\,ds$ (только x).

### Жёсткий край при нулевом угле входа (тело это не содержит)
Поле $B_y = B_0\Theta(s)$ с плоской границей ⊥ орбите. Из $\nabla\times B = 0$: $\partial B_s/\partial y = \partial B_y/\partial s$ ⇒ $B_s = yB_0\delta_D(s)$.
Сила $F_y = q(v_sB_x - v_xB_s) = -q v_x y B_0\delta_D(s)$ ⇒ на входе $\Delta p_y = -h\,x'\,y$, на выходе $\Delta p_y = +h\,x'\,y$
($x' = p_x/p_s$). Это обычная краевая фокусировка $-h\tan e\,y$ с эффективным углом $e = x'$: частица на дисперсионной
орбите ($x' = D'\delta$) входит в магнит под углом. Тонкая линза $k = \pm hD'\delta$ ⇒
**F** $= +\frac{h}{4\pi}\sum_{BM}\big[\beta_yD'\big|_{in} - \beta_yD'\big|_{out}\big]$. По x: $F_x = q(v_yB_s - \dots)$ ∝ $y'y$ — не хроматично, поэтому ξx совпадает везде.

### Численно (оптика из линейных матриц, 400 срезов на диполь; D — на dp/p)""")
code(r'''import subprocess, sys
src("analysis/chrom_terms.py", "python")
print(subprocess.run([sys.executable, str(R/"analysis/chrom_terms.py")], capture_output=True, text=True, encoding="utf-8").stdout)''')
md(r"""### Численные проверки гипотез
* Python RK4 точного гамильтониана (`exact_hamiltonian_check.py`): только тело → ξy = −2.225485; + краевые удары $\Delta p_y=\mp h x' y$ → **−1.817063** (= MAD-X до 6 знаков). С обратным знаком ударов −2.634 (не совпадает).
* MAD-X (`../MADX/edge_tests.madx`): TWISS толстый SBEND −1.817060; PTC exact model 2, `fringe=false` −1.817063 (жёсткий край PTC всегда учитывает этот кинематический член);
  PTC `exact=false` −2.932 (без всех членов h: только A); MAKETHIN (тонкие диполи) −2.932 / −2.813 — тонкие удары теряют и тело C/E, и края F.""")
code(r'''print(open(R/"MADX/output/edge_tests.txt").read()); print(open(R/"analysis/exact_hamiltonian_check.out").read())''')
md(r"""### Вывод
| член | ξx | ξy | BELL / COSY FR 0 | OptiM / MAD-X / PTC |
|---|---|---|---|---|
| A: −(1/4π)∮γ ds | −2.8134 | −2.9322 | да | да |
| C: (1/4π)∫hDγ ds | +0.6615 | +0.7068 | да | да |
| E: −(1/2π)∫hD'α_x ds | +0.6905 | — | да | да |
| F: жёсткий край, Δp_y = ∓h x' y | 0 | +0.4084 | **нет** | да |
| сумма | −1.4615 | −2.2254 / −1.8170 | | |

* **Доказано:** разница +0.408 — ровно член F (аналитически 0.40842, численно в RK4 −2.225485 → −1.817063, совпадает с MAD-X/PTC до 1e-6).
  Он следует из уравнений Максвелла для поля с резкой границей (B_s на краю) и присутствует в любом реальном магните с краем,
  даже при нулевом угле входа.
* **Физически правильно** для секторного магнита с жёстким краем: **ξy = −1.817** (OptiM, MAD-X, PTC).
  BELL и COSY с `FR 0` считают только тело магнита (поле без края, не удовлетворяющее Максвеллу на границе), поэтому −2.2255.
  BELL и COSY действительно согласуются между собой, но в модели «без краёв».
* **Не проверено:** реализация полной матрицы TRANSPORT 2-го порядка (Brown) не делалась — член F проверен прямым выводом и RK4;
  BELL с добавленным краевым отображением не запускался (ядро BELL не меняли); COSY FR 1/2 не запускается (нет SYSCA.DAT),
  поэтому что именно COSY включает в FR 1–3 — не проверено. С конечным fringe (FINT, HGAP ≠ 0) добавятся ещё члены.""")


md(r"""## 9. Краевые поля в COSY: FR 3 (и FR 1, FR 2) — проверка гипотезы «жёсткого края»

**SYSCA.DAT.** На машине его не было (поиск по всем дискам: только `C:\Program Files\COSY 10.2` без него).
По Beam Physics Manual 10.2 (разд. 3.3.7) SYSCA.DAT нужен только для FR 2 (FR 3 и FR 2.5 интегрируют напрямую),
а создаётся штатной процедурой `CRSYSCA` из `cosy.fox` («This procedure produces the file SYSCA.DAT, which is shipped with the code»).
Запустили её как чёрный ящик: `COSY/sysca/mk_sysca.fox` → `SYSCA.DAT` (295 кБ, ~10 мин; предупреждения RK «minimum stepsize» и POTI error norm 3e-4).
Практически COSY 10.2 читает SYSCA.DAT и при FR 1 (без файла FR 1 падал), поэтому FR 1 тоже стал доступен только после этого.
Альтернатива — официальный SYSCA.DAT с сайта загрузки COSY INFINITY (bmtdynamics.org).

**Как устроено.** `lattice_md_fr.fox` — то же кольцо, но `RINGFR FD FQ GD`: `FR FD` перед каждым `MS`, `FR FQ` для квадруполей/секступолей,
`GD` — апертура диполя (масштаб функции Энге). (a) FR 3 только на диполях; (b) FR 3 на всём кольце; плюс FR 1, FR 2 и FR 3 только на квадруполях.
Каждый случай считается двумя способами: OV 3 2 1 + TP (`ov321_fr.fox`) и OV 3 3 0 → карта → замкнутая орбита при ±δK → след (`ov330_fr.fox` + `fr_analysis.py`).""")
code(r'''src("COSY/lattice_md_fr.fox", "pascal"); src("COSY/ov321_fr.fox", "pascal"); src("COSY/ov330_fr.fox", "pascal")
src("COSY/fr_analysis.py", "python"); src("COSY/sysca/mk_sysca.fox", "pascal")''')
code(r'''sys.path.insert(0, str(R/"COSY")); import importlib, fr_analysis; importlib.reload(fr_analysis)
subprocess.run([sys.executable, str(R/"COSY/fr_analysis.py")], check=True, capture_output=True)
fr = pd.DataFrame(json.load(open(R/"COSY/output/fr_analysis.json")))
fr["dxi_x(330-321)"] = fr.xi_x_CO330-fr.xi_x_TP321; fr["dxi_y(330-321)"] = fr.xi_y_CO330-fr.xi_y_TP321
pd.set_option("display.float_format", "{:.6f}".format); pd.set_option("display.width", 250)
fr[["FR_dip","FR_quad","gap_m","Qx_TP","Qy_TP","xi_x_TP321","xi_x_CO330","xi_y_TP321","xi_y_CO330","dxi_x(330-321)","dxi_y(330-321)"]]''')
md(r"""### 9a. Сходимость FR 3 (только диполи) к жёсткому краю при апертуре → 0
Функция Энге масштабируется апертурой D; при D→0 край становится жёстким, и ξy должен стремиться к MAD-X/PTC (−1.817063).""")
code(r'''g = fr[(fr.FR_dip==3)&(fr.FR_quad==0)].sort_values("gap_m")
k = np.polyfit(g.gap_m, g.xi_y_TP321, 1); print(f"linear fit xi_y(gap) = {k[1]:.6f} + {k[0]:.4f}*gap  -> extrapolated hard edge xi_y(0) = {k[1]:.6f};  MAD-X = -1.817063;  body only (FR 0) = -2.225485")
k2 = np.polyfit(g.gap_m, g.Qy_TP, 1); print(f"Qy(gap) -> {k2[1]:.6f} at gap=0;  FR 0: 0.461337;  MAD-X Qy_frac = {mp['Q2']-2:.6f}")
fig, ax = plt.subplots(1, 2, figsize=(12, 4))
ax[0].plot(g.gap_m*1e3, g.xi_y_TP321, "o-", label="COSY FR 3 dipoles, OV 3 2 1 TP")
ax[0].plot(g.gap_m*1e3, g.xi_y_CO330, "x", ms=10, label="COSY FR 3, OV 3 3 0 closed orbit")
ax[0].axhline(-1.817063, c="r", ls="--", label="MAD-X / PTC hard edge"); ax[0].axhline(-2.225485, c="k", ls=":", label="body only (FR 0)")
ax[0].set_xlabel("dipole aperture D, mm"); ax[0].set_ylabel("xi_y"); ax[0].legend()
ax[1].plot(g.gap_m*1e3, g.Qy_TP, "o-"); ax[1].axhline(0.461337, c="k", ls=":"); ax[1].set_xlabel("D, mm"); ax[1].set_ylabel("Qy frac")
plt.tight_layout(); plt.show()''')
md(r"""### 9b. Выводы по краевым полям
* FR 3 на диполях: ξy = −2.2255 → **−1.8077** (D = 50 мм), ξx почти не меняется (−1.46150). При D → 0 ξy → −1.8171 = MAD-X/PTC.
  Гипотеза подтверждена: разница BELL/COSY-FR0 и OptiM/MAD-X — краевой член диполя Δp_y = ∓h x′ y (в FR 0 COSY его нет), а не ошибка программ.
* Конечный край (D = 50 мм) ещё и уменьшает Qy (0.46134 → 0.45573: вертикальная фокусировка краевого поля ~ D·h²) — у жёсткого края её нет; отсюда же малая разница ξy с MAD-X (0.009 при D=50 мм, линейно по D).
* FR 3 на квадруполях: Qx,y падают на ~0.003, ξ меняется на ~0.009 (x и y) — эффект размазанного градиента, к гипотезе не относится.
* FR 1 ≈ FR 2 на диполях (одинаково, т.к. для `MS` FR 2 берёт ту же процедуру) совпадают с FR 3 до 1e-5; FR 1/2 на всём кольце — до 4e-5.
* **OV 3 3 0 и OV 3 2 1 дают одну и ту же хроматичность** (FR 0 и FR 3): разница ≤ 6e-6 — это ошибка конечной разности ±δK в методе замкнутой орбиты (3-й порядок карты), не физика.
  OV 3 3 0 сам по себе считается без проблем (и с FR 3); не работает только `TP` на 6D-карте без ВЧ (жорданов блок l–δ, нет продольной нормальной формы).""")
md(r"""## 10. 6D-карта с ВЧ (OV 3 3 0): можно ли достать линейную ξ? TS и метод неподвижной точки
1. **COSY `TS`** (нормальная форма) на OV 3 3 0 + тонкий резонатор: сдвиги частот от амплитуд (r_x², r_s², r_y²).
2. **OV 3 3 1 c `PART 270*PARA(1)`**: в 6D энергия уже переменная, параметр в карту не попадает (все коэффициенты по PARA = 0) → TP/TS дают ξ = 0. Не работает.
3. **Метод неподвижной точки** (Python, `rf_ts_analysis.py`): кольцевая карта COSY OV 3 3 0 + свой тонкий ВЧ-удар + постоянный сдвиг фазы за оборот Δl
   (= смещение частоты ВЧ) → неподвижная точка с δ* ≠ 0; Q — собственные частоты 6×6 якобиана; ξ = dQ/dδ*.""")
code(r'''src("COSY/rf_ts.fox", "pascal"); src("COSY/rf_ts_analysis.py", "python")
subprocess.run([sys.executable, str(R/"COSY/rf_ts_analysis.py")], check=True, capture_output=True, cwd=R/"COSY")
rr = json.load(open(R/"COSY/output/rf_ts_analysis.json"))
print(open(R/"COSY/output/rf_ts_ov331.txt").read()[:600])''')
code(r'''rows = []
for nm in ("FR0", "FR3dip"):
    for r in rr["fixed_point_"+nm]:
        q, d = r["Q_sorted"], r["dQ_ddelta_p_sorted"]
        rows.append(dict(lattice=nm, V_kV=r["V_kV"], Qs=q[0], Qy=q[1], Qx=q[2], xi_x_fixedpoint=d[2], xi_y_fixedpoint=d[1]))
fp = pd.DataFrame(rows)
ref = fr.set_index(["FR_dip","FR_quad","gap_m"])
fp["xi_x_OV321"] = np.where(fp.lattice=="FR0", ref.loc[(0,0,0.05),"xi_x_TP321"], ref.loc[(3,0,0.05),"xi_x_TP321"])
fp["xi_y_OV321"] = np.where(fp.lattice=="FR0", ref.loc[(0,0,0.05),"xi_y_TP321"], ref.loc[(3,0,0.05),"xi_y_TP321"])
fp["diff_x"] = fp.xi_x_fixedpoint-fp.xi_x_OV321; fp["diff_y"] = fp.xi_y_fixedpoint-fp.xi_y_OV321
pd.set_option("display.float_format", "{:.9g}".format); fp''')
code(r'''ts = pd.DataFrame(rr["TS"]); Qs = {100: 0.004097138126850692, 10: 0.001295596640668595, 1: 0.0004097026125950379, 0.1: 0.0001295593139383740, 0.01: 4.097025139526789e-05, 0.001: 1.295593108310074e-05}
ts["Qs_COSY"] = ts.V_kV.map(Qs)
# second-order chromaticity d2Qx/ddK2 from the FR0 OV 3 3 0 map (closed orbit, central 2nd difference)
b = fr_analysis.parse_maps("map_ov330_fr.txt")[0]; F = 1+1/b["g"]
def Qx_of(d): return fr_analysis.xi_closed_orbit(dict(b), h=1e-9)["Qx"] if d == 0 else None
maps = b["maps"]
def T(z, d):
    v = [z[0], z[1], z[2], z[3], 0.0, d]; return np.array([sum(c*np.prod([v[k]**e[k] for k in range(6)]) for e, c in maps[i].items()) for i in range(4)])
def Jm(z, d, e=1e-7): return np.array([(T(z+e*np.eye(4)[k], d)-T(z-e*np.eye(4)[k], d))/(2*e) for k in range(4)]).T
def Qx(d):
    z = np.zeros(4)
    for _ in range(40): z = z+np.linalg.solve(np.eye(4)-Jm(z, d), T(z, d)-z)
    m = Jm(z, d); c = (m[0,0]+m[1,1])/2; return np.arccos(c)/(2*np.pi)
h = 2e-3; q2 = (Qx(h)-2*Qx(0)+Qx(-h))/h**2
ts["beta_s=(l|dK)/sin2piQs"] = 17.836/np.sin(2*np.pi*ts.Qs_COSY)
ts["ratio dQx/dr_s2 : (Q''/2)/beta_s"] = ts.dQx_drs2/(0.5*q2/ts["beta_s=(l|dK)/sin2piQs"])
print(f"d2Qx/d(dK)^2 = {q2:.6f}  (second-order chromaticity in dK units)"); ts''')
md(r"""### 10a. Что работает и что нет
* **TS (нормальная форма 6D с ВЧ) линейную ξ не даёт.** Перекрёстный член dQx/d(r_s²) ∝ Q_s (≈3.07e-4 при 100 кВ → 9.7e-7 при 1 В):
  отношение к (Q″/2)/β_s постоянно (0.1540) на 5 порядках V, т.е. это **вторая** хроматичность, усреднённая по синхротронному колебанию;
  линейная ξ·⟨δ⟩ = 0 усредняется и не входит. При V → 0 этот член → 0, а не к ξ. При 10 кВ и 10 В TS вернул только константы (сбой нормальной формы — малые знаменатели).
  dQx/d(r_x²), dQx/d(r_y²) от V не зависят — обычные амплитудные сдвиги.
* **OV 3 3 1 с энергией как параметром не работает**: в 6D δ — фазовая переменная, `RP EK*PARA(1)` в карту не проходит (все коэффициенты по PARA = 0, TP → ξ = 0).
* **Работает метод неподвижной точки**: смещаем неподвижную точку по δ (сдвигом фазы ВЧ за оборот, т.е. частотой), берём тюны 6×6 якобиана.
  ξ совпадает с OV 3 2 1 до ~1e-8 **при любом V** (от 100 кВ до 1 В) — синхротронная связь линейную ξ не меняет; для FR 3 на диполях тоже (ξy = −1.80767).
  Q_s у нас = Q_s(COSY)/√(1+1/γ) — другая нормировка амплитуды удара (наш удар по δK = V/K); на ξ не влияет.""")


md(r"""## 11. Вывод: краевое поле диполя и FR 3 в COSY

Итог разделов 8–10.

* **Причина разницы ξy.** Жёсткий край диполя даёт краевой удар Δp_y = ∓h x′ y на входе и выходе, и он добавляет **+0.408** к ξy (разд. 8).
  OptiM, MAD-X и PTC его учитывают (−1.817). BELL и COSY FR 0 считают только тело магнита (−2.2255).
* **COSY FR 3 на диполях** (разд. 9): ξy = **−1.807670** при апертуре 50 мм. При апертуре → 0 линейная экстраполяция даёт **−1.817037**, у MAD-X −1.817063. Гипотеза жёсткого края подтверждена.
* FR 1 и FR 2 совпадают с FR 3 в пределах 4e-5. Краевые поля квадруполей (≈0.009 в ξ, ≈0.003 в Q) к эффекту не относятся.
* **OV 3 3 0 = OV 3 2 1:** разница ≤ 6e-6 при FR 0 и FR 3 (это ошибка конечной разности). Сам OV 3 3 0 считается без проблем; ВЧ нужна только для `TP`.
* **SYSCA.DAT** создан штатной процедурой `CRSYSCA` (`COSY/sysca/mk_sysca.fox`) с предупреждениями RK и POTI. Для эталона лучше взять официальный файл с сайта COSY INFINITY.
* **6D с ВЧ** (разд. 10): `TS` даёт только усреднённую 2-ю хроматичность (∝ Q_s). При OV 3 3 1 параметр энергии в карту не входит (ξ = 0).
  Метод сдвига неподвижной точки совпадает с OV 3 2 1 до 1e-8 при V = 100 кВ … 1 В.

| Модель | Qy | ξx | ξy |
|---|---|---|---|
| COSY FR 0 (только тело), OV 3 2 1 / OV 3 3 0 | 0.461337 | −1.461536 / −1.461536 | −2.225485 / −2.225486 |
| COSY FR 3 диполи, D = 50 мм, OV 3 2 1 / OV 3 3 0 | 0.455729 | −1.461498 / −1.461499 | −1.807670 / −1.807671 |
| COSY FR 1 / FR 2 диполи, 50 мм | 0.455728 | −1.461505 | −1.807680 |
| COSY FR 3 диполи, D → 0 (экстраполяция) | 0.461330 | ≈ −1.461536 | −1.817037 |
| MAD-X / PTC (жёсткий край) | 0.461337 | −1.461536 | −1.817063 |
| COSY + ВЧ, сдвиг неподвижной точки (FR 0 / FR 3) | — | −1.461536 / −1.461498 | −2.225485 / −1.807670 |""")

nb["cells"] = c; nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nbf.write(nb, str(pathlib.Path(__file__).parent/"analysis.ipynb"))