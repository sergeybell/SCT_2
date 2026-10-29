import nbformat as nbf
nb = nbf.v4.new_notebook(); c = []
md = lambda s: c.append(nbf.v4.new_markdown_cell(s)); code = lambda s: c.append(nbf.v4.new_code_cell(s))
md("""# Хроматичность magnetic_2_md: OptiM, BELL, COSY, MAD-X
Общее определение: **ξ = dQ/dδ, δ = dp/p**, натуральная (секступоли = 0), первый порядок.
Пересчёты: COSY δ_K = ΔK/K = (1+1/γ)δ ⇒ ξ = (1+1/γ)·dQ/dδ_K; MAD-X TWISS PT = β0δ ⇒ ξ = β0·DQ.
eta = ΔT/T/δ = alpha_c − 1/γ².""")
code(r'''import re, pathlib, pandas as pd
R = pathlib.Path(r"C:\GIT_REPOS\SCT_2\chromaticity")
def kv(p):
    d = {}
    for l in open(p, encoding="utf-8-sig"):
        m = re.match(r"\s*([\w]+)\s*=\s*([-+0-9.Ee]+)\s*$", l)
        if m: d[m[1]] = float(m[2])
    return d
def tfs(p):
    lines = open(p).read().splitlines()
    hdr = next(l for l in lines if l.startswith("*")).split()[1:]
    rows = [l.split() for l in lines if not l.startswith(("@","*","$"))]
    return hdr, rows
opt = open(R/"OptiM/output/optim_reference.txt", encoding="utf-8-sig").read()
f = lambda pat: float(re.search(pat, opt)[1])
optim = dict(Qx=f(r"Qx =\s*([-\d.]+)"), Qy=f(r"Qy =\s*([-\d.]+)"), ac=f(r"compaction :\s*([-\d.]+)"),
             eta=f(r"Slip Factor :\s*([-\d.]+)"), xx=f(r"Chromaticities :\s*([-\d.]+)"), xy=f(r"\(hor\)\s*([-\d.]+)"))
bell = kv(R/"BELL/output/bell_results.txt"); c1 = kv(R/"COSY/output/cosy_ov321.txt"); c3 = kv(R/"COSY/output/cosy_ov330.txt")
mp = kv(R/"MADX/output/madx_plain_results.txt")
import json; oa = json.load(open(R/"COSY/output/ov330_map_analysis.json")); F = c3["factor_1p1g"]
h, r = tfs(R/"MADX/output/ptc_twiss_summary.tfs"); ps = dict(zip(h, map(float, r[0])))
h, r = tfs(R/"MADX/output/ptc_normal.tfs"); pn = {(x[0].strip('"'), int(x[1])): float(x[-1]) for x in r}
ex = open(R/"analysis/exact_hamiltonian_check.out").read()
exq = [float(v) for v in re.findall(r"[-\d.]+", ex.splitlines()[0].split("[")[1])]
NaN = float("nan")
rows = [
 ("OptiM", "GUI Twiss, linear formulas", optim["Qx"], optim["Qy"], optim["ac"], optim["eta"], optim["xx"], optim["xy"]),
 ("BELL", "one-turn matrix (linear)", bell["linear_Qx"], bell["linear_Qy"], bell["linear_alpha_c"], bell["linear_eta"], NaN, NaN),
 ("BELL", "Taylor map o3, Q(δ) central diff", 2+bell["taylor_Qx_frac"], 2+bell["taylor_Qy_frac"], bell["taylor_alpha0"], bell["taylor_eta0"], bell["taylor_xi_x"], bell["taylor_xi_y"]),
 ("BELL", "Lie factor map, Q(δ) central diff", 2+bell["lie_Qx"], 2+bell["lie_Qy"], bell["lie_alpha0"], bell["lie_eta0"], bell["lie_xi_x"], bell["lie_xi_y"]),
 ("COSY", "OV 3 2 1, TP, ×(1+1/γ)", 2+c1["Qx_frac"], 2+c1["Qy_frac"], NaN, NaN, c1["xi_x_p"], c1["xi_y_p"]),
 ("COSY", "OV 3 3 0, dTr/dδK from ME on closed orbit (in .fox), ×(1+1/γ)", 2+c3["Qx_frac"], 2+c3["Qy_frac"], c3["alpha_c"], c3["eta"], c3["xi_x_p"], c3["xi_y_p"]),
 ("COSY", "OV 3 3 0 map dump → closed orbit(δK) + 2x2 Jacobian tunes, FD ±1e-4, ×(1+1/γ)", 2+c3["Qx_frac"], 2+c3["Qy_frac"], c3["alpha_c"], c3["eta"], oa["closed_orbit_FD"][0]*F, oa["closed_orbit_FD"][1]*F),
 ("COSY", "OV 3 3 0 map dump → 4D map, δK as parameter (drop l), ×(1+1/γ)", 2+c3["Qx_frac"], 2+c3["Qy_frac"], c3["alpha_c"], c3["eta"], oa["param4D_x"]*F, oa["param4D_y"]*F),
 ("Python", "exact H, RK4, body + hard-edge kick Δp_y=∓h x' y, FD", 2+exq[0], 2+exq[1], NaN, NaN, *[float(v) for v in re.findall(r"[-\d.]+", ex.splitlines()[2].split("xi")[1].strip(" []"))]),
 ("MAD-X", "TWISS, β0·DQ(PT)", mp["Q1"], mp["Q2"], mp["alfa"], mp["eta"], mp["DQ1_dp"], mp["DQ2_dp"]),
 ("MAD-X", "TWISS CHROM, β0·DQ(PT)", mp["Q1"], mp["Q2"], mp["alfa"], mp["eta"], mp["DQ1_chrom_dp"], mp["DQ2_chrom_dp"]),
 ("MAD-X", "TWISS DELTAP=±1e-4 (dp/p) FD", mp["Q1"], mp["Q2"], mp["alfa"], mp["eta"], mp["DQ1_fd_dp"], mp["DQ2_fd_dp"]),
 ("MAD-X PTC", "ptc_twiss icase=56 no=5 exact", 2+ps["Q1"], 2+ps["Q2"], ps["ALPHA_C"], ps["ETA_C"], ps["DQ1"], ps["DQ2"]),
 ("MAD-X PTC", "ptc_normal icase=5 no=5", 2+pn[("Q1",0)], 2+pn[("Q2",0)], NaN, NaN, pn[("DQ1",1)], pn[("DQ2",1)]),
 ("Python", "exact H, RK4, body only (no edge kick), FD", 2+exq[0], 2+exq[1], NaN, NaN, *[float(v) for v in re.findall(r"[-\d.]+", ex.splitlines()[1].split("[")[1])]),
]
df = pd.DataFrame(rows, columns=["program", "model / method", "Qx", "Qy", "alpha_c", "eta", "xi_x (dQ/dδp)", "xi_y (dQ/dδp)"])
pd.set_option("display.precision", 6); pd.set_option("display.width", 200)
df''')
md("## Сырые величины до пересчёта")
code(r'''pd.DataFrame([
 ("COSY OV321 dQ/dδ_K", c1["dQx_ddK_raw"], c1["dQy_ddK_raw"], c1["factor_1p1g"]),
 ("COSY OV330 dQ/dδ_K", c3["dQx_ddK_raw"], c3["dQy_ddK_raw"], c3["factor_1p1g"]),
 ("MAD-X DQ (dQ/dPT)", mp["DQ1_pt"], mp["DQ2_pt"], mp["beta0"]),
], columns=["quantity", "x", "y", "factor -> dp/p"])''')
md("## Высшие порядки (δ = dp/p)")
code(r'''pd.DataFrame([
 ("BELL Taylor d2Q/dδ2", bell["taylor_xi2_x"], bell["taylor_xi2_y"]),
 ("PTC normal DQ order 2", pn[("DQ1",2)], pn[("DQ2",2)]),
 ("PTC normal DQ order 3", pn[("DQ1",3)], pn[("DQ2",3)]),
 ("BELL alpha1 (dC/C=a0δ+a1δ²)", bell["taylor_alpha1"], float("nan")),
 ("PTC ALPHA_C_P (=dα_c/dδ=2α1)", ps["ALPHA_C_P"], float("nan")),
], columns=["quantity", "x", "y"])''')
md("""## Выводы
* **Qx, Qy, alpha_c, eta, ξx** совпадают во всех программах (ξx = −1.46154 ± 2e-6) после приведения к δ = dp/p.
  Это подтверждает пересчёты: COSY ×(1+1/γ)=1.874162, MAD-X TWISS ×β0=0.485635.
* **ξy расходится на две группы:**
  * OptiM, MAD-X TWISS (все 3 способа), PTC exact (twiss и normal): **ξy = −1.8171**;
  * BELL (Taylor и Ли), COSY (OV 3 2 1 и OV 3 3 0) и независимый Python-интегратор точного гамильтониана
    секторного магнита (фиксированное поле, жёсткий край, без fringe): **ξy = −2.2255**.
  * Разница +0.408 — краевой член жёсткого края (Δp_y = ∓h x' y на входе/выходе), который есть в OptiM/MAD-X/PTC и отсутствует в моделях «только тело» (BELL, COSY FR 0). Вывод и проверка — в разделе «Почему ξy…» ниже.
* COSY OV 3 3 0: TP не работает, т.к. без ВЧ продольный блок [[1,(l|δK)],[0,1]] — жорданова клетка с двойным собственным числом 1 (не вращение); нормальная форма 6D невозможна. Из дампа 6D-отображения все способы (замкнутая орбита + след, 4D с δK как параметром) дают тот же ξ, что OV 3 2 1 TP. С тонким ВЧ (фаза 0) TP работает и даёт Qx, Qy, Qs∝√V, но линейной хроматичности в нём нет: δ не инвариант, ξ1 усредняется синхротронными колебаниями.
* Определение δ (dp/p vs ΔK/K vs PT) даёт множители 1.87 и 0.486 — это основная ловушка при сравнении «сырых» чисел.
* Порядок: всё — первый порядок ξ = dQ/dδ; ξ2 горизонтальный совпадает у BELL и PTC (8.555), вертикальный — нет (8.10 vs 7.37).""")
md(r"""## COSY OV 3 3 0: коэффициенты отображения и почему TP не работает
Файл `../COSY/output/ov330_map_analysis.txt` (скрипт `../COSY/ov330_map_analysis.py`, дамп `map_ov330_full.txt`).
Методы (все ×(1+1/γ) к dp/p): (1) замкнутая орбита поперечного отображения при фиксированном δK + якобиан 2×2,
Q = acos(Tr/2)/2π, центральная разность; (2) 4D-отображение с δK как параметром: dTr/dδK = Σ(i|iδK) + подача с орбиты,
dQ/dδK = −dTr/(4π sin μ). Оба = OV 3 2 1 TP. TP в OV 3 3 0 падает: продольный блок [[1,(l|δK)=17.836 м],[0,1]] —
жорданова клетка, собственные числа 1,1, нет третьей пары e^{±iμ} ⇒ нет 6D нормальной формы. С тонким ВЧ (фаза 0)
TP работает (Qs ∝ √V), но линейной хроматичности в TP/TS нет: δ уже не параметр, а координата синхротронных колебаний.""")
code(r'''print(open(R/"COSY/output/ov330_map_analysis.txt", encoding="utf-8").read())
print("\n".join(l for l in open(R/"COSY/output/rf_tp_ov330.txt") if l.startswith(("CASE","TP_Q"))))''')
md(r"""## Почему ξy: OptiM/MAD-X = −1.817, а BELL/COSY = −2.2255

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
nb["cells"] = c; nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nbf.write(nb, "analysis.ipynb")
