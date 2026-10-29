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
 ("MAD-X", "TWISS, β0·DQ(PT)", mp["Q1"], mp["Q2"], mp["alfa"], mp["eta"], mp["DQ1_dp"], mp["DQ2_dp"]),
 ("MAD-X", "TWISS CHROM, β0·DQ(PT)", mp["Q1"], mp["Q2"], mp["alfa"], mp["eta"], mp["DQ1_chrom_dp"], mp["DQ2_chrom_dp"]),
 ("MAD-X", "TWISS DELTAP=±1e-4 (dp/p) FD", mp["Q1"], mp["Q2"], mp["alfa"], mp["eta"], mp["DQ1_fd_dp"], mp["DQ2_fd_dp"]),
 ("MAD-X PTC", "ptc_twiss icase=56 no=5 exact", 2+ps["Q1"], 2+ps["Q2"], ps["ALPHA_C"], ps["ETA_C"], ps["DQ1"], ps["DQ2"]),
 ("MAD-X PTC", "ptc_normal icase=5 no=5", 2+pn[("Q1",0)], 2+pn[("Q2",0)], NaN, NaN, pn[("DQ1",1)], pn[("DQ2",1)]),
 ("Python", "exact H, RK4, fixed field, hard edge, FD", 2+exq[0], 2+exq[1], NaN, NaN, *[float(v) for v in re.findall(r"[-\d.]+", ex.splitlines()[1].split("[")[1])]),
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
  * Разница 0.408. PTC с `exact=false` даёт −2.932 (без кривизны); вклад члена (1+hx)·p_y²/2 на дисперсионной
    орбите, (1/4π)∮hDγ_y ds = +0.707, объясняет −2.932 → −2.2255. Остаток 0.408 у MAD-X/OptiM, вероятно,
    — вклад краёв диполя (кинематика/fringe жёсткого края при нулевом угле), которого нет в моделях «только тело»
    (COSY FR 0, BELL, Python). Проверить fringe в COSY (FR 1/2) не удалось — нет SYSCA.DAT.
* COSY OV 3 3 0: TP не работает, т.к. без ВЧ продольный блок [[1,(l|δK)],[0,1]] — жорданова клетка с двойным собственным числом 1 (не вращение); нормальная форма 6D невозможна. Из дампа 6D-отображения все способы (замкнутая орбита + след, 4D с δK как параметром) дают тот же ξ, что OV 3 2 1 TP. С тонким ВЧ (фаза 0) TP работает и даёт Qx, Qy, Qs∝√V, но линейной хроматичности в нём нет: δ не инвариант, ξ1 усредняется синхротронными колебаниями.
* Определение δ (dp/p vs ΔK/K vs PT) даёт множители 1.87 и 0.486 — это основная ловушка при сравнении «сырых» чисел.
* Порядок: всё — первый порядок ξ = dQ/dδ; ξ2 горизонтальный совпадает у BELL и PTC (8.555), вертикальный — нет (8.10 vs 7.37).""")
nb["cells"] = c; nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nbf.write(nb, "analysis.ipynb")
