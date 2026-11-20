# make_notebook.py - builds magnetic_ring.ipynb from scripts/*.jl and output/*; then
#   jupyter nbconvert --to notebook --execute --inplace magnetic_ring.ipynb
# RUN_JULIA=1 in the environment makes the notebook re-run every step script (short defaults).
import nbformat as nbf
nb = nbf.v4.new_notebook(); C = nb.cells
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))

md(r"""# Равновесная энергия в магнитном кольце `magnetic_2`: разбор шагов 1–7
Исследование `studies/magnetic_ring`: кольцо `magnetic_2` (дейтрон, 270 МэВ, $h=1$) на текущем API BELL, физика — главы 1–3
`notes/equilibrium-energy-v3`. Главная величина — сдвиг равновесного импульса
$$\delta_{eq}=-\frac{\eta_1\langle\delta^2\rangle+\Lambda}{\eta_0},\qquad \langle\delta^2\rangle=\frac{\delta_m^2}{2},\qquad
\Lambda=-\frac{2\pi}{C}\left(J_x\xi_x+J_y\xi_y\right),$$
и производные от него $\gamma_{eff}=\gamma+\beta^2\gamma\,\delta_{eq}+\beta^2\langle\delta^2\rangle/(2\gamma)$, $\Delta E/E$.
Оптика, хроматичность, slip и расстройка читаются из факторов Ли (гамильтониан → Драгт–Финн); трекинг и DA — только перекрёстная проверка.

Каждый шаг: физическая преамбула → скрипт шага и его отчёт → таблицы и рисунки из `output/` → выводы по рисункам.
По умолчанию тетрадь читает сохранённые результаты прогонов; `RUN_JULIA=1` перед запуском Jupyter перезапускает скрипты
(короткие значения по умолчанию, шаг 7 ≈ 13 мин) и печатает их вывод. Ни одно число здесь не подгоняется: всё прочитано из прогонов.""")
code(r"""import os, json, subprocess
from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from IPython.display import Markdown, SVG, display
OUT = Path('output'); RUN_JULIA = os.environ.get('RUN_JULIA') == '1'
pd.set_option('display.float_format', '{:.6g}'.format); pd.set_option('display.max_rows', 200)

def header(script):
    lines = Path('scripts', script).read_text(encoding='utf-8').splitlines()
    return '\n'.join(l for l in lines[:next(i for i, l in enumerate(lines) if not l.startswith('#'))])

def run_step(script, report):
    print(header(script), end='\n\n')
    if RUN_JULIA:
        cmd = ['julia', '--project=.', f'scripts/{script}']; print('$', ' '.join(cmd), flush=True)
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding='utf-8', errors='replace') as p:
            for line in p.stdout: print(line, end='', flush=True)
        if p.returncode: raise RuntimeError(f'{script} exited with {p.returncode}')
    else:
        print(f'RUN_JULIA не задан: показан сохранённый отчёт output/{report}.')
    text = (OUT / report).read_text(encoding='utf-8')
    display(Markdown(f'<details><summary><b>Отчёт <code>output/{report}</code></b></summary>\n\n{text}\n\n</details>'))

def kv(name):
    df = pd.read_csv(OUT / name); return dict(zip(df.quantity, df.value))

seq = pd.read_csv(OUT / 'sequence.csv'); seq['s1'] = seq.length_m.cumsum(); seq['s0'] = seq.s1 - seq.length_m; seq['sc'] = (seq.s0 + seq.s1) / 2
opt = pd.read_csv(OUT / '02_optics.csv'); at = lambda col, s: np.interp(s, opt.s_m, opt[col])
CIRC = seq.s1.iloc[-1]
FAM = {'SXD': seq[seq.name == 'SF1'], 'SXB': seq[seq.name == 'SF2'],
       'SYB': seq[seq.name.isin(['SD', 'OS1']) & (at('Dx_m', seq.sc) > 0.5)]}""")

# ---------------------------------------------------------------- step 1
md(r"""## 1. Решётка: замыкание поворота до $2\pi$
Замкнутая орбита существует, только если суммарный угол поворота равен $2\pi$. В исходном файле угол BM записан как
`0.3926990753258494`, что на $6.4\cdot10^{-9}$ рад меньше $\pi/8$; шестнадцать магнитов недогибают на $\Delta\theta\approx1.0\cdot10^{-7}$ рад.
Однородный недогиб эквивалентен сдвигу опорного импульса порядка $\Delta\theta/2\pi\approx1.6\cdot10^{-8}$ — того же порядка, что искомый
$\delta_{eq}\sim10^{-8}$. Поэтому угол каждого BM ставится ровно $\pi/8$ (длина дуги не меняется, $C$ тоже), и замкнутая решётка
сохраняется в `lattices/magnetic_2_closed.bell`.""")
code("run_step('01_lattice.jl', '01_lattice.md')")
code(r"""cl = pd.read_csv(OUT / 'closure.csv', index_col=0); display(cl)
dth = -cl.loc['total_bend_minus_2pi_rad', 'before']
print(f'недогиб Δθ = {dth:.3e} рад, Δθ/2π = {dth / (2 * np.pi):.3e}')
display(pd.read_csv(OUT / 'unique_elements.csv'))
display(SVG(filename=str(OUT / 'bend_deficit.svg')))""")
md(r"""**Рисунок 1.** Суммарный угол поворота минус $2\pi$ до и после замыкания: $-1.02\cdot10^{-7}$ рад → $1.8\cdot10^{-15}$ рад (округление Float64).

**Выводы.** После замыкания поворот равен $2\pi$ до машинной точности, $\Delta C=0$, $Q_x$ сдвинулся на $1.1\cdot10^{-8}$, $Q_y$ не изменился.
Без этой правки ошибка геометрии была бы сравнима с самим $\delta_{eq}$. Побочная находка: `write_structure` не сохраняет угол бит в бит
($-5.6\cdot10^{-17}$ рад на магнит), поэтому записан отредактированный объект, а не текст `write_structure` (BELL issue #64; угол — #63).""")

# ---------------------------------------------------------------- step 2
md(r"""## 2. Линейная оптика и выбор режимов ВЧ
Из фактора Ли и матрицы одного оборота: функции Твисса, тюны, дисперсия и коэффициенты slip
$$\alpha_0=\frac1C\oint\frac{D_x}{\rho}\,ds,\qquad \eta_0=\alpha_0-\frac1{\gamma^2},\qquad
\eta_1=\alpha_1-\frac{\alpha_0}{\gamma^2}+\frac{3\gamma^2-1}{2\gamma^4},\qquad \gamma_t=\alpha_0^{-1/2}.$$
При $\gamma=1.144<\gamma_t$ кольцо работает ниже критической энергии, $\eta_0<0$. В обозначениях BELL поздний приход — $l>0$,
поэтому стационарная фаза полости в BELL $\phi_s^B=\pi$ (в заметках $\phi_s^N=0$, словарь $\phi^B=-\phi^N$).
Синхротронный тюн и амплитуда фазы малых колебаний:
$$Q_s=\sqrt{\frac{h\,eV\,|\eta_0\cos\phi_s|}{2\pi\beta^2E}},\qquad \hat\phi=\frac{h|\eta_0|}{Q_s}\,\delta_m.$$
Два режима при $h=1$, $\delta_m=10^{-4}$: «линейный» $\hat\phi=0.05$ и «кубический» $\hat\phi=0.3$, где поправка
$1-\hat\phi^2/32$ к $\langle\delta^2\rangle$ (глава 3) различима. Оба удовлетворяют неравенствам заметок:
$|\eta_1\delta_m/\eta_0|$, $|\Lambda/(\eta_0\delta_m)|$, $\omega_sT_s\ll1$.""")
code("run_step('02_linear_optics.jl', '02_linear_optics.md')")
code(r"""n2 = kv('02_numbers.csv')
display(pd.Series({k: n2[k] for k in ['C_m', 'Qx', 'Qy', 'alpha0', 'alpha1', 'eta0', 'eta1', 'gamma_t']}, name='значение').to_frame())
fig, (ax, axe) = plt.subplots(2, 1, figsize=(13, 5.6), sharex=True, gridspec_kw=dict(height_ratios=[5, 1]))
ax.plot(opt.s_m, opt.betx_m, label=r'$\beta_x$'); ax.plot(opt.s_m, opt.bety_m, label=r'$\beta_y$'); ax.plot(opt.s_m, opt.Dx_m, 'k', label=r'$D_x$')
for (f, d), c in zip(FAM.items(), ['C3', 'C4', 'C2']):
    ax.plot(d.sc, at('Dx_m' if f == 'SXD' else 'betx_m' if f == 'SXB' else 'bety_m', d.sc), 'v', ms=9, color=c, label=f'слоты {f}')
ax.set_ylabel(r'$\beta$, $D_x$ [м]'); ax.grid(alpha=.3); ax.legend(ncol=6, fontsize=8, loc='upper center')
for _, e in seq.iterrows():
    if e.type == 'SBend': axe.axvspan(e.s0, e.s1, .3, .7, color='0.7')
    elif e.type == 'Quadrupole': axe.axvspan(e.s0, e.s1, *((.5, 1) if e.k1 > 0 else (0, .5)), color='C3' if e.k1 > 0 else 'C0')
    elif e.type == 'Sextupole': axe.axvline(e.sc, color='C2', lw=2)
axe.set_yticks([]); axe.set_xlabel('s [м]'); axe.set_xlim(0, CIRC); plt.tight_layout(); plt.show()""")
md(r"""**Рисунок 2.** Функции $\beta_x$, $\beta_y$ и дисперсия $D_x$ вдоль кольца (шаг 0.02 м, $k_2=0$). Внизу — схема: серые — BM, красные
вверх — QF, синие вниз — QD, зелёные — секступоли. Треугольники — слоты семейств шага 5: SXD в максимуме $D_x$, SXB в максимумах
$\beta_x$, SYB в максимумах $\beta_y$.""")
code(r"""cells = pd.read_csv(OUT / '02_cells.csv'); cells = cells[cells.kind == 'arc_cell']
fig, ax = plt.subplots(figsize=(9, 3.4)); w = .38
ax.bar(cells['index'] - w / 2, cells.mu_x_deg - 90, w, label=r'$\mu_x$'); ax.bar(cells['index'] + w / 2, cells.mu_y_deg - 90, w, label=r'$\mu_y$')
ax.axhline(np.degrees(cells.mu_x_periodic_rad.iloc[0]) - 90, color='C0', ls='--', lw=.8, label=r'периодическая $\mu_x$')
ax.axhline(np.degrees(cells.mu_y_periodic_rad.iloc[0]) - 90, color='C1', ls='--', lw=.8, label=r'периодическая $\mu_y$')
ax.set_xlabel('ячейка арки'); ax.set_ylabel(r'$\mu-90^\circ$ [град]'); ax.grid(alpha=.3); ax.legend(fontsize=8, ncol=2); plt.tight_layout(); plt.show()
j = json.load(open(OUT / '02_rf_regimes.json', encoding='utf-8')); reg = pd.DataFrame(j['regimes']).T
display(reg)
for name, r in j['regimes'].items():
    print(f"{name:12s} φ̂ = h|η0|δm/Qs = {1 * abs(n2['eta0']) * r['delta_m'] / r['Qs']:.4f}")""")
md(r"""**Рисунок 3.** Набег фазы на ячейку арки в кольце минус $90^\circ$ (столбцы) и периодический набег ячейки, взятой отдельно (пунктир).

**Выводы.** Оптика симметрична: две суперпериодические копии по 42 элемента, $\beta_{max}\approx19.4$ м в обеих плоскостях,
$D_{x,max}=11.70$ м в слоте SF1. Ячейки арки — практически $90^\circ$ (разброс $0.08^\circ$ по $x$ и $0.4^\circ$ по $y$, его создают прямые
промежутки). $Q_x=2.46172$ и $Q_y=2.46134$ почти равны: $|Q_x-Q_y|=3.8\cdot10^{-4}$ — это важно для шага 7. Ниже критической энергии
$\eta_0=-0.534$, а $\eta_1=+1.017$: $\delta_{eq}$ без секступолей положителен. Матрица, фактор Ли и `path_expansion` совпадают до $10^{-8}$.
Замечание: JSON режимов хранит фазу заметок $\phi_s^N=0$, а полость в BELL на шагах 3–7 стоит на $\pi$.""")

# ---------------------------------------------------------------- step 2b
md(r"""## 2b. Хроматичность из генераторов Ли
Фактор одного оборота без ВЧ $\mathcal M=e^{:f_2:}e^{:f_3:}e^{:f_4:}$ переводится в нормальные координаты матрицей $A$:
симплектический дисперсионный сдвиг ($x\to x+D\delta$, $p_x\to p_x+D'\delta$, $l\to l-D'x+Dp_x$), затем нормировка Куранта–Снайдера.
В этих координатах $f_2=-\mu_xJ_x-\mu_yJ_y+\dots$, а кубический генератор содержит $c\,J\delta$. Так как поворот на $+\mu$ — это генератор $-\mu J$,
$$\mu(\delta)=\mu-c\,\delta\;\Rightarrow\;\xi_1=\frac{dQ}{d\delta}=-\frac{c}{2\pi}.$$
Второй порядок даёт член $c_2J\delta^2$ в $f_4-\tfrac12[f_3,F_3]$, где $F_3$ решает гомологическое уравнение $F_3-F_3\circ r=-f_3$:
$d^2Q/d\delta^2=-c_2/\pi$. Независимая проверка — отклик одного секступоля $\Delta\xi_{x,y}=\pm k_2L\beta_{x,y}D_x/(4\pi)$.
`symplectic_eigen_A` тут не годится: продольный блок магнитной карты — сдвиг, а не вращение (BELL issue #65).""")
code("run_step('02b_generator_chromaticity.jl', '02b_chromaticity.md')")
code(r"""ch = pd.read_csv(OUT / '02b_chromaticity.csv'); display(ch)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.2))
d = np.linspace(-5e-3, 5e-3, 201); Q0 = {'x': n2['Qx'], 'y': n2['Qy']}
for case, ls in (('k2=0', '-'), ('one_sextupole', '--')):
    for p, c in (('x', 'C0'), ('y', 'C1')):
        g = ch[(ch.case == case) & (ch.plane == p)].set_index('quantity').generator
        a1.plot(d * 1e3, Q0[p] + g.xi1 * d + g.xi2 * d**2 / 2, ls, color=c, label=f'$Q_{p}$, ' + ('$k_2=0$' if case == 'k2=0' else 'один SF1, $k_2=1$'))
a1.set_xlabel(r'$\delta\cdot10^3$'); a1.set_ylabel(r'$Q_0+\xi_1\delta+\xi_2\delta^2/2$'); a1.grid(alpha=.3); a1.legend(fontsize=8)
steps = [1e-4, 1e-5, 1e-6]
for _, r in ch[ch.quantity.isin(['xi1', 'xi2'])].iterrows():
    err = [abs(r.generator - r[f'orders_dd_{h:.0e}'.replace('e-0', 'e-')]) for h in steps]
    a2.loglog(steps, np.maximum(err, 1e-12), 'o-' if r.quantity == 'xi1' else 's--', label=f'{r.case}, {r.plane}, {r.quantity}')
a2.set_xlabel('шаг центральной разности $h$'); a2.set_ylabel('|генератор − chromaticity_orders|'); a2.grid(alpha=.3, which='both'); a2.legend(fontsize=7)
plt.tight_layout(); plt.show()""")
md(r"""**Рисунок 4.** Слева: $Q(\delta)$ из генераторных $\xi_1$, $\xi_2$ для голого кольца и для одного SF1 с $k_2=1$ м$^{-3}$ (линейная матрица
одна и та же). Справа: расхождение генераторного значения с центральной разностью `chromaticity_orders` в зависимости от шага.

**Выводы.** Генератор даёт $\xi_1=(-1.4615,\,-2.2255)$ и $\xi_2=(8.555,\,8.104)$, совпадая с `chromaticity_series` до всех знаков, а с
`chromaticity_orders` — с расхождением, которое для $\xi_1$ падает как $h^2$, а для $\xi_2$ минимально при $h=10^{-5}$: при меньшем шаге
вторую разность портит округление $\sim\varepsilon/h^2$. Генератор свободен от обеих ошибок; шаг по умолчанию $10^{-4}$ смещает $\xi_2$ на $\sim10^{-4}$. Отклик одного секступоля $\Delta\xi=(+2.1642,\,-0.4531)$ совпадает с $k_2\!\int\!\beta D\,ds/(4\pi)$ до $10^{-8}$.
Знак скобки $-\tfrac12[f_3,F_3]$ подтверждён: с плюсом $\xi_2$ одного секступоля меняет знак. Этот способ чтения $\xi$ используется на шагах 3–7.""")

# ---------------------------------------------------------------- step 3
md(r"""## 3. Инварианты Ли: $\eta_1$, нормальная форма с полостью, действия против трека
Проверяется формула заметки для $\eta_1$ (тождество скорости $-c_2+c_1^2=(3\gamma^2-1)/(2\gamma^4)$ объясняет, почему она совпадает
с `path_expansion`). Затем к решётке добавляется одна тонкая полость (линейный режим) и строится 6D нормальная форма
`dragt_finn_nf6d` фактора `ring_lie_factor`:
$$h=-\sum_k\mu_kJ_k+\sum_{i\le j}c_{ij}J_iJ_j,\qquad \delta_{eq}=z^*_6+d\cdot J.$$
Её действия $J_x,J_y,J_s$ должны сохраняться на треке, а действие Куранта–Снайдера — нет: оно не учитывает синхротронную модуляцию
дисперсионной орбиты $D\delta(t)$. Полость уже входит в `ring_lie_factor`, поэтому `one_turn_factor_rf` на такой решётке считает её
дважды: $Q_s\to\sqrt2\,Q_s$ (BELL issue #66).""")
code("run_step('03_lie_invariants.jl', '03_lie_invariants.md')")
code(r"""n3 = kv('03_numbers.csv'); display(pd.Series(n3, name='значение').to_frame())
act = pd.read_csv(OUT / '03_actions.csv'); rel = lambda c: act[c] / act[c].mean() - 1
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
ax[0].plot(act.turn, rel('Jx_cs'), lw=.6, label='Курант–Снайдер'); ax[0].plot(act.turn, rel('Jx_cs_dispersion_free'), lw=.6, label=r'К–С без $D\delta$')
ax[0].plot(act.turn, rel('Jx_nf'), 'k', lw=1, label='нормальная форма'); ax[0].set_title(r'$J_x/\langle J_x\rangle-1$')
ax[1].plot(act.turn, rel('Jy_cs'), lw=.6, label='Курант–Снайдер'); ax[1].plot(act.turn, rel('Jy_nf'), 'k', lw=1, label='нормальная форма')
ax[1].set_title(r'$J_y/\langle J_y\rangle-1$')
ax[2].plot(act.l * 1e3, act.delta * 1e4, '.', ms=1.5); ax[2].set_xlabel('l [мм]'); ax[2].set_ylabel(r'$\delta\cdot10^4$')
ax[2].set_title(rf"продольная плоскость, $J_s$: разброс {act.Js_nf.max() / act.Js_nf.min() - 1:.1e}")
for a in ax[:2]: a.set_xlabel('оборот'); a.legend(fontsize=8)
for a in ax: a.grid(alpha=.3)
plt.tight_layout(); plt.show()""")
md(r"""**Рисунок 5.** Трек `exact_track_element` + `rf_kick`, 1874 оборота (два синхротронных периода). Слева и в центре — относительное
изменение действий: нормальная форма (чёрный) против Куранта–Снайдера. Справа — продольный фазовый портрет.

**Выводы.** $\eta_1-\eta_1^{формула}=0$ точно, $Q_s$ нормальной формы совпадает с формулой до $2\cdot10^{-9}$. Действия нормальной формы
сохраняются: $J_x$ в пределах 1.2 %, $J_y$ — $5\cdot10^{-4}$, $J_s$ — $4\cdot10^{-7}$. Горизонтальное действие Куранта–Снайдера меняется
в 2.3 раза (вычитание $D\delta$ уменьшает это до 23 %): его искажает синхротронная модуляция дисперсионной орбиты, которую нормальная
форма с полостью учитывает. По $y$ дисперсии нет, и оба действия сохраняются на уровне $10^{-4}$; у нормальной формы $J_y$ модулирован
с синхротронным периодом (амплитуда $5\cdot10^{-4}$), у Куранта–Снайдера — биениями $2.7\cdot10^{-4}$. Это различие здесь не разбиралось. Двойной учёт полости подтверждён: $Q_s$ после второго удара в $1.41421\approx\sqrt2$ раз больше.""")

# ---------------------------------------------------------------- step 4
md(r"""## 4. Равновесный импульс без секступолей (главы 1–3)
Три случая в одной решётке ($k_2=0$), четыре независимых чтения в каждом:
**формула** заметки, **тор** нормальной формы ($z^*_6+d\cdot J$), итерация **карты степени 3** фактора `ring_lie_factor` и **точный трек**.
Из трека $\delta_{eq}$ — постоянный член гармонического фита (гармоники 1–3, тюн уточняется).

- **A**: линейный режим, $J=0$ — глава 1, $\delta_{eq}=-\eta_1\delta_m^2/(2\eta_0)$.
- **B**: линейный режим, конечные $J_x,J_y$ — глава 2, добавляется $\Lambda$ от генераторной $\xi$.
- **C**: кубический режим $\hat\phi=0.3$ — глава 3: $\langle\delta^2\rangle=\tfrac{\delta_m^2}{2}(1-\hat\phi^2/32)$, $\Omega/\omega_s=1-\hat\phi^2/16$.

Тор нормальной формы — квадратичный полином по $J$, поэтому поправку $\hat\phi^2$ он не содержит и сравнивается с главами 1–2.""")
code("run_step('04_equilibrium.jl', '04_equilibrium.md')")
code(r"""n4 = pd.read_csv(OUT / '04_numbers.csv')
cols = ['formula', 'nf', 'map', 'exact']
tab = n4[n4.column.isin(cols) & n4.quantity.isin(['delta_eq', 'unc_delta_eq', 'delta2', 'Q', 'gamma_eff', 'delta_E_over_E'])]
display(tab.pivot_table(index=['case', 'quantity'], columns='column', values='value')[cols])
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for a, (case, g) in zip(ax, n4.groupby('case')):
    v = g.set_index(['column', 'quantity']).value; ref = v['target', 'delta_eq_ch12']
    y = [(v[c, 'delta_eq'] / ref - 1) * 1e3 for c in cols]; e = [v[c, 'unc_delta_eq'] / ref * 1e3 for c in cols]
    a.errorbar(range(4), y, e, fmt='o', capsize=4, color='k')
    for t, c, lab in (('delta_eq_ch12', 'C0', 'главы 1–2'), ('delta_eq_ch3', 'C3', r'глава 3 ($1-\hat\phi^2/32$)'), ('delta_eq_pendulum', 'C2', 'маятник')):
        a.axhline((v['target', t] / ref - 1) * 1e3, color=c, ls='--', lw=1, label=lab)
    a.set_xticks(range(4), ['формула', 'тор НФ', 'карта ст. 3', 'трек']); a.set_title(f'{case}: δ_eq(гл. 1–2) = {ref:.4e}')
    a.set_ylabel(r'$(\delta_{eq}/\delta_{eq}^{1-2}-1)\cdot10^3$'); a.grid(alpha=.3); a.legend(fontsize=7)
plt.tight_layout(); plt.show()""")
md(r"""**Рисунок 6.** $\delta_{eq}$ в четырёх чтениях относительно значения глав 1–2 (в промилле), с погрешностью фита. Линии — значения формул.

**Выводы.**
- **A и B.** Все четыре чтения совпадают в пределах $1.5\cdot10^{-4}$ относительно: $\delta_{eq}=9.53\cdot10^{-9}$ при $J=0$ и
  $1.908\cdot10^{-8}$ при конечных $J$. Вклад $\Lambda=5.1\cdot10^{-9}$ (хроматичность × бетатронное действие) удваивает сдвиг —
  при этих действиях он сравним с синхротронным членом $\eta_1\langle\delta^2\rangle$.
- **C.** Трек ($9.5013\cdot10^{-9}$) подтверждает главу 3 ($9.5027\cdot10^{-9}$) в пределах погрешности фита и отличается от глав 1–2
  на $-3\cdot10^{-3}$, т. е. на $\hat\phi^2/32$. Тюн $\Omega/\omega_s$ из карты и трека тоже следует $1-\hat\phi^2/16$.
- **STOP скрипта** — карта степени 3 даёт $9.545\cdot10^{-9}$, ближе к главам 1–2, чем к главе 3. Это усечение: поправка $\hat\phi^2$ к
  среднему приходит из члена $\phi^4$ гамильтониана полости, которого нет в карте степени 3. Тюн и третья гармоника $\hat\phi^2$
  несут, среднее — нет. Принято как свойство усечения, числа не правились.
- $\Delta E/E\approx2.7\cdot10^{-9}$ (A) и $4.9\cdot10^{-9}$ (B); $\gamma_{eff}-\gamma\sim10^{-8}$.""")

# ---------------------------------------------------------------- step 5
md(r"""## 5. Семейства секступолей: $\xi_x=\xi_y=\eta_1=0$
Три семейства по картине `place_sextupole_families`: SXD — в максимуме $D_x$ (2 магнита в SF1), SXB — в максимумах $\beta_x$ (4 в SF2),
SYB — в максимумах $\beta_y$ (8: SD и короткие дрейфы OS1 арок). В приближении тонкого магнита
$$\Delta\xi_x=\frac{1}{4\pi}\sum k_2L\beta_xD_x,\qquad \Delta\xi_y=-\frac{1}{4\pi}\sum k_2L\beta_yD_x,\qquad
\Delta\alpha_1=\Delta\eta_1=-\frac{1}{2C}\sum k_2LD_x^3 .$$
$\alpha_0$ секступоли не меняют, поэтому $\eta_1$ обнуляется через $\alpha_1=\alpha_0/\gamma^2-(3\gamma^2-1)/(2\gamma^4)=-0.678$.
Свой отклик (генераторная $\xi$ + `path_expansion`) решается на $\xi=\eta_1=0$. Сосед — `solve_A`: он обнуляет $\alpha_1$, а не $\eta_1$.""")
code("run_step('05_sextupoles.jl', '05_sextupoles.md')")
code(r"""n5 = pd.read_csv(OUT / '05_numbers.csv')
resp = n5[n5.setting == 'response'].pivot(index='quantity', columns='source', values='value').loc[['dxi_x', 'dxi_y', 'deta1'], ['SXD', 'SXB', 'SYB']]
L = 0.12; thin = pd.DataFrame({f: [(L * at('betx_m', d.sc) * at('Dx_m', d.sc)).sum() / (4 * np.pi),
                                   -(L * at('bety_m', d.sc) * at('Dx_m', d.sc)).sum() / (4 * np.pi),
                                   -(L * at('Dx_m', d.sc)**3).sum() / (2 * CIRC)] for f, d in FAM.items()}, index=resp.index)
display(pd.concat({'отклик BELL/генератор': resp, 'тонкая формула': thin, 'отношение': resp / thin}, axis=1))
sr = n5[n5.source == 'study_reading'].pivot(index='quantity', columns='setting', values='value')[['k2=0', 'solve_A', 'own']]
display(sr.loc[['k2_SXD', 'k2_SXB', 'k2_SYB', 'xi_x', 'xi_y', 'alpha1', 'eta1', 'Lambda', 'formula_delta_eq', 'exact_delta_eq', 'exact_unc_delta_eq']])
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4)); w = .38; x = np.arange(3)
a1.bar(x - w / 2, sr.loc[['k2_SXD', 'k2_SXB', 'k2_SYB'], 'own'], w, label=r'свой: $\xi=\eta_1=0$')
a1.bar(x + w / 2, sr.loc[['k2_SXD', 'k2_SXB', 'k2_SYB'], 'solve_A'], w, label=r'solve_A: $\xi=\alpha_1=0$')
a1.set_xticks(x, ['SXD', 'SXB', 'SYB']); a1.set_ylabel(r'$k_2$ [м$^{-3}$]'); a1.axhline(0, color='k', lw=.6); a1.grid(alpha=.3); a1.legend(fontsize=8)
S = ['k2=0', 'solve_A', 'own']
a2.bar(x - w / 2, np.abs(sr.loc['formula_delta_eq', S]).clip(1e-16), w, label='формула')
a2.bar(x + w / 2, np.abs(sr.loc['exact_delta_eq', S]), w, yerr=sr.loc['exact_unc_delta_eq', S], capsize=4, label='точный трек')
a2.set_yscale('log'); a2.set_xticks(x, [r'$k_2=0$', 'solve_A', 'свой']); a2.set_ylabel(r'$|\delta_{eq}|$'); a2.grid(alpha=.3, which='both'); a2.legend(fontsize=8)
plt.tight_layout(); plt.show()""")
md(r"""**Рисунок 7.** Слева: силы семейств для своего решения и для `solve_A`. Справа: $|\delta_{eq}|$ в линейном режиме при конечных $J$
для трёх настроек (формула — нижняя граница $10^{-16}$ для визуализации нуля).

**Выводы.**
- Отклик BELL совпадает с тонкой формулой (отношение ≈ 1 с точностью до толщины 0.12 м и до интерполяции $\beta$, $D$ в центре слота);
  в частности $\Delta\eta_1=-\sum k_2LD^3/(2C)$ — поэтому SXD в максимуме дисперсии главный рычаг по $\eta_1$. Обусловленность 8.8.
- Своё решение: $k_2=(0.818,\,-0.425,\,-0.611)$ м$^{-3}$, $\xi_x=\xi_y=\eta_1=0$ до $10^{-15}$, $\alpha_1=-0.678$. Формула даёт
  $\delta_{eq}=0$; трек — $-7.1\cdot10^{-13}\pm1.8\cdot10^{-14}$, т. е. в $2.7\cdot10^4$ раз меньше голого кольца ($1.91\cdot10^{-8}$).
  Остаток статистически значим и главами 1–3 не описывается; это порог следующего порядка, а не расхождение двух чтений одной величины.
- `solve_A` убирает $\alpha_1$, но оставляет кинематический $\eta_1=0.678$ и $\delta_{eq}=6.35\cdot10^{-9}$ (формула и трек совпадают).
  Для задачи равновесной энергии целевой должна быть $\eta_1$, а не $\alpha_1$ (BELL issue #67 — также про маску $D$ и общие слоты).""")

# ---------------------------------------------------------------- step 6
md(r"""## 6. Краевые поля диполей
Силы шага 5 не меняются; к диполям применяются три модели края: `NoEdge`, `HardEdge()` (порядок 2, $e_1=e_2=0$, зазор 0) и
`SoftEdge(LinearRamp(0.05))` (ширина 0.05 м — допущение шага, $\int F(1-F)\,dz=L_f/6=8.33\cdot10^{-3}$ м).
Жёсткий край с нулевым углом не меняет линейную карту, но его члены второго порядка дают хроматическую вертикальную фокусировку:
$\xi_y\neq0$ при прежних тюнах. Мягкий край добавляет линейную вертикальную фокусировку, по порядку величины
$|\Delta Q_y|\approx\frac1{4\pi}\sum_{\text{края}}\beta_y h^2\!\int\!F(1-F)\,dz$, $h=\theta/L$. Через $\Lambda$ любое $\xi\ne0$ снова сдвигает $\delta_{eq}$.""")
code("run_step('06_fringe.jl', '06_fringe.md')")
code(r"""n6 = pd.read_csv(OUT / '06_numbers.csv'); n6 = n6[n6.model != 'setup'].pivot(index='quantity', columns='model', values='value')
M = ['NoEdge', 'HardEdge', 'SoftEdge']
display(n6.loc[['Qx', 'Qy', 'dQy', 'xi_x', 'xi_y', 'xi_orders_y', 'eta1', 'formula_delta_eq', 'track_delta_eq', 'track_unc_delta_eq'], M])
bm = seq[seq.name == 'BM']; h = (bm.angle_closed_rad / bm.length_m).to_numpy()
est = (at('bety_m', np.r_[bm.s0, bm.s1]) * np.r_[h, h]**2).sum() * (0.05 / 6) / (4 * np.pi)
print(f"оценка |ΔQy| = Σ β_y h² ∫F(1−F)dz / 4π = {est:.4e};  BELL SoftEdge: ΔQy = {n6.loc['dQy', 'SoftEdge']:.4e}")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4)); w = .38; x = np.arange(3)
a1.bar(x - w / 2, n6.loc['xi_y', M], w, label=r'$\xi_y$ генератор'); a1.bar(x + w / 2, n6.loc['xi_orders_y', M], w, label=r'$\xi_y$ chromaticity_orders')
a1.set_xticks(x, M); a1.set_ylabel(r'$\xi_y$'); a1.grid(alpha=.3); a1.legend(fontsize=8)
a2.bar(x - w / 2, n6.loc['formula_delta_eq', M] * 1e9, w, label='формула')
a2.bar(x + w / 2, n6.loc['track_delta_eq', M] * 1e9, w, yerr=n6.loc['track_unc_delta_eq', M] * 1e9, capsize=4, label='трек')
a2.set_xticks(x, M); a2.set_ylabel(r'$\delta_{eq}\cdot10^9$'); a2.axhline(0, color='k', lw=.6); a2.grid(alpha=.3); a2.legend(fontsize=8)
plt.tight_layout(); plt.show()""")
md(r"""**Рисунок 8.** Слева: вертикальная хроматичность при трёх моделях края (генератор против `chromaticity_orders`). Справа: $\delta_{eq}$
(формула против трека) на силах своего решения шага 5.

**Выводы.**
- Край ломает коррекцию шага 5 только по $\xi_y$: жёсткий край даёт $\xi_y=0.408$ при неизменных тюнах и $\eta_1$, мягкий — $0.411$.
  Генератор и `chromaticity_orders` совпадают до $10^{-6}$.
- $\delta_{eq}$ становится $-1.06\cdot10^{-9}$ — 5.5 % голого значения; формула и трек совпадают до $10^{-3}$ относительно. Значит, для
  $\delta_{eq}\approx0$ секступоли надо решать на факторе **с краем** (достаточно `HardEdge`).
- Мягкий край: $\Delta Q_y=-1.27\cdot10^{-3}$ при $\Delta Q_x\approx-10^{-6}$; простая оценка через $\beta_yh^2\!\int\!F(1-F)$ даёт ту же
  величину. Разница мягкого и жёсткого края в $\xi_y$ (0.7 %) и $\eta_1=-1.0\cdot10^{-4}$ на $\delta_{eq}$ почти не влияет.
- `exact_track_element` не принимает `EdgeMap`; в треке край — фактор Ли порядка 3 (BELL issue #68).""")

# ---------------------------------------------------------------- step 7
md(r"""## 7. Нелинейность и динамическая апертура (перекрёстная проверка)
Своё решение шага 5, `NoEdge`. Аналитическая граница — критерий Гиньяра–Хагеля: дробный тюн приближается $p/q$ ($q\le100$, здесь
$\{Q_x\}\approx6/13$), вдоль линейной орбиты $x_n^{(0)}=R^n(x_0,0,0,0,0,\delta)$ перемножаются якобианы карты, и устойчивость требует
$$\Bigl|\operatorname{Tr}\prod_{n=1}^{q}J\bigl(x_n^{(0)}\bigr)\Bigr|<2 .$$
При $x_0=0$ след равен $2\cos(q\mu_x)$. Карта берётся усечённой до степени 2 (только секступоли) и до степени 3.
Трек `exact_track_turn` (300 оборотов, апертура 5 см) — перекрёстная проверка; оба результата показаны рядом, без усреднения.""")
code("run_step('07_dynamic_aperture.jl', '07_dynamic_aperture.md')")
code(r"""n7 = kv('07_numbers.csv')
display(pd.Series({k: n7[k] for k in ['Qx_frac', 'Qy_frac', 'p_over_q_p', 'p_over_q_q', 'p_over_q_err', 'linear_Tr', 'max_abs_f3', 'max_abs_f4',
                                     'n_launches', 'n_survived']}, name='значение').to_frame())
print(f"запас следа при x0 = 0: 2 − 2cos(qμx) = {2 - n7['linear_Tr']:.2e};  |Qx − Qy| = {abs(n7['Qx_frac'] - n7['Qy_frac']):.2e}")
gu = pd.read_csv(OUT / '07_guignard.csv'); gu = gu[gu.delta == 0]
sv = pd.read_csv(OUT / '07_survival.csv'); s0 = sv[(sv.y0 == 0) & (sv.delta == 0)]
fig, ax = plt.subplots(1, 2, figsize=(14, 4.4))
for a, lim, scale in ((ax[0], 0.05, 1e4), (ax[1], 1.0, None)):
    for deg, c in ((2, 'C0'), (3, 'C3')):
        g = gu[(gu.degree == deg) & (gu.kind == 'sample') & (gu.x0.abs() <= lim)].sort_values('x0')
        a.plot(g.x0 * 1e3, (g.tr_horizontal - 2) * scale if scale else g.tr_horizontal, '.-', ms=3, color=c, label=f'след, карта степени {deg}')
        for xb in gu[(gu.degree == deg) & (gu.kind == 'boundary')].x0:
            if abs(xb) <= lim: a.axvline(xb * 1e3, color=c, ls=':', lw=1)
    a.set_xlabel(r'$x_0$ [мм]'); a.grid(alpha=.3)
ax[0].axhline(0, color='g', lw=1, label=r'$\mathrm{Tr}=2$'); ax[0].set_ylabel(r'$(\mathrm{Tr}\prod J-2)\cdot10^4$')
ax[0].plot(s0.x0 * 1e3, np.full(len(s0), .04), 'k^', ms=7, transform=ax[0].get_xaxis_transform(), label='трек: выжил 300 оборотов')
ax[1].axhspan(-2, 2, color='g', alpha=.07); ax[1].set_yscale('symlog', linthresh=2); ax[1].set_ylabel(r'$\mathrm{Tr}\prod J$')
ax[0].set_title('окно апертуры ±50 мм'); ax[1].set_title('всё окно поиска ±1 м'); ax[0].legend(fontsize=8)
plt.tight_layout(); plt.show()""")
md(r"""**Рисунок 9.** След произведения якобианов за $q=13$ оборотов вдоль линейной орбиты при $\delta=0$. Слева — окно апертуры,
отклонение следа от 2 (выше нуля — неустойчиво); треугольники — запуски точного трека, все выжили. Справа — всё окно поиска,
зелёная полоса $|\mathrm{Tr}|<2$. Пунктир — первый выход из полосы.""")
code(r"""ph = pd.read_csv(OUT / '07_phase.csv')
fig, ax = plt.subplots(1, 4, figsize=(16, 3.8))
for a, (x0, g) in zip(ax, ph.groupby('x0')):
    a.plot(g.x * 1e3, g.px * 1e3, '.', ms=2); a.set_title(f'$x_0$ = {x0 * 1e3:.0f} мм'); a.set_xlabel('x [мм]'); a.grid(alpha=.3)
ax[0].set_ylabel(r'$p_x$ [мрад]'); plt.tight_layout(); plt.show()""")
md(r"""**Рисунок 10.** Фазовые портреты точного трека в начале кольца, $y=0$, $\delta=0$, 300 оборотов.

**Выводы.**
- Резонансы низких порядков далеко: до $1/3$ — 0.128, до $1/4$ — 0.038 (по ближайшему кратному). Портреты — эллипсы без островов
  вплоть до 40 мм; все 83 запуска выжили, потерь внутри апертуры 5 см нет. Штриховой вид портретов при 5–20 мм — это близость к $6/13$:
  за 300 оборотов точки ложатся 13 медленно дрейфующими группами. При 40 мм эллипс заполнен — расстройка от амплитуды уводит тюн от $6/13$.
- Граница Гиньяра по карте степени 2 — $\pm8.8$ мм, по степени 3 — $\pm0.46$ м. **Это расхождение не сглажено.** Рисунок 9 показывает
  механизм: $\{Q_x\}$ лежит в $1.8\cdot10^{-4}$ от $6/13$, поэтому при $x_0=0$ след уже $1.99978$ — запас до границы $2.2\cdot10^{-4}$.
  След степени 2 растёт с амплитудой квадратично и пересекает 2 уже при 9 мм; след степени 3 уходит от 2 вниз, т. е. расстройка
  уводит тюн от $6/13$ — так же, как на портрете 40 мм. В окне ±50 мм с треком согласуется степень 3, а не степень 2; степень 2 ещё и
  возвращается в полосу при $|x_0|>0.5$ м (re-entrant). Вывод: на этом тюне критерий плохо обусловлен и чувствителен к усечению карты,
  а 300 оборотов не разрешают резонанс 13-го порядка. Это гипотеза о механизме, а не выбор победителя.
- Отдельно: $|Q_x-Q_y|=3.8\cdot10^{-4}$ — кольцо сидит на разностном резонансе связи. Без косых элементов связи нет, но любая косая
  ошибка смешает плоскости полностью.""")

# ---------------------------------------------------------------- summary
md(r"""## Итоги
| шаг | результат |
| --- | --- |
| 1 | BM $=\pi/8$, поворот $2\pi$ до $2\cdot10^{-15}$; недогиб исходного файла эквивалентен $\delta\sim1.6\cdot10^{-8}$ |
| 2 | $Q=(2.4617,\,2.4613)$, $\eta_0=-0.534$, $\eta_1=1.017$, $\gamma_t=2.08$; режимы $\hat\phi=0.05$ и $0.3$ |
| 2b | $\xi_1=(-1.46,\,-2.23)$, $\xi_2=(8.55,\,8.10)$ из генераторов, совпадение с series |
| 3 | $\eta_1$ = формула; действия НФ сохраняются, К–С по $x$ нет; двойной учёт полости $\sqrt2$ |
| 4 | $\delta_{eq}=9.53\cdot10^{-9}$ ($J=0$), $1.91\cdot10^{-8}$ (конечные $J$); трек подтверждает $1-\hat\phi^2/32$ |
| 5 | $\xi=\eta_1=0$ тремя семействами; $\delta_{eq}$: $1.9\cdot10^{-8}\to-7\cdot10^{-13}$ (трек); `solve_A` оставляет $6.4\cdot10^{-9}$ |
| 6 | край диполя: $\xi_y=0.41$, $\delta_{eq}=-1.06\cdot10^{-9}$ |
| 7 | трек: 83/83 до 4.5 см; Гиньяр: 8.8 мм (ст. 2) / 0.46 м (ст. 3) |

**Открытые вопросы** (не сглажены, числа не правились):
1. Шаг 4: карта степени 3 не несёт $\hat\phi^2$ в среднем $\delta_{eq}$ — усечение; для кубического режима нужна карта степени ≥ 4.
2. Шаг 5: остаток трека $-7.1\cdot10^{-13}$ после коррекции значим и не описывается главами 1–3.
3. Шаг 7: граница Гиньяра степени 2 (8.8 мм) против трека (без потерь до 45 мм); степень 3 с треком согласуется, критерий на тюне
   в $1.8\cdot10^{-4}$ от $6/13$ плохо обусловлен.

**Следующие шаги.** Решать секступоли на факторе с `HardEdge` (цель $\xi=\eta_1=0$ с краем); длинный трек DA (`LONG_DA=1`) и
чувствительность к выбору $q$; разнести $Q_x$ и $Q_y$ хотя бы на $10^{-2}$, прежде чем вводить ошибки и спин.""")

nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
nbf.write(nb, 'magnetic_ring.ipynb')
