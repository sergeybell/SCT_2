# benchmarks/lie_flow — `Bell.LieFlow` (оператор потока Ли) против COSY INFINITY 10.2

Карты s-независимых элементов (дрейф, квадруполь, секторный диполь, секступоль) из полного гамильтониана рядом Ли
$e^{-L:H:}z$ до машинной точности, порядки 3/5/7/10: совпадение с COSY, скорость (против COSY, «силовой» формы DA-потока,
DOP853/GBS из `../da_integrators`), симплектичность (Bell и экспортированная карта COSY), кольцо: Драгт–Финн против потока.

| файл | что делает |
|---|---|
| `cosy/elem6_template.fox`, `cosy/run_cosy.ps1` | `DL`, `MQ`, `MS`, `MH` (протон 200 МэВ, масса 938.27208816 МэВ как в Bell; COSY `RP 200 938.27208816/AMUMEV 1` вместо `RPP`), все 6 строк карты и время, N = 3, 5, 7, 10 |
| `julia/lie_flow_bench.jl` | всё сравнение → `output/flow_vs_cosy.csv`, `flow_speed.csv`, `flow_symp.csv`, `flow_symp_L.csv`, `ring_df_vs_flow.csv` |
| `make_notebook.py` → `lie_flow.ipynb` | блокнот (русский, выполнен), рисунки `fig_*.png` |

Окружение Julia — соседнее `../da_integrators/julia` (Bell через `Pkg.develop`, BenchmarkTools); там же используются `cosy_io.jl`, `da_fast.jl`.
```
powershell -File cosy\run_cosy.ps1                                   # ~1 мин
julia --project=..\da_integrators\julia julia\lie_flow_bench.jl      # ~1.5 мин
python make_notebook.py; jupyter nbconvert --to notebook --execute --inplace lie_flow.ipynb
```
Итог: коэффициенты = COSY до ≤ 1.5e-15 отн. (диполь N=10: 4e-14); J^T S J − S на уровне округления у обоих;
время Bell/COSY: диполь 0.6–1.2, секступоль 0.9–1.5, квадруполь 1.5–2.1, дрейф 3–15; поток в 5–30 раз быстрее DOP853/GBS;
кольцо порядка 3: поток = Драгт–Финн до 1e-14 отн.