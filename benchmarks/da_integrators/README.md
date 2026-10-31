# benchmarks/da_integrators — Bell.DiffAlgebra / Bell.Integrators против COSY INFINITY 10.2

Точность и скорость DA-алгебры Bell и интеграторов ОДУ с DA-начальными условиями в сравнении с COSY.
Ядро BELL не изменялось (только `Pkg.develop` в отдельном окружении `julia/`). Устроено так, чтобы позже стать тестом Bell:
задачи, допуски и эталонные значения целиком в `julia/*.jl`, COSY-результаты — в `cosy/output/`.

| файл | что делает |
|---|---|
| `cosy/da_ops_template.fox`, `cosy/run_cosy_da.ps1` | DA-задачи nv = 6, N = 3/5/7/10: коэффициенты (`acc_o<N>.txt`) и время операций (`time_o<N>.txt`) |
| `cosy/hh_template.fox` | Хенона–Хейлеса, RK4 на FOX (тот же алгоритм, что в Bell), карты и время |
| `cosy/elem_template.fox`, `cosy/run_cosy_ode.ps1` | одиночные MQ / MS (протон 200 МэВ, масса 938.27208816 МэВ как в Bell; COSY `RP 200 938.27208816/AMUMEV 1` вместо `RPP`), карты COSY и время |
| `julia/setup.jl` | однократная настройка окружения (Bell dev, BenchmarkTools) |
| `julia/cosy_io.jl` | чтение DA из вывода COSY `WRITE` |
| `julia/da_compare.jl` | (1) совпадение коэффициентов, (2) скорость DA-операций → `output/da_accuracy.csv`, `da_timing.csv`, `da_compile.csv` |
| `julia/ode_compare.jl` | (3) все интеграторы Bell + DA-поток vs COSY → `output/ode_workprecision.csv`, `ode_cosy.csv` |
| `julia/hh_rk4_same.jl` | одинаковый RK4 Bell vs COSY → `output/hh_rk4_same.csv` |
| `julia/profile.jl` | (4) профилирование DOP853 по квадруполю → `output/profile.csv`, `profile_flat_o7.txt` |
| `julia/da_fast.jl`, `da_fast_bench.jl` | прототипы: in-place умножение, `compose_fast`, `inverse_fast`, `daflow` → `output/da_fast.csv` |
| `julia/rk_fast.jl` | прототип шага DOP853 с in-place стадиями → `output/rk_fast.csv` |
| `make_notebook.py` → `da_integrators.ipynb` | блокнот с результатами (русский, выполнен), рисунки `fig_*.png` |

Порядок запуска (≈ 15 мин):
```
powershell -File cosy\run_cosy_da.ps1      # ~3.5 мин
powershell -File cosy\run_cosy_ode.ps1     # ~1.5 мин
julia julia\setup.jl                       # один раз
julia --project=julia julia\da_compare.jl  # ~2 мин
julia --project=julia julia\ode_compare.jl # ~6 мин
julia --project=julia julia\hh_rk4_same.jl
julia --project=julia julia\profile.jl
julia --project=julia julia\da_fast_bench.jl
julia --project=julia julia\rk_fast.jl
python make_notebook.py; jupyter nbconvert --to notebook --execute --inplace da_integrators.ipynb
```

Главное: коэффициенты Bell = COSY до округления (≤ 3e-15 отн.); умножение Bell в ~2 раза быстрее COSY, элементарные функции
в ~1.5 раза, композиция/обращение в 2–9 раз медленнее (прототипы ускоряют в 3–10 раз). Для s-независимых элементов COSY
использует DA-оператор потока, а не RK; аналогичный прототип в Bell по скорости на уровне COSY и в 5–55 раз быстрее DOP853/GBS.
Подробности и рекомендации — в блокноте, раздел 6.