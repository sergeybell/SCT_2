# EquilibriumEnergy — δ_eq по шагам

Зависит от Bell в `C:\GIT_REPOS\BELL`. Физика — в `note/main.tex`.
OptiM-файл `C:\GIT_REPOS\SCT_2\OptiM\magnetic\magnetic_2.opt` при загрузке
переписывается в `lattices/magnetic_2.bell`.

Тонкая полость — элемент кольца (`RFCavity` в `lat.elements`), напряжение
включено с первого шага. Фаза в отчёте: `φ_note > 0` значит частица опаздывает.

```
julia --project=. setup.jl
julia --project=. scripts/run_v08.jl     # девять шагов, рисунки и gen/*.tex
cd note; latexmk -xelatex main.tex
```

Тетрадь одна: `notebooks/delta_eq.ipynb`. Ядро **Julia EquilibriumEnergy**.
Короткие прогоны считает сама; полные рисунки на 12 000 оборотов лежат в `note/figs`.

Четыре модели, полость в списке элементов:

- M1 — формула из времени пролёта
- M2 — среднее по тору нормальной формы произведения Ли (`build_factors`, `dragt_finn_nf6d`)
- M3 — итерация полинома степени 3 того же произведения (`factor_to_map`, `evaluate_map`)
- M4 — точный интегратор (`exact_track_element` + `rf_kick` на полости)

Числа в заметке пишет `scripts/run_v08.jl` в `note/gen/`.
