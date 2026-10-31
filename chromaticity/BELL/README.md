# BELL (Julia) — хроматичность magnetic_2_md

Запуск: `run_bell.bat` (или `julia --project=C:\GIT_REPOS\BELL chrom_bell.jl`).
Ядро BELL не менялось, используется как есть (v0.10.0).

* Решётка читается OptiM-ридером BELL `convert_optim` → `magnetic_2_md.bell` → `load_structure`.
* Частица: `DEUTERON` из BELL (m = 1875.61294257 МэВ, CODATA 2018; та же масса задана и в COSY `PART`; заголовок OptiM 1875.612945 не используется), угол BM ровно π/8 (`convert_optim(...; close_ring = true)`: поле B из .opt даёт угол на 1.6e-8 меньше),
  W = 270 МэВ (из заголовка .opt).
* Что считается (`output/bell_results.txt`):
  * `linear_*` — одно-оборотная матрица (`ring_matrix`, `periodic_optics`, `optics_along`): Qx, Qy, alpha_c, eta.
  * `taylor_*` — одно-оборотное отображение Тейлора 3-го порядка (`ring_taylor_map`), `path_expansion`
    (alpha0, alpha1, eta0, eta1) и `chromaticity_orders`.
  * `lie_*` — то же из одно-оборотного Ли-фактора (`observables`: `ring_lie_factor` → `factor_to_map`).
* Определение хроматичности: δ = dp/p (канонический импульс BELL), ξ1 = dQ/dδ — центральная разность
  (шаг 1e-4) дробных тюнов якобиана отображения на замкнутой орбите при фиксированном δ;
  ξ2 = d²Q/dδ². Натуральная (секступоли 0). Поля магнитов фиксированы, диполь — секторный, без краёв и fringe.
* Нормальная форма `dragt_finn_nf6d` в BELL рассчитана на кольцо с ВЧ (6D) и тюны/хроматичность в
  данном скрипте не даёт; хроматичность берётся из Ли-отображения через `observables`.

* Twiss вдоль кольца: `output/bell_twiss.csv` (s, β_x, β_y, D_x на dp/p) из `optics_along`, шаг 2 см.
