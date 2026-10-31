# Step 01: one particle, no betatron, RF on. Raw δ_n over a few synchrotron periods.
function step_01!(S, ctx)
    section("01. Ряд δ при включённом ВЧ, без бетатрона")
    dm = 1e-3
    N = 400
    z0 = init_canonical(S, dm, 0.0, 0.0)
    Z = run_turns(z -> step_exact(S, z), z0, N)
    d = deltas_of(Z)
    ctx.series[:Z4_short] = Z
    ctx.series[:d4_short] = d
    ctx.series[:z0_dm] = z0
    ctx.series[:dm] = dm
    m1 = formula_delta_eq(S.pe, S.xi, S.C, dm, 0.0, 0.0)
    p = plot(xlabel = "номер оборота n", ylabel = "δ_n · 10³", legend = :topright, size = (760, 380))
    plot!(p, 1:N, d .* 1e3, lw = 1.3, marker = :circle, ms = 2, label = "точный интегратор, δ_m = 10⁻³")
    hline!(p, [0.0], color = :black, lw = 0.8, label = "δ = 0")
    hline!(p, [m1 * 1e3], ls = :dash, lw = 2, label = "формула −η₁δ_m²/(2η₀)")
    savefig_step(p, "fig_01_delta")
    half = (maximum(d) - minimum(d)) / 2
    mac!(ctx, "dm", "10^{-3}")
    mac!(ctx, "Nshort", string(N))
    mac!(ctx, "dHalf", texnum(half; d = 2))
    mac!(ctx, "MoneSimple", texnum(m1; d = 3))
    mac!(ctx, "Qs", @sprintf("%.6f", S.Qs))
    mac!(ctx, "Ts", @sprintf("%.1f", 1 / S.Qs))
    conc = @sprintf("Из-за включённого ВЧ δ колеблется с периодом ≈%.0f оборотов, размах %.2e. Центр %.2e глазом не виден: он в тысячу раз меньше размаха.",
                    1 / S.Qs, half, m1)
    write_conclusion("step01", conc)
    write_step_tex("step01", """
\\section{Ряд \$\\delta\$ при включённом ВЧ}
\\textbf{Вопрос.} Видно ли равновесное \$\\delta\$ на сыром ряде оборотов, если полость стоит в кольце и напряжение включено?

Частица стартует с \$\\delta_m=\\vdm\$, без бетатронной амплитуды. Тонкая полость --- элемент кольца, удар есть на каждом обороте. Без этого удара ряд был бы прямой: \$\\delta\$ не меняется. Здесь удар есть, поэтому есть синхротронное колебание с тюном \$Q_s=\\vQs\$ (период \$T_s\\approx\\vTs\$ оборота).

\\begin{figure}[H]\\centering
\\includegraphics[width=0.92\\textwidth]{figs/fig_01_delta.png}
\\caption{Точный интегратор, \$N=\\vNshort\$ оборотов. Размах \$\\vdHalf\$. Горизонталь формулы \$\\vMoneSimple\$ на этом масштабе совпадает с нулём.}
\\end{figure}

\\textbf{Вывод.} \$\\delta\$ колеблется из-за ВЧ. Центр колебания глазом не виден: формула даёт \$\\vMoneSimple\$, это в тысячу раз меньше размаха.
""")
    return ctx
end
