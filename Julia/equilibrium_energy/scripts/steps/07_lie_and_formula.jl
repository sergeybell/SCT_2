# Step 07: formula vs Lie normal form, in words, on the running-mean plot.
function step_07!(S, ctx)
    section("07. Формула рядом с генератором Ли")
    z0 = ctx.series[:z0_dm]
    d = ctx.series[:d4]
    J0 = nf6d_actions(S.nf, z0)
    m2 = nf6d_delta_eq(S.nf, J0)
    dm = ctx.series[:dm]
    m1 = formula_delta_eq(S.pe, S.xi, S.C, dm, 0.0, 0.0)
    # the same shift written two ways, for the synchrotron action only
    ds_formula = -S.pe.eta1 * S.nf.sigma_delta_s^2 / S.pe.eta0
    ds_nf = S.nf.d[3]
    dx_formula = 2pi / S.C * S.xi.x[2] / S.pe.eta0
    dx_nf = S.nf.d[1]
    dy_formula = 2pi / S.C * S.xi.y[2] / S.pe.eta0
    dy_nf = S.nf.d[2]
    rm = running_mean(d)
    p = plot(xlabel = "номер оборота N", ylabel = "⟨δ⟩_N · 10⁶", legend = :topright,
             size = (760, 440), xscale = :log10)
    plot!(p, 1:length(d), rm .* 1e6, lw = 1.6, label = "точный интегратор")
    hline!(p, [m1 * 1e6], lw = 2, ls = :dash, label = "формула из η, ξ")
    hline!(p, [m2 * 1e6], lw = 2, ls = :dot, label = "среднее нормальной формы")
    savefig_step(p, "fig_07_lie")
    mac!(ctx, "dxForm", texnum(dx_formula; d = 5))
    mac!(ctx, "dxNF", texnum(dx_nf; d = 5))
    mac!(ctx, "dxRel", texnum(dx_nf / dx_formula - 1; d = 1))
    mac!(ctx, "dyForm", texnum(dy_formula; d = 5))
    mac!(ctx, "dyNF", texnum(dy_nf; d = 5))
    mac!(ctx, "dsForm", texnum(ds_formula; d = 5))
    mac!(ctx, "dsNF", texnum(ds_nf; d = 5))
    mac!(ctx, "dsRel", texnum(ds_nf / ds_formula - 1; d = 1))
    mac!(ctx, "etazero", @sprintf("%.6f", S.pe.eta0))
    mac!(ctx, "etaone", @sprintf("%.6f", S.pe.eta1))
    mac!(ctx, "xix", @sprintf("%.4f", S.xi.x[2]))
    mac!(ctx, "xiy", @sprintf("%.4f", S.xi.y[2]))
    conc = @sprintf("Один и тот же сдвиг записан дважды. Формула: ∂δ_eq/∂J_x = 2πξ_x/(Cη₀) = %.5e. Нормальная форма, без этих имён: %.5e (отн. разница %.1e). По синхротронному действию: формула %.5e, NF %.5e.",
                    dx_formula, dx_nf, dx_nf / dx_formula - 1, ds_formula, ds_nf)
    write_conclusion("step07", conc)
    write_step_tex("step07", """
\\section{Формула и произведение Ли}
\\textbf{Вопрос.} Собирает ли формула те же числа, которые уже лежат в произведении Ли?

Формула говорит: равновесный \$\\delta\$ сдвигается на \$2\\pi\\xi_x/(C\\eta_0)\$ на единицу горизонтального действия, и на \$-\\eta_1\\sigma_\\delta^2/\\eta_0\$ на единицу синхротронного. Нормальная форма этих имён не использует: она усредняет лабораторный \$\\delta\$ по фазам генераторов.

\\begin{center}
\\begin{tabular}{lcc}
\\toprule
сдвиг \$\\delta_{\\mathrm{eq}}\$ на единицу действия & формула & нормальная форма \\\\
\\midrule
горизонтальное действие & \$\\vdxForm\$~м\$^{-1}\$ & \$\\vdxNF\$ \\\\
вертикальное действие & \$\\vdyForm\$~м\$^{-1}\$ & \$\\vdyNF\$ \\\\
синхротронное действие & \$\\vdsForm\$~м\$^{-1}\$ & \$\\vdsNF\$ \\\\
\\bottomrule
\\end{tabular}
\\end{center}

Относительная разница по горизонтали \$\\vdxRel\$, по синхротрону \$\\vdsRel\$. Это не таблица «насколько теория верна в процентах к трекингу»: это две записи одного сдвига.

\\begin{figure}[H]\\centering
\\includegraphics[width=0.92\\textwidth]{figs/fig_07_lie.png}
\\caption{Скользящее среднее интегратора и две горизонтали: формула из \$\\eta,\\xi\$ и среднее нормальной формы. Линия M2 посчитана из генераторов, без подстановки формулы.}
\\end{figure}

\\textbf{Вывод.} Формула собрала те же числа, которые алгебра уже содержит. Дальше среднее по оборотам садится на этот общий уровень.
""")
    return ctx
end
