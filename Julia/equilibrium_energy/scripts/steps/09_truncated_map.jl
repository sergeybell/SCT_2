# Step 09: truncated map M3, Lie-product iteration, exact M4. Mean δ and action.
function step_09!(S, ctx)
    section("09. Усечённая карта, произведение Ли, точный интегратор")
    dm = ctx.series[:dm]
    N = length(ctx.series[:d4])
    z0 = ctx.series[:z0_dm]
    d4 = ctx.series[:d4]
    println("  truncated map $N turns …")
    t3 = @elapsed Z3 = run_turns(z -> step_map(S, z), z0, N)
    d3 = deltas_of(Z3)
    @printf("  map %.1f s; Lie product …\n", t3)
    tL = @elapsed ZL = run_turns(z -> step_lie(S, z), z0, N)
    dL = deltas_of(ZL)
    @printf("  Lie %.1f s\n", tL)
    m1 = formula_delta_eq(S.pe, S.xi, S.C, dm, 0.0, 0.0)
    m2 = nf6d_delta_eq(S.nf, nf6d_actions(S.nf, z0))
    rm4 = running_mean(d4); rm3 = running_mean(d3); rmL = running_mean(dL)
    p = plot(xlabel = "номер оборота N", ylabel = "⟨δ⟩_N · 10⁶", legend = :topright,
             size = (760, 440), xscale = :log10)
    plot!(p, 1:N, rm4 .* 1e6, lw = 1.6, label = "точный интегратор (M4)")
    plot!(p, 1:N, rmL .* 1e6, lw = 1.6, ls = :dash, label = "произведение Ли")
    plot!(p, 1:N, rm3 .* 1e6, lw = 1.6, ls = :dot, label = "усечённая карта степени 3 (M3)")
    hline!(p, [m1 * 1e6], lw = 1.5, color = :black, label = "формула")
    hline!(p, [m2 * 1e6], lw = 1.5, color = :gray, ls = :dash, label = "M2")
    savefig_step(p, "fig_09_mean")

    Ax = 10e-3
    zB = init_canonical(S, dm, Ax, 0.0)
    Nb = 4000
    println("  action drift $Nb turns …")
    Z4b = run_turns(z -> step_exact(S, z), zB, Nb)
    Z3b = run_turns(z -> step_map(S, z), zB, Nb)
    ZLb = run_turns(z -> step_lie(S, z), zB, Nb)
    J4 = [cs_J_turn(Z4b[:, n], S.po.x) for n in 1:Nb]
    J3 = [cs_J_turn(Z3b[:, n], S.po.x) for n in 1:Nb]
    JL = [cs_J_turn(ZLb[:, n], S.po.x) for n in 1:Nb]
    pJ = plot(xlabel = "номер оборота n", ylabel = "(J − J_{точн}) / J_{точн}", legend = :topright, size = (760, 420))
    plot!(pJ, 1:Nb, (JL .- J4) ./ J4, lw = 1.4, label = "произведение Ли")
    plot!(pJ, 1:Nb, (J3 .- J4) ./ J4, lw = 1.4, label = "усечённая карта")
    savefig_step(pJ, "fig_09_action")
    amps = [1e-3, 3e-3, 1e-2, 2e-2, 3e-2]
    err_map = [symplectic_error(jacobian_of_map(S.M, init_canonical(S, dm, A, 0.0))) for A in amps]
    pS = plot(xlabel = "A_x, мм", ylabel = "||Jᵀ S J − S||", legend = :topleft,
              size = (760, 400), yscale = :log10)
    plot!(pS, amps .* 1e3, max.(err_map, 1e-18), marker = :circle, lw = 1.6, label = "усечённая карта")
    savefig_step(pS, "fig_09_symplectic")
    mac!(ctx, "SympMap", texnum(err_map[end]; d = 2))
    mac!(ctx, "SympLie", "машинная")
    mac!(ctx, "SympSmall", texnum(err_map[1]; d = 2))
    mac!(ctx, "JxRelMap", texnum(mean((J3[end-199:end] .- J4[end-199:end]) ./ J4[end-199:end]); d = 2))
    mac!(ctx, "JxRelLie", texnum(mean((JL[end-199:end] .- J4[end-199:end]) ./ J4[end-199:end]); d = 2))
    mac!(ctx, "MeanMap", texnum(rm3[end]; d = 4))
    mac!(ctx, "MeanLie", texnum(rmL[end]; d = 4))
    mac!(ctx, "MeanExact", texnum(rm4[end]; d = 4))
    mac!(ctx, "DiffMapExact", texnum(rm3[end] - rm4[end]; d = 2))
    mac!(ctx, "DiffLieExact", texnum(rmL[end] - rm4[end]; d = 2))
    mac!(ctx, "DiffMapMtwo", texnum(rm3[end] - m2; d = 2))
    mac!(ctx, "JxExact", texnum(mean(J4); d = 3))
    mac!(ctx, "JxLie", texnum(mean(JL); d = 3))
    mac!(ctx, "JxMap", texnum(mean(J3); d = 3))
    mac!(ctx, "JxMapEnd", texnum(mean(J3[end-199:end]); d = 3))
    mac!(ctx, "JxExactEnd", texnum(mean(J4[end-199:end]); d = 3))
    mac!(ctx, "Nbact", string(Nb))
    conc = @sprintf("Среднее δ: интегратор %.3e, Ли %.3e (разница %.2e), карта %.3e (разница %.2e). Относительный уход Jx карты за %d оборотов %.2e, Ли %.2e. Ошибка симплектичности карты: %.2e при 1 мм, %.2e при 30 мм.",
                    rm4[end], rmL[end], rmL[end] - rm4[end], rm3[end], rm3[end] - rm4[end],
                    Nb, mean((J3[end-199:end] .- J4[end-199:end]) ./ J4[end-199:end]),
                    mean((JL[end-199:end] .- J4[end-199:end]) ./ J4[end-199:end]),
                    err_map[1], err_map[end])
    write_conclusion("step09", conc)
    write_step_tex("step09", """
\\section{Усечённая карта и произведение Ли}
\\textbf{Вопрос.} Что портит полином степени 3: среднее \$\\delta\$ или само движение?

M3 --- итерация полинома степени 3 того же произведения Ли (это \$\\mathtt{MAP\\_ORDER}\$ Bell). Отдельно итерируется само произведение экспонент (\$\\mathtt{track}\$ / \$\\mathtt{apply\\_factor}\$). M4 --- точный интегратор. Формула и M2 --- горизонтали, не эталон.

\\begin{figure}[H]\\centering
\\includegraphics[width=0.92\\textwidth]{figs/fig_09_mean.png}
\\caption{Скользящее среднее \$\\delta\$. В конце: интегратор \$\\vMeanExact\$, произведение Ли \$\\vMeanLie\$ (разница \$\\vDiffLieExact\$), карта \$\\vMeanMap\$ (разница к интегратору \$\\vDiffMapExact\$).}
\\end{figure}

\\begin{figure}[H]\\centering
\\includegraphics[width=0.92\\textwidth]{figs/fig_09_action.png}
\\caption{Относительный уход \$J_x\$ от точного интегратора, \$A_x=10\$~мм, \$N=\\vNbact\$. Карта: \$\\vJxRelMap\$, произведение Ли: \$\\vJxRelLie\$.}
\\end{figure}

\\begin{figure}[H]\\centering
\\includegraphics[width=0.92\\textwidth]{figs/fig_09_symplectic.png}
\\caption{Ошибка симплектичности однооборотного шага усечённой карты: \$\\vSympSmall\$ при 1~мм, \$\\vSympMap\$ при 30~мм. Произведение Ли симплектично по построению (поток генератора).}
\\end{figure}

\\textbf{Вывод.} Среднее \$\\delta\$ у трёх кривых рядом. На 10~мм за \$\\vNbact\$ оборотов уход действия карты ещё мал (\$\\vJxRelMap\$). Симплектичность карты портится с амплитудой; произведение Ли остаётся симплектичным.
""")
    return ctx
end
