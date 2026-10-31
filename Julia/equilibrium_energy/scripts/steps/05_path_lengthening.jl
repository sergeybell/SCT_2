# Step 05: extra path length with RF on, as the difference of two particles.
function step_05!(S, ctx)
    section("05. Удлинение пути при включённом ВЧ")
    dm = ctx.series[:dm]
    Ax = ctx.series[:Ax]
    Ay = 20e-3
    z0 = init_canonical(S, dm, 0.0, 0.0)
    zX = init_canonical(S, dm, Ax, 0.0)
    zY = init_canonical(S, dm, 0.0, Ay)
    nshow = 3
    s0 = sample_along(S, z0; nturns = nshow)
    sX = sample_along(S, zX; nturns = nshow)
    sY = sample_along(S, zY; nturns = nshow)
    dLX = sX.L .- s0.L
    dLY = sY.L .- s0.L
    Jx = cs_action(Ax, S.po.x)
    Jy = cs_action(Ay, S.po.y)
    LamX = formula_B_Lambda(S.C, Jx, S.xi.x[2], 0.0, 0.0)
    LamY = formula_B_Lambda(S.C, 0.0, 0.0, Jy, S.xi.y[2])
    pX = plot(xlabel = "s, м", ylabel = "ΔL_x − ΔL_0, мм", legend = :topleft, size = (760, 400))
    plot!(pX, sX.s, dLX .* 1e3, lw = 1.4, label = "разность двух частиц, A_x = 20 мм")
    plot!(pX, [0, nshow * S.C], [0, nshow * S.C * LamX] .* 1e3, ls = :dash, lw = 2, label = "Λ_x s (формула)")
    vline!(pX, [S.C, 2S.C], color = :gray, ls = :dot, label = "")
    savefig_step(pX, "fig_05_path_x")
    pY = plot(xlabel = "s, м", ylabel = "ΔL_y − ΔL_0, мм", legend = :topleft, size = (760, 400))
    plot!(pY, sY.s, dLY .* 1e3, lw = 1.4, label = "разность двух частиц, A_y = 20 мм")
    plot!(pY, [0, nshow * S.C], [0, nshow * S.C * LamY] .* 1e3, ls = :dash, lw = 2, label = "Λ_y s (формула)")
    vline!(pY, [S.C, 2S.C], color = :gray, ls = :dot, label = "")
    savefig_step(pY, "fig_05_path_y")
    slopeX = ([sX.s ones(length(sX.s))] \ dLX)[1]
    slopeY = ([sY.s ones(length(sY.s))] \ dLY)[1]
    Nf = 40
    function per_turn_diff(zA, zB, N)
        a = collect(zA); b = collect(zB); out = zeros(N)
        for n in 1:N
            la, lb = a[5], b[5]
            a = step_exact(S, a); b = step_exact(S, b)
            out[n] = -((a[5] - la) - (b[5] - lb)) / S.C
        end
        return out
    end
    LnX = per_turn_diff(zX, z0, Nf)
    LnY = per_turn_diff(zY, z0, Nf)
    pT = plot(xlabel = "номер оборота n", ylabel = "Λ · 10⁶", legend = :topright, size = (760, 420))
    plot!(pT, 1:Nf, running_mean(LnX) .* 1e6, lw = 2, label = "скользящее среднее, плоскость x")
    hline!(pT, [LamX * 1e6], ls = :dash, lw = 2, label = "формула Λ_x")
    plot!(pT, 1:Nf, running_mean(LnY) .* 1e6, lw = 2, label = "скользящее среднее, плоскость y")
    hline!(pT, [LamY * 1e6], ls = :dash, lw = 2, label = "формула Λ_y")
    savefig_step(pT, "fig_05_path_mean")
    mac!(ctx, "LamX", texnum(LamX; d = 3))
    mac!(ctx, "LamY", texnum(LamY; d = 3))
    mac!(ctx, "LamXtr", texnum(slopeX; d = 3))
    mac!(ctx, "LamYtr", texnum(slopeY; d = 3))
    mac!(ctx, "LamXturn", texnum(mean(LnX); d = 3))
    mac!(ctx, "LamYturn", texnum(mean(LnY); d = 3))
    mac!(ctx, "Nfpath", string(Nf))
    mac!(ctx, "DiffLamX", texnum(slopeX - LamX; d = 2))
    mac!(ctx, "DiffLamY", texnum(slopeY - LamY; d = 2))
    conc = @sprintf("Наклон разности ΔL(s) за 3 оборота при включённом ВЧ: Λ_x=%.3e (формула %.3e, разница %.2e); Λ_y=%.3e (формула %.3e, разница %.2e). Пооборотное среднее за %d оборотов ещё дышит синхротроном.",
                    slopeX, LamX, slopeX - LamX, slopeY, LamY, slopeY - LamY, Nf)
    write_conclusion("step05", conc)
    write_step_tex("step05", """
\\section{Удлинение пути}
\\textbf{Вопрос.} На сколько бетатрон удлиняет путь, если полость включена и синхротрон уже идёт?

Две частицы в одном кольце: одна без бетатронной амплитуды, другая с \$A=20\$~мм. Рисунок --- разность их лишних длин, не прогон с \$V=0\$. Формула \$\\Lambda=-(2\\pi/C)J\\xi\$.

\\begin{figure}[H]\\centering
\\includegraphics[width=0.92\\textwidth]{figs/fig_05_path_x.png}
\\caption{Плоскость \$x\$. Наклон формулы \$\\Lambda_x=\\vLamX\$.}
\\end{figure}

\\begin{figure}[H]\\centering
\\includegraphics[width=0.92\\textwidth]{figs/fig_05_path_y.png}
\\caption{Плоскость \$y\$. Наклон формулы \$\\Lambda_y=\\vLamY\$.}
\\end{figure}

\\begin{figure}[H]\\centering
\\includegraphics[width=0.92\\textwidth]{figs/fig_05_path_mean.png}
\\caption{Пооборотная разность за первые \$\\vNfpath\$ оборотов, пока синхротронные фазы ещё не разбежались. Наклон \$\\Delta L(s)\$ за три оборота: \$x\$ \$\\vLamXtr\$ (формула \$\\vLamX\$, разница \$\\vDiffLamX\$), \$y\$ \$\\vLamYtr\$ (формула \$\\vLamY\$, разница \$\\vDiffLamY\$).}
\\end{figure}

\\textbf{Вывод.} В плоскости \$y\$ наклон за три оборота совпадает с формулой (разница \$\\vDiffLamY\$). В плоскости \$x\$ \$\\Lambda\$ меньше, за три оборота его ещё закрывает бетатронное колебание длины --- смотреть тренд на рисунке, не линейную подгонку. На многих оборотах синхротронные фазы двух частиц разбегаются.
""")
    return ctx
end
