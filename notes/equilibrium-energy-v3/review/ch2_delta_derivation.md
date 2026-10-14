# Вывод главы 2: бетатронное запаздывание в переменной δ

Все формулы получены заново из определения периода и из главы 1. Старый текст `sections/2_mag_beta.tex` не использовался. Проверки с фиксированными числами — в конце.

Обозначения главы 1: $\delta=(p-p_s)/p_s$, $\tau=(T-T_s)/T_s$, $T_s=C_s/v_s$, $\eta_0=\alpha_0-1/\gamma_s^2$, $\eta_1=\alpha_1-\alpha_0/\gamma_s^2+(3\gamma_s^2-1)/(2\gamma_s^4)$, $\omega_s^2=K\omega_{\mathrm{rf}}|\eta_0|$, $K\cos\varphi_s$ входит так, что $\kappa\omega_{\mathrm{rf}}=-\omega_s^2/\eta_0$ при $\kappa=K\cos\varphi_s$. Множитель зазора $1-\delta/\gamma_s^2$ сохранён.

## 1. Счёт порядков

Поперечная амплитуда и $\delta$ — один малый параметр. Действие $J\sim(\text{амплитуда})^2$. В пути удерживаются $\delta^2$ и члены, линейные по $J$. Отбрасываются $\delta^3$, $\delta J$ и $J^2$.

Отсюда же: функции Твисса берутся при $\delta=0$. Поправка $\gamma_x(\delta)-\gamma_x(0)$, умноженная на $J$, есть $\delta J$.

## 2. Разложение элемента пути

Метрика плоского кольца без кручения:

$$
\mathrm{d}l^2=(1+x/\rho)^2\mathrm{d}s^2+\mathrm{d}x^2+\mathrm{d}y^2.
$$

$$
\frac{\mathrm{d}l}{\mathrm{d}s}=\sqrt{(1+x/\rho)^2+x'^2+y'^2},\qquad x'=\mathrm{d}x/\mathrm{d}s.
$$

Подкоренное выражение равно $1+u$,

$$
u=2\frac{x}{\rho}+\Bigl(\frac{x}{\rho}\Bigr)^2+x'^2+y'^2.
$$

$$
(1+u)^{1/2}=1+\frac{u}{2}-\frac{u^2}{8}+O(u^3).
$$

$$
\frac{u}{2}=\frac{x}{\rho}+\frac{1}{2}\Bigl(\frac{x}{\rho}\Bigr)^2+\frac{x'^2+y'^2}{2}.
$$

Ведущий кусок $u^2/8$ есть $(2x/\rho)^2/8=x^2/(2\rho^2)$. Он сокращает $(x/\rho)^2/2$. Следующие куски $u^2$ содержат $(x/\rho)(x'^2+y'^2)$ — третий порядок по амплитуде — и $(x'^2)^2$ — четвёртый. До второго порядка

$$
\mathrm{d}l=\Bigl(1+\frac{x}{\rho}+\frac{x'^2+y'^2}{2}\Bigr)\mathrm{d}s.
$$

Члена $x^2/\rho^2$ не остаётся, поэтому нет ни $\langle x_\beta^2\rangle/\rho^2$, ни перекрёстного $x_D x_\beta/\rho^2$.

## 3. Почему время совпадает с путём и почему коэффициент при $\Lambda$ равен 1

$$
T=\oint\frac{\mathrm{d}l}{v}.
$$

Магнитная сила работы не совершает, $|p|$ между зазорами постоянен, $v=v(\delta)$ одна и та же на траектории с бетатронным движением и на замкнутой орбите $J=0$ при том же $\delta$. Поэтому $T=L/v$, $T_{J=0}=C(\delta)/v$ и

$$
\frac{T-T_{J=0}}{T_{J=0}}=\frac{L-C(\delta)}{C(\delta)}.
$$

Определение $\tau$ разбивается точно:

$$
\tau=\frac{T-T_s}{T_s}
=\frac{T_{J=0}-T_s}{T_s}+\frac{T-T_{J=0}}{T_s}.
$$

Первое слагаемое есть глава 1: орбита $J=0$ имеет длину $C(\delta)$ и скорость $v(\delta)$, откуда $\eta_0\delta+\eta_1\delta^2$. Второе после среднего по быстрым фазам есть

$$
\Lambda=\frac{\langle T-T_{J=0}\rangle_\psi}{T_s}.
$$

Коэффициент равен 1, потому что оба члена — части одного отношения $(T-T_s)/T_s$.

Перевод в длину:

$$
\frac{T-T_{J=0}}{T_s}=\frac{L-C(\delta)}{vT_s}=\frac{L-C(\delta)}{C_s}\frac{v_s}{v},
$$

поскольку $vT_s=v C_s/v_s$. По главе 1 $v_s/v=1+O(\delta)$. Избыток длины есть $O(J)$, произведение $O(\delta J)$ отбрасывается. На удержанном порядке

$$
\Lambda=\frac{\langle L-C(\delta)\rangle_\psi}{C_s},
\qquad
\tau=\eta_0\delta+\eta_1\delta^2+\Lambda.
$$

Часть $x_D$ уже сидит в $C(\delta)$, то есть в $\alpha_0$ и $\alpha_1$. Новые вклады — $x_\beta$ и $x_{\mathrm{nl}}$.

## 4. Перекрёстные члены $x_D$ и $x_\beta$

$$
x=x_D+x_\beta+x_{\mathrm{nl}},\qquad x_D=D\delta+D_2\delta^2,\qquad x_\beta=\sqrt{2J_x\beta_x}\cos\psi_x,
$$

$$
x'=x_D'+x_\beta'+x_{\mathrm{nl}}'.
$$

Избыток над орбитой $x=x_D$:

$$
\begin{aligned}
L-C(\delta)=\oint\Biggl[&
\frac{x_\beta+x_{\mathrm{nl}}}{\rho}
+\frac{x_\beta'^2+y_\beta'^2}{2}
+x_D'x_\beta'+x_D'x_{\mathrm{nl}}'\\
&+x_\beta'x_{\mathrm{nl}}'+\frac{x_{\mathrm{nl}}'^2}{2}+\cdots
\Biggr]\mathrm{d}s.
\end{aligned}
$$

- $x_\beta/\rho$ линейно по амплитуде, $\langle\cos\psi_x\rangle_\psi=0$.
- $x_D'x_\beta'$ имеет второй порядок (оба множителя первого), но $x_D'$ от фазы не зависит и $\langle x_\beta'\rangle_\psi=0$, поэтому среднее равно нулю.
- $x_D'x_{\mathrm{nl}}'$: $x_{\mathrm{nl}}\sim J\sim(\text{амплитуда})^2$, $x_D'\sim\delta$, произведение есть $\delta J$, третий порядок.
- $x_\beta'x_{\mathrm{nl}}'$ и $x_{\mathrm{nl}}'^2$ выше второго порядка.
- $y_{\mathrm{nl}}$ в плоском кольце не входит в член $y/\rho$; $y_{\mathrm{nl}}'^2$ четвёртого порядка.

Остаётся

$$
\Lambda=\frac{1}{C_s}\oint\Biggl(
\frac{\langle x_{\mathrm{nl}}\rangle_\psi}{\rho}
+\frac{\langle x_\beta'^2\rangle_\psi}{2}
+\frac{\langle y_\beta'^2\rangle_\psi}{2}
\Biggr)\mathrm{d}s.
$$

## 5. Средний квадрат наклона

$$
\alpha_x=-\beta_x'/2,\qquad \psi_x'=1/\beta_x,\qquad \gamma_x=(1+\alpha_x^2)/\beta_x.
$$

$$
x_\beta=\sqrt{2J_x\beta_x}\cos\psi_x.
$$

$$
\begin{aligned}
x_\beta'
&=\sqrt{2J_x}\Bigl(\frac{\beta_x'}{2\sqrt{\beta_x}}\cos\psi_x-\sqrt{\beta_x}\cdot\frac{1}{\beta_x}\sin\psi_x\Bigr)\\
&=\sqrt{2J_x}\Bigl(\frac{\beta_x'}{2\sqrt{\beta_x}}\cos\psi_x-\frac{\sin\psi_x}{\sqrt{\beta_x}}\Bigr).
\end{aligned}
$$

$\beta_x'=-2\alpha_x$, поэтому $\beta_x'/(2\sqrt{\beta_x})=-\alpha_x/\sqrt{\beta_x}$ и

$$
x_\beta'=-\sqrt{\frac{2J_x}{\beta_x}}(\sin\psi_x+\alpha_x\cos\psi_x).
$$

$$
x_\beta'^2=\frac{2J_x}{\beta_x}(\sin\psi_x+\alpha_x\cos\psi_x)^2.
$$

$$
\langle\sin^2\rangle_\psi=\langle\cos^2\rangle_\psi=\frac12,\qquad\langle\sin\cos\rangle_\psi=0,
$$

$$
\langle x_\beta'^2\rangle_\psi=\frac{2J_x}{\beta_x}\cdot\frac{1+\alpha_x^2}{2}=\frac{J_x}{\beta_x}(1+\alpha_x^2)=J_x\gamma_x.
$$

В пути стоит половина этого среднего:

$$
\oint\frac{\langle x_\beta'^2\rangle_\psi}{2}\mathrm{d}s=\frac{J_x}{2}\oint\gamma_x\,\mathrm{d}s.
$$

Итого

$$
\Lambda=\frac{1}{C_s}\Biggl[
\frac{J_x}{2}\oint\gamma_x\,\mathrm{d}s
+\frac{J_y}{2}\oint\gamma_y\,\mathrm{d}s
+\oint\frac{\langle x_{\mathrm{nl}}\rangle_\psi}{\rho}\,\mathrm{d}s
\Biggr].
$$

При постоянных $\beta$ и $\alpha=0$ выполнено $\gamma=1/\beta$ и $\oint\mathrm{d}s/\beta=2\pi\nu$. Без сдвига орбиты

$$
\frac{J}{2}\oint\gamma\,\mathrm{d}s=\frac{J}{2}\cdot 2\pi\nu=\pi J\nu,
\qquad
\Lambda=\frac{\pi}{C_s}(J_x\nu_x+J_y\nu_y).
$$

При $\alpha\neq 0$ равенство $\oint\gamma\,\mathrm{d}s=2\pi\nu$ ложно: $\gamma=(1+\alpha^2)/\beta\ge 1/\beta$. Подстановка $\pi J\nu/C_s$ вместо интеграла $\gamma$ тогда неверна. Это проверено на ячейке FODO (числа ниже).

## 6. Вторая производная и тождество $\oint\gamma=\oint\beta K$

Дифференцирование $x_\beta'$ при $\mathrm{d}(\beta^{-1/2})/\mathrm{d}s=\alpha\beta^{-3/2}$ (из $\beta'=-2\alpha$) и $\mathrm{d}(\sin\psi)/\mathrm{d}s=\cos\psi/\beta$. Символьно (SymPy, подстановки $\psi'=1/\beta$, $\beta'=-2\alpha$):

$$
x_\beta''=-\sqrt{2J_x}\;\beta_x^{-3/2}(1+\alpha_x^2+\alpha_x'\beta_x)\cos\psi_x.
$$

Уравнение $x_\beta''+K_x x_\beta=0$ с $x_\beta=\sqrt{2J_x\beta_x}\cos\psi_x$:

$$
-\beta_x^{-3/2}(1+\alpha_x^2+\alpha_x'\beta_x)+K_x\beta_x^{1/2}=0,
$$

$$
K_x\beta_x^2=1+\alpha_x^2+\alpha_x'\beta_x,
\qquad
\alpha_x'=K_x\beta_x-\gamma_x.
$$

Интеграл $\alpha_x'$ по периоду равен нулю, поэтому

$$
\oint\gamma_x\,\mathrm{d}s=\oint\beta_x K_x\,\mathrm{d}s.
$$

$K_x$ — коэффициент в $x''+K_x x=0$ на синхронном импульсе, включая $1/\rho^2$ диполей. Для вертикали $K_y$ берётся из $y''+K_y y=0$, со своим знаком.

## 7. Хроматичность линейной решётки

Уравнение вокруг замкнутой орбиты: $x_\beta''+K_x x_\beta/(1+\delta)=0$. Жёсткость $K_x/(1+\delta)=K_x-K_x\delta+O(\delta^2)$.

Гамильтониан при $\delta=0$: $p_x^2/2+K_x x^2/2$, добавка $\Delta K_x\, x^2/2$. Среднее $\langle x^2\rangle_\psi=J_x\beta_x$, средняя добавка $\Delta K_x\beta_x J_x/2$. Производная по действию — добавочный набег фазы на единицу длины $\Delta K_x\beta_x/2$. За оборот

$$
\Delta\nu_x=\frac{1}{4\pi}\oint\beta_x\Delta K_x\,\mathrm{d}s.
$$

Проверка множителя на постоянной $\beta$: $K=1/\beta^2$, $\alpha=0$, $x=\sqrt{2J\beta}\cos\psi$, $p=-\sqrt{2J/\beta}\sin\psi$,

$$
\frac{p^2}{2}+\frac{Kx^2}{2}=\frac{J}{\beta}(\sin^2+\cos^2)=\frac{J}{\beta},
\qquad
\frac{\partial H}{\partial J}=\frac{1}{\beta}=\psi'.
$$

Для хроматичности $\Delta K_x=-K_x\delta$,

$$
\xi_x=-\frac{1}{4\pi}\oint\beta_x K_x\,\mathrm{d}s.
$$

Вместе с тождеством пункта 6:

$$
\frac{J_x}{2}\oint\gamma_x\,\mathrm{d}s=-2\pi\xi_x J_x.
$$

При $x_{\mathrm{nl}}=0$ геометрическая формула совпадает с

$$
\Lambda=-\frac{2\pi}{C_s}(J_x\xi_x+J_y\xi_y).
$$

### Однородный канал, проверка знака

$K=\mathrm{const}>0$, секступолей нет, $\beta=1/\sqrt{K}$, $\alpha=0$, $\gamma=\sqrt{K}$, $\nu=C_s\sqrt{K}/(2\pi)$, $x_{\mathrm{nl}}=0$.

Из $x''+Kx/(1+\delta)=0$ волновое число $\sqrt{K/(1+\delta)}$,

$$
\nu(\delta)=\nu(1+\delta)^{-1/2}= \nu\Bigl(1-\frac{\delta}{2}+O(\delta^2)\Bigr),
\qquad
\xi=-\frac{\nu}{2}.
$$

Та же величина из интеграла: $\oint\beta K\,\mathrm{d}s=\sqrt{K}\,C_s=2\pi\nu$, $\xi=-(2\pi\nu)/(4\pi)=-\nu/2$.

Путь: $\Lambda=\pi J\nu/C_s>0$.

Хроматическая запись:

$$
-\frac{2\pi}{C_s}\xi J=-\frac{2\pi}{C_s}\Bigl(-\frac{\nu}{2}\Bigr)J=\frac{\pi J\nu}{C_s}.
$$

Совпали. Знак минус перед $\xi$ обязателен: путь от наклонов длиннее, натуральная хроматичность отрицательна. Конечная разность $\nu(\delta)$ при $\delta=\pm 10^{-4}$, $K=0{,}25$, $C=40$ даёт $\xi$, согласованное с $-\nu/2$ до $10^{-8}$ (остаток порядка $\delta^2$).

Вертикаль с постоянной положительной $K_y$ имеет тот же $\xi_y=-\nu_y/2$. В FODO с равными силами линз $\beta_y$ в фокусирующей по $x$ линзе равно $\beta_x$ в дефокусирующей, и $\oint\beta_y K_y\,\mathrm{d}s=\oint\beta_x K_x\,\mathrm{d}s$, так что $\xi_y=\xi_x<0$ и тождество с $\gamma_y$ сохраняет тот же знак.

## 8. Полная хроматичность, симплектичность

Гамильтониан ($s$ — независимая переменная, импульсы в единицах $p_s$):

$$
H=\frac{p_x^2+p_y^2}{2(1+\delta)}+\frac{k_1}{2}(x^2-y^2)-\frac{x\delta}{\rho}
$$

плюс члены только от $\delta$. Отсюда $x'=p_x/(1+\delta)$, $p_x'=-k_1 x+\delta/\rho$, то есть $(1+\delta)x''+k_1 x=\delta/\rho$. Однородное уравнение для отклонения от замкнутой орбиты есть $x_\beta''+k_1 x_\beta/(1+\delta)=0$.

$$
\frac{\partial H}{\partial\delta}=-\frac{p_x^2+p_y^2}{2(1+\delta)^2}-\frac{x}{\rho}=-\frac{x'^2+y'^2}{2}-\frac{x}{\rho},
$$

потому что $p_x=(1+\delta)x'$. Это $z'=-({\mathrm{d}l}/{\mathrm{d}s}-1)$ в геометрической части. Член потенциала $\propto x^2$, не зависящий от $\delta$, меняет $K$ и в $z'$ не входит. Члены только от $\delta$ — скорость главы 1, без $J$. Поэтому бетатронная часть

$$
\Delta z_\beta=-\Delta L_\beta.
$$

Отображение за оборот в первом порядке по $J$: $\psi\to\psi+2\pi\nu(\delta)$, $J\to J$, $\delta\to\delta$, $z\to z+\Delta z(J,\delta)$, $\nu=\nu+\xi\delta$ при $J=0$. Форма $\mathrm{d}\psi\wedge\mathrm{d}J+\mathrm{d}z\wedge\mathrm{d}\delta$.

$$
\mathrm{d}\psi_2\wedge\mathrm{d}J_2=\mathrm{d}\psi\wedge\mathrm{d}J+2\pi\xi\,\mathrm{d}\delta\wedge\mathrm{d}J,
$$

$$
\mathrm{d}z_2\wedge\mathrm{d}\delta_2=\mathrm{d}z\wedge\mathrm{d}\delta+\frac{\partial\Delta z}{\partial J}\,\mathrm{d}J\wedge\mathrm{d}\delta.
$$

$\mathrm{d}J\wedge\mathrm{d}\delta=-\mathrm{d}\delta\wedge\mathrm{d}J$, а $(\partial\Delta z/\partial\delta)\,\mathrm{d}\delta\wedge\mathrm{d}\delta=0$. Лишний член формы есть

$$
\Bigl(2\pi\xi-\frac{\partial\Delta z}{\partial J}\Bigr)\mathrm{d}\delta\wedge\mathrm{d}J.
$$

Он равен нулю при $\partial\Delta z/\partial J=2\pi\xi$. Зависимость $\nu$ от $J$ даёт $\mathrm{d}J\wedge\mathrm{d}J=0$ и условия не меняет. Если $\xi=\xi(J)$, интегрирование производит $J^2$; в первом порядке $\xi$ берётся при $J=0$.

$$
\Delta z_\beta=2\pi(\xi_x J_x+\xi_y J_y),
\qquad
\Delta L_\beta=-2\pi(\xi_x J_x+\xi_y J_y),
$$

$$
\Lambda=-\frac{2\pi}{C_s}(J_x\xi_x+J_y\xi_y).
$$

Это первый порядок по $J$. Полная $\xi$ уже включает секступоли, поэтому $x_{\mathrm{nl}}$ не прибавляется к этой формуле второй раз. Геометрическая запись с $\gamma$ и $x_{\mathrm{nl}}$ равна ей: интеграл $\gamma$ даёт часть $\xi$ из $K$, интеграл $x_{\mathrm{nl}}/\rho$ — остаток.

### Секступоль, согласование двух записей

Тонкий секступоль, $\Delta x'=-(m/2)(x^2-y^2)$. Средний пинок $\theta=-(m/2)(J_x\beta_x-J_y\beta_y)$.

Замкнутая орбита от пина: $\oint x/\rho\,\mathrm{d}s=D(s_0)\theta$. Вывод: $\oint x/\rho\,\mathrm{d}s=\oint x(D''+KD)\,\mathrm{d}s=\oint D(x''+Kx)\,\mathrm{d}s$, а $x''+Kx$ равно нулю всюду, кроме пина, где интеграл равен $\theta$.

Подстановка feed-down $x=D\delta+x_\beta$ в пинок даёт $\Delta x'=-(m D\delta) x_\beta$, то есть тонкая линза силы $q=m D\delta$. Сдвиг тюна $\Delta\nu=\beta q/(4\pi)$, откуда вклад в хроматичность $\xi=\beta m D/(4\pi)$ (горизонталь, без вертикального среднего). Тогда

$$
-2\pi\xi J_x=-\frac{1}{2}m\beta_x D J_x=D\theta
$$

при $\theta=-(m/2)J_x\beta_x$. Путь от $x_{\mathrm{nl}}$ равен $-2\pi\xi_{\mathrm{sext}}J$. Численно тождество $\oint x/\rho=D\theta$ на однородном канале с $1/\rho$ (трапеции, 4000 шагов) сошлось с относительной ошибкой $2\cdot 10^{-6}$. Формула $\Delta\nu=\beta q/(4\pi)$ при $q=10^{-4}$ сошлась с матрицей оборота с относительной ошибкой $2\cdot 10^{-5}$; ошибка падает пропорционально $q$.

## 9. Продольная система

$\Lambda$ постоянна: $J$ за синхротронный период не меняется, а $\Lambda(\delta)-\Lambda(0)$ есть $O(\delta J)$.

В счёте этой главы поперечная амплитуда того же порядка, что $\delta_m$, поэтому $\Lambda=O(\delta_m^2)$, как $\eta_1\delta^2$.

$$
\frac{\mathrm{d}\varphi}{\mathrm{d}t}=\omega_{\mathrm{rf}}(\eta_0\delta+\eta_1\delta^2+\Lambda),
\qquad
\frac{\mathrm{d}\delta}{\mathrm{d}t}=K\cos\varphi_s\,\varphi\Bigl(1-\frac{\delta}{\gamma_s^2}\Bigr).
$$

$\langle\tau\rangle=0$ точно, из возвращения фазы, как в главе 1. Поэтому

$$
\delta_{\mathrm{eq}}=-\frac{\eta_1\langle\delta^2\rangle+\Lambda}{\eta_0}.
$$

Дифференцирование второго уравнения, $\kappa\varphi=\dot\delta/(1-\delta/\gamma_s^2)$:

$$
\ddot\delta=\kappa\omega_{\mathrm{rf}}(\eta_0\delta+\eta_1\delta^2+\Lambda)\Bigl(1-\frac{\delta}{\gamma_s^2}\Bigr)
-\frac{(\dot\delta)^2}{\gamma_s^2(1-\delta/\gamma_s^2)}.
$$

$$
(\eta_0\delta+\eta_1\delta^2+\Lambda)\Bigl(1-\frac{\delta}{\gamma_s^2}\Bigr)
=\eta_0\delta+\Bigl(\eta_1-\frac{\eta_0}{\gamma_s^2}\Bigr)\delta^2+\Lambda-\frac{\Lambda\delta}{\gamma_s^2}+O(\delta^3).
$$

$\Lambda\delta=O(\delta_m^3)$ отбрасывается. $(\dot\delta)^2\cdot O(\delta)$ тоже. С $\kappa\omega_{\mathrm{rf}}=-\omega_s^2/\eta_0$:

$$
\ddot\delta+\omega_s^2\delta
=-\frac{\omega_s^2}{\eta_0}\Bigl(\eta_1-\frac{\eta_0}{\gamma_s^2}\Bigr)\delta^2
-\frac{1}{\gamma_s^2}(\dot\delta)^2
-\frac{\omega_s^2\Lambda}{\eta_0}.
$$

Это уравнение приложения с $c_0=-\omega_s^2\Lambda/\eta_0$. Решение

$$
\delta=-\frac{\Lambda}{\eta_0}-\frac{\eta_1 a^2}{2\eta_0}
+a\cos\psi
+\frac{a^2}{6}\Bigl(\frac{\eta_1}{\eta_0}-\frac{2}{\gamma_s^2}\Bigr)\cos 2\psi,
\qquad\psi=\omega_s t.
$$

$c_0$ входит только в постоянную. Подстановка решения в уравнение, разложенная по амплитуде до порядка $a^2$ и $\Lambda$ (порядок $a^2$), даёт тождественный нуль (SymPy, коэффициент при $\varepsilon^1$ и $\varepsilon^2$). Остаток численно есть отброшенный член $\Lambda\cdot a$: при $a=1/500$, $\eta_0=1/50$, $\eta_1=1/20$, $\gamma_s^2=9$, $\Lambda=10^{-6}$ коэффициент $\cos\psi$ в невязке равен $-5{,}18\cdot 10^{-7}$, а оценка $2a(\Lambda/\eta_0)$ через множитель $(\omega^2/\eta_0)(\eta_1-\eta_0/\gamma_s^2)$ даёт $-5{,}3\cdot 10^{-7}$.

Разность $\delta(0)-\delta(\pi)=2a$: постоянная и $\cos 2\psi=1$ сокращаются. Производная содержит $\sin\psi$ и от $\Lambda$ не зависит, так что $a=\delta_m$ при тех же условиях, что в главе 1.

При $\Lambda=O(\delta_m^2)$ квадрат постоянной есть $O(\delta_m^4)$, поэтому $\langle\delta^2\rangle=\delta_m^2/2$ до $O(\delta_m^4)$ и

$$
\delta_{\mathrm{eq}}=-\frac{1}{\eta_0}\Biggl(\frac{\eta_1\delta_m^2}{2}+\Lambda\Biggr).
$$

Усреднение уравнения второго порядка: $\langle\ddot\delta\rangle=0$, на ведущей гармонике $\langle\dot\delta^2\rangle=\omega_s^2\langle\delta^2\rangle$, члены $\pm 1/\gamma_s^2$ сокращаются, $\Lambda$ остаётся. Множитель зазора в $\delta_{\mathrm{eq}}$ не входит.

Поправка к $\omega_s$ порядка $\delta_m^2$ не выписывается: её не фиксирует обрыв $\tau$ на $\delta^2$. Линейный по $\Lambda$ сдвиг частоты того же порядка.

При $|\Lambda|\sim|\eta_0|\delta_m$ сдвиг $-\Lambda/\eta_0$ сравним с $\delta_m$, и $\langle\delta^2\rangle$ содержит $\delta_{\mathrm{eq}}^2$. Условие $|\Lambda/(\eta_0\delta_m)|\ll 1$ это исключает.

Колебательная часть $\varphi$ совпадает с формулой главы 1: постоянная $\eta_0\delta_{\mathrm{eq}}+\eta_1\delta_m^2/2+\Lambda=0$, а $\Lambda$ рядом с гармоникой $\delta$ даёт $O(\delta_m^3)$.

## 10. Энергия

$\gamma$ зависит только от модуля $p$, то есть от $\delta$. Усреднение ряда главы 1:

$$
\gamma_{\mathrm{eff}}=\gamma_s+\beta_s^2\gamma_s\delta_{\mathrm{eq}}+\frac{\beta_s^2\delta_m^2}{4\gamma_s}+O(\delta_m^4),
$$

$$
\Delta E_{\mathrm{eq}}=\beta_s^2 E_s\Biggl(\delta_{\mathrm{eq}}+\frac{\delta_m^2}{4\gamma_s^2}\Biggr)+O(\delta_m^4 E_s).
$$

Вклад $\Lambda$: $-\beta_s^2\gamma_s\Lambda/\eta_0$ в $\gamma_{\mathrm{eff}}$ и $-\beta_s^2 E_s\Lambda/\eta_0$ в энергии.

## 11. Проверки фиксированными числами

Скрипт `review/_check_ch2_delta.py`.

1. Путь при $x/\rho=0{,}02$, $x'=0{,}03$, $y'=-0{,}04$: точный корень минус $1+x/\rho+(x'^2+y'^2)/2$ равен $-2{,}52\cdot 10^{-5}$. Кусок третьего порядка $(x/\rho)(x'^2+y'^2)/2=2{,}5\cdot 10^{-5}$. Остаток после его учёта — четвёртый порядок.

2. $\langle x'^2\rangle=J\gamma$ при $\beta=2{,}5$, $\alpha=0{,}4$, $J=1{,}7\cdot 10^{-6}$: разность нулевая.

3. Канал $K=0{,}25$, $C=40$: $\nu=3{,}18309886$, $\xi=-\nu/2=-1{,}59154943$, конечная разность $-1{,}59154944$. Обе формулы $\Lambda$ при $J=2{,}5\cdot 10^{-6}$ дают $6{,}25\cdot 10^{-7}$.

4. FODO, дрейф $L=2$, сила линзы $0{,}4\,\mathrm{м}^{-1}$: $\nu=0{,}13099$, $\oint\gamma\,\mathrm{d}s=1{,}74574312$, $\oint\beta K\,\mathrm{d}s$ то же до $2\cdot 10^{-13}$, $-4\pi\xi$ то же. При этом $2\pi\nu=0{,}82303$, а $\pi\nu=0{,}41152\neq \tfrac12\oint\gamma=0{,}87287$. Формула с $\nu$ вместо $\gamma$ на этой ячейке неверна; формула с $\xi$ верна.

5. Подстановка решения в уравнение до порядка $a^2$ и $\Lambda$: коэффициенты при $\varepsilon$ и $\varepsilon^2$ равны нулю.

6. Секступоль: $m=0{,}3$, $J=1{,}7\cdot 10^{-6}$, $\beta=2$, $D=0{,}4$ дают $\Delta L=-2{,}04\cdot 10^{-7}$ и из $-2\pi\xi J$, и из $D\theta$.

## 12. Сомнения

- Кольцо плоское, без кручения и без связи $x$–$y$. Иначе метрика и раздельные $J_x$, $J_y$ меняются.
- Линейная хроматичность выведена из $K\to K/(1+\delta)$. Если какой-то элемент масштабируется иначе, интеграл $\xi$ из $K$ меняется. Запись через полную $\mathrm{d}\nu/\mathrm{d}\delta$ при этом остаётся: она опирается на симплектичность и на $z'=-(\mathrm{d}l/\mathrm{d}s-1)$, а не на этот масштаб.
- Краевые поля не разложены по отдельности. Они входят в $x_{\mathrm{nl}}$ и в полную $\xi$.
- Поправка к синхротронной частоте порядка $\delta_m^2$, в том числе член, линейный по $\Lambda$, не приводится: обрыв $\tau$ на $\delta^2$ её не фиксирует.
- Формула Шоджи $\Delta C=-2\pi(\xi_x J_x+\xi_y J_y)$ при действии, равном половине инварианта Куранта–Снайдера, совпадает с выведенной здесь. Вывод в тексте самостоятельный; статья указана как источник симплектического соотношения.
