# Вывод главы 1 в переменной δ

Все формулы ниже получены заново. Проверка рациональной арифметикой — в конце.

Обозначения: $p=p_s(1+\delta)$, $E_s=\gamma_s mc^2$, $v_s=\beta_s c$, $\beta_s^2=1-1/\gamma_s^2$, $T_s=C_s/v_s$, $\omega_{\mathrm{rf}}=2\pi h/T_s$, $t=kT_s$. Ускорения нет, $\sin\varphi_s=0$.

## 1. Скорость

Точно $v=pc^2/E$, поэтому

$$
\frac{v_s}{v}=\frac{p_s}{p}\frac{E}{E_s}=\frac{\gamma}{\gamma_s(1+\delta)}.
$$

Точно $\gamma^2-1=(p/(mc))^2$ и $\gamma_s^2-1=(p_s/(mc))^2$, откуда

$$
\gamma^2=1+(\gamma_s^2-1)(1+\delta)^2.
$$

Раскрытие $(1+\delta)^2=1+2\delta+\delta^2$ и деление на $\gamma_s^2$:

$$
\frac{\gamma^2}{\gamma_s^2}=1+2b\delta+b\delta^2,
\qquad
b=1-\frac{1}{\gamma_s^2}.
$$

Это тождество, не ряд. Корень: $x=2b\delta+b\delta^2$, $\gamma/\gamma_s=(1+x)^{1/2}$.

$$
(1+x)^{1/2}=1+\frac{x}{2}-\frac{x^2}{8}+O(x^3),
\qquad
x^2=4b^2\delta^2+O(\delta^3).
$$

$$
\frac{x}{2}=b\delta+\frac{b}{2}\delta^2,
\qquad
\frac{x^2}{8}=\frac{b^2}{2}\delta^2+O(\delta^3).
$$

Коэффициент при $\delta^2$:

$$
\frac{b}{2}-\frac{b^2}{2}=\frac{b}{2}(1-b)=\frac{b}{2\gamma_s^2}=\frac{\gamma_s^2-1}{2\gamma_s^4},
$$

поскольку $1-b=1/\gamma_s^2$. Итак

$$
\frac{\gamma}{\gamma_s}=1+b\delta+\frac{\gamma_s^2-1}{2\gamma_s^4}\delta^2+O(\delta^3).
$$

$$
(1+\delta)^{-1}=1-\delta+\delta^2+O(\delta^3).
$$

Произведение:

- $\delta$: $b-1=-1/\gamma_s^2$;
- $\delta^2$: $1$ (из $\delta^2$ в обратном) $+\,(-b)$ (перекрёстный) $+\,(\gamma_s^2-1)/(2\gamma_s^4)$.

$$
1-b+\frac{\gamma_s^2-1}{2\gamma_s^4}
=\frac{1}{\gamma_s^2}+\frac{\gamma_s^2-1}{2\gamma_s^4}
=\frac{2\gamma_s^2+\gamma_s^2-1}{2\gamma_s^4}
=\frac{3\gamma_s^2-1}{2\gamma_s^4}.
$$

$$
\frac{v_s}{v}=1-\frac{\delta}{\gamma_s^2}+\frac{3\gamma_s^2-1}{2\gamma_s^4}\delta^2+O(\delta^3).
$$

Предел $\gamma_s\to 1$: $v\propto p$, $v_s/v=1/(1+\delta)=1-\delta+\delta^2+O(\delta^3)$. Формула даёт $1-\delta+(3-1)/2\,\delta^2=1-\delta+\delta^2$.

## 2. Запаздывание

$T=C/v$, поэтому $T/T_s=(C/C_s)(v_s/v)$ и $\tau=T/T_s-1$.

$$
\frac{C}{C_s}=1+\alpha_0\delta+\alpha_1\delta^2+O(\delta^3).
$$

Произведение до $\delta^2$:

$$
\Bigl(1+\alpha_0\delta+\alpha_1\delta^2\Bigr)
\Bigl(1-\frac{\delta}{\gamma_s^2}+\frac{3\gamma_s^2-1}{2\gamma_s^4}\delta^2\Bigr)
=1+\Bigl(\alpha_0-\frac{1}{\gamma_s^2}\Bigr)\delta
+\Bigl(\alpha_1-\frac{\alpha_0}{\gamma_s^2}+\frac{3\gamma_s^2-1}{2\gamma_s^4}\Bigr)\delta^2+O(\delta^3).
$$

Перекрёстный вклад в $\delta^2$ есть $\alpha_0\delta\cdot(-\delta/\gamma_s^2)$. После вычитания единицы и отбрасывания $O(\delta^3)$

$$
\tau=\eta_0\delta+\eta_1\delta^2,
\qquad
\eta_0=\alpha_0-\frac{1}{\gamma_s^2},
\qquad
\eta_1=\alpha_1-\frac{\alpha_0}{\gamma_s^2}+\frac{3\gamma_s^2-1}{2\gamma_s^4}.
$$

## 3. Отображение и непрерывная система

Опоздание $\tau T_s$ сдвигает фазу на $\omega_{\mathrm{rf}}\tau T_s=2\pi h\tau$. При $\tau>0$ частица позже, $\varphi$ растёт:

$$
\varphi_{k+1}-\varphi_k=2\pi h\,\tau_k.
$$

Это точно для точного $\tau=(T-T_s)/T_s$. Дальше $\tau$ — усечённый ряд.

Зазор: $\Delta E=qV\sin(\varphi_s+\varphi)=qV\cos\varphi_s\sin\varphi$ точно, потому что $\sin\varphi_s=0$.

Из $E\,dE=pc^2\,dp$ следует $dE=v\,dp$. При $dp=p_s\,d\delta$

$$
d\delta=\frac{dE}{vp_s}.
$$

$v_sp_s=\beta_sc\cdot(E_s\beta_s/c)=\beta_s^2E_s$, поэтому

$$
d\delta=\frac{v_s}{v}\frac{dE}{\beta_s^2E_s}.
$$

Скачок за оборот мал по сравнению с $E_s$ (много оборотов на синхротронный период), и тот же множитель берётся на $\delta$ данного оборота. В $v_s/v$ нужен член, линейный по $\delta$: член $\delta^2$, умноженный на $\sin\varphi=O(\delta_m)$, даёт $O(\delta_m^3)$.

$$
\delta_{k+1}-\delta_k
=\frac{qV}{\beta_s^2E_s}\cos\varphi_s\sin\varphi_k\Bigl(1-\frac{\delta_k}{\gamma_s^2}\Bigr).
$$

Деление на $T_s$, $K=qV/(\beta_s^2E_sT_s)$:

$$
\frac{d\varphi}{dt}=\omega_{\mathrm{rf}}\tau,
\qquad
\frac{d\delta}{dt}=K\cos\varphi_s\sin\varphi\Bigl(1-\frac{\delta}{\gamma_s^2}\Bigr).
$$

Линейный предел $\tau=\eta_0\delta$, $\sin\varphi=\varphi$, множитель скорости равен 1:

$$
\ddot\delta=K\omega_{\mathrm{rf}}\eta_0\cos\varphi_s\,\delta.
$$

Колебания ограничены при $K\omega_{\mathrm{rf}}\eta_0\cos\varphi_s<0$. При $K>0$ и $\omega_{\mathrm{rf}}>0$ это $\eta_0\cos\varphi_s<0$, то есть

$$
\cos\varphi_s=-\operatorname{sgn}\eta_0,
\qquad
\omega_s^2=K\omega_{\mathrm{rf}}|\eta_0|,
\qquad
K\omega_{\mathrm{rf}}\cos\varphi_s=-\frac{\omega_s^2}{\eta_0}.
$$

Второе равенство для $\omega_s^2$: $\omega_s^2=-K\omega_{\mathrm{rf}}\eta_0\cos\varphi_s=K\omega_{\mathrm{rf}}\eta_0\operatorname{sgn}\eta_0=K\omega_{\mathrm{rf}}|\eta_0|$.

## 4. Усреднение до решения

На замкнутой траектории внутри сепаратрисы сумма приращений фазы за период равна нулю. Она равна $2\pi h\sum\tau_k$. Отсюда $\sum\tau_k=0$ и $\langle\tau\rangle=0$ точно. В усечённом ряде

$$
\langle\tau\rangle=\eta_0\langle\delta\rangle+\eta_1\langle\delta^2\rangle=0,
\qquad
\delta_{\mathrm{eq}}=-\frac{\eta_1}{\eta_0}\langle\delta^2\rangle.
$$

$\delta_{\mathrm{eq}}$ уже найден. Решение нужно только для $\langle\delta^2\rangle$.

## 5. Уравнение для δ

$\sin\varphi\approx\varphi$, $\kappa=K\cos\varphi_s$, $\kappa\omega_{\mathrm{rf}}=-\omega_s^2/\eta_0$.

$$
\dot\varphi=\omega_{\mathrm{rf}}(\eta_0\delta+\eta_1\delta^2),
\qquad
\dot\delta=\kappa\varphi\Bigl(1-\frac{\delta}{\gamma_s^2}\Bigr).
$$

Дифференцирование второго:

$$
\ddot\delta=\kappa\dot\varphi\Bigl(1-\frac{\delta}{\gamma_s^2}\Bigr)-\kappa\varphi\frac{\dot\delta}{\gamma_s^2}.
$$

Из второго уравнения $\kappa\varphi=\dot\delta/(1-\delta/\gamma_s^2)$, и последнее слагаемое равно $(\dot\delta)^2/[\gamma_s^2(1-\delta/\gamma_s^2)]$. Подстановка $\dot\varphi$:

$$
\ddot\delta
=\kappa\omega_{\mathrm{rf}}(\eta_0\delta+\eta_1\delta^2)\Bigl(1-\frac{\delta}{\gamma_s^2}\Bigr)
-\frac{(\dot\delta)^2}{\gamma_s^2(1-\delta/\gamma_s^2)}.
$$

$$
(\eta_0\delta+\eta_1\delta^2)\Bigl(1-\frac{\delta}{\gamma_s^2}\Bigr)
=\eta_0\delta+\Bigl(\eta_1-\frac{\eta_0}{\gamma_s^2}\Bigr)\delta^2+O(\delta^3).
$$

$(\dot\delta)^2\cdot O(\delta)$ есть третий порядок по амплитуде. С $\kappa\omega_{\mathrm{rf}}=-\omega_s^2/\eta_0$:

$$
\ddot\delta+\omega_s^2\delta
=-\frac{\omega_s^2}{\eta_0}\Bigl(\eta_1-\frac{\eta_0}{\gamma_s^2}\Bigr)\delta^2
-\frac{1}{\gamma_s^2}(\dot\delta)^2.
$$

Если к $\tau$ добавить $\Lambda=O(\delta^2)$, в скобке появляется $\Lambda$. Умножение на $\kappa\omega_{\mathrm{rf}}$ даёт постоянную $-\omega_s^2\Lambda/\eta_0$. Член $\Lambda\delta$ есть $O(\delta_m^3)$ и отбрасывается. В приложении это $c_0=-\omega_s^2\Lambda/\eta_0$.

## 6. Метод Боголюбова–Митропольского

$$
\ddot\delta+\omega^2\delta=\varepsilon\bigl(c_0+c_2\delta^2+c_d(\dot\delta)^2\bigr),
\quad
\omega=\omega_s,
$$

$$
c_2=-\frac{\omega^2}{\eta_0}\Bigl(\eta_1-\frac{\eta_0}{\gamma_s^2}\Bigr),
\qquad
c_d=-\frac{1}{\gamma_s^2}.
$$

Подстановка $\delta=a\cos\psi+\varepsilon u_1$, $\dot a=\varepsilon A_1$, $\dot\psi=\omega+\varepsilon B_1$.

Производная ведущего члена:

$$
\frac{d}{dt}(a\cos\psi)=-a\omega\sin\psi+\varepsilon(A_1\cos\psi-a B_1\sin\psi)+O(\varepsilon^2).
$$

Производная $\varepsilon u_1$ на этом порядке есть $\varepsilon\omega\partial_\psi u_1$. Обозначим

$$
v_1=A_1\cos\psi-a B_1\sin\psi+\omega\partial_\psi u_1.
$$

Вторая производная $-a\omega\sin\psi$ даёт на порядке $\varepsilon$ вклад $-A_1\omega\sin\psi-a\omega B_1\cos\psi$. Вторая производная $\varepsilon v_1$ даёт $\varepsilon\omega\partial_\psi v_1$, и

$$
\partial_\psi v_1=-A_1\sin\psi-a B_1\cos\psi+\omega\partial_{\psi\psi}u_1
$$

снова даёт $-A_1\omega\sin\psi-a\omega B_1\cos\psi$. Два вклада складываются. С $\omega^2\varepsilon u_1$ левая часть порядка $\varepsilon$:

$$
\omega^2(\partial_{\psi\psi}u_1+u_1)-2\omega A_1\sin\psi-2a\omega B_1\cos\psi=f_1.
$$

Та же левая часть на порядке $\varepsilon^k$, потому что $A_k$, $B_k$, $u_k$ входят линейно так же, а произведения младших поправок уходят в $f_k$.

Правило: $\partial_{\psi\psi}\cos(m\psi)=-m^2\cos(m\psi)$, значит

$$
(\partial_{\psi\psi}+1)\cos(m\psi)=(1-m^2)\cos(m\psi).
$$

При $m\neq 1$ коэффициент $F_m$ гармоники в правой части $\omega^2(\partial_{\psi\psi}+1)u$ переходит в $u$ как $F_m/[\omega^2(1-m^2)]$. При $m=1$ оператор даёт нуль, и гармоника определяет $A_k$, $B_k$.

На ведущем решении $\delta=a\cos\psi$, $\dot\delta=-a\omega\sin\psi$:

$$
\cos^2\psi=\frac{1+\cos 2\psi}{2},
\qquad
\sin^2\psi=\frac{1-\cos 2\psi}{2},
$$

$$
f_1=c_0+\frac{a^2}{2}(c_2+c_d\omega^2)+\frac{a^2}{2}(c_2-c_d\omega^2)\cos 2\psi.
$$

Первой гармоники нет, поэтому $A_1=0$, $B_1=0$: амплитуда постоянна, частота на этом порядке равна $\omega$.

Постоянная: $(\partial_{\psi\psi}+1)$ на константе равен 1, 

$$
\omega^2 U_0=c_0+\frac{a^2}{2}(c_2+c_d\omega^2).
$$

$\cos 2\psi$: $1-4=-3$,

$$
U_2=-\frac{a^2(c_2-c_d\omega^2)}{6\omega^2}.
$$

Подстановка коэффициентов:

$$
c_2+c_d\omega^2
=-\frac{\omega^2}{\eta_0}\eta_1+\frac{\omega^2}{\gamma_s^2}-\frac{\omega^2}{\gamma_s^2}
=-\frac{\omega^2\eta_1}{\eta_0}.
$$

Два члена $\omega^2/\gamma_s^2$ гасят друг друга. Отсюда

$$
U_0=\frac{c_0}{\omega^2}-\frac{\eta_1 a^2}{2\eta_0}.
$$

$$
c_2-c_d\omega^2=-\frac{\omega^2\eta_1}{\eta_0}+\frac{2\omega^2}{\gamma_s^2},
$$

потому что $-c_d\omega^2=+\omega^2/\gamma_s^2$ добавляется к $\omega^2/\gamma_s^2$ из раскрытия $c_2$. Тогда

$$
U_2=\frac{a^2}{6}\Bigl(\frac{\eta_1}{\eta_0}-\frac{2}{\gamma_s^2}\Bigr).
$$

При $\varepsilon=1$ и $\psi=\omega_s t$

$$
\delta=\frac{c_0}{\omega_s^2}-\frac{\eta_1 a^2}{2\eta_0}
+a\cos\psi
+\frac{a^2}{6}\Bigl(\frac{\eta_1}{\eta_0}-\frac{2}{\gamma_s^2}\Bigr)\cos 2\psi.
$$

В главе 1 $c_0=0$.

Проекция на $\cos\psi$, дающая поправку частоты порядка $a^2$, появляется на следующем порядке метода. Её физическое значение не фиксируется обрывом $\tau$ на $\delta^2$: член $\eta_2\delta^3$ в $\tau$ и $\varphi^3/6$ в $\sin\varphi$ дают вклад того же порядка в уравнение для частоты. Поправка не выписывается. В решении $\psi=\omega_s t$.

## 7. Полуразмах и среднее

$$
\dot\delta=-\omega_s\sin\psi\,(a+4D_2\cos\psi),
$$

где $D_2$ — коэффициент при $\cos 2\psi$. Множитель 4: производная $\cos 2\psi$ даёт $-2\omega_s\sin 2\psi$, а $\sin 2\psi=2\sin\psi\cos\psi$. Нули при $\sin\psi=0$, то есть $\psi=0,\pi$. В обеих точках $\cos 2\psi=1$, разность $\delta(0)-\delta(\pi)=2a$. При $|\eta_1 a/\eta_0|\ll 1$ и $a\ll 1$ других нулей нет. При $a>0$

$$
\delta_m=\frac{\delta_{\max}-\delta_{\min}}{2}=a.
$$

$$
\langle\delta^2\rangle=\frac{\delta_m^2}{2}+O(\delta_m^4):
$$

$\langle\cos^2\psi\rangle=1/2$, средние разных гармоник равны нулю, $\delta_{\mathrm{eq}}^2$ и $D_2^2$ порядка $\delta_m^4$.

Постоянный член $-\eta_1\delta_m^2/(2\eta_0)=-(\eta_1/\eta_0)\langle\delta^2\rangle$ совпадает с $\delta_{\mathrm{eq}}$ из фазового условия. Проверка усреднением уравнения: $\langle\ddot\delta\rangle=0$, потому что $\dot\delta$ периодична,

$$
\omega_s^2\delta_{\mathrm{eq}}
=-\frac{\omega_s^2}{\eta_0}\Bigl(\eta_1-\frac{\eta_0}{\gamma_s^2}\Bigr)\langle\delta^2\rangle
-\frac{1}{\gamma_s^2}\langle(\dot\delta)^2\rangle.
$$

На ведущей гармонике $\langle(\dot\delta)^2\rangle=\omega_s^2\delta_m^2/2=\omega_s^2\langle\delta^2\rangle$, члены $\pm 1/\gamma_s^2$ сокращаются.

$$
\delta(t)=\delta_{\mathrm{eq}}+\delta_m\cos(\omega_s t)
+\frac{\delta_m^2}{6}\Bigl(\frac{\eta_1}{\eta_0}-\frac{2}{\gamma_s^2}\Bigr)\cos(2\omega_s t),
\qquad
\delta_{\mathrm{eq}}=-\frac{\eta_1\delta_m^2}{2\eta_0}.
$$

## 8. Фаза

$\delta^2=(\delta_m^2/2)(1+\cos 2\psi)+O(\delta_m^3)$. Постоянная $\eta_0\delta_{\mathrm{eq}}+\eta_1\delta_m^2/2=0$. Коэффициент при $\cos\psi$ равен $\eta_0\delta_m$. Коэффициент при $\cos 2\psi$:

$$
\eta_0\cdot\frac{\delta_m^2}{6}\Bigl(\frac{\eta_1}{\eta_0}-\frac{2}{\gamma_s^2}\Bigr)+\frac{\eta_1\delta_m^2}{2}
=\frac{\eta_1\delta_m^2}{6}-\frac{\eta_0\delta_m^2}{3\gamma_s^2}+\frac{\eta_1\delta_m^2}{2}
=\frac{2\eta_1}{3}\delta_m^2-\frac{\eta_0}{3\gamma_s^2}\delta_m^2.
$$

$\int\cos\psi\,dt=\sin\psi/\omega_s$, $\int\cos 2\psi\,dt=\sin 2\psi/(2\omega_s)$. Постоянная интегрирования равна нулю: среднее $\tau$ равно нулю, так что $\varphi$ за период не нарастает, а при $\psi=0$ выполнено $\dot\delta=0$ и из $\dot\delta=\kappa\varphi(1-\delta/\gamma_s^2)$ следует $\varphi=0$. Умножение на $\omega_{\mathrm{rf}}$:

$$
\varphi=\frac{\omega_{\mathrm{rf}}}{\omega_s}
\Biggl[
\eta_0\delta_m\sin\psi
+\delta_m^2\Bigl(\frac{\eta_1}{3}-\frac{\eta_0}{6\gamma_s^2}\Bigr)\sin 2\psi
\Biggr],
\qquad
\psi=\omega_s t.
$$

Половина от $2\eta_1/3$ есть $\eta_1/3$, половина от $\eta_0/(3\gamma_s^2)$ есть $\eta_0/(6\gamma_s^2)$.

## 9. Энергия и γ

$$
\gamma=\gamma_s\Bigl[1+b\delta+\frac{\gamma_s^2-1}{2\gamma_s^4}\delta^2\Bigr]+O(\delta^3)
=\gamma_s+\beta_s^2\gamma_s\,\delta+\frac{\beta_s^2}{2\gamma_s}\delta^2+O(\delta^3).
$$

Проверка коэффициентов: $\gamma_s b=\gamma_s(1-1/\gamma_s^2)=\beta_s^2\gamma_s$ и

$$
\gamma_s\cdot\frac{\gamma_s^2-1}{2\gamma_s^4}=\frac{\gamma_s^2-1}{2\gamma_s^3}=\frac{\beta_s^2\gamma_s^2}{2\gamma_s^3}=\frac{\beta_s^2}{2\gamma_s}.
$$

Среднее, $\langle\delta^2\rangle=\delta_m^2/2$:

$$
\gamma_{\mathrm{eff}}=\langle\gamma\rangle
=\gamma_s+\beta_s^2\gamma_s\,\delta_{\mathrm{eq}}+\frac{\beta_s^2}{2\gamma_s}\cdot\frac{\delta_m^2}{2}+O(\delta_m^4)
=\gamma_s+\beta_s^2\gamma_s\,\delta_{\mathrm{eq}}+\frac{\beta_s^2\delta_m^2}{4\gamma_s}+O(\delta_m^4).
$$

$\langle\delta^3\rangle=O(\delta_m^4)$, потому что $\langle\cos^3\psi\rangle=0$, а постоянная и вторая гармоника порядка $\delta_m^2$.

Запись $\gamma_{\mathrm{eff}}=\gamma_s+\beta_s^2\gamma_s\delta_{\mathrm{eq}}$ удерживает первые два члена. Добавка полного среднего есть $\beta_s^2\delta_m^2/(4\gamma_s)$. Оба слагаемых сдвига — второго порядка по $\delta_m$. Отношение добавки к члену с $\delta_{\mathrm{eq}}$:

$$
\frac{\beta_s^2\delta_m^2/(4\gamma_s)}{\beta_s^2\gamma_s|\delta_{\mathrm{eq}}|}
=\frac{\delta_m^2}{4\gamma_s^2\cdot|\eta_1|\delta_m^2/(2|\eta_0|)}
=\frac{|\eta_0|}{2|\eta_1|\gamma_s^2},
$$

где подставлено $|\delta_{\mathrm{eq}}|=|\eta_1|\delta_m^2/(2|\eta_0|)$.

$$
\Delta E_{\mathrm{eq}}=\langle E\rangle-E_s=E_s\frac{\gamma_{\mathrm{eff}}-\gamma_s}{\gamma_s}
=\beta_s^2 E_s\Bigl(\delta_{\mathrm{eq}}+\frac{\delta_m^2}{4\gamma_s^2}\Bigr)+O(\delta_m^4 E_s).
$$

## 10. Независимая проверка через отклонение энергии

Пусть $w=(E-E_s)/(\beta_s^2 E_s)$. Тейлор $E(p)$ с $dE/dp=v$ и $d^2E/dp^2=c^2/(\gamma^2 E)$ даёт $w=\delta+\delta^2/(2\gamma_s^2)+O(\delta^3)$. Зазор меняет энергию точно, поэтому $\dot w=\kappa\varphi$ без множителя скорости. Обращение ряда: $\delta=w-w^2/(2\gamma_s^2)+O(w^3)$, и

$$
\tau=\eta_0 w+\tilde\eta_1 w^2,
\qquad
\tilde\eta_1=\eta_1-\frac{\eta_0}{2\gamma_s^2}.
$$

Тогда $\ddot w+\omega_s^2 w=-(\omega_s^2\tilde\eta_1/\eta_0)w^2$. Гармонический баланс без $(\dot w)^2$: постоянная $-\tilde\eta_1 a^2/(2\eta_0)$, вторая гармоника $\tilde\eta_1 a^2/(6\eta_0)$.

Обратный переход $w=\delta+\delta^2/(2\gamma_s^2)$, $\delta^2=(a^2/2)(1+\cos 2\psi)+O(a^3)$:

$$
\text{постоянная }w=D_0+\frac{a^2}{4\gamma_s^2}
=-\frac{\eta_1 a^2}{2\eta_0}+\frac{a^2}{4\gamma_s^2},
$$

откуда $D_0=-\eta_1 a^2/(2\eta_0)$. Вторая гармоника:

$$
D_2+\frac{a^2}{4\gamma_s^2}
=\frac{\tilde\eta_1 a^2}{6\eta_0}
=\frac{\eta_1 a^2}{6\eta_0}-\frac{a^2}{12\gamma_s^2},
$$

$$
D_2=\frac{\eta_1 a^2}{6\eta_0}-\frac{a^2}{12\gamma_s^2}-\frac{a^2}{4\gamma_s^2}
=\frac{a^2}{6}\Bigl(\frac{\eta_1}{\eta_0}-\frac{2}{\gamma_s^2}\Bigr).
$$

Совпадает с прямым решением в $\delta$. Та же подстановка воспроизводит коэффициент при $\sin 2\psi$ в $\varphi$.

## 11. Что проверено числом

Скрипт на рациональных числах (SymPy, точные дроби):

- ряд $v_s/v$ и $\gamma(\delta)$ при $\gamma_s^2=9/4$ совпал с формулами;
- предел $\gamma_s^2=1$ совпал с $1/(1+\delta)$ до $\delta^2$;
- при $\eta_0=2$, $\eta_1=3$, $\omega^2=5$, $\gamma_s^2=4$, $a=1/5$ получено $c_2+c_d\omega^2=-\omega^2\eta_1/\eta_0$, $D_0=-3/100$, $D_2=1/150$;
- переход к $w$ на тех же дробях воспроизвёл постоянную и вторую гармонику;
- половина коэффициента $\cos 2\psi$ в $\tau$ совпала с коэффициентом $\sin 2\psi$ в $\varphi$;
- подстановка решения в уравнение для $\delta$ даёт тождественный нуль в невязке до порядка $a^2$ включительно;
- фаза, восстановленная из $\dot\delta=\kappa\varphi(1-\delta/\gamma_s^2)$, совпадает с проинтегрированным $\tau$ до порядка $a^2$.

Численное интегрирование усечённой системы ($\sin\varphi=\varphi$, $\tau$ до $\delta^2$, множитель $1-\delta/\gamma_s^2$) при $\eta_0=-0{,}53$, $\eta_1=1{,}02$, $\gamma_s=1{,}1439$, $\omega_\rf/\omega_s=3$ воспроизводит $\delta_\eq$ и коэффициент $\cos 2\psi$. Сдвиг частоты на этих же траекториях не совпадает с формулой, полученной только из квадратичного уравнения для $w$: поправка порядка $\delta_m^2$ чувствительна к отброшенным членам, и в тексте она не выписывается.

## 12. Оговорки

Поправка частоты порядка $\delta_m^2$ не выписана: она не замкнута без $\delta^3$ в $\tau$ и без $\varphi^3$.

Если в главе 2 постоянная $\Lambda$ не порядка $\delta^2$, а порядка самой амплитуды, в уравнении для $\delta$ появляется линейный член $(\omega_s^2\Lambda/(\eta_0\gamma_s^2))\delta$. Он сдвигает частоту на относительную величину порядка $\Lambda/(\eta_0\gamma_s^2)$. При счёте $\Lambda=O(\delta_m^2)$ этот член есть $O(\delta_m^3)$ и отброшен законно. Решение приложения даёт постоянную $c_0/\omega_s^2-\eta_1 a^2/(2\eta_0)$ до членов, линейных по $c_0$ и квадратичных по $a$; произведения $c_0 a$ и $c_0^2$ в него не входят.

Конечный скачок энергии заменён на $dE=v\,dp$ с $v$ на текущем $\delta$. Это отвечает малости прироста за один оборот.
