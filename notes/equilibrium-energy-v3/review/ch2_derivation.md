# Вывод формул главы 2 и поправок приложения A

Обозначения совпадают с запиской. Скобки \(\langle\,\rangle_\psi\) --- среднее по бетатронной фазе; скобки без индекса --- среднее за синхротронный период, как в главе 1. Индекс Твисса у \(\beta,\alpha,\gamma\) пишется явно, лоренц-фактор --- только \(\gamma_s\).

Удерживаемый порядок: члены \(A^2\) и \(\Lambda\), перекрёстные \(\Lambda A\) в первой гармонике \(\tau\) и \(\tilde\eta_1\Lambda\) в частоте. Отбрасываются \(A^3\), \(\Lambda^2\), \(\Lambda A^2\).

## 1. Элемент пути

Положение относительно плоской опорной орбиты: продольное смещение вдоль касательной плюс \(x\) по главной нормали и \(y\) по бинормали. Квадрат элемента длины

\[
\mathrm{d}l^2=\bigl((1+x/\rho)^2+x'^2+y'^2\bigr)\mathrm{d}s^2,
\]

где штрих есть \(\mathrm{d}/\mathrm{d}s\), кручение на этом порядке не удерживается. Тогда

\[
\frac{\mathrm{d}l}{\mathrm{d}s}=\sqrt{1+u},
\qquad
u=2\frac{x}{\rho}+\Bigl(\frac{x}{\rho}\Bigr)^2+x'^2+y'^2.
\]

Ряд \(\sqrt{1+u}=1+u/2-u^2/8+O(u^3)\). Величины \(x/\rho\) и углы \(x',y'\) одного порядка малости, поэтому \(O(u^3)\) есть третий порядок.

\[
\frac{u}{2}=\frac{x}{\rho}+\frac12\Bigl(\frac{x}{\rho}\Bigr)^2+\frac{x'^2+y'^2}{2}.
\]

Квадрат \(u^2=(2x/\rho)^2+O(3)=4(x/\rho)^2+O(3)\),

\[
\frac{u^2}{8}=\frac12\Bigl(\frac{x}{\rho}\Bigr)^2+O(3).
\]

Разность

\[
\frac{u}{2}-\frac{u^2}{8}
=\frac{x}{\rho}+\frac{x'^2+y'^2}{2}+O(3):
\]

члены \((x/\rho)^2/2\) и \(-u^2/8\) сокращаются. Итог

\[
\mathrm{d}l=\Bigl(1+\frac{x}{\rho}+\frac{x'^2+y'^2}{2}\Bigr)\mathrm{d}s.
\]

На орбите \(J=0\), \(w=0\) выполнено \(x=x'=y'=0\), поэтому \(\mathrm{d}l=\mathrm{d}s\) и длина равна \(C_s\). Избыток длины траектории

\[
\Delta L_\beta=\oint\Bigl(\frac{x}{\rho}+\frac{x'^2+y'^2}{2}\Bigr)\mathrm{d}s.
\]

## 2. Время и путь

В магнитном поле модуль скорости постоянен вдоль траектории и определяется энергией. При \(w=0\) энергия равна \(E_s\), скорость на траектории с данным \(J\) и на траектории \(J=0\) одна и та же, \(v=v_s\). Периоды

\[
T=\frac{C_s+\Delta L_\beta}{v_s},
\qquad
T(J=0)=\frac{C_s}{v_s}=T_s,
\]

\[
\frac{T-T(J=0)}{T_s}=\frac{\Delta L_\beta}{C_s}.
\]

После среднего по \(\psi\) это определение \(\Lambda\). Если скорости различаются, отношение времён не равно отношению длин: множитель \(v_s/v\) из главы 1 здесь равен единице только на обеих сравниваемых траекториях сразу.

## 3. Производная \(x_\beta\) и \(\langle x'^2\rangle\)

\[
x_\beta=\sqrt{2J_x\beta_x}\cos\psi_x,
\qquad
\frac{\mathrm{d}\psi_x}{\mathrm{d}s}=\frac{1}{\beta_x},
\qquad
\alpha_x=-\frac12\beta_x'.
\]

Действие \(J_x\) вдоль оборота постоянно. Дифференцирование:

\begin{align*}
x_\beta'
&=\sqrt{2J_x}\Bigl(\frac{\beta_x'}{2\sqrt{\beta_x}}\cos\psi_x-\sqrt{\beta_x}\sin\psi_x\cdot\frac{1}{\beta_x}\Bigr)\\
&=\sqrt{\frac{2J_x}{\beta_x}}\Bigl(\frac{\beta_x'}{2}\cos\psi_x-\sin\psi_x\Bigr)\\
&=-\sqrt{\frac{2J_x}{\beta_x}}\bigl(\sin\psi_x+\alpha_x\cos\psi_x\bigr),
\end{align*}

где \(\beta_x'/2=-\alpha_x\). Квадрат

\[
x_\beta'^2=\frac{2J_x}{\beta_x}\bigl(\sin\psi_x+\alpha_x\cos\psi_x\bigr)^2
=\frac{2J_x}{\beta_x}\bigl(\sin^2\psi_x+2\alpha_x\sin\psi_x\cos\psi_x+\alpha_x^2\cos^2\psi_x\bigr).
\]

Средние \(\langle\sin^2\rangle_\psi=\langle\cos^2\rangle_\psi=1/2\), \(\langle\sin\cos\rangle_\psi=0\):

\[
\langle x_\beta'^2\rangle_\psi
=\frac{2J_x}{\beta_x}\cdot\frac{1+\alpha_x^2}{2}
=J_x\frac{1+\alpha_x^2}{\beta_x}
=J_x\gamma_x.
\]

То же для \(y\): \(\langle y_\beta'^2\rangle_\psi=J_y\gamma_y\).

Разложение \(x=x_\beta+x_{\mathrm{nl}}\) даёт \(x'=x_\beta'+x_{\mathrm{nl}}'\). Среднее \(\langle x_\beta'\rangle_\psi=0\), потому что \(\langle\sin\psi_x\rangle_\psi=\langle\cos\psi_x\rangle_\psi=0\), а \(x_{\mathrm{nl}}'\) от быстрой фазы не зависит. Перекрёстный член \(2\langle x_\beta'x_{\mathrm{nl}}'\rangle_\psi\) равен нулю. Квадрат \(x_{\mathrm{nl}}'^2\) квадратичен по сдвигу, а сам сдвиг квадратичен по амплитуде, то есть линеен по \(J\); его квадрат есть \(O(J^2)\) и в первом порядке по \(J\) отбрасывается. Поэтому

\[
\Bigl\langle\frac{x'^2}{2}\Bigr\rangle_\psi=\frac{J_x\gamma_x}{2},
\qquad
\Bigl\langle\frac{y'^2}{2}\Bigr\rangle_\psi=\frac{J_y\gamma_y}{2}.
\]

В члене \(x/\rho\) среднее \(\langle x_\beta\rangle_\psi=0\), остаётся \(\langle x_{\mathrm{nl}}\rangle_\psi\). Интеграл по кольцу и деление на \(C_s\):

\[
\Lambda=\frac{1}{C_s}\Biggl[
\frac{J_x}{2}\oint\gamma_x\,\mathrm{d}s
+\frac{J_y}{2}\oint\gamma_y\,\mathrm{d}s
+\oint\frac{\langle x_{\mathrm{nl}}\rangle_\psi}{\rho}\,\mathrm{d}s
\Biggr].
\]

\(J_{x,y}\) вынесены: на одном обороте действие линейного движения постоянно.

## 4. Постоянные бета-функции

При \(\beta_x=\mathrm{const}\) выполнено \(\alpha_x=-\beta_x'/2=0\) и \(\gamma_x=1/\beta_x\). Набег фазы за оборот

\[
\oint\frac{\mathrm{d}s}{\beta_x}=2\pi\nu_x
\]

и есть определение \(\nu_x\) (при переменной \(\beta_x\) то же равенство верно для \(\oint\mathrm{d}s/\beta_x\), но \(\gamma_x\neq 1/\beta_x\)). Без \(x_{\mathrm{nl}}\)

\begin{align*}
\Lambda
&=\frac{1}{C_s}\Bigl(\frac{J_x}{2}\cdot 2\pi\nu_x+\frac{J_y}{2}\cdot 2\pi\nu_y\Bigr)
=\frac{\pi}{C_s}\bigl(J_x\nu_x+J_y\nu_y\bigr).
\end{align*}

## 5. Хроматичность

Производящая функция однооборотного отображения \(G(J,\delta)\) в переменных действие--фаза и \((\text{продольная координата},\delta)\) устроена так, что набег бетатронной фазы и избыток пути --- смешанные производные одной \(G\):

\[
\Delta\psi=\frac{\partial G}{\partial J}=2\pi\nu,
\qquad
\Delta L=-\frac{\partial G}{\partial\delta}.
\]

Равенство смешанных производных

\[
\frac{\partial\Delta L}{\partial J}
=-\frac{\partial^2 G}{\partial\delta\,\partial J}
=-\frac{\partial\Delta\psi}{\partial\delta}
=-2\pi\frac{\partial\nu}{\partial\delta}
=-2\pi\xi.
\]

Деление на \(C_s\) и определение \(\Lambda=\langle\Delta L_\beta\rangle_\psi/C_s\) дают

\[
\frac{\partial\Lambda}{\partial J_{x,y}}=-\frac{2\pi}{C_s}\,\frac{\partial\nu_{x,y}}{\partial\delta}.
\]

В первом порядке по \(J\) хроматичность берётся при \(J=0\), и интегрирование по действию (аддитивная по \(J\) часть пути, не зависящая от амплитуды, уже включена в \(\alpha_0,\alpha_1\)) даёт

\[
\Lambda=-\frac{2\pi}{C_s}\bigl(J_x\xi_x+J_y\xi_y\bigr).
\]

Секступоли входят в полные \(\xi_{x,y}=\mathrm{d}\nu_{x,y}/\mathrm{d}\delta\). Сдвиг \(x_{\mathrm{nl}}\) меняет и длину, и набег фазы; симплектичность связывает эти две производные, поэтому отдельный интеграл \(\oint\langle x_{\mathrm{nl}}\rangle/\rho\,\mathrm{d}s\) в формулу через \(\xi\) писать не нужно.

Проверка знаком. В чистом фокусирующем канале с постоянной силой жёсткость квадруполя \(K\propto 1/p\propto 1/(1+\delta)\), а \(\nu\propto\sqrt{K}\), поэтому

\[
\nu(\delta)=\frac{\nu}{\sqrt{1+\delta}}=\nu\Bigl(1-\frac{\delta}{2}+O(\delta^2)\Bigr),
\qquad
\xi=\frac{\mathrm{d}\nu}{\mathrm{d}\delta}=-\frac{\nu}{2}.
\]

Подстановка в формулу через \(\xi\):

\[
\Lambda=-\frac{2\pi}{C_s}J\Bigl(-\frac{\nu}{2}\Bigr)=\frac{\pi J\nu}{C_s},
\]

что совпадает с п. 4. Знак фиксирован этим совпадением: при \(\xi<0\) путь удлиняется, как того требует положительный член \(x'^2/2\).

## 6. Вход \(\Lambda\) в систему

Ряд главы 1 есть запаздывание орбиты \(J=0\). Добавка среднего бетатронного запаздывания при \(w=0\) есть \(\Lambda\), поэтому

\[
\tau=\eta_0w+\tilde\eta_1w^2+\Lambda.
\]

Эквивалентность с рядом по \(\delta\) переносится из главы 1 без изменений: там показано

\[
\eta_0\delta+\eta_1\delta^2=\eta_0w+\tilde\eta_1w^2
\]

до отброшенного \(O(w^3)\), а \(\Lambda\) добавляется к обеим записям одним слагаемым. Зависимость действий и бета-функций от \(w\) дала бы \(\Lambda(w)-\Lambda(0)=O(wJ)\); это следующий порядок, и \(\Lambda\) оставляется постоянной. \(J_x,J_y\) за синхротронный период не меняются: излучение и медленный обмен с продольным движением в эту главу не входят.

При \(\sin\varphi\approx\varphi\) первое уравнение \eqref{eq:smooth} с этим \(\tau\) и второе уравнение \eqref{eq:system} дают систему главы.

## 7. Уравнение для \(w\)

Из второго уравнения системы \(\ddot w=K\cos\varphi_s\,\mathrm{d}\varphi/\mathrm{d}t\). Первое даёт \(\mathrm{d}\varphi/\mathrm{d}t=\omega_{\mathrm{rf}}\tau\), поэтому

\[
\ddot w=(K\omega_{\mathrm{rf}}\cos\varphi_s)\,(\eta_0w+\tilde\eta_1w^2+\Lambda).
\]

Подстановка \(K\omega_{\mathrm{rf}}\cos\varphi_s=-\omega_s^2/\eta_0\) из \eqref{eq:ws}:

\begin{align*}
\ddot w
&=\Bigl(-\frac{\omega_s^2}{\eta_0}\Bigr)(\eta_0w+\tilde\eta_1w^2+\Lambda)\\
&=-\omega_s^2w-\frac{\omega_s^2\tilde\eta_1}{\eta_0}w^2-\frac{\omega_s^2\Lambda}{\eta_0}.
\end{align*}

Перенос линейного члена:

\[
\ddot w+\omega_s^2w=c_0+c_2w^2,
\qquad
c_0=-\frac{\omega_s^2\Lambda}{\eta_0},
\qquad
c_2=-\frac{\omega_s^2\tilde\eta_1}{\eta_0}.
\]

## 8. Метод Боголюбова--Митропольского при \(c_0\neq 0\)

Решается \(\ddot x+\omega^2 x=\varepsilon(c_0+c_2 x^2)\). Анзац и уравнение порядка \(\varepsilon^k\) те же, что в приложении. Множитель 2 при \(A_k\) и \(B_k\) выведен в записке главы 1 и здесь не повторяется.

### 8.1. Первый порядок

\[
f_1=c_0+c_2a^2\cos^2\psi
=c_0+\frac{c_2a^2}{2}+\frac{c_2a^2}{2}\cos 2\psi,
\]

поскольку \(\cos^2\psi=(1+\cos 2\psi)/2\). Проекций на \(\sin\psi\) и \(\cos\psi\) нет, поэтому \(A_1=B_1=0\).

Постоянная: \(\omega^2 U_0=c_0+c_2a^2/2\),

\[
U_0=\frac{c_0+c_2a^2/2}{\omega^2}.
\]

Гармоника \(\cos 2\psi\): множитель оператора \(1-4=-3\),

\[
\omega^2(-3)U_2=\frac{c_2a^2}{2},
\qquad
U_2=-\frac{c_2a^2}{6\omega^2}.
\]

\[
u_1=U_0+U_2\cos 2\psi
=\frac{c_0+c_2a^2/2}{\omega^2}-\frac{c_2a^2}{6\omega^2}\cos 2\psi.
\]

При \(c_0=0\) это формула главы 1.

### 8.2. Второй порядок

Правая часть \(\varepsilon(c_0+c_2 x^2)\) при \(x=a\cos\psi+\varepsilon u_1+O(\varepsilon^2)\) в порядке \(\varepsilon^2\) есть \(2c_2 a\cos\psi\, u_1\). Постоянная \(c_0\) умножена на \(\varepsilon\) один раз и вклада \(\varepsilon^2\) не даёт. Подстановка \(u_1\):

\begin{align*}
f_2
&=2c_2a\cos\psi\,(U_0+U_2\cos 2\psi)\\
&=2c_2a\Bigl(U_0\cos\psi+U_2\cdot\frac{\cos 3\psi+\cos\psi}{2}\Bigr)\\
&=2c_2a\Bigl(U_0+\frac{U_2}{2}\Bigr)\cos\psi+c_2a U_2\cos 3\psi.
\end{align*}

Синуса нет, проекция на \(\sin\psi\) даёт \(A_2=0\). Сумма

\begin{align*}
U_0+\frac{U_2}{2}
&=\frac{c_0}{\omega^2}+\frac{c_2a^2}{2\omega^2}+\frac12\Bigl(-\frac{c_2a^2}{6\omega^2}\Bigr)\\
&=\frac{c_0}{\omega^2}+\frac{c_2a^2}{\omega^2}\Bigl(\frac{6}{12}-\frac{1}{12}\Bigr)
=\frac{c_0}{\omega^2}+\frac{5}{12}\frac{c_2a^2}{\omega^2}.
\end{align*}

Слагаемое \(1/2\) --- постоянная часть \(u_1\), порождённая средним \(x^2\) и силой \(c_0\). Слагаемое \(-1/12\) --- вторая гармоника. Проекция на \(\cos\psi\): левая часть уравнения порядка равна \(-2a\omega B_2\), правая равна \(2c_2a(U_0+U_2/2)\),

\begin{align*}
B_2
&=-\frac{c_2}{\omega}\Bigl(U_0+\frac{U_2}{2}\Bigr)
=-\frac{c_2c_0}{\omega^3}-\frac{5}{12}\frac{c_2^2 a^2}{\omega^3}.
\end{align*}

Коэффициент при \(\cos 3\psi\) равен \(c_2 a U_2=-c_2^2 a^3/(6\omega^2)\) и от \(c_0\) не зависит. Множитель оператора \(1-9=-8\):

\[
\omega^2(-8)U_3=-\frac{c_2^2 a^3}{6\omega^2},
\qquad
U_3=\frac{c_2^2 a^3}{48\omega^4}.
\]

Эта гармоника порядка \(c_2^2 a^3\) и не удерживается. При \(\varepsilon=1\)

\begin{align*}
x&=\frac{c_0+c_2a^2/2}{\omega^2}+a\cos\psi-\frac{c_2a^2}{6\omega^2}\cos 2\psi,\\
\Omega&=\omega+B_2
=\omega-\frac{c_2c_0}{\omega^3}-\frac{5}{12}\frac{c_2^2 a^2}{\omega^3}.
\end{align*}

### 8.3. Подстановка коэффициентов главы 2

\[
\frac{c_0}{\omega_s^2}=-\frac{\Lambda}{\eta_0},
\qquad
\frac{c_2}{\omega_s^2}=-\frac{\tilde\eta_1}{\eta_0},
\qquad
a=A.
\]

Постоянный член:

\[
\frac{c_0}{\omega_s^2}+\frac{c_2}{\omega_s^2}\frac{A^2}{2}
=-\frac{\Lambda}{\eta_0}-\frac{\tilde\eta_1 A^2}{2\eta_0}
=-\frac{1}{\eta_0}\Bigl(\frac{\tilde\eta_1 A^2}{2}+\Lambda\Bigr)
=w_{\mathrm{eq}}.
\]

Коэффициент при \(\cos 2\psi\):

\[
-\frac{c_2 A^2}{6\omega_s^2}
=-\frac{A^2}{6}\Bigl(-\frac{\tilde\eta_1}{\eta_0}\Bigr)
=\frac{\tilde\eta_1 A^2}{6\eta_0}.
\]

Частота. Произведение

\[
\frac{c_2}{\omega_s^2}\cdot\frac{c_0}{\omega_s^2}
=\Bigl(-\frac{\tilde\eta_1}{\eta_0}\Bigr)\Bigl(-\frac{\Lambda}{\eta_0}\Bigr)
=\frac{\tilde\eta_1\Lambda}{\eta_0^2},
\]

\[
\frac{B_2}{\omega_s}
=-\frac{c_2 c_0}{\omega_s^4}-\frac{5}{12}\frac{c_2^2 A^2}{\omega_s^4}
=-\frac{\tilde\eta_1\Lambda}{\eta_0^2}-\frac{5}{12}\Bigl(\frac{\tilde\eta_1 A}{\eta_0}\Bigr)^2.
\]

\[
\frac{\Omega}{\omega_s}
=1-\frac{\tilde\eta_1\Lambda}{\eta_0^2}-\frac{5}{12}\Bigl(\frac{\tilde\eta_1 A}{\eta_0}\Bigr)^2.
\]

При \(c_0=0\) член \(\tilde\eta_1\Lambda/\eta_0^2\) выпадает, и остаются формулы главы 1.

Тот же перекрёстный член частоты виден из наклона \(\tau\) в сдвинутой точке. При \(A=0\) и отбрасывании \(\tilde\eta_1 w^2\) корень \(\tau=0\) есть \(w=-\Lambda/\eta_0\). Производная

\[
\frac{\mathrm{d}\tau}{\mathrm{d}w}=\eta_0+2\tilde\eta_1 w,
\qquad
\Bigl(\frac{\mathrm{d}\tau}{\mathrm{d}w}\Bigr)_{w=-\Lambda/\eta_0}
=\eta_0\Bigl(1-\frac{2\tilde\eta_1\Lambda}{\eta_0^2}\Bigr).
\]

Линейная частота пропорциональна корню из этого наклона, делённого на \(\eta_0\):

\[
\frac{\Omega}{\omega_s}
=\sqrt{1-\frac{2\tilde\eta_1\Lambda}{\eta_0^2}}
=1-\frac{\tilde\eta_1\Lambda}{\eta_0^2}+O(\Lambda^2).
\]

Двойка в наклоне становится единицей под корнем. Это совпадает с членом \(-c_2 c_0/\omega^3\).

## 9. Перекрёстные члены в \(\varphi\)

\[
w=w_{\mathrm{eq}}+A\cos\psi+D\cos 2\psi,
\qquad
D=\frac{\tilde\eta_1 A^2}{6\eta_0},
\qquad
w_{\mathrm{eq}}=-\frac{1}{\eta_0}\Bigl(\frac{\tilde\eta_1 A^2}{2}+\Lambda\Bigr).
\]

Квадрат до удерживаемого порядка. Произведения \(w_{\mathrm{eq}}^2\), \(D^2\), \(2w_{\mathrm{eq}}D\) и \(2AD\cos\psi\cos 2\psi\) не удерживаются: первые три дают \(\Lambda^2\), \(A^4\) или \(\Lambda A^2\), последнее --- \(A^3\) в первой гармонике и \(A^3\) в третьей. Остаётся

\[
w^2=\frac{A^2}{2}(1+\cos 2\psi)+2w_{\mathrm{eq}}A\cos\psi.
\]

Подстановка в \(\tau=\eta_0 w+\tilde\eta_1 w^2+\Lambda\):

\begin{align*}
\tau
&=\eta_0(w_{\mathrm{eq}}+A\cos\psi+D\cos 2\psi)
+\tilde\eta_1\Bigl(\frac{A^2}{2}+\frac{A^2}{2}\cos 2\psi+2w_{\mathrm{eq}}A\cos\psi\Bigr)
+\Lambda.
\end{align*}

Постоянный член:

\[
\eta_0 w_{\mathrm{eq}}+\tilde\eta_1\frac{A^2}{2}+\Lambda
=\eta_0\Biggl(-\frac{1}{\eta_0}\Bigl(\frac{\tilde\eta_1 A^2}{2}+\Lambda\Bigr)\Biggr)
+\frac{\tilde\eta_1 A^2}{2}+\Lambda
=0.
\]

Коэффициент при \(\cos\psi\):

\[
\eta_0 A+2\tilde\eta_1 A w_{\mathrm{eq}}.
\]

Полный \(w_{\mathrm{eq}}\) содержит \(-\tilde\eta_1 A^2/(2\eta_0)\). Его вклад

\[
2\tilde\eta_1 A\Bigl(-\frac{\tilde\eta_1 A^2}{2\eta_0}\Bigr)
=-\frac{\tilde\eta_1^2 A^3}{\eta_0}
\]

есть \(O(A^3)\) и отбрасывается, как в главе 1. В перекрёстном члене остаётся \(w_{\mathrm{eq}}=-\Lambda/\eta_0\):

\[
\eta_0 A+2\tilde\eta_1 A\Bigl(-\frac{\Lambda}{\eta_0}\Bigr)
=\eta_0 A\Bigl(1-\frac{2\tilde\eta_1\Lambda}{\eta_0^2}\Bigr).
\]

Коэффициент при \(\cos 2\psi\):

\[
\eta_0 D+\tilde\eta_1\frac{A^2}{2}
=\tilde\eta_1 A^2\Bigl(\frac16+\frac12\Bigr)
=\frac{2}{3}\tilde\eta_1 A^2.
\]

Члена \(\Lambda\) здесь нет. Итог

\[
\tau
=\eta_0 A\Bigl(1-\frac{2\tilde\eta_1\Lambda}{\eta_0^2}\Bigr)\cos\psi
+\frac{2}{3}\tilde\eta_1 A^2\cos 2\psi.
\]

Интегрирование \(\mathrm{d}\varphi/\mathrm{d}t=\omega_{\mathrm{rf}}\tau\) при \(\mathrm{d}t=\mathrm{d}\psi/\Omega\):

\begin{align*}
\int\cos\psi\,\mathrm{d}\psi&=\sin\psi,\\
\int\frac{2}{3}\tilde\eta_1 A^2\cos 2\psi\,\mathrm{d}\psi
&=\frac{2}{3}\tilde\eta_1 A^2\cdot\frac12\sin 2\psi
=\frac{\tilde\eta_1 A^2}{3}\sin 2\psi.
\end{align*}

Постоянная интегрирования равна нулю: \(\varphi\) отсчитывается от устойчивой точки, правая часть нечётна по \(\psi\). Поэтому

\[
\varphi
=\frac{\omega_{\mathrm{rf}}}{\Omega}\Biggl[
\eta_0 A\Bigl(1-\frac{2\tilde\eta_1\Lambda}{\eta_0^2}\Bigr)\sin\psi
+\frac{\tilde\eta_1 A^2}{3}\sin 2\psi
\Biggr].
\]

При \(\Lambda=0\) множитель первой гармоники равен \(\eta_0 A\), и формула совпадает с \eqref{eq:phi}.

Порядок отброшенного куска второй гармоники. Произведение \(2w_{\mathrm{eq}}D\) в \(w^2\) дало бы в \(\tau\) добавку \(\tilde\eta_1\cdot 2w_{\mathrm{eq}}D\). При \(w_{\mathrm{eq}}=-\Lambda/\eta_0\) её модуль есть \(|\tilde\eta_1\Lambda/\eta_0^2|\) от уже удержанной второй гармоники \(\tilde\eta_1 A^2\), то есть лишний малый множитель сверх порядка самой второй гармоники. В первой гармонике удержанный перекрёстный член имеет тот же малый множитель относительно ведущего \(\eta_0 A\), а не относительно члена, который сам уже мал. Поэтому первая гармоника поправляется, вторая --- нет.

## 10. Среднее

Из \eqref{eq:tau-avg} и ряда по \(\delta\)

\[
\eta_0\delta_{\mathrm{eq}}+\eta_1\langle\delta^2\rangle+\Lambda=0,
\qquad
\delta_{\mathrm{eq}}=-\frac{\eta_1\langle\delta^2\rangle+\Lambda}{\eta_0}.
\]

На решении ведущая часть \(w-\,w_{\mathrm{eq}}\) есть \(A\cos\psi\). В том же порядке, что и в главе 1, \(\langle\cos^2\psi\rangle=1/2\), а \(w_{\mathrm{eq}}^2\) и \(D^2\) не удерживаются, так что \(\langle w^2\rangle=A^2/2\). Из связи \(w\) и \(\delta\) до второго порядка \(\langle\delta^2\rangle=\langle w^2\rangle\) на этом порядке:

\[
\langle\delta^2\rangle=\frac{A^2}{2},
\qquad
\delta_{\mathrm{eq}}=-\frac{1}{\eta_0}\Bigl(\frac{\eta_1 A^2}{2}+\Lambda\Bigr).
\]

Условие \(|\Lambda/(\eta_0 A)|\ll 1\) как раз делает \((\Lambda/\eta_0)^2\) малым рядом с \(A^2/2\): их отношение есть \(2(\Lambda/(\eta_0 A))^2\).

Из \eqref{eq:w-series}

\[
w_{\mathrm{eq}}
=\delta_{\mathrm{eq}}+\frac{\langle\delta^2\rangle}{2\gamma_s^2}
=\delta_{\mathrm{eq}}+\frac{A^2}{4\gamma_s^2}.
\]

Подстановка \(\delta_{\mathrm{eq}}\) и \(\tilde\eta_1=\eta_1-\eta_0/(2\gamma_s^2)\):

\begin{align*}
\delta_{\mathrm{eq}}+\frac{A^2}{4\gamma_s^2}
&=-\frac{\eta_1 A^2}{2\eta_0}-\frac{\Lambda}{\eta_0}+\frac{A^2}{4\gamma_s^2}\\
&=-\frac{1}{\eta_0}\Biggl(\Bigl(\eta_1-\frac{\eta_0}{2\gamma_s^2}\Bigr)\frac{A^2}{2}+\Lambda\Biggr)\\
&=-\frac{1}{\eta_0}\Bigl(\frac{\tilde\eta_1 A^2}{2}+\Lambda\Bigr).
\end{align*}

Это постоянный член решения. Сдвиг энергии \(\Delta E_{\mathrm{eq}}=\beta_s^2 E_s w_{\mathrm{eq}}\) есть определение \(w\), усреднённое по периоду.

При \(A\to 0\) ряд \(\tau=\eta_0 w+\Lambda\) (квадрат \(w\) на этом шаге уже отброшен) обращается в нуль при \(w=-\Lambda/\eta_0\). Второе уравнение системы при \(\varphi=0\) даёт \(\dot w=0\). Точка \((\varphi,w)=(0,-\Lambda/\eta_0)\) неподвижна, и решение для \(\varphi\) содержит только \(\sin\psi\) и \(\sin 2\psi\), так что центр фазовых колебаний остаётся в \(\varphi=0\).

Невязка полного квадратичного ряда в этой точке есть \(\tilde\eta_1(\Lambda/\eta_0)^2\). Она порядка \(\Lambda^2\) и в удержанные члены не входит. Точный корень \(\eta_0 w+\tilde\eta_1 w^2+\Lambda=0\) отличается от \(-\Lambda/\eta_0\) на величину порядка \(\tilde\eta_1\Lambda^2/\eta_0^3\).
