# Вывод главы 3 в переменной δ

Все формулы получены заново из главы 1 и подраздела A.1. Старый текст главы 3 (переменная \(w\)) не использовался. Проверки на фиксированных числах — в конце.

Обозначения главы 1: \(\delta=(p-p_s)/p_s\), \(\tau=(T-T_s)/T_s=\eta_0\delta+\eta_1\delta^2\), \(\omega_s^2=K\omega_{\mathrm{rf}}|\eta_0|\), \(K\omega_{\mathrm{rf}}\cos\varphi_s=-\omega_s^2/\eta_0\), \(\delta_m=(\delta_{\max}-\delta_{\min})/2\), \(\delta_{\mathrm{eq}}=\langle\delta\rangle=-(\eta_1/\eta_0)\langle\delta^2\rangle\) из \(\langle\tau\rangle=0\).

## 1. Система

\(\sin\varphi=\varphi-\varphi^3/6+O(\varphi^5)\). Производная удержанного куска по \(\varphi\) равна \(1-\varphi^2/2\), то есть началу ряда \(\cos\varphi\). Множитель зазора \(1-\delta/\gamma_s^2\) сохранён.

\[
\dot\varphi=\omega_{\mathrm{rf}}(\eta_0\delta+\eta_1\delta^2),
\qquad
\dot\delta=K\cos\varphi_s\Bigl(\varphi-\frac{\varphi^3}{6}\Bigr)\Bigl(1-\frac{\delta}{\gamma_s^2}\Bigr).
\]

## 2. Только нелинейность синуса

При \(\eta_1=0\) и \(\Lambda=0\) выполнено \(\tau=\eta_0\delta\), поэтому \(\langle\tau\rangle=0\) даёт \(\delta_{\mathrm{eq}}=0\). Второе уравнение не входит. Кубический член равновесие не сдвигает.

Симметрия при множителе скорости, заменённом единицей: правые части нечётны, замена \((\varphi,\delta)\to(-\varphi,-\delta)\) — симметрия. Тождество \(\cos^3\psi=(3\cos\psi+\cos 3\psi)/4\) не содержит постоянной и чётных гармоник. Нечётная сила не смещает центр.

Множитель \(1-\delta/\gamma_s^2\) симметрию нарушает, но в главе 1 его вклад в постоянную уничтожается членом \((\dot\delta)^2\), а вклад в \(\langle\delta^2\rangle\) при фиксированном полуразмахе есть \(O(\delta_m^4)\). Кубическая поправка есть \(O(\hat\varphi^2\delta_m^2)\). Перекрёстный член отброшен, дальше множитель заменён единицей. Вторая гармоника главы 1 возвращается при сборке с \(\eta_1\).

Что меняет кубический член: частоту (сила слабее линейной) и третью гармонику. Что не меняет: \(\langle\delta\rangle\), чётные гармоники.

Частота не входит в \(\delta_{\mathrm{eq}}\). Среднее по периоду есть среднее по \(\psi\). В \(\delta=\dot\varphi/(\omega_{\mathrm{rf}}\eta_0)\) частота — общий множитель; при фиксированном \(\delta_m\) он сокращается. При \(\eta_1=0\) сдвиг нулевой при любой частоте. При \(\eta_1\neq 0\) кубический член меняет \(\delta_{\mathrm{eq}}\) только через \(\langle\delta^2\rangle\).

Уравнение фазы. \(\dot\varphi=\omega_{\mathrm{rf}}\eta_0\delta\), дифференцирование и второе уравнение:

\[
\ddot\varphi=\omega_{\mathrm{rf}}\eta_0 K\cos\varphi_s\Bigl(\varphi-\frac{\varphi^3}{6}\Bigr)=-\omega_s^2\Bigl(\varphi-\frac{\varphi^3}{6}\Bigr),
\]

потому что \(\omega_{\mathrm{rf}}\eta_0 K\cos\varphi_s=-\omega_s^2\). Это

\[
\ddot\varphi+\omega_s^2\Bigl(\varphi-\frac{\varphi^3}{6}\Bigr)=0.
\]

Полный синус даёт маятник \(\ddot\varphi+\omega_s^2\sin\varphi=0\). Полуразмах \(\hat\varphi=(\varphi_{\max}-\varphi_{\min})/2\).

Сдвиг частоты \(\Omega/\omega_s=1-\hat\varphi^2/16\) есть вклад этого кубического члена. Член \(\tau\) порядка \(\delta^3\) даёт вклад того же формального порядка по амплитуде и при обрыве \(\tau\) на \(\delta^2\) не определён. На \(\delta_{\mathrm{eq}}\) не влияет ни тот, ни другой.

## 3. Метод в приложении A.2

Схема A.1 применяется к \(\varphi\) вместо \(\delta\). Уравнение:

\[
\ddot\varphi+\omega^2\varphi=\varepsilon c_3\varphi^3,
\qquad
\omega=\omega_s,\quad c_3=\omega^2/6.
\]

На \(\varphi=a\cos\psi\):

\[
\cos^3\psi=\frac{3\cos\psi+\cos 3\psi}{4},
\]

из \(\cos 3\psi=4\cos^3\psi-3\cos\psi\). Правая часть первого порядка

\[
f_1=\frac{3c_3 a^3}{4}\cos\psi+\frac{c_3 a^3}{4}\cos 3\psi.
\]

Синуса нет, в левой части схемы коэффициент при \(\sin\psi\) равен \(-2\omega A_1\), значит \(A_1=0\).

При \(m=1\) оператор \(\partial_{\psi\psi}+1\) аннулирует \(\cos\psi\), первая гармоника в \(u_1\) запрещена, коэффициент при \(\cos\psi\) слева есть \(-2a\omega B_1\):

\[
B_1=-\frac{3c_3 a^2}{8\omega}.
\]

При \(m=3\) множитель \(1-9=-8\):

\[
-8\omega^2 U_3=\frac{c_3 a^3}{4},
\qquad
U_3=-\frac{c_3 a^3}{32\omega^2}.
\]

Подстановка \(c_3=\omega^2/6\):

\[
B_1=-\frac{\omega a^2}{16},
\qquad
U_3=-\frac{a^3}{192}.
\]

Проверка: \(3(\omega^2/6)a^2/(8\omega)=\omega a^2/16\), \((\omega^2/6)a^3/(32\omega^2)=a^3/192\).

\[
\varphi=a\cos\psi-\frac{a^3}{192}\cos 3\psi,
\qquad
\frac{\Omega}{\omega}=1-\frac{a^2}{16}.
\]

\(a\) — коэффициент первой гармоники. Главе нужен полуразмах. При \(\psi=0\) оба косинуса равны 1, производная \(-a\sin\psi+(a^3/64)\sin 3\psi\) имеет нули только при \(\sin\psi=0\), если \(a^2\ll 1\) (отношение \(\sin 3\psi/\sin\psi=3-4\sin^2\psi\) ограничено). Поэтому

\[
\hat\varphi=a\Bigl(1-\frac{a^2}{192}\Bigr),
\qquad
a=\hat\varphi+\frac{\hat\varphi^3}{192}+O(\hat\varphi^5).
\]

Подстановка:

\[
\varphi=\hat\varphi\cos\psi+\frac{\hat\varphi^3}{192}(\cos\psi-\cos 3\psi),
\qquad
\frac{\Omega}{\omega_s}=1-\frac{\hat\varphi^2}{16}.
\]

В частоте \(a^2=\hat\varphi^2+O(\hat\varphi^4)\).

Сдвиг \(\psi\to\psi-\pi/2\): \(\cos(\psi-\pi/2)=\sin\psi\), \(\cos(3\psi-3\pi/2)=-\sin 3\psi\), потому что \(\cos(\beta-\pi/2)=\sin\beta\) при \(\beta=3\psi-\pi\) и \(\sin(3\psi-\pi)=-\sin 3\psi\). Множитель \(\operatorname{sgn}\eta_0\) уравнивает ведущий член с главой 1 при \(\hat\varphi\to 0\):

\[
\varphi=\operatorname{sgn}\eta_0\Biggl[\hat\varphi\sin\psi+\frac{\hat\varphi^3}{192}(\sin\psi+\sin 3\psi)\Biggr].
\]

В максимуме \(\psi=\pi/2\) выполнено \(\sin 3\psi=-1\), два члена \(\hat\varphi^3/192\) вычитаются, и пик равен \(\hat\varphi\).

## 4. Импульс и связь полуразмахов

\(\delta=\dot\varphi/(\omega_{\mathrm{rf}}\eta_0)\). Дифференцирование по \(\psi\) даёт множители \(\hat\varphi+\hat\varphi^3/192\) и \(3\hat\varphi^3/192=\hat\varphi^3/64\). Отношение \(\operatorname{sgn}\eta_0/\eta_0=1/|\eta_0|\). Умножение на \(\Omega\):

\[
\delta=\frac{\Omega}{\omega_{\mathrm{rf}}|\eta_0|}
\Biggl[\Bigl(\hat\varphi+\frac{\hat\varphi^3}{192}\Bigr)\cos\psi+\frac{\hat\varphi^3}{64}\cos 3\psi\Biggr].
\]

В точках \(\psi=0,\pi\) косинусы равны \(\pm 1\) и складываются. Сумма коэффициентов, делённая на \(\hat\varphi\), есть \(1+\hat\varphi^2/48\), потому что \(1/64=3/192\). Вместе с \(\Omega/\omega_s=1-\hat\varphi^2/16\):

\[
\delta_m=\frac{\omega_s\hat\varphi}{\omega_{\mathrm{rf}}|\eta_0|}
\Bigl(1-\frac{\hat\varphi^2}{16}\Bigr)\Bigl(1+\frac{\hat\varphi^2}{48}\Bigr)
=\frac{\omega_s\hat\varphi}{\omega_{\mathrm{rf}}|\eta_0|}\Bigl(1-\frac{\hat\varphi^2}{24}\Bigr).
\]

\(1/16-1/48=1/24\). Первая скобка — частота, вторая — форма \(\mathrm{d}\varphi/\mathrm{d}\psi\).

Деление на \(\delta_m\) сокращает \(\Omega\). Коэффициент первой гармоники:

\[
\frac{1+\hat\varphi^2/192}{1+\hat\varphi^2/48}=1-\frac{\hat\varphi^2}{64},
\]

потому что \(1/48-1/192=1/64\). Третья: \(\hat\varphi^2/64\). Итак

\[
\delta=\delta_m\Biggl[\Bigl(1-\frac{\hat\varphi^2}{64}\Bigr)\cos\psi+\frac{\hat\varphi^2}{64}\cos 3\psi\Biggr].
\]

Точный маятник. Умножение на \(\dot\varphi\):

\[
\frac{\mathrm{d}}{\mathrm{d}t}\Bigl(\frac{\dot\varphi^2}{2}-\omega_s^2\cos\varphi\Bigr)=0,
\]

\[
\dot\varphi^2=2\omega_s^2(\cos\varphi-\cos\hat\varphi)=4\omega_s^2\Bigl(\sin^2\frac{\hat\varphi}{2}-\sin^2\frac{\varphi}{2}\Bigr).
\]

Максимум модуля скорости \(2\omega_s|\sin(\hat\varphi/2)|\) при \(\varphi=0\). Симметрия \(\delta_{\min}=-\delta_{\max}\) даёт

\[
\delta_m=\frac{2\omega_s}{\omega_{\mathrm{rf}}|\eta_0|}\sin\frac{\hat\varphi}{2}.
\]

Ряд \(2\sin(\hat\varphi/2)=\hat\varphi(1-\hat\varphi^2/24)+O(\hat\varphi^5)\) совпадает с кубическим результатом.

## 5. Средний квадрат

Ортогональность \(\cos\psi\) и \(\cos 3\psi\), \(\langle\cos^2\rangle=1/2\), квадрат третьей гармоники есть \(O(\hat\varphi^4)\):

\[
\langle\delta^2\rangle=\frac{\delta_m^2}{2}\Bigl(1-\frac{\hat\varphi^2}{64}\Bigr)^2
=\frac{\delta_m^2}{2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr),
\]

остаток \(O(\hat\varphi^4\delta_m^2)\). Удвоение \(\hat\varphi^2/64\) даёт \(\hat\varphi^2/32\).

Поправка есть только потому, что \(\delta_m\) — полуразмах. В максимуме третья гармоника складывается с первой, поэтому коэффициент первой гармоники равен \(\delta_m(1-\hat\varphi^2/64)\), а не \(\delta_m\). Средний квадрат несёт первая гармоника. Если бы \(\delta_m\) было её коэффициентом, поправка порядка \(\hat\varphi^2\) отсутствовала бы. Частота сократилась при нормировке на \(\delta_m\).

\[
\delta_{\mathrm{eq}}=-\frac{\eta_1\delta_m^2}{2\eta_0}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr).
\]

При \(\eta_1=0\) правая часть нулевая.

## 6. Эллиптическое среднее

Для \(\ddot\varphi+\omega_s^2\sin\varphi=0\) внутри сепаратрисы \(0\le\hat\varphi<\pi\) (неустойчивая точка \(\varphi=\pi\): \(\varphi=\pi+\varepsilon\) даёт \(\ddot\varepsilon-\omega_s^2\varepsilon=0\)).

\(k=\sin(\hat\varphi/2)\), \(\sin(\varphi/2)=k\sin\theta\). Тогда \(\mathrm{d}\varphi=2k\cos\theta\,\mathrm{d}\theta/\sqrt{1-k^2\sin^2\theta}\) и \(\dot\varphi=2\omega_s k\cos\theta\) на четверти, где скорость положительна. Отношение

\[
\frac{\mathrm{d}\varphi}{\dot\varphi}=\frac{\mathrm{d}\theta}{\omega_s\sqrt{1-k^2\sin^2\theta}}.
\]

\(T/4=\mathrm{K}(k)/\omega_s\). На подстановке \(\delta^2/\delta_m^2=\cos^2\theta\). Частота \(\omega_s\) сокращается в среднем:

\[
\langle\delta^2\rangle=\frac{\delta_m^2}{\mathrm{K}}\int_0^{\pi/2}\frac{\cos^2\theta\,\mathrm{d}\theta}{\sqrt{1-k^2\sin^2\theta}}.
\]

Тождество \(\cos^2\theta=1-1/k^2+(1-k^2\sin^2\theta)/k^2\) даёт интеграл \([\mathrm{E}-(1-k^2)\mathrm{K}]/k^2\). Деление на \(\mathrm{K}\):

\[
\langle\delta^2\rangle=\delta_m^2\frac{\mathrm{E}(k)/\mathrm{K}(k)-1+k^2}{k^2}.
\]

Это точное среднее для полного синуса. Ряд \((1-\hat\varphi^2/32)\) есть начало его разложения.

Бином: \((1-z)^{-1/2}=1+z/2+3z^2/8+O(z^3)\), \((1-z)^{1/2}=1-z/2-z^2/8+O(z^3)\), \(z=k^2\sin^2\theta\).

\(\int_0^{\pi/2}\sin^2=\pi/4\), \(\int_0^{\pi/2}\sin^4=3\pi/16\).

\[
\mathrm{K}=\frac{\pi}{2}+\frac{\pi k^2}{8}+\frac{9\pi k^4}{128}+O(k^6)
=\frac{\pi}{2}\Bigl(1+\frac{k^2}{4}+\frac{9k^4}{64}+O(k^6)\Bigr),
\]

\[
\mathrm{E}=\frac{\pi}{2}-\frac{\pi k^2}{8}-\frac{3\pi k^4}{128}+O(k^6)
=\frac{\pi}{2}\Bigl(1-\frac{k^2}{4}-\frac{3k^4}{64}+O(k^6)\Bigr).
\]

Обратная скобка \(\mathrm{K}\): \(1-k^2/4-5k^4/64+O(k^6)\), потому что \((k^2/4)^2=4k^4/64\) и \(-9+4=-5\). Произведение со скобкой \(\mathrm{E}\): при \(k^4\) члены \(-5/64-3/64+4/64=-1/16\).

\[
\frac{\mathrm{E}}{\mathrm{K}}=1-\frac{k^2}{2}-\frac{k^4}{16}+O(k^6),
\]

\[
\frac{\mathrm{E}/\mathrm{K}-1+k^2}{k^2}=\frac12-\frac{k^2}{16}+O(k^4).
\]

\(k^2=\hat\varphi^2/4+O(\hat\varphi^4)\), \(k^2/16=\hat\varphi^2/64+O(\hat\varphi^4)\), отсюда множитель \(1/2-\hat\varphi^2/64=(1/2)(1-\hat\varphi^2/32)\).

## 7. Член \(\eta_1\)

В первом порядке \(\langle\delta^2\rangle\) берётся на решении при \(\eta_1=0\): добавка порядка \(\eta_1\delta_m^2\) (постоянная и \(\cos 2\psi\)) даёт в средний квадрат либо второй порядок по \(\eta_1\), либо нулевое среднее с нечётными гармониками.

Вторая гармоника главы 1 на ведущем члене \(\delta_m\cos\psi\) сохраняется; замена ведущего члена на форму с третьей гармоникой меняет её на долю \(\hat\varphi^2\), произведение отброшено. Постоянная не копируется из главы 1: там \(\langle\delta^2\rangle=\delta_m^2/2\). Здесь она равна новому \(\delta_{\mathrm{eq}}\).

Косинусы имеют нулевую производную при \(\psi=0,\pi\), постоянная и вторая гармоника из полуразмаха выпадают. Сдвиг полуразмаха фазы второй гармоникой есть второй порядок по \(\eta_1\delta_m/\eta_0\).

\[
\delta=\delta_{\mathrm{eq}}
+\delta_m\Biggl[\Bigl(1-\frac{\hat\varphi^2}{64}\Bigr)\cos\psi+\frac{\hat\varphi^2}{64}\cos 3\psi\Biggr]
+\frac{\delta_m^2}{6}\Bigl(\frac{\eta_1}{\eta_0}-\frac{2}{\gamma_s^2}\Bigr)\cos 2\psi.
\]

Поправка частоты порядка \((\eta_1\delta_m/\eta_0)^2\) не фиксируется. Удержан только \(1-\hat\varphi^2/16\).

## 8. С \(\Lambda\)

\(\tau=\eta_0\delta+\eta_1\delta^2+\Lambda\) из главы 2. Отсюда \(\delta_{\mathrm{eq}}=-(\eta_1\langle\delta^2\rangle+\Lambda)/\eta_0\).

При \(\eta_1=0\) переменная \(u=\delta+\Lambda/\eta_0\) удовлетворяет тому же маятнику. Постоянная \(\Lambda/(\eta_0\gamma_s^2)\) в множителе зазора меняет только частоту. \(\langle u^2\rangle\) есть среднее п. 5, \(\delta_{\mathrm{eq}}=-\Lambda/\eta_0\).

\(\Lambda\) порядка \(\delta_m^2\). Квадрат постоянной в \(\langle\delta^2\rangle\) отброшен (условие \(|\Lambda/(\eta_0\delta_m)|\ll 1\)). Подстановка среднего:

\[
\delta_{\mathrm{eq}}
=-\frac{1}{\eta_0}\Biggl[\frac{\eta_1\delta_m^2}{2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr)+\Lambda\Biggr].
\]

В \(\gamma=\gamma_s+\beta_s^2\gamma_s\delta+(\beta_s^2/(2\gamma_s))\delta^2+O(\delta^3)\) член \(\langle\delta^3\rangle\) на этом решении есть \(O(\delta_m^4)\), если \(\Lambda\) порядка \(\delta_m^2\). Множитель \(1-\hat\varphi^2/32\) стоит при каждом \(\langle\delta^2\rangle\), но не при \(\Lambda\):

\[
\gamma_{\mathrm{eff}}
=\gamma_s+\beta_s^2\gamma_s\delta_{\mathrm{eq}}
+\frac{\beta_s^2\delta_m^2}{4\gamma_s}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr),
\]

\[
\Delta E_{\mathrm{eq}}
=\beta_s^2 E_s\Biggl[\delta_{\mathrm{eq}}+\frac{\delta_m^2}{4\gamma_s^2}\Bigl(1-\frac{\hat\varphi^2}{32}\Bigr)\Biggr].
\]

Проверка связи: \((\beta_s^2/(2\gamma_s))\langle\delta^2\rangle\) при \(\langle\delta^2\rangle=\delta_m^2(1-\hat\varphi^2/32)/2\) равно \(\beta_s^2\delta_m^2(1-\hat\varphi^2/32)/(4\gamma_s)\). Деление \((\gamma_{\mathrm{eff}}-\gamma_s)/\gamma_s\) и умножение на \(E_s\) дают энергию: \(\langle\delta^2\rangle/(2\gamma_s^2)=\delta_m^2(1-\hat\varphi^2/32)/(4\gamma_s^2)\).

При \(\Lambda=0\) первая строка есть \(\delta_{\mathrm{eq}}\) без бетатронного члена. При \(\hat\varphi\to 0\) формулы совпадают с главой 2.

Условие обрыва синуса: \(\hat\varphi^2\ll 1\). Остальные условия — из главы 2.

## 9. Численные проверки

Скрипт `review/_check_ch3_delta.py`. Фиксированные числа, без подгонки. `scipy.special.ellipk` и `ellipe` принимают параметр \(m=k^2\).

Дроби: \(1/16-1/48=1/24\), \(1/48-1/192=1/64\), \(2/64=1/32\), \(11/192-3/192=1/24\).

Ряд эллиптического отношения. Величина \((2\langle\delta^2\rangle/\delta_m^2-1)/\hat\varphi^2\) должна стремиться к \(-1/32=-0.03125\).

| \(\hat\varphi\) | отношение | \((2r-1)/\hat\varphi^2\) | разность с \(-1/32\) |
|---|---|---|---|
| 0.4 | 0.497483206510 | −0.03145992 | −2.1·10⁻⁴ |
| 0.2 | 0.499373956365 | −0.03130218 | −5.2·10⁻⁵ |
| 0.1 | 0.499843684865 | −0.03126303 | −1.3·10⁻⁵ |
| 0.05 | 0.499960933430 | −0.03125326 | −3.3·10⁻⁶ |
| 0.025 | 0.499990234121 | −0.03125081 | −8.1·10⁻⁷ |

Остаток падает как \(\hat\varphi^2\), то есть следующий член ряда есть \(O(\hat\varphi^4)\).

Квадратура по \(\theta\) (подынтегральное \(\cos^2\theta/\sqrt{1-k^2\sin^2\theta}\), без вызова `ellipk` в числителе) против замкнутой формулы: при \(\hat\varphi=0.3,1,2,3\) расхождение не больше \(10^{-14}\).

Квадратура по \(\varphi\) с заменой \(\varphi=\hat\varphi\sin\alpha\), тоже без эллиптических функций в ответе: при \(\hat\varphi=0.5,1.2,2.2\) расхождение около \(3\cdot 10^{-11}\).

Интегрирование маятника \(\ddot\varphi+\sin\varphi=0\) (`DOP853`, допуски \(10^{-12}\)), \(\delta=\dot\varphi\), \(\omega_s=\omega_{\mathrm{rf}}=|\eta_0|=1\).

- \(\hat\varphi=0.2\): \(\delta_m\) совпал с \(2\sin(\hat\varphi/2)\) до \(6\cdot 10^{-13}\); \(\langle\delta^2\rangle\) с эллиптической формулой до \(10^{-13}\); частота с \(\pi/(2\mathrm{K})\) совпала.
- \(\hat\varphi=1.2\) и \(2.5\): те же совпадения до \(10^{-12}\)–\(10^{-13}\). Ряд \(1-\hat\varphi^2/32\) при \(2.5\) уже груб, как и должно быть вне области \(\hat\varphi^2\ll 1\).

Ряды \(\mathrm{K}\) и \(\mathrm{E}\) до \(k^4\) и отношение \(1/2-k^2/16\): при \(m=10^{-4},10^{-3},10^{-2}\) остаток отношения есть \(O(m^2)\).

Кубическое уравнение \(\ddot\varphi+(\varphi-\varphi^3/6)=0\), старт из интеграла энергии \(v_0=\hat\varphi\sqrt{1-\hat\varphi^2/12}\).

При \(\hat\varphi=0.15\):

- полуразмах фазы вернулся равным \(0.15\);
- среднее \(\delta\) есть \(-5\cdot 10^{-16}\) (сдвига нет);
- коэффициент \(\sin\psi\): \(0.1500176256\) против \(0.1500175781\);
- коэффициент \(\sin 3\psi\): \(1.763\cdot 10^{-5}\) против \(\hat\varphi^3/192=1.758\cdot 10^{-5}\);
- \(\delta_m\): \(0.149859309\) против \(\hat\varphi(1-\hat\varphi^2/24)=0.149859375\);
- \(\langle\delta^2\rangle\): расхождение \(1.7\cdot 10^{-8}\);
- \((2\langle\delta^2\rangle/\delta_m^2-1)/\hat\varphi^2=-0.031319\) против \(-0.03125\);
- \(\Omega=0.99859259\) против \(1-\hat\varphi^2/16=0.99859375\);
- \(\langle\delta^2\rangle/(c_1^2/2)-1=1.2\cdot 10^{-7}\): если делить на коэффициент первой гармоники, а не на полуразмах, поправки порядка \(\hat\varphi^2\) нет.

При \(\hat\varphi=0.3\) расхождения больше примерно в \((0.3/0.15)^4=16\) раз для членов \(O(\hat\varphi^4)\) и \(O(\hat\varphi^5)\) у гармоник, как ждёт оборванный ряд.

## 10. Сомнения

1. Сдвиг \(\Omega/\omega_s=1-\hat\varphi^2/16\) — только вклад \(\varphi^3\). Член \(\delta^3\) в \(\tau\) даёт другой вклад того же формального порядка и здесь не вычислен. На \(\delta_{\mathrm{eq}}\) оба не влияют. В текст это внесено одной фразой.
2. Множитель \(1-\delta/\gamma_s^2\) в системе оставлен, но при построении маятника заменён единицей. Его вторая гармоника возвращена из главы 1. Перекрёстные члены порядка \(\hat\varphi^2\delta_m^2/\gamma_s^2\) и \(\eta_1\hat\varphi^2\delta_m^2\) отброшены.
3. В главе 4 ещё стоит ссылка на старую метку `eq:w2-ell`. Новая метка точной формулы — `eq:d2-ell`. Метки `eq:deq-ph` и `eq:deq-ph-L` сохранены.
