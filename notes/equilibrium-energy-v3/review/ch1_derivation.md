# Вывод формул главы 1 и приложения A

Обозначения совпадают с запиской. Индекс $s$ опущен у $\beta$, $\gamma$, $E$, $p$, $v$ только в промежуточных строках, где величина явно берётся на синхронной частице; в окончательных формулах индекс восстановлен.

## 1. Обозначения

### 1.1. Производные энергии

$$E=\sqrt{p^2 c^2+m^2 c^4},\qquad E^2=p^2 c^2+m^2 c^4.$$

Дифференциал: $2E\,dE=2p c^2\,dp$, поэтому

$$\frac{dE}{dp}=\frac{p c^2}{E}.$$

Скорость $v=p c^2/E$, значит $dE/dp=v$.

Вторая производная:

$$\frac{d^2E}{dp^2}
=\partial_p\Bigl(\frac{p c^2}{E}\Bigr)
=\frac{c^2}{E}+p c^2\cdot\Bigl(-\frac{1}{E^2}\Bigr)\frac{dE}{dp}
=\frac{c^2}{E}-\frac{v^2}{E}
=\frac{c^2-v^2}{E}.$$

Так как $c^2-v^2=c^2(1-\beta^2)=c^2/\gamma^2$,

$$\frac{d^2E}{dp^2}=\frac{c^2}{\gamma^2 E}.$$

### 1.2. Ряд для $w$

Тейлор в точке $p_s$:

$$E=E_s+(dE/dp)_s(p-p_s)+\frac12(d^2E/dp^2)_s(p-p_s)^2+O((p-p_s)^3).$$

Подстановка $(dE/dp)_s=v_s$, $(d^2E/dp^2)_s=c^2/(\gamma_s^2 E_s)$, $p-p_s=\delta p_s$:

$$E-E_s=v_s p_s\,\delta+\frac12\frac{c^2}{\gamma_s^2 E_s}p_s^2\delta^2+O(\delta^3).$$

Связь $p_s=E_s v_s/c^2$ даёт

$$v_s p_s=\beta_s c\cdot\frac{E_s\beta_s c}{c^2}=\beta_s^2 E_s,$$

$$p_s^2=\frac{E_s^2 v_s^2}{c^4}=\frac{E_s^2\beta_s^2}{c^2}.$$

Деление на $\beta_s^2 E_s$:

$$\frac{v_s p_s\,\delta}{\beta_s^2 E_s}=\delta.$$

Второй член:

$$\frac{1}{\beta_s^2 E_s}\cdot\frac12\cdot\frac{c^2}{\gamma_s^2 E_s}\cdot\frac{E_s^2\beta_s^2}{c^2}\,\delta^2
=\frac{1}{2\gamma_s^2}\delta^2.$$

($\beta_s^2$ и $c^2$ сокращаются, $E_s^2/E_s^2=1$.)

Итог:

$$w=\frac{E-E_s}{\beta_s^2 E_s}=\delta+\frac{\delta^2}{2\gamma_s^2}+O(\delta^3).$$

### 1.3. Фаза и запаздывание

Опоздание на время $\Delta t=T-T_s=\tau T_s$ сдвигает фазу ВЧ на

$$\omega_{\mathrm{rf}}\Delta t=\omega_{\mathrm{rf}}\tau T_s.$$

По определению $\omega_{\mathrm{rf}}=2\pi h/T_s$, поэтому

$$\omega_{\mathrm{rf}}\tau T_s=2\pi h\tau.$$

При $\tau>0$ частица опаздывает, и $\varphi$ растёт. Это даёт отображение фазы в п. 3.

## 2. Запаздывание за оборот

### 2.1. Точное отношение периодов

Скорость в магнитном поле постоянна вдоль траектории, $T=C/v$, $T_s=C_s/v_s$, поэтому

$$\frac{T}{T_s}=\frac{C}{C_s}\frac{v_s}{v}$$

точно, и $\tau=T/T_s-1$.

### 2.2. От $p$ к $v_s/v$

$$p=\beta\gamma m c,\qquad p_s=\beta_s\gamma_s m c,\qquad \frac{p}{p_s}=1+\delta,$$

значит $\beta\gamma=\beta_s\gamma_s(1+\delta)$. Тогда

$$\beta=\frac{\beta\gamma}{\gamma}=\frac{\beta_s\gamma_s(1+\delta)}{\gamma},$$

$$\frac{v_s}{v}=\frac{\beta_s}{\beta}=\frac{\gamma}{\gamma_s(1+\delta)}.$$

### 2.3. Ряд для $\gamma$

$$(\beta\gamma)^2=\beta_s^2\gamma_s^2(1+\delta)^2,\qquad \gamma^2=1+(\beta\gamma)^2.$$

Тождество $\beta_s^2\gamma_s^2=\gamma_s^2-1$:

$$\gamma^2=1+(\gamma_s^2-1)(1+2\delta+\delta^2)
=1+(\gamma_s^2-1)+2(\gamma_s^2-1)\delta+(\gamma_s^2-1)\delta^2$$

$$=\gamma_s^2+2(\gamma_s^2-1)\delta+(\gamma_s^2-1)\delta^2
=\gamma_s^2\bigl(1+2\beta_s^2\delta+\beta_s^2\delta^2\bigr),$$

поскольку $(\gamma_s^2-1)/\gamma_s^2=\beta_s^2$. Отсюда

$$\frac{\gamma}{\gamma_s}=\sqrt{1+x},\qquad x=2\beta_s^2\delta+\beta_s^2\delta^2.$$

Ряд $\sqrt{1+x}=1+x/2-x^2/8+O(x^3)$, и $x=O(\delta)$, так что $O(x^3)=O(\delta^3)$.

$$x/2=\beta_s^2\delta+(\beta_s^2/2)\delta^2.$$

$$x^2=(2\beta_s^2\delta+\beta_s^2\delta^2)^2=4\beta_s^4\delta^2+O(\delta^3),$$

$$x^2/8=\beta_s^4\delta^2/2+O(\delta^3).$$

$$x/2-x^2/8
=\beta_s^2\delta+\bigl(\beta_s^2/2-\beta_s^4/2\bigr)\delta^2
=\beta_s^2\delta+(\beta_s^2/2)(1-\beta_s^2)\delta^2.$$

$1-\beta_s^2=1/\gamma_s^2$, поэтому

$$\frac{\gamma}{\gamma_s}
=1+\beta_s^2\delta+\frac{\beta_s^2}{2\gamma_s^2}\delta^2+O(\delta^3).$$

### 2.4. Деление на $1+\delta$

$$(1+\delta)^{-1}=1-\delta+\delta^2+O(\delta^3).$$

$$\frac{v_s}{v}
=\Bigl(1+\beta_s^2\delta+\frac{\beta_s^2}{2\gamma_s^2}\delta^2\Bigr)(1-\delta+\delta^2)+O(\delta^3).$$

Раскрытие:

- $1\cdot 1=1$
- $1\cdot(-\delta)=-\delta$
- $1\cdot\delta^2=\delta^2$
- $\beta_s^2\delta\cdot 1=\beta_s^2\delta$
- $\beta_s^2\delta\cdot(-\delta)=-\beta_s^2\delta^2$
- $(\beta_s^2/(2\gamma_s^2))\delta^2\cdot 1=\beta_s^2\delta^2/(2\gamma_s^2)$
- остальные произведения имеют порядок $\delta^3$

Коэффициент при $\delta$: $\beta_s^2-1=-1/\gamma_s^2$.

Коэффициент при $\delta^2$:

$$1-\beta_s^2+\frac{\beta_s^2}{2\gamma_s^2}
=\frac{1}{\gamma_s^2}+\frac{\beta_s^2}{2\gamma_s^2}
=\Bigl(1+\frac{\beta_s^2}{2}\Bigr)\frac{1}{\gamma_s^2}.$$

Итог:

$$\frac{v_s}{v}
=1-\frac{\delta}{\gamma_s^2}
+\Bigl(1+\frac{\beta_s^2}{2}\Bigr)\frac{\delta^2}{\gamma_s^2}
+O(\delta^3).$$

### 2.5. Произведение с длиной орбиты

$$\frac{C}{C_s}=1+\alpha_0\delta+\alpha_1\delta^2+O(\delta^3).$$

$$\frac{T}{T_s}
=(1+\alpha_0\delta+\alpha_1\delta^2)
\Bigl(1-\frac{\delta}{\gamma_s^2}
+\Bigl(1+\frac{\beta_s^2}{2}\Bigr)\frac{\delta^2}{\gamma_s^2}\Bigr)
+O(\delta^3).$$

Обозначим $b_1=-1/\gamma_s^2$, $b_2=(1+\beta_s^2/2)/\gamma_s^2$. Тогда

$$(1+\alpha_0\delta+\alpha_1\delta^2)(1+b_1\delta+b_2\delta^2)
=1+(\alpha_0+b_1)\delta+(\alpha_1+b_2+\alpha_0 b_1)\delta^2+O(\delta^3).$$

$$\alpha_0+b_1=\alpha_0-\frac{1}{\gamma_s^2}=\eta_0.$$

$$\alpha_1+b_2+\alpha_0 b_1
=\alpha_1+\frac{1+\beta_s^2/2}{\gamma_s^2}+\alpha_0\Bigl(-\frac{1}{\gamma_s^2}\Bigr)
=\alpha_1-\frac{\alpha_0}{\gamma_s^2}+\frac{1+\beta_s^2/2}{\gamma_s^2}
=\eta_1.$$

$$\tau=\frac{T}{T_s}-1=\eta_0\delta+\eta_1\delta^2+O(\delta^3).$$

В главе удерживается

$$\tau=\eta_0\delta+\eta_1\delta^2.$$

## 3. Система продольных уравнений

### 3.1. Прирост энергии на зазоре

$$\sin(\varphi_s+\varphi)=\sin\varphi_s\cos\varphi+\cos\varphi_s\sin\varphi.$$

При $\varphi_s=0$: $\sin\varphi_s=0$, $\cos\varphi_s=1$, произведение равно $\sin\varphi$.

При $\varphi_s=\pi$: $\sin\varphi_s=0$, $\cos\varphi_s=-1$, произведение равно $-\sin\varphi$, и $\sin(\pi+\varphi)=-\sin\varphi$.

В обоих случаях

$$\sin(\varphi_s+\varphi)=\cos\varphi_s\sin\varphi,$$

и прирост энергии равен $qV\cos\varphi_s\sin\varphi$.

На синхронной энергии

$$\Bigl(\frac{dE}{d\delta}\Bigr)_s
=\Bigl(\frac{dE}{dp}\Bigr)_s p_s
=v_s p_s
=\beta_s^2 E_s,$$

поэтому в удерживаемом порядке $dE=\beta_s^2 E_s\,d\delta$ и

$$\delta_{n+1}-\delta_n
=\frac{qV\cos\varphi_s\sin\varphi_n}{\beta_s^2 E_s}.$$

Фаза за оборот, п. 1.3:

$$\varphi_{n+1}-\varphi_n=2\pi h\,\tau(\delta_{n+1}).$$

### 3.2. Сглаживание

$$t=n T_s,\qquad
\frac{\varphi_{n+1}-\varphi_n}{T_s}
=\frac{2\pi h}{T_s}\tau
=\omega_{\mathrm{rf}}\tau,$$

$$\frac{\delta_{n+1}-\delta_n}{T_s}
=\frac{qV\cos\varphi_s\sin\varphi}{\beta_s^2 E_s T_s}
=K\cos\varphi_s\sin\varphi,$$

$$K=\frac{qV}{\beta_s^2 E_s T_s}.$$

Сглаженная система:

$$\frac{d\varphi}{dt}=\omega_{\mathrm{rf}}\tau,
\qquad
\frac{d\delta}{dt}=K\cos\varphi_s\sin\varphi.$$

### 3.3. Устойчивая точка и $\omega_s$

Линеаризация $\sin\varphi\approx\varphi$, $\tau\approx\eta_0\delta$:

$$\frac{d\varphi}{dt}=\omega_{\mathrm{rf}}\eta_0\delta,
\qquad
\frac{d\delta}{dt}=K\cos\varphi_s\,\varphi.$$

Дифференцирование второго уравнения и подстановка первого:

$$\ddot\delta
=K\cos\varphi_s\,\frac{d\varphi}{dt}
=K\cos\varphi_s\,\omega_{\mathrm{rf}}\eta_0\delta
=K\omega_{\mathrm{rf}}\eta_0\cos\varphi_s\,\delta.$$

Гармонический режим $\ddot\delta=-\omega_s^2\delta$ требует

$$K\omega_{\mathrm{rf}}\eta_0\cos\varphi_s=-\omega_s^2<0.$$

При $K>0$ и $\omega_{\mathrm{rf}}>0$ нужно $\eta_0\cos\varphi_s<0$, то есть

$$\cos\varphi_s=-\operatorname{sgn}\eta_0.$$

Тогда $\eta_0\cos\varphi_s=\eta_0(-\operatorname{sgn}\eta_0)=-|\eta_0|$ и

$$\omega_s^2
=-K\omega_{\mathrm{rf}}\eta_0\cos\varphi_s
=K\omega_{\mathrm{rf}}|\eta_0|.$$

Дальше

$$K\omega_{\mathrm{rf}}\cos\varphi_s
=K\omega_{\mathrm{rf}}(-\operatorname{sgn}\eta_0)
=-\frac{K\omega_{\mathrm{rf}}|\eta_0|}{\eta_0}
=-\frac{\omega_s^2}{\eta_0}.$$

Проверка знака при $\eta_0<0$: $\operatorname{sgn}\eta_0=-1$, $\cos\varphi_s=1$, $-\omega_s^2/\eta_0=-K\omega_{\mathrm{rf}}|\eta_0|/\eta_0=K\omega_{\mathrm{rf}}$. Совпадает с $K\omega_{\mathrm{rf}}\cos\varphi_s$.

### 3.4. Система главы

$\sin\varphi\approx\varphi$ и $\tau=\eta_0\delta+\eta_1\delta^2$ подставляются в сглаженную систему:

$$\frac{d\varphi}{dt}=\omega_{\mathrm{rf}}(\eta_0\delta+\eta_1\delta^2),
\qquad
\frac{d\delta}{dt}=K\cos\varphi_s\,\varphi.$$

## 4. Решение

### 4.1. Уравнение для $\delta$

Из второго уравнения системы главы $\ddot\delta=K\cos\varphi_s\,d\varphi/dt$. Первое даёт $d\varphi/dt=\omega_{\mathrm{rf}}\tau$, поэтому

$$\ddot\delta=K\cos\varphi_s\,\omega_{\mathrm{rf}}\,\tau
=(K\omega_{\mathrm{rf}}\cos\varphi_s)(\eta_0\delta+\eta_1\delta^2).$$

Подстановка $K\omega_{\mathrm{rf}}\cos\varphi_s=-\omega_s^2/\eta_0$:

$$\ddot\delta
=\Bigl(-\frac{\omega_s^2}{\eta_0}\Bigr)(\eta_0\delta+\eta_1\delta^2)
=-\omega_s^2\delta-\frac{\omega_s^2\eta_1}{\eta_0}\delta^2.$$

Перенос линейного члена:

$$\ddot\delta+\omega_s^2\delta=-\frac{\omega_s^2\eta_1}{\eta_0}\delta^2.$$

### 4.2. Полуразмах

Решение (вывод --- приложение A, подстановка ниже):

$$\delta=\delta_{\mathrm{eq}}+A\cos\psi+D\cos 2\psi,
\qquad
D=\frac{\eta_1 A^2}{6\eta_0},
\qquad
\delta_{\mathrm{eq}}=-\frac{\eta_1 A^2}{2\eta_0},
\qquad
\psi=\Omega t.$$

В точках экстремума $\psi=0$ и $\psi=\pi$ (внутри области применимости $|\eta_1 A/\eta_0|\ll 1$ других нулей $d\delta/d\psi$ нет: $d\delta/d\psi=-\sin\psi\,(A+4D\cos\psi)$ и $|4D/A|=|\eta_1 A/(6\eta_0)|\cdot 4\ll 1$):

$$\cos 0=\cos 2\cdot 0=1,\qquad \cos\pi=-1,\qquad \cos 2\pi=1.$$

$$\delta(0)=\delta_{\mathrm{eq}}+A+D,
\qquad
\delta(\pi)=\delta_{\mathrm{eq}}-A+D.$$

Обе точки сдвинуты на одну и ту же величину $\delta_{\mathrm{eq}}+D$. Разность

$$\delta(0)-\delta(\pi)=2A,
\qquad
\frac{\delta_{\max}-\delta_{\min}}{2}=A.$$

Коэффициент при $\cos\psi$ совпадает с полуразмахом. Частота $\Omega$ входит только в закон $\psi=\Omega t$ и в $\varphi$; постоянный член $\delta_{\mathrm{eq}}$ от $\Omega$ не зависит.

### 4.3. Квадратура фазы

Подробная подстановка $\delta$ в $\tau$ выполнена в п. A.4; здесь проверка дифференцированием готовой формулы.

$$\varphi=\frac{\omega_{\mathrm{rf}}}{\Omega}
\Bigl(\eta_0 A\sin\psi+\frac{\eta_1 A^2}{3}\sin 2\psi\Bigr),
\qquad
\frac{d\psi}{dt}=\Omega.$$

$$\frac{d\varphi}{dt}
=\omega_{\mathrm{rf}}
\Bigl(\eta_0 A\cos\psi+\frac{\eta_1 A^2}{3}\cdot 2\cos 2\psi\Bigr)
=\omega_{\mathrm{rf}}
\Bigl(\eta_0 A\cos\psi+\frac{2}{3}\eta_1 A^2\cos 2\psi\Bigr).$$

Это совпадает с $\omega_{\mathrm{rf}}\tau$ из п. A.4 до порядка $A^2$.

## 5. Усреднение

### 5.1. Среднее запаздывание

На замкнутой траектории $\varphi(t+T_{\mathrm{syn}})=\varphi(t)$, поэтому

$$\Bigl\langle\frac{d\varphi}{dt}\Bigr\rangle
=\frac{1}{T_{\mathrm{syn}}}\int_0^{T_{\mathrm{syn}}}\frac{d\varphi}{dt}\,dt
=\frac{\varphi(T_{\mathrm{syn}})-\varphi(0)}{T_{\mathrm{syn}}}
=0.$$

Из $d\varphi/dt=\omega_{\mathrm{rf}}\tau$ и $\omega_{\mathrm{rf}}\neq 0$ следует $\langle\tau\rangle=0$ точно (до обрезания ряда $\tau$). Средний период равен $T_s$.

Для удержанного ряда $\tau=\eta_0\delta+\eta_1\delta^2$

$$\langle\tau\rangle=\eta_0\langle\delta\rangle+\eta_1\langle\delta^2\rangle=0,$$

$$\langle\delta\rangle=-\frac{\eta_1}{\eta_0}\langle\delta^2\rangle.$$

По определению $\delta_{\mathrm{eq}}:=\langle\delta\rangle$.

### 5.2. Среднее $\delta^2$ на решении

$$\delta=\delta_{\mathrm{eq}}+A\cos\psi+D\cos 2\psi.$$

$$\delta^2
=\delta_{\mathrm{eq}}^2
+A^2\cos^2\psi
+D^2\cos^2 2\psi
+2\delta_{\mathrm{eq}} A\cos\psi
+2\delta_{\mathrm{eq}} D\cos 2\psi
+2AD\cos\psi\cos 2\psi.$$

Средние по периоду $\psi$:

$$\langle\cos\psi\rangle=0,
\qquad
\langle\cos 2\psi\rangle=0,$$

$$\langle\cos\psi\cos 2\psi\rangle
=\Bigl\langle\frac{\cos 3\psi+\cos\psi}{2}\Bigr\rangle=0,$$

$$\langle\cos^2\psi\rangle=\Bigl\langle\frac{1+\cos 2\psi}{2}\Bigr\rangle=\frac12,
\qquad
\langle\cos^2 2\psi\rangle=\frac12.$$

$$\langle\delta^2\rangle
=\delta_{\mathrm{eq}}^2+\frac{A^2}{2}+\frac{D^2}{2}.$$

$\delta_{\mathrm{eq}}=-\eta_1 A^2/(2\eta_0)=O(A^2)$, $D=\eta_1 A^2/(6\eta_0)=O(A^2)$, поэтому $\delta_{\mathrm{eq}}^2=O(A^4)$ и $D^2=O(A^4)$:

$$\langle\delta^2\rangle=\frac{A^2}{2}+O(A^4).$$

$$\delta_{\mathrm{eq}}
=-\frac{\eta_1}{\eta_0}\cdot\frac{A^2}{2}
=-\frac{\eta_1 A^2}{2\eta_0}$$

в удерживаемом порядке. Это постоянный член решения.

### 5.3. Энергия

$$w=\delta+\frac{\delta^2}{2\gamma_s^2}+O(\delta^3),$$

$$\langle w\rangle
=\langle\delta\rangle+\frac{\langle\delta^2\rangle}{2\gamma_s^2}
=\delta_{\mathrm{eq}}+\frac{A^2/2}{2\gamma_s^2}
=-\frac{\eta_1 A^2}{2\eta_0}+\frac{A^2}{4\gamma_s^2}.$$

Общий знаменатель $2\eta_0$:

$$\frac{A^2}{4\gamma_s^2}
=\frac{\eta_0 A^2}{2\eta_0\cdot 2\gamma_s^2}
=\frac{\eta_0}{2\gamma_s^2}\cdot\frac{A^2}{2\eta_0}.$$

$$w_{\mathrm{eq}}
=-\frac{\eta_1 A^2}{2\eta_0}
+\frac{\eta_0}{2\gamma_s^2}\cdot\frac{A^2}{2\eta_0}
=-\Bigl(\eta_1-\frac{\eta_0}{2\gamma_s^2}\Bigr)\frac{A^2}{2\eta_0}.$$

$$\Delta E_{\mathrm{eq}}=\beta_s^2 E_s w_{\mathrm{eq}}=\langle E\rangle-E_s,$$

где второе равенство есть определение $w$, усреднённое по периоду.

### 5.4. Полуамплитуда фазы

Ведущий член $\varphi$ есть $(\omega_{\mathrm{rf}}/\Omega)\eta_0 A\sin\psi$. При $\Omega=\omega_s(1+O((\eta_1 A/\eta_0)^2))$ модуль коэффициента в линейном приближении равен

$$\hat\varphi=\frac{\omega_{\mathrm{rf}}|\eta_0|A}{\omega_s}.$$

Условие $\hat\varphi\ll 1$ есть условие $\sin\varphi\approx\varphi$. Условие $|\eta_1 A/\eta_0|\ll 1$ есть малость квадратичной поправки относительно линейного члена $\eta_0\delta$.

## Приложение A. Метод Боголюбова--Митропольского

Книжный параметр порядка $\varepsilon$ вводится так: решается

$$\ddot x+\omega^2 x=\varepsilon c_2 x^2,$$

а в уравнении главы $\varepsilon=1$ и $c_2=-\omega_s^2\eta_1/\eta_0$. При $\varepsilon=1$ анзац совпадает с записью $x=a\cos\psi+u_1+u_2$.

$$\dot a=\varepsilon A_1+\varepsilon^2 A_2,
\qquad
\dot\psi=\omega+\varepsilon B_1+\varepsilon^2 B_2,$$

$$x=a\cos\psi+\varepsilon u_1(a,\psi)+\varepsilon^2 u_2(a,\psi).$$

$u_k$ не содержат $\cos\psi$ и $\sin\psi$.

### A.1. Откуда множитель 2 при $A_k$ и $B_k$

Ниже $A_1$, $B_1$ ещё не положены нулю; члены выписаны до порядка $\varepsilon$. Пусть $\nu=\dot\psi=\omega+\varepsilon B_1+O(\varepsilon^2)$, $\dot a=\varepsilon A_1+O(\varepsilon^2)$.

$$\frac{d}{dt}(a\cos\psi)=\dot a\cos\psi-a\nu\sin\psi.$$

$$\frac{d^2}{dt^2}(a\cos\psi)
=\ddot a\cos\psi
-2\dot a\,\nu\sin\psi
-a\dot\nu\sin\psi
-a\nu^2\cos\psi.$$

$\ddot a=O(\varepsilon^2)$, $\dot\nu=O(\varepsilon^2)$ (производная $B_1(a)$ умножается на $\dot a=O(\varepsilon)$). В порядке $\varepsilon$:

$$\nu^2=\omega^2+2\omega\varepsilon B_1+O(\varepsilon^2),$$

$$-a\nu^2\cos\psi=-a\omega^2\cos\psi-2a\omega\varepsilon B_1\cos\psi+O(\varepsilon^2),$$

$$-2\dot a\,\nu\sin\psi=-2\omega\varepsilon A_1\sin\psi+O(\varepsilon^2).$$

Сложение с $\omega^2\cdot a\cos\psi$ уничтожает $-a\omega^2\cos\psi$. Вклад $a\cos\psi$ в $\ddot x+\omega^2 x$ в порядке $\varepsilon$:

$$\varepsilon(-2\omega A_1\sin\psi-2a\omega B_1\cos\psi).$$

Вклад $\varepsilon u_1$: вторая производная по быстрому углу даёт $\varepsilon\omega^2\partial_{\psi\psi}u_1$, а $\omega^2\cdot\varepsilon u_1$ добавляет $\varepsilon\omega^2 u_1$. Поправки из $\dot a$ и из $B_1$, действующие на $u_1$, начинаются с $\varepsilon\cdot\varepsilon$.

Коэффициент при $\varepsilon$ в $\ddot x+\omega^2 x-\varepsilon c_2 x^2=0$, с правой частью $c_2(a\cos\psi)^2$ (следующий порядок $\varepsilon^2$ от произведения $2a\cos\psi\cdot u_1$):

$$\omega^2(\partial_{\psi\psi}u_1+u_1)-2\omega A_1\sin\psi-2a\omega B_1\cos\psi
=c_2 a^2\cos^2\psi.$$

### A.2. Почему первая гармоника уходит в $A_k$, $B_k$

Оператор $\partial_{\psi\psi}+1$ на гармонике $e^{in\psi}$ даёт множитель $1-n^2$. При $n=\pm 1$ множитель равен нулю: периодическая $u_k$ имеет нулевую проекцию на $\cos\psi$ и $\sin\psi$.

Проекция уравнения на $\sin\psi$ определяет $A_k$ (член $-2\omega A_k$), проекция на $\cos\psi$ определяет $B_k$ (член $-2a\omega B_k$).

### A.3. Первый порядок

$$\cos^2\psi=\frac{1+\cos 2\psi}{2},$$

$$c_2 a^2\cos^2\psi=\frac{c_2 a^2}{2}+\frac{c_2 a^2}{2}\cos 2\psi.$$

Проекции на $\sin\psi$ и $\cos\psi$ справа равны нулю, слева равны $-2\omega A_1$ и $-2a\omega B_1$. Значит $A_1=0$, $B_1=0$.

Постоянная: $\omega^2 u_{10}=c_2 a^2/2$,

$$u_{10}=\frac{c_2 a^2}{2\omega^2}.$$

$\cos 2\psi$: множитель $1-4=-3$,

$$\omega^2(-3)u_{12}=\frac{c_2 a^2}{2},
\qquad
u_{12}=-\frac{c_2 a^2}{6\omega^2}.$$

$$u_1=\frac{c_2 a^2}{2\omega^2}-\frac{c_2 a^2}{6\omega^2}\cos 2\psi
=\frac{c_2 a^2}{2\omega^2}\Bigl(1-\frac13\cos 2\psi\Bigr).$$

### A.4. Второй порядок, проекции

При $A_1=B_1=0$ имеем $\dot a=\varepsilon^2 A_2+O(\varepsilon^3)$, $\nu=\omega+\varepsilon^2 B_2+O(\varepsilon^3)$.

Скорость до порядка $\varepsilon^2$:

$$\dot x
=-a\omega\sin\psi
+\varepsilon\omega\partial_\psi u_1
+\varepsilon^2(A_2\cos\psi-a B_2\sin\psi+\omega\partial_\psi u_2)
+O(\varepsilon^3).$$

(Член $\varepsilon\cdot\varepsilon^2$ из $B_2\partial_\psi u_1$ и член $\dot a\cdot\varepsilon\partial_a u_1$ имеют порядок $\varepsilon^3$.)

Вторая производная, коэффициент при $\varepsilon^2$, после добавления $\omega^2 x$:

$$\omega^2(\partial_{\psi\psi}u_2+u_2)-2\omega A_2\sin\psi-2a\omega B_2\cos\psi.$$

Двойка собирается из двух равных вкладов: дифференцирования ведущего $-a\omega\sin\psi$ при медленных $\dot a$ и $\nu-\omega$ и дифференцирования самой поправки $\varepsilon^2(A_2\cos\psi-a B_2\sin\psi)$ быстрым $\nu\approx\omega$.

Правая часть $\varepsilon c_2(a\cos\psi+\varepsilon u_1)^2$ в порядке $\varepsilon^2$:

$$2 c_2 a u_1\cos\psi.$$

Подстановка $u_1$:

$$2a u_1\cos\psi
=2a\cdot\frac{c_2 a^2}{2\omega^2}\cos\psi\Bigl(1-\frac13\cos 2\psi\Bigr)
=\frac{c_2 a^3}{\omega^2}
\Bigl(\cos\psi-\frac13\cos\psi\cos 2\psi\Bigr).$$

$$\cos\psi\cos 2\psi=\frac12(\cos 3\psi+\cos\psi),$$

$$\frac13\cos\psi\cos 2\psi=\frac16\cos 3\psi+\frac16\cos\psi.$$

$$\cos\psi-\frac13\cos\psi\cos 2\psi
=\Bigl(1-\frac16\Bigr)\cos\psi-\frac16\cos 3\psi
=\frac56\cos\psi-\frac16\cos 3\psi.$$

Умножение на $c_2$:

$$2 c_2 a u_1\cos\psi
=\frac{c_2^2 a^3}{\omega^2}
\Bigl(\frac56\cos\psi-\frac16\cos 3\psi\Bigr).$$

Проекция на $\sin\psi$ справа равна нулю, слева равна $-2\omega A_2$. Значит $A_2=0$, амплитуда $a$ постоянна.

Проекция на $\cos\psi$:

$$-2a\omega B_2=\frac56\frac{c_2^2 a^3}{\omega^2},$$

$$B_2
=-\frac{5}{6}\cdot\frac{c_2^2 a^3}{2a\omega\cdot\omega^2}
=-\frac{5}{12}\frac{c_2^2 a^2}{\omega^3}.$$

Проекция на $\cos 3\psi$. Множитель оператора $1-9=-8$:

$$\omega^2(-8)U_3=-\frac16\frac{c_2^2 a^3}{\omega^2},$$

$$U_3
=\frac{1}{8}\cdot\frac{1}{6}\frac{c_2^2 a^3}{\omega^4}
=\frac{c_2^2 a^3}{48\omega^4}.$$

$$u_2=\frac{c_2^2 a^3}{48\omega^4}\cos 3\psi.$$

Эта гармоника не удерживается в решении главы: она порядка $c_2^2 a^3$, то есть порядка $A^3$ после подстановки $c_2\sim\omega^2\eta_1/\eta_0$.

### A.5. Подстановка в уравнение главы

$$\varepsilon=1,\qquad \omega=\omega_s,\qquad a=A,\qquad c_2=-\frac{\omega_s^2\eta_1}{\eta_0}.$$

Постоянный член:

$$\frac{c_2 A^2}{2\omega_s^2}
=\frac{1}{2\omega_s^2}\Bigl(-\frac{\omega_s^2\eta_1}{\eta_0}\Bigr)A^2
=-\frac{\eta_1 A^2}{2\eta_0}
=\delta_{\mathrm{eq}}.$$

Коэффициент при $\cos 2\psi$:

$$-\frac{c_2 A^2}{6\omega_s^2}
=-\frac{1}{6\omega_s^2}\Bigl(-\frac{\omega_s^2\eta_1}{\eta_0}\Bigr)A^2
=\frac{\eta_1 A^2}{6\eta_0}.$$

$$c_2^2=\omega_s^4\Bigl(\frac{\eta_1}{\eta_0}\Bigr)^2,
\qquad
\frac{c_2^2}{\omega_s^4}=\Bigl(\frac{\eta_1}{\eta_0}\Bigr)^2.$$

При $\varepsilon=1$

$$\Omega=\omega+B_2
=\omega_s-\frac{5}{12}\frac{c_2^2 A^2}{\omega_s^3}
=\omega_s\Bigl(1-\frac{5}{12}\frac{c_2^2 A^2}{\omega_s^4}\Bigr),$$

$$\frac{\Omega}{\omega_s}
=1-\frac{5}{12}\Bigl(\frac{\eta_1 A}{\eta_0}\Bigr)^2.$$

### A.6. Квадратура $\varphi$

$$\delta=\delta_{\mathrm{eq}}+A\cos\psi+D\cos 2\psi,
\qquad
D=\frac{\eta_1 A^2}{6\eta_0},
\qquad
\delta_{\mathrm{eq}}=-\frac{\eta_1 A^2}{2\eta_0}=-3D.$$

До порядка $A^2$ (слагаемые $2\delta_{\mathrm{eq}} A\cos\psi$, $2AD\cos\psi\cos 2\psi$, $\delta_{\mathrm{eq}}^2$, $D^2$ суть $O(A^3)$ или $O(A^4)$):

$$\delta^2=A^2\cos^2\psi+O(A^3)
=\frac{A^2}{2}(1+\cos 2\psi)+O(A^3).$$

$$\tau=\eta_0(\delta_{\mathrm{eq}}+A\cos\psi+D\cos 2\psi)
+\eta_1\frac{A^2}{2}(1+\cos 2\psi)+O(A^3).$$

Постоянный член:

$$\eta_0\delta_{\mathrm{eq}}+\eta_1\frac{A^2}{2}
=\eta_0\Bigl(-\frac{\eta_1 A^2}{2\eta_0}\Bigr)+\frac{\eta_1 A^2}{2}
=0.$$

Коэффициент при $\cos\psi$: $\eta_0 A$.

Коэффициент при $\cos 2\psi$:

$$\eta_0 D+\eta_1\frac{A^2}{2}
=\eta_0\cdot\frac{\eta_1 A^2}{6\eta_0}+\frac{\eta_1 A^2}{2}
=\eta_1 A^2\Bigl(\frac16+\frac12\Bigr)
=\eta_1 A^2\Bigl(\frac16+\frac36\Bigr)
=\frac{2}{3}\eta_1 A^2.$$

$$\tau=\eta_0 A\cos\psi+\frac{2}{3}\eta_1 A^2\cos 2\psi+O(A^3).$$

$$d\varphi/dt=\omega_{\mathrm{rf}}\tau,
\qquad dt=d\psi/\Omega,$$

$$\varphi
=\frac{\omega_{\mathrm{rf}}}{\Omega}
\int\Bigl(\eta_0 A\cos\psi+\frac{2}{3}\eta_1 A^2\cos 2\psi\Bigr)d\psi.$$

$$\int\cos\psi\,d\psi=\sin\psi,
\qquad
\int\cos 2\psi\,d\psi=\frac12\sin 2\psi,$$

$$\int\frac{2}{3}\eta_1 A^2\cos 2\psi\,d\psi
=\frac{2}{3}\eta_1 A^2\cdot\frac12\sin 2\psi
=\frac{\eta_1 A^2}{3}\sin 2\psi.$$

Постоянная интегрирования равна нулю: $\varphi$ есть отклонение от устойчивой точки, и правая часть нечётна по $\psi$. Итог:

$$\varphi=\frac{\omega_{\mathrm{rf}}}{\Omega}
\Bigl(\eta_0 A\sin\psi+\frac{\eta_1 A^2}{3}\sin 2\psi\Bigr).$$
