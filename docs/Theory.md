# Automatic-Partial-Wave-Decomposition of Three-Body Force (aPWD3)

The method is from https://doi.org/10.1140/epja/i2009-10903-6 and https://doi.org/10.1140/epja/i2011-11048-9. We specially thank professor Kacper Topolnicki for discussions.

## spin-coupled 3N states

starting from uncoupled 3N states in the (231) particle order

$$
\ket{\dfrac{1}{2}m_2}\otimes \ket{\dfrac{1}{2}m_3}\otimes \ket{\dfrac{1}{2}m_1}
$$

in (23)1 representation, particle-1 is the spectator.

first couple 2 and 3:

$$
\begin{aligned}
\left| (\frac{1}{2}\frac{1}{2})s m_s \right\rangle=&
\sum_{m_2=\pm 1/2}\sum_{m_3=\pm 1/2}
C(\frac{1}{2},\frac{1}{2},s;m_2,m_3,m_s)
\left| \frac{1}{2} m_2 \right\rangle \otimes
\left| \frac{1}{2} m_3 \right\rangle\\
=&
\sum_{m_2=\pm 1/2}
C(\frac{1}{2},\frac{1}{2},s;m_2,m_s-m_2,m_s)
\left| \frac{1}{2} m_2 \right\rangle \otimes
\left| \frac{1}{2}, m_s-m_2 \right\rangle\\
\end{aligned}
$$

where CG coefficients: C(j1,j2,j3;m1,m2,m3)

and then couple with particle-1:

$$
\begin{aligned}
\left| (s\frac{1}{2})S M_S \right\rangle=&
\sum_{m_1=\pm 1/2}\sum_{m_s=0,1}
C(\frac{1}{2},s,S;m_1,m_s,M_S)
\left| s m_s \right\rangle \otimes
\left| \frac{1}{2}m_1 \right\rangle\\
=&
\sum_{m_1=\pm 1/2}
C(\frac{1}{2},s,S;m_1,M_S-m_1,M_S)
\left| s ,M_S-m_1 \right\rangle \otimes
\left| \frac{1}{2}m_1 \right\rangle\\
=&
\sum_{m_1=\pm 1/2}\sum_{m_2=\pm 1/2}
C(\frac{1}{2},s,S;m_1,M_S-m_1,M_S)
C(\frac{1}{2},\frac{1}{2},s;m_2,M_S-m_1-m_2,M_S-m_1)\\
&\left| \frac{1}{2} m_2 \right\rangle \otimes
\left| \frac{1}{2}, M_S-m_1-m_2 \right\rangle \otimes
\left| \frac{1}{2}m_1 \right\rangle
\end{aligned}
$$

in Mathematica, th spin-coupled 3N states are constructed to be vectors of length 8,

```mathematica
spin1N[m_] := Switch[m, 1/2, {1, 0}, -1/2, {0, 1}, _, {0, 0}];
spin3N[s_, S_, MS_] := Sum[CG[s, 1/2, S, MS - m1, m1, MS]* CG[1/2, 1/2, s, m2, MS - m1 - m2, MS - m1]* ArrayFlatten[
 KroneckerProduct[spin1N[m2], spin1N[MS - m1 - m2], spin1N[m1]], 1], {m1, -1/2, 1/2}, {m2, -1/2, 1/2}];
```

## chiral 3NF at N2LO

define $\boldsymbol{p}_i(\boldsymbol{p}_i')$ the initial(final) momentum of the nucleon i, $\boldsymbol{q}_i=\boldsymbol{p}_i'-\boldsymbol{p}_i$ the momentum change.

three topologies: TPE, OPE, Contacts.

### TPE

$$
V_{\text{3N}}^{2\pi}=V_{\text{TPE1}}+V_{\text{TPE2}}
$$

where

$$
V_{\text{TPE1}}=\dfrac{1}{2}\left(\dfrac{g_A}{2f_\pi}\right)^2\sum_{(i,j,k)\in \mathcal{I}}
\dfrac{\boldsymbol{\sigma}_i\cdot \boldsymbol{q}_i \;\boldsymbol{\sigma}_j\cdot \boldsymbol{q}_j}{(q_i^2+m_\pi^2)(q_j^2+m_\pi^2)}
\left( -\dfrac{4 c_1 m_\pi^2}{f_\pi^2}+\dfrac{2 c_3}{f_\pi^2} \boldsymbol{q}_i \cdot \boldsymbol{q}_j \right)\boldsymbol{\tau}_i \cdot \boldsymbol{\tau}_j
$$

$$
V_{\text{TPE2}}=\dfrac{1}{2}\left(\dfrac{g_A}{2f_\pi}\right)^2\sum_{(i,j,k)\in \mathcal{I}}
\dfrac{\boldsymbol{\sigma}_i\cdot \boldsymbol{q}_i \;\boldsymbol{\sigma}_j\cdot \boldsymbol{q}_j}{(q_i^2+m_\pi^2)(q_j^2+m_\pi^2)}
\left( \dfrac{c_4}{f_\pi^2} \right)
\boldsymbol{\sigma}_k\cdot (\boldsymbol{q}_i\times \boldsymbol{q}_j)
\boldsymbol{\tau}_k\cdot (\boldsymbol{\tau}_i\times \boldsymbol{\tau}_j)
$$

6 possible summations for $i\neq j \neq k$ :

$$
\mathcal{I}=\set{(123),(132),(213),(312),(231),(321)}
$$

### OPE

$$
V_{\text{3N}}^{1\pi}=-\dfrac{g_A}{8f_\pi^2}\dfrac{c_D}{f_\pi^2 \Lambda_\chi}
\sum_{(i,j,k)\in \mathcal{I}} \dfrac{\boldsymbol{\sigma}_j\cdot \boldsymbol{q}_j}{q_j^2+m_\pi^2}(\boldsymbol{\sigma}_i\cdot \boldsymbol{q}_j)(\boldsymbol{\tau}_i \cdot \boldsymbol{\tau}_j)
$$

### Contact

$$
V_{\text{3N}}^{\text{contact}}=\dfrac{c_E}{2f_\pi^4 \Lambda_\chi}\sum_{(i,j,k)\in \mathcal{I}}
(\boldsymbol{\tau}_j \cdot \boldsymbol{\tau}_k)
$$

### Sum

total 3NF at N2LO:

$$
\begin{aligned}
V_{\text{3N}}=\sum_{(i,j,k)\in \mathcal{I}} \cdot \Big[ &
F_{\text{TPE1}}^{ij}(\boldsymbol{\sigma}_i\cdot \boldsymbol{q}_i) \;(\boldsymbol{\sigma}_j\cdot \boldsymbol{q}_j)\;(\boldsymbol{\tau}_i \cdot \boldsymbol{\tau}_j)\\
+&F_{\text{TPE2}}^{ij}(\boldsymbol{\sigma}_i\cdot \boldsymbol{q}_i) \;(\boldsymbol{\sigma}_j\cdot \boldsymbol{q}_j)\;\boldsymbol{\sigma}_k\cdot (\boldsymbol{q}_i\times \boldsymbol{q}_j)
\;\boldsymbol{\tau}_k\cdot (\boldsymbol{\tau}_i\times \boldsymbol{\tau}_j)\\
+&F_{\text{OPE}}^{j}(\boldsymbol{\sigma}_j\cdot \boldsymbol{q}_j)\;(\boldsymbol{\sigma}_i\cdot \boldsymbol{q}_j)\;(\boldsymbol{\tau}_i \cdot \boldsymbol{\tau}_j)\\
+&F_{\text{contact}}(\boldsymbol{\tau}_j \cdot \boldsymbol{\tau}_k)
\Big]
\end{aligned}
$$

where the spin-independent functions:

$$
F_{\text{TPE1}}^{ij}=\dfrac{1}{2}\left(\dfrac{g_A}{2f_\pi}\right)^2
\dfrac{1}{(q_i^2+m_\pi^2)(q_j^2+m_\pi^2)}\left( -\dfrac{4 c_1 m_\pi^2}{f_\pi^2}+\dfrac{2 c_3}{f_\pi^2} \boldsymbol{q}_i \cdot \boldsymbol{q}_j \right)
$$

$$
F_{\text{TPE2}}^{ij}=\dfrac{1}{2}\left(\dfrac{g_A}{2f_\pi}\right)^2
\dfrac{1}{(q_i^2+m_\pi^2)(q_j^2+m_\pi^2)}\left( \dfrac{c_4}{f_\pi^2} \right)
$$

$$
F_{\text{OPE}}^{j}=-\dfrac{g_A}{8f_\pi^2}\dfrac{c_D}{f_\pi^2 \Lambda_\chi}\dfrac{1}{q_j^2+m_\pi^2}
$$

$$
F_{\text{contact}}=\dfrac{c_E}{2f_\pi^4 \Lambda_\chi}
$$

using permutation, we only need to evaluate the Faddeev component of $V_{\text{3N}}$, $V_{\mathrm{3N}}^{(1)}$, which is invariant under switching particle 2 and 3:

$$
V_{3 \mathrm{N}}=V_{3 \mathrm{N}}^{(1)}+P_{123}^{-1} V_{3 \mathrm{N}}^{(1)} P_{123}+P_{132}^{-1} V_{3 \mathrm{N}}^{(1)} P_{132}
$$

which is:

$$
V_{\text{3N}}^{(1)}=\; c_1\cdot V_{c_1} + c_3\cdot V_{c_3} + c_4\cdot V_{c_4} + \frac{c_D}{\Lambda_\chi}\cdot V_{c_D} + \frac{c_E}{\Lambda_\chi}\cdot V_{c_E}
$$

where:

$$
V_{c_1} =F_{\text{TPE1}}(\boldsymbol{\sigma}_2\cdot \boldsymbol{q}_2) \;(\boldsymbol{\sigma}_3\cdot \boldsymbol{q}_3)\;(\boldsymbol{\tau}_2 \cdot \boldsymbol{\tau}_3)
$$

$$
V_{c_3} =F_{\text{TPE2}}(\boldsymbol{\sigma}_2\cdot \boldsymbol{q}_2) \;(\boldsymbol{\sigma}_3\cdot \boldsymbol{q}_3)\;(\boldsymbol{\tau}_2 \cdot \boldsymbol{\tau}_3)
$$

$$
V_{c_4} =F_{\text{TPE3}}\,\boldsymbol{\sigma}_1\cdot (\boldsymbol{q}_2\times \boldsymbol{q}_3)\;(\boldsymbol{\sigma}_2\cdot \boldsymbol{q}_2) \;(\boldsymbol{\sigma}_3\cdot \boldsymbol{q}_3)
\;\boldsymbol{\tau}_1\cdot (\boldsymbol{\tau}_2\times \boldsymbol{\tau}_3)
$$

$$
V_{c_D}=F_{\text{OPE1}}(\boldsymbol{\sigma}_3\cdot \boldsymbol{q}_3)\;(\boldsymbol{\sigma}_2\cdot \boldsymbol{q}_3)\;(\boldsymbol{\tau}_2 \cdot \boldsymbol{\tau}_3)+F_{\text{OPE2}}(\boldsymbol{\sigma}_2\cdot \boldsymbol{q}_2)\;(\boldsymbol{\sigma}_3\cdot \boldsymbol{q}_2)\;(\boldsymbol{\tau}_2 \cdot \boldsymbol{\tau}_3)
$$

$$
V_{c_E}=F_{\text{ct}}(\boldsymbol{\tau}_2 \cdot \boldsymbol{\tau}_3)
$$

with scalar functions:

$$
F_{\text{TPE1}}=\left(\dfrac{g_A}{2f_\pi}\right)^2
\dfrac{1}{(q_2^2+m_\pi^2)(q_3^2+m_\pi^2)}\left( -\dfrac{4  m_\pi^2}{f_\pi^2}\right)
$$

$$
F_{\text{TPE2}}=\left(\dfrac{g_A}{2f_\pi}\right)^2
\dfrac{1}{(q_2^2+m_\pi^2)(q_3^2+m_\pi^2)}\left( \dfrac{2}{f_\pi^2} \boldsymbol{q}_2 \cdot \boldsymbol{q}_3 \right)
$$

$$
F_{\text{TPE3}}=\left(\dfrac{g_A}{2f_\pi}\right)^2
\dfrac{1}{(q_2^2+m_\pi^2)(q_3^2+m_\pi^2)}\left( \dfrac{1}{f_\pi^2} \right)
$$

$$
F_{\text{OPE1}}=-\dfrac{g_A}{8f_\pi^2}\dfrac{1}{f_\pi^2 }\dfrac{1}{q_3^2+m_\pi^2}
$$

$$
F_{\text{OPE2}}=-\dfrac{g_A}{8f_\pi^2}\dfrac{1}{f_\pi^2 }\dfrac{1}{q_2^2+m_\pi^2}
$$

$$
F_{\text{ct}}=\dfrac{1}{f_\pi^4}
$$

with Mathematica, the spin structures of 3NF operators can be efficiently represented using KroneckerProduct, such as TPE1 term:

```mathematica
Um = IdentityMatrix[2];
Sigma1m = PauliMatrix[1];
Sigma2m = PauliMatrix[2];
Sigma3m = PauliMatrix[3];
Sigmam = {Sigma1m, Sigma2m, Sigma3m};
q1 = {q1x, q1y, q1z};
q2 = {q2x, q2y, q2z};
q3 = {q3x, q3y, q3z};
Vtpe1 = Ftpe1* KroneckerProduct[q2.Sigmam, q3.Sigmam, Um];
```

## LS scheme

partial-wave 3N states in LS coupling scheme:

$$
|\boldsymbol{p} \boldsymbol{q} \beta\rangle \equiv\left|\boldsymbol{p} \boldsymbol{q}(l \lambda) L\left(s \frac{1}{2}\right) S(L S) J M_J\right\rangle
 \otimes
\left|\left(t \frac{1}{2}\right) T m_T\right\rangle
$$

where total spin $s$ and relative angular momentum $l$ of the pair (23); $\lambda$ the angular momentum between pair (23) and 1; $l,\lambda$ coupled to $L$ the total 3N angular momentum; $s,1/2$ coupled to $S$ the total 3N spin; finally $L,S$ coupled to $J$.

and:

$$
\boldsymbol{q}_1=\boldsymbol{q}'-\boldsymbol{q}
$$

$$
\boldsymbol{q}_2=\boldsymbol{p}'-\dfrac{1}{2}\boldsymbol{q}'-\Big(\boldsymbol{p}-\dfrac{1}{2}\boldsymbol{q}\Big)
$$

$$
\boldsymbol{q}_3=-\boldsymbol{p}'-\dfrac{1}{2}\boldsymbol{q}'-\Big(-\boldsymbol{p}-\dfrac{1}{2}\boldsymbol{q}\Big)
$$

abbreviation:

$$
|\boldsymbol{p} \boldsymbol{q} \beta\rangle \equiv\left|\boldsymbol{p} \boldsymbol{q} l \lambda L s S t T J \right\rangle
$$

parity:

$$
P=(-1)^{l+\lambda}
$$

$T,J,P$ are good quantum numbers for 3NF.

(23) pair need to satisfy:

$$
(-1)^{l+s+t}=-1
$$

LS scheme matrix elements:

$$
\left\langle
\boldsymbol{p}' \boldsymbol{q}' l' \lambda' L' s' S' t' T J \right|V^{\mathrm{3N}}\left|
\boldsymbol{p} \boldsymbol{q} l \lambda L s S t T J
\right\rangle
\left\langle t' T\right|\hat{I}\left|t T\right\rangle
$$

define:

$$
\left\langle t' T\right|\hat{I}\left|t T\right\rangle=I(t',t,T)
$$

and:

$$
G(\beta', \beta)=I(t',t,T)\times
\left\langle
\boldsymbol{p}' \boldsymbol{q}' l' \lambda' L' s' S' t' T J \right|V^{\mathrm{3N}}\left|
\boldsymbol{p} \boldsymbol{q} l \lambda L s S t T J
\right\rangle
$$

using angular momentum coupling we have:

$$
\begin{aligned}
G(\beta',\beta)&=I(t',t,T)\cdot\dfrac{1}{2J+1}\sum_{M_J=-J}^{J}\sum_{m_{L'}=-L'}^{L'}\sum_{m_L=-L}^{L}\\
&\cdot C\left(L', S', J ; m_{L'}, M_J-m_{L'}, M_J\right) C\left(L, S, J ; m_L, M_J-m_L, M_J\right)\\
&\cdot \int d \hat{p}^{\prime} \int d \hat{q}^{\prime} \int d \hat{p} \int d \hat{q}\;
\mathcal{Y}_{l',\lambda'}^{*L',m_{L'}}(\hat{p}',\hat{q}') \mathcal{Y}_{l,\lambda}^{L,m_{L}}(\hat{p},\hat{q})\\
&\cdot
\left\langle s', S^{\prime}, M_J-m_{L^{\prime}}\right|V_{\mathrm{3N}}\left(\boldsymbol{p}^{\prime}, \boldsymbol{q}^{\prime}, \boldsymbol{p}, \boldsymbol{q}\right)\left|s, S, M_J-m_{L}\right\rangle
\end{aligned}
$$

where CG coefficients $C(j_1,j_2,j_3;m_1,m_2,m_3)$ and coupled spherical harmonics:

$$
\mathcal{Y}_{l,\lambda}^{L,m_{L}}(\hat{p},\hat{q})=\sum_{m_l=-l}^{l}C(l,\lambda,L,m_l,m_L-m_l,m_L)Y_{l,m_l}(\hat{p})Y_{\lambda,m_L-m_l}(\hat{q})
$$

we write into the following form:

$$
\begin{aligned}
\tilde{G}(\beta',\beta)&=I(t',t,T)\cdot\dfrac{1}{2J+1}\sum_{M_J=-J}^{J}\sum_{m_{L'}=-L'}^{L'}\sum_{m_L=-L}^{L}\\
&\cdot C\left(L', S', J ; m_{L'}, M_J-m_{L'}, M_J\right) C\left(L, S, J ; m_L, M_J-m_L, M_J\right)\\
&\cdot \mathcal{Y}_{l',\lambda'}^{*L',m_{L'}}(\hat{p}',\hat{q}') \mathcal{Y}_{l,\lambda}^{L,m_{L}}(\hat{p},\hat{q})
\left\langle s', S^{\prime}, M_J-m_{L^{\prime}}\right|V_{\mathrm{3N}}\left(\boldsymbol{p}^{\prime}, \boldsymbol{q}^{\prime}, \boldsymbol{p}, \boldsymbol{q}\right)\left|s, S, M_J-m_{L}\right\rangle
\end{aligned}
$$

which is calculated analytically in Mathematica:

```mathematica
Gt[lp_, lamp_, Lp_, sp_, twoSp_, l_, lam_, L_, s_, twoS_, twoJ_] := 1/(twoJ + 1)* ParallelSum[(CG[Lp, twoSp/2, twoJ/2, mLp, twoMJ/2 - mLp, twoMJ/2]* CG[L, twoS/2, twoJ/2, mL, twoMJ/2 - mL, twoMJ/2])* Ybra[Lp, mLp, lp, lamp, thetapp, phipp, thetaqp, phiqp, wigner]* Yket[L, mL, l, lam, thetaq, wigner]* (spin3N[sp, twoSp/2, twoMJ/2 - mLp].V.spin3N[s, twoS/2, twoMJ/2 - mL]), {twoMJ, -twoJ, twoJ, 2}, {mLp, -Lp, Lp}, {mL, -L, L}];
```

the isospin matrix elements can be calculated easily. There are two isospin operators:

$$
\hat{I}_1=(\boldsymbol{\tau}_2 \cdot \boldsymbol{\tau}_3)
$$

$$
\hat{I}_2=\boldsymbol{\tau}_1\cdot (\boldsymbol{\tau}_2\times \boldsymbol{\tau}_3)
$$

which can be easily checked using Mathematica:

$$
I_1(t',t,T)=(2t(t+1)-3)\,\delta_{t,t'}
$$

$$
I_2(t',t,T)=2\sqrt{3}i(-1)^{t}\,\delta_{t+t',1}\delta_{T,1/2}
$$

the final 3NF matrix elements in LS scheme:

$$
G(\beta',\beta)=\int d \hat{p}^{\prime} \int d \hat{q}^{\prime} \int d \hat{p} \int d \hat{q}\;\tilde{G}(\beta',\beta)
$$

due to rotational symmetry, we take:

$$
\hat{p}=\hat{z}
$$

$$
\phi_q=0
$$

so that:

$$
G(\beta',\beta)=8\pi^2 \int \sin \theta_{\hat{q}} \sin \theta_{\hat{p}'} \sin \theta_{\hat{q}'} d\theta_{\hat{q}} d\theta_{\hat{p}'} d\theta_{\hat{q}'} d\phi_{\hat{p}'} d\phi_{\hat{q}'} \;\tilde{G}(\beta',\beta)
$$

## JJ scheme

partial-wave 3N states in JJ coupling scheme:

$$
|\boldsymbol{p} \boldsymbol{q} \alpha\rangle \equiv\left|\boldsymbol{p} \boldsymbol{q} (ls)j (\lambda\frac{1}{2})j_3 (j j_3)J M_J\right\rangle
 \otimes
\left|\left(t \frac{1}{2}\right) T m_T\right\rangle
$$

abbreviation:

$$
|\boldsymbol{p} \boldsymbol{q} \alpha\rangle =\left|\boldsymbol{p} \boldsymbol{q} l s j \lambda j_3 t T J \right\rangle
$$

related to LS coupling by:

$$
|\boldsymbol{p} \boldsymbol{q} \alpha\rangle=\sum_{L,S}\sqrt{\hat{L}\hat{S}\hat{j}\hat{j}_3}
\left\{\begin{array}{ccc}
l & s & j \\
\lambda & \frac{1}{2} & j_3 \\
L & S & J
\end{array}\right\}
|\boldsymbol{p} \boldsymbol{q} \beta\rangle
$$

denote 3NF matrix elements under JJ scheme by $H(\alpha',\alpha)$, then:

$$
H(\alpha',\alpha)=
\sum_{S=1/2}^{3/2}\sum_{L=|J-S|}^{J+S}
\sum_{S'=1/2}^{3/2}\sum_{L'=|J-S'|}^{J+S'}
\sqrt{\hat{L}'\hat{L}\hat{S}'\hat{S}\hat{j}'\hat{j}\hat{j}_3'\hat{j}_3}
\left\{\begin{array}{ccc}
l' & s' & j' \\
\lambda' & \frac{1}{2} & j_3' \\
L' & S' & J
\end{array}\right\}
\left\{\begin{array}{ccc}
l & s & j \\
\lambda & \frac{1}{2} & j_3 \\
L & S & J
\end{array}\right\}
G(\beta',\beta)
$$
