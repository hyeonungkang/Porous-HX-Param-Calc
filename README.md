# Annular Fin → Porous Media 근사: 수식 유도

## 목표

**Annular finned tube에서 핀이 존재하는 Annular 영역만을 등가 Porous Media로 치환**

```
              Annular Zone (Porous로 근사)
                    ↓
            ┌──────────────┐
            │   ████████   │ ← 핀 영역
    ────────┤      ○       ├────────  유동 방향 →
            │   ████████   │
            └──────────────┘
                 ↑
              튜브 (Solid)
```

**입력:** $F_s$ (핀 간격), $h_f$ (핀 높이), 환경조건 $(T, v)$
**출력:** Porous 파라미터 $(1/K, C_2)$

---

## 1. 기하학적 관계식

### 1.1 고정 파라미터

| 기호 | 의미 | 값 |
|------|------|-----|
| $D_c$ | 튜브 외경 | 24 mm |
| $\delta_f$ | 핀 두께 | 0.5 mm |
| $S_1$ | 횡방향 피치 | 55.333 mm |
| $N$ | 튜브 열 수 | 4 |

### 1.2 종속 파라미터 유도

**핀 외경** (Annular Zone의 외경):
$$D_o = D_c + 2h_f$$

**핀 피치** (핀 중심 간 거리):
$$F_p = F_s + \delta_f$$

**공극률** (Annular Zone 내 유체가 차지하는 체적 비율):
$$\varepsilon = 1 - \frac{\delta_f}{F_p} = \frac{F_s}{F_s + \delta_f}$$

**최소 자유유동 면적비** (유동이 통과하는 최소 단면적 / 전체 정면 면적):
$$\sigma = \frac{S_1 - D_c - 2h_f \cdot (\delta_f / F_p)}{S_1}$$

---

## 2. 면적비 (Area Ratio) 계산

Nir 상관식의 핵심 인자인 **면적비**를 계산한다.

### 2.1 1 피치 구간의 표면적

**핀 표면적** (양면 + 외주 팁):
$$A_{fin} = 2 \cdot \frac{\pi}{4}(D_o^2 - D_c^2) + \pi D_o \delta_f$$

**튜브 노출 면적** (핀 사이 간격):
$$A_{base} = \pi D_c F_s$$

**총 열전달 표면적:**
$$A_{total} = A_{fin} + A_{base}$$

**맨 튜브 표면적** (핀이 없을 때):
$$A_{bare} = \pi D_c F_p$$

### 2.2 면적비

$$\boxed{AR = \frac{A_{total}}{A_{bare}}}$$

> **물리적 의미:** 핀으로 인해 표면적이 몇 배 증가했는지를 나타냄.
> $F_s \downarrow$ (핀 조밀) → $AR \uparrow$ → 저항 증가

---

## 3. Nir (1991) 마찰계수 상관식

### 3.1 경험식

$$\boxed{f_N = 1.1 \cdot Re_{D_c}^{-0.25} \cdot \left(\frac{S_1}{D_c}\right)^{-0.4} \cdot AR^{0.15}}$$

### 3.2 Reynolds 수 정의

$$Re_{D_c} = \frac{\rho \cdot v_{max} \cdot D_c}{\mu}$$

여기서 **최대 속도** (최소 단면적 통과 시):
$$v_{max} = \frac{v_{inlet}}{\sigma}$$

### 3.3 각 항의 역할

| 항 | 지수 | 효과 |
|----|------|------|
| $Re^{-0.25}$ | 음수 | 속도 증가 → 마찰계수 감소 (난류화) |
| $(S_1/D_c)^{-0.4}$ | 음수 | 피치 좁음 → 마찰 증가 |
| $AR^{0.15}$ | 양수 | 면적비 큼 → 마찰 증가 |

---

## 4. 압력강하 계산

### 4.1 총 압력강하

$$\boxed{\Delta P_{total} = N \cdot f_N \cdot \frac{\rho \cdot v_{max}^2}{2}}$$

### 4.2 단위 길이당 압력강하

유동 방향 깊이: $L = S_2 \cdot N$

$$Y(v) \equiv \frac{\Delta P}{L} = \frac{f_N \cdot \rho \cdot v_{max}^2}{2 \cdot S_2}$$

### 4.3 속도 의존성

$f_N \propto v_{max}^{-0.25}$ 이므로:

$$Y(v) \propto v^{-0.25} \cdot v^2 = v^{1.75}$$

> **문제:** Darcy-Forchheimer 형태 $(Av + Bv^2)$와 불일치
> **해결:** 2-Point Fitting으로 근사

---

## 5. Darcy-Forchheimer 방정식

### 5.1 Porous Media 기본 모델

CFD에서 porous zone의 압력강하:

$$\boxed{\frac{\Delta P}{L} = \underbrace{\frac{\mu}{K} \cdot v}_{\text{점성 저항}} + \underbrace{\frac{C_2 \rho}{2} \cdot v^2}_{\text{관성 저항}}}$$

### 5.2 물리적 의미

| 항 | 지배 현상 | 지배 영역 |
|----|-----------|-----------|
| $\frac{\mu}{K}v$ | 벽면 전단, 점성 마찰 | 저속 (Re 작음) |
| $\frac{C_2 \rho}{2}v^2$ | 형상 저항, 와류 손실 | 고속 (Re 큼) |

---

## 6. 2-Point Fitting: 핵심 변환

### 6.1 목표

Nir 식의 $Y(v) \propto v^{1.75}$를 Darcy-Forchheimer 형태로 근사:

$$Y(v) \approx A \cdot v + B \cdot v^2$$

### 6.2 두 속도점 선택

| 점 | 속도 | 용도 |
|----|------|------|
| 고속점 | $v_1 = v_{design}$ | 설계 조건 |
| 저속점 | $v_2 = 0.3 \cdot v_{design}$ | 저속 영역 보정 |

### 6.3 연립방정식

$$\begin{cases}
Y_1 = A \cdot v_1 + B \cdot v_1^2 \\[6pt]
Y_2 = A \cdot v_2 + B \cdot v_2^2
\end{cases}$$

### 6.4 해법

**Step 1:** 양변을 속도로 나누어 선형화

$$Z_1 = \frac{Y_1}{v_1} = A + B \cdot v_1$$
$$Z_2 = \frac{Y_2}{v_2} = A + B \cdot v_2$$

**Step 2:** $B$ 도출

$$\boxed{B = \frac{Z_1 - Z_2}{v_1 - v_2} = \frac{Y_1/v_1 - Y_2/v_2}{v_1 - v_2}}$$

**Step 3:** $A$ 도출

$$\boxed{A = Z_1 - B \cdot v_1 = \frac{Y_1}{v_1} - B \cdot v_1}$$

---

## 7. Porous 파라미터 최종 도출

### 7.1 계수 대응

Darcy-Forchheimer: $\frac{\Delta P}{L} = \frac{\mu}{K}v + \frac{C_2 \rho}{2}v^2$

Fitting 결과: $\frac{\Delta P}{L} = Av + Bv^2$

비교하면:

$$\frac{\mu}{K} = A \quad \Rightarrow \quad \boxed{\frac{1}{K} = \frac{A}{\mu}}$$

$$\frac{C_2 \rho}{2} = B \quad \Rightarrow \quad \boxed{C_2 = \frac{2B}{\rho}}$$

### 7.2 최종 공식 정리

| Porous 파라미터 | 공식 | 단위 | CFD 입력명 |
|-----------------|------|------|------------|
| **점성 저항** | $\displaystyle\frac{1}{K} = \frac{A}{\mu}$ | 1/m² | Viscous Resistance |
| **관성 저항** | $\displaystyle C_2 = \frac{2B}{\rho}$ | 1/m | Inertial Resistance |
| 투과도 | $K = \mu / A$ | m² | Permeability |

---

## 8. 전체 계산 흐름

```
┌─────────────────────────────────────────────────────────────┐
│  INPUT: Fs, hf, T, v                                        │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  Step 1: 기하 파라미터                                       │
│    Do = Dc + 2hf                                            │
│    Fp = Fs + δf                                             │
│    σ = (S1 - Dc - 2hf·δf/Fp) / S1                          │
│    AR = A_total / A_bare                                    │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  Step 2: Nir 마찰계수 (두 속도점에서)                        │
│    v_max = v / σ                                            │
│    Re = ρ·v_max·Dc / μ                                      │
│    f_N = 1.1·Re^(-0.25)·(S1/Dc)^(-0.4)·AR^(0.15)           │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  Step 3: 압력강하 (두 속도점에서)                            │
│    Y = f_N·ρ·v_max² / (2·S2)                                │
│    → Y₁ @ v₁,  Y₂ @ v₂                                      │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  Step 4: 2-Point Fitting                                    │
│    Z₁ = Y₁/v₁,  Z₂ = Y₂/v₂                                  │
│    B = (Z₁ - Z₂) / (v₁ - v₂)                                │
│    A = Z₁ - B·v₁                                            │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  Step 5: Porous 파라미터                                    │
│    1/K = A / μ      [1/m²]                                  │
│    C₂ = 2B / ρ      [1/m]                                   │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  OUTPUT: 1/K, C₂  → CFD Porous Zone 입력                    │
└─────────────────────────────────────────────────────────────┘
```

---

## 9. 적용 예시

### 9.1 조건

- 온도: $T = 14.80177$ °C
- 풍속: $v = 2.019723$ m/s
- 공기 물성: $\rho = 1.2258$ kg/m³, $\mu = 1.788 \times 10^{-5}$ Pa·s

### 9.2 케이스: $F_s = 4$ mm, $h_f = 4$ mm

**Step 1: 기하 파라미터**

| 계산 | 결과 |
|------|------|
| $D_o = 24 + 2(4) = 32$ mm | |
| $F_p = 4 + 0.5 = 4.5$ mm | |
| $\varepsilon = 4/4.5 = 0.889$ | |
| $\sigma = 0.550$ | (수치 계산) |
| $AR = 3.11$ | (면적 계산) |

**Step 2: 두 속도점에서 $f_N$ 계산**

고속점 ($v_1 = 2.02$ m/s):
- $v_{max,1} = 2.02 / 0.550 = 3.67$ m/s
- $Re_1 = 1.2258 \times 3.67 \times 0.024 / 1.788 \times 10^{-5} = 6039$
- $f_{N,1} = 1.1 \times 6039^{-0.25} \times 2.306^{-0.4} \times 3.11^{0.15} = 0.106$

저속점 ($v_2 = 0.606$ m/s):
- $Re_2 = 1812$
- $f_{N,2} = 0.146$

**Step 3: 압력강하**

$$Y_1 = \frac{0.106 \times 1.2258 \times 3.67^2}{2 \times 0.0553} = 15.8 \text{ Pa/m}$$

$$Y_2 = \frac{0.146 \times 1.2258 \times 1.10^2}{2 \times 0.0553} = 1.96 \text{ Pa/m}$$

**Step 4: 2-Point Fitting**

$$Z_1 = 15.8 / 2.02 = 7.82$$
$$Z_2 = 1.96 / 0.606 = 3.23$$

$$B = \frac{7.82 - 3.23}{2.02 - 0.606} = 3.24 \text{ Pa·s}^2/\text{m}^3$$

$$A = 7.82 - 3.24 \times 2.02 = 1.18 \text{ Pa·s/m}^2$$

**Step 5: Porous 파라미터**

$$\frac{1}{K} = \frac{1.18}{1.788 \times 10^{-5}} = 6.59 \times 10^4 \text{ 1/m}^2$$

$$C_2 = \frac{2 \times 3.24}{1.2258} = 5.29 \text{ 1/m}$$

### 9.3 결과 요약

| $F_s$ | $h_f$ | $\varepsilon$ | $AR$ | $1/K$ [1/m²] | $C_2$ [1/m] |
|-------|-------|---------------|------|--------------|-------------|
| 2 | 4 | 0.800 | 4.80 | 7.33×10⁴ | 5.97 |
| **4** | **4** | **0.889** | **3.11** | **6.59×10⁴** | **5.37** |
| 6 | 4 | 0.923 | 2.46 | 6.26×10⁴ | 5.11 |
| 8 | 4 | 0.941 | 2.12 | 6.07×10⁴ | 4.95 |
| 4 | 6 | 0.889 | 4.39 | 7.12×10⁴ | 5.80 |
| 4 | 8 | 0.889 | 5.81 | 7.62×10⁴ | 6.21 |

---

## 10. 핵심 공식 요약

### 입력 → 출력 변환의 핵심 3단계

**① Nir 경험식으로 압력강하 계산**
$$\frac{\Delta P}{L} = \frac{N \cdot f_N \cdot \rho \cdot v_{max}^2}{2L}$$

**② 2-Point Fitting으로 계수 추출**
$$A = \frac{Y_1}{v_1} - \frac{Y_1/v_1 - Y_2/v_2}{v_1 - v_2} \cdot v_1$$
$$B = \frac{Y_1/v_1 - Y_2/v_2}{v_1 - v_2}$$

**③ Porous 파라미터로 변환**
$$\boxed{\frac{1}{K} = \frac{A}{\mu}} \quad \boxed{C_2 = \frac{2B}{\rho}}$$

---

## 11. 사용 방법

### 11.1 코드 실행

```bash
python cfd_porous_calc.py
```

실행하면 다음과 같은 결과를 얻을 수 있습니다:
- 단일 케이스 상세 결과 (Fs=4mm, hf=4mm)
- 여러 케이스 비교표
- CFD 소프트웨어 입력 가이드
- JSON 형식 출력

### 11.2 커스텀 계산

Python 코드 내에서 직접 함수 호출:

```python
from cfd_porous_calc import calculate_porous_parameters

result = calculate_porous_parameters(
    Fs_mm=4.0,      # 핀 간격 [mm]
    hf_mm=4.0,      # 핀 높이 [mm]
    T_celsius=15.0, # 온도 [°C]
    v_design=2.0    # 설계 풍속 [m/s]
)

print(f"점성 저항 (1/K): {result['porous']['inv_K']:.4e} [1/m²]")
print(f"관성 저항 (C2): {result['porous']['C2']:.4f} [1/m]")
```

---

## 12. CFD 소프트웨어 입력 가이드

### ANSYS Fluent

**Porous Zone 설정:**
1. `Cell Zone Conditions` → Porous Zone 체크
2. `Viscous Resistance (1/m²)`: 계산 결과의 `1/K` 값 입력
3. `Inertial Resistance (1/m)`: 계산 결과의 `C2` 값 입력
4. **방향성 설정:**
   - Direction-1 (유동 방향): 계산값 그대로
   - Direction-2, 3 (횡방향): 계산값 × 1000 (유동 차단)

**LTNE (Local Thermal Non-Equilibrium) 설정:**
1. `Non-Equilibrium Thermal Model` 활성화
2. `Surface Area Density`: 계산 결과의 `a_fs` 값 입력
3. `Fluid Solid Heat Transfer Coefficient`: 계산 결과의 `h_fs` 값 입력

### SimScale

1. `Advanced Concepts` → `Porous Media`
2. `Permeability K (m²)`: $K = 1 / (1/K)$ 계산하여 입력
3. `Forchheimer coefficient C_F (1/m)`: `C2` 값 직접 입력

### COMSOL Multiphysics

1. `Porous Media Flow` → `Darcy-Forchheimer`
2. `Permeability (m²)`: `K` 값 입력
3. `Forchheimer coefficient (1/m)`: `C2` 값 입력

---

## 13. 참고문헌

1. **Nir, A.** (1991). "Heat Transfer and Friction Factor Correlations for Crossflow over Staggered Finned Tube Banks", *Heat Transfer Engineering*, Vol.12, No.1, pp.43-58.

2. **Briggs, D.E., & Young, E.H.** (1963). "Convection heat transfer and pressure drop of air flowing across triangular pitch banks of finned tubes."

---

## 14. 라이센스 및 기여

이 프로젝트는 교육 및 연구 목적으로 제공됩니다. 개선 사항이나 버그가 있으면 이슈를 등록해주세요.

**작성자:** Claude (Anthropic)
**최종 업데이트:** 2026-01-17
**버전:** 2.0 (Nir 상관식 기반)
