# CFD Porous Zone 파라미터 계산 - 빠른 사용 가이드

## 개요

이 도구는 **Annular Fin 열교환기**를 **Porous Media**로 근사하여 CFD 시뮬레이션에 필요한 파라미터를 계산합니다.

**핵심 특징:**
- ✅ **Nir (1991) 상관식 기반** - Annular Finned Tube 전용 마찰계수 상관식
- ✅ **2-Point Fitting** - Darcy-Forchheimer 계수 정확 도출
- ✅ **온도 의존 물성치** - Sutherland's Law로 정밀 계산
- ✅ **LTNE 지원** - 열전달 계수 및 비표면적 출력

---

## 빠른 시작

### 1. 코드 실행

```bash
python cfd_porous_calc.py
```

### 2. 출력 확인

코드를 실행하면 다음 정보가 출력됩니다:

#### (1) 핵심 CFD 입력값

```
============================================================
  ★★★ CFD Porous Media 입력값 (핵심 출력) ★★★
============================================================
  점성 저항 (1/K)     : 6.5900e+04  [1/m²]
  관성 저항 (C2)      : 5.3700  [1/m]
  투과도 (K)          : 2.7100e-10  [m²]
============================================================
```

#### (2) 여러 케이스 비교표

```
====================================================================================================
  Porous 파라미터 비교표 (T = 14.80°C, v = 2.0197 m/s)
====================================================================================================
  Fs   hf |       ε       σ     AR |   1/K [1/m²]  C2 [1/m] |  ΔP [Pa]  h [W/m²K]
----------------------------------------------------------------------------------------------------
   2    4 |  0.8000  0.4778   4.80 |   7.3300e+04     5.9700 |     16.52     103.45
   4    4 |  0.8889  0.5502   3.11 |   6.5900e+04     5.3700 |     11.95      94.82
   6    4 |  0.9231  0.5923   2.46 |   6.2600e+04     5.1100 |     10.20      89.63
   8    4 |  0.9412  0.6204   2.12 |   6.0700e+04     4.9500 |      9.19      86.05
====================================================================================================
```

---

## CFD 소프트웨어 입력

### ANSYS Fluent

#### Porous Zone 설정

1. **Cell Zone Conditions** → Porous Zone 체크

2. **Viscous Resistance (1/m²)**
   ```
   Direction-1 (유동 방향):  6.59e+04  [계산값]
   Direction-2 (횡방향):     6.59e+07  [계산값 × 1000]
   Direction-3 (횡방향):     6.59e+07  [계산값 × 1000]
   ```

3. **Inertial Resistance (1/m)**
   ```
   Direction-1 (유동 방향):  5.37
   Direction-2 (횡방향):     5370
   Direction-3 (횡방향):     5370
   ```

4. **Porosity**
   ```
   모든 방향: 0.8889
   ```

#### LTNE (비평형 열전달) 설정

1. **Non-Equilibrium Thermal Model** 활성화

2. **Surface Area Density (a_fs)**: `698.85` [1/m]

3. **Fluid Solid Heat Transfer Coefficient (h_fs)**: `94.82` [W/(m²·K)]

---

### SimScale

1. **Porous Media** 활성화

2. **Permeability (K)**:
   ```
   K = 1 / (1/K) = 1 / 6.59e+04 = 1.52e-05 [m²]
   ```

3. **Forchheimer Coefficient (C2)**: `5.37` [1/m]

---

### COMSOL Multiphysics

1. **Porous Media Flow** → **Darcy-Forchheimer**

2. **Permeability (κ)**: `1.52e-05` [m²]

3. **Forchheimer Coefficient (β)**: `5.37` [1/m]

---

## 입력 변수 커스터마이징

코드 내에서 다음 변수를 수정할 수 있습니다:

```python
# cfd_porous_calc.py 파일 하단 (542행부터)

if __name__ == "__main__":
    # 환경 조건 수정
    T_ambient = 14.80177    # 대기 온도 [°C]
    v_wind = 2.019723       # 설계 풍속 [m/s]

    # 단일 케이스 계산
    result = calculate_porous_parameters(
        Fs_mm=4.0,      # 핀 간격 [mm]
        hf_mm=4.0,      # 핀 높이 [mm]
        T_celsius=T_ambient,
        v_design=v_wind
    )
```

또는 함수를 직접 호출:

```python
from cfd_porous_calc import calculate_porous_parameters

result = calculate_porous_parameters(
    Fs_mm=3.0,           # 핀 간격 [mm]
    hf_mm=5.0,           # 핀 높이 [mm]
    T_celsius=20.0,      # 온도 [°C]
    v_design=1.5,        # 설계 풍속 [m/s]
    Dc_mm=24.0,          # 튜브 외경 [mm] (선택사항)
    delta_f_mm=0.5,      # 핀 두께 [mm] (선택사항)
    s1_mm=55.333,        # 횡방향 피치 [mm] (선택사항)
    pitch_ratio=1.0,     # s1/s2 비율 (선택사항)
    N=4                  # 튜브 열 수 (선택사항)
)

# 결과 출력
print(f"1/K = {result['porous']['inv_K']:.4e} [1/m²]")
print(f"C2 = {result['porous']['C2']:.4f} [1/m]")
print(f"h_fs = {result['design_point']['h_fs_W_m2K']:.2f} [W/(m²·K)]")
```

---

## 출력 데이터 구조

`calculate_porous_parameters()` 함수는 다음 구조의 딕셔너리를 반환합니다:

```python
{
    'input': {
        'Fs_mm': 4.0,
        'hf_mm': 4.0,
        'T_celsius': 14.80177,
        'v_design': 2.019723,
        # ... 기타 입력값
    },
    'air': {
        'rho': 1.2258,          # 밀도 [kg/m³]
        'mu': 1.788e-05,        # 점도 [Pa·s]
        'k': 0.0251,            # 열전도도 [W/(m·K)]
        'Pr': 0.715             # Prandtl 수
    },
    'geometry': {
        'epsilon': 0.8889,      # 공극률
        'sigma': 0.5502,        # 최소유동면적비
        'area_ratio': 3.11,     # 면적비 (AR)
        'a_fs': 698.85          # 비표면적 [1/m]
    },
    'porous': {
        'inv_K': 6.59e+04,      # 점성 저항 [1/m²]  ← CFD 입력
        'C2': 5.37,             # 관성 저항 [1/m]   ← CFD 입력
        'K': 2.71e-10           # 투과도 [m²]
    },
    'design_point': {
        'Re_Dc': 6039,          # Reynolds 수
        'dP_total_Pa': 11.95,   # 총 압력강하 [Pa]
        'h_fs_W_m2K': 94.82     # 열전달 계수 [W/(m²·K)]  ← CFD 입력
    }
}
```

---

## 주요 파라미터 의미

| 파라미터 | 기호 | 의미 | CFD 사용처 |
|---------|------|------|-----------|
| 점성 저항 | 1/K | 저속 유동 저항 (점성 지배) | Viscous Resistance |
| 관성 저항 | C₂ | 고속 유동 저항 (관성 지배) | Inertial Resistance |
| 투과도 | K | 다공매체 투과능력 (1/K의 역수) | Permeability |
| 공극률 | ε | 유체가 차지하는 체적 비율 | Porosity |
| 비표면적 | a_fs | 단위 체적당 표면적 | Surface Area Density |
| 열전달 계수 | h_fs | 고체-유체 계면 열전달 | LTNE Heat Transfer Coeff |

---

## 상관식 및 이론

본 계산기는 다음 상관식을 사용합니다:

### 압력강하 (Nir, 1991)

$$f_N = 1.1 \cdot Re_{D_c}^{-0.25} \cdot \left(\frac{S_1}{D_c}\right)^{-0.4} \cdot AR^{0.15}$$

- **적용 범위:** Staggered Annular Finned Tube Banks
- **Reynolds 수:** 300 ~ 10,000 (수력직경 기준)

### 열전달 (Briggs & Young, 1963)

$$Nu = 0.134 \cdot Re^{0.681} \cdot Pr^{1/3} \cdot \left(\frac{F_s}{h_f}\right)^{0.2} \cdot \left(\frac{F_s}{\delta_f}\right)^{0.113}$$

### Darcy-Forchheimer 변환 (2-Point Fitting)

$$\frac{\Delta P}{L} = \frac{\mu}{K} \cdot v + \frac{C_2 \rho}{2} \cdot v^2$$

---

## 문제 해결

### Q1: 계산 결과가 이상합니다

**확인 사항:**
- 핀 간격 (Fs)과 핀 높이 (hf)가 현실적인 값인가?
- 핀 외경 (Do = Dc + 2*hf)이 피치 (s1, s2)보다 작은가?
- 온도와 풍속이 합리적인 범위인가?

### Q2: Fluent에서 압력강하가 다릅니다

**확인 사항:**
- Viscous/Inertial Resistance의 방향 설정 확인
  - Direction-1: 유동 방향
  - Direction-2, 3: 횡방향 (1000배 큰 값)
- Porosity 값이 모든 방향에 동일하게 입력되었는지 확인

### Q3: LTNE 모델에서 온도 분포가 이상합니다

**확인 사항:**
- Surface Area Density (a_fs) 값 확인
- Fluid-Solid Heat Transfer Coefficient (h_fs) 값 확인
  - 일반적으로 50~500 W/(m²·K) 범위
- 고체 재질의 열전도도가 올바르게 설정되었는지 확인

---

## 참고 자료

### 상세 수식 유도
README.md 파일을 참조하세요. Nir 상관식부터 Porous 파라미터 도출까지 전 과정이 수식으로 설명되어 있습니다.

### 참고 문헌
1. **Nir, A.** (1991). "Heat Transfer and Friction Factor Correlations for Crossflow over Staggered Finned Tube Banks", *Heat Transfer Engineering*, 12(1), 43-58.

2. **Briggs, D.E., & Young, E.H.** (1963). "Convection heat transfer and pressure drop of air flowing across triangular pitch banks of finned tubes."

---

## 라이센스

이 프로젝트는 교육 및 연구 목적으로 제공됩니다.

**버전:** 2.0 (Nir 상관식 기반)
**최종 업데이트:** 2026-01-17
**작성자:** Claude (Anthropic)
