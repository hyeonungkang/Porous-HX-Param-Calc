# CFD Porous Zone 파라미터 계산 가이드

## 개요

이 문서는 **핀-튜브 열교환기**의 Porous Zone 설정을 위한 CFD 파라미터 계산 방법과 Fluent 입력 가이드를 제공합니다.

**적용 조건:**
- **유체:** 15°C 공기, 1 atm
- **설계 유속:** 1.5 m/s (Frontal Velocity)
- **열전달 모델:** LTNE (Local Thermal Non-Equilibrium, 비평형 열전달)

---

## 1. 입력 변수 및 고정 파라미터

### 1.1 설계 변수 (Design Variables)

| 파라미터 | 기호 | 단위 | 설명 |
|---------|------|------|------|
| 핀 간격 | $F_s$ | mm | Fin Spacing |
| 핀 높이 | $h_f$ | mm | Fin Height |
| 피치 비율 | $Ratio = s_1/s_2$ | - | 종방향/횡방향 피치 비율 (0.5 < Ratio < 2.0) |

### 1.2 고정 파라미터 (Fixed Parameters)

| 파라미터 | 기호 | 값 | 단위 | 비고 |
|---------|------|-----|------|------|
| 튜브 외경 | $D_c$ | 24 | mm | Tube OD |
| 핀 두께 | $\delta_f$ | 0.5 | mm | - |
| 횡방향 피치 | $s_1$ | 55.333 | mm | Spanwise Pitch |
| 튜브 열 수 | $N$ | 4 | - | 유동 방향 튜브 열 수 (단위 크기 0.2m 기준) |
| 배열 형태 | - | Staggered | - | 엇갈림 배열 (열전달 효율 우수) |

### 1.3 종속 변수 (Derived Variables)

| 파라미터 | 기호 | 수식 | 단위 |
|---------|------|------|------|
| 핀 외경 | $D_o$ | $D_c + 2h_f$ | mm |
| 핀 피치 | $F_p$ | $F_s + \delta_f$ | mm | Center-to-Center |
| 종방향 피치 | $s_2$ | $s_1 / Ratio$ | mm |

### 1.4 공기 물성치 (15°C, 1 atm 기준)

**명료화된 가정:**

| 물성치 | 기호 | 값 | 단위 |
|--------|------|-----|------|
| 밀도 | $\rho$ | 1.225 | kg/m³ |
| 점성 | $\mu$ | 1.789 × 10⁻⁵ | Pa·s |
| 열전도율 | $k$ | 0.0253 | W/m·K |
| Prandtl 수 | $Pr$ | 0.71 | - |

**중요:** 코드 내부에서 이 값들이 **고정값으로 하드코딩**되어 있습니다. 다른 온도 조건을 사용할 경우 코드 수정이 필요합니다.

---

## 2. 계산 결과값 (Output Parameters)

코드 실행 시 다음 5가지 주요 파라미터가 출력됩니다:

### 2.1 공극률 (Porosity)

**기호:** $\epsilon$  
**단위:** - (무차원)  
**수식:**
$$\epsilon = 1 - \frac{\delta_f}{F_p}$$

**의미:** Annular Zone 내에서 공기가 흐를 수 있는 체적 비율 (핀 고체 부피 제외)

### 2.2 비표면적 (Surface Area Density)

**기호:** $a_{fs}$  
**단위:** m⁻¹  
**수식:**
$$a_{fs} = \frac{A_{surface}}{V_{annular}} = \frac{0.5\pi(D_o^2 - D_c^2) + \pi D_o \delta_f + \pi D_c F_s}{0.25\pi(D_o^2 - D_c^2) \cdot F_p}$$

**의미:** 단위 체적당 열전달 표면적

### 2.3 점성 저항 (Viscous Resistance)

**기호:** $1/K$  
**단위:** m⁻²  
**계산 방법:** 2-Point Fitting Method를 통한 Darcy-Forchheimer 계수 변환

**의미:** 유동 저항의 점성 성분

### 2.4 관성 저항 (Inertial Resistance)

**기호:** $C_2$  
**단위:** m⁻¹  
**계산 방법:** 2-Point Fitting Method를 통한 Darcy-Forchheimer 계수 변환

**의미:** 유동 저항의 관성 성분

### 2.5 계면 열전달 계수 (Interfacial Heat Transfer Coefficient)

**기호:** $h_{fs}$  
**단위:** W/m²·K  
**상관식:** Briggs & Young Correlation (Annular Finned Tubes 전용)

$$Nu = 0.134 \cdot Re_D^{0.681} \cdot Pr^{0.33} \cdot \left(\frac{s}{l}\right)^{0.2} \cdot \left(\frac{s}{t}\right)^{0.113}$$

$$h_{fs} = \frac{Nu \cdot k}{D_c}$$

**의미:** 핀 표면에서 공기로의 열전달 계수 (LTNE 모델 필수)

---

## 3. 계산 프로세스 (Calculation Flow)

### 3.1 기하학적 파라미터 계산

1. **종속 변수 계산**
   - $F_p = F_s + \delta_f$
   - $D_o = D_c + 2h_f$
   - $s_2 = s_1 / Ratio$

2. **안전성 검사 (Overlap Check)**
   - `if D_o > s_1`: 횡방향 충돌 경고
   - `if D_o > s_2`: 종방향 충돌 경고
   - 충돌 시 오류 메시지 출력 및 계산 중단

3. **Porous Zone 파라미터**
   - $\epsilon$, $a_{fs}$, $D_h$ 계산

### 3.2 유동 저항 계산 (Wang Correlation)

**Step 1: 최소 유동 단면적 비율**
$$\sigma = \frac{A_{min}}{A_{face}} = \frac{s_1 - D_c - 2h_f \delta_f / F_p}{s_1}$$

**Step 2: 최대 유속 (Physical Velocity)**
$$v_{max} = \frac{v_{inlet}}{\sigma}$$

**Step 3: 레이놀즈 수**
$$Re_{Dc} = \frac{\rho \cdot v_{max} \cdot D_c}{\mu}$$

**Step 4: Wang 마찰 계수**
$$f = 0.0267 \cdot Re_{Dc}^{F_1} \cdot \left(\frac{s_1}{s_2}\right)^{F_2} \cdot \left(\frac{F_s}{D_c}\right)^{F_3}$$

여기서:
- $F_1 = -0.764 + 0.739(s_1/s_2) + 0.177(F_p/D_c) - \frac{0.00758}{N}$
- $F_2 = -15.689 + \frac{64.021}{\ln(Re_{Dc})}$
- $F_3 = 1.696 - \frac{15.695}{\ln(Re_{Dc})}$

**Step 5: 압력강하**
$$\frac{\Delta P}{L} = \frac{f}{D_h} \cdot \frac{\rho v_{max}^2}{2}$$

**Step 6: 2-Point Fitting (Darcy-Forchheimer 변환)**

두 가지 속도 포인트 ($v_1 = 1.5$ m/s, $v_2 = 0.5$ m/s)에서 계산된 압력강하값을 사용하여:

$$\frac{\Delta P}{L} = A \cdot v + B \cdot v^2$$

형태로 변환하고, 다음 관계식으로 CFD 계수를 도출:

- $1/K = \frac{A}{\mu}$ (Viscous Resistance)
- $C_2 = \frac{2B}{\rho}$ (Inertial Resistance)

### 3.3 열전달 계수 계산 (Briggs & Young)

- 설계점 레이놀즈 수 ($Re_{design}$) 사용
- Nusselt 수 계산 후 $h_{fs}$ 도출

---

## 4. Fluent 입력 가이드

### 4.1 Porous Zone 설정

**위치:** `Cell Zone Conditions` > `Porous Zone`

#### 4.1.1 Porosity

- **입력 위치:** `Porosity`
- **값:** 계산 결과의 `Porosity (ε)`
- **단위:** 무차원
- **입력 방법:** 
  ```
  Direction-1: [계산값]
  Direction-2: [계산값] (동일)
  Direction-3: [계산값] (동일)
  ```

#### 4.1.2 Viscous Resistance (1/K)

- **입력 위치:** `Viscous Resistance`
- **값:** 계산 결과의 `Viscous Resistance (1/K)`
- **단위:** m⁻²
- **입력 방법:**
  ```
  Direction-1: [계산값] (주유동 방향)
  Direction-2: [계산값 × 1000] (큰 값을 넣어 유동 차단)
  Direction-3: [계산값 × 1000] (큰 값을 넣어 유동 차단)
  ```
  
  **주의:** Direction-2, Direction-3에는 **주유동 방향 값의 1,000배**를 입력하여 유동이 옆으로 새지 않게 막아야 합니다.

#### 4.1.3 Inertial Resistance (C2)

- **입력 위치:** `Inertial Resistance`
- **값:** 계산 결과의 `Inertial Resistance (C2)`
- **단위:** m⁻¹
- **입력 방법:**
  ```
  Direction-1: [계산값] (주유동 방향)
  Direction-2: [계산값 × 1000] (큰 값을 넣어 유동 차단)
  Direction-3: [계산값 × 1000] (큰 값을 넣어 유동 차단)
  ```

### 4.2 LTNE (Local Thermal Non-Equilibrium) 설정

**위치:** `Cell Zone Conditions` > `Porous Zone` > `Non-Equilibrium Thermal Model`

#### 4.2.1 Non-Equilibrium Thermal Model 활성화

- `Enable Non-Equilibrium Thermal Model` 체크

#### 4.2.2 Surface Area Density

- **입력 위치:** `Surface Area Density`
- **값:** 계산 결과의 `Surface Area Density (a_fs)`
- **단위:** m⁻¹
- **입력 방법:** 계산 결과값 그대로 입력

#### 4.2.3 Fluid-Solid Heat Transfer Coefficient

- **입력 위치:** `Fluid Solid Heat Transfer Coefficient`
- **값:** 계산 결과의 `Interfacial Heat Transfer Coeff (h_fs)`
- **단위:** W/m²·K
- **의미:** 핀 표면에서 공기로의 열전달 계수
- **중요:** 이 값이 있어야 핀 효율(Fin Efficiency)과 핀 온도 분포가 정확히 계산됩니다.

---

## 5. 사용 예시

### 5.1 단일 케이스 계산

```python
# cfd_porous_calc.py 파일에서 다음 변수들을 수정:
Fs_input = 2.0      # 핀 간격 [mm]
hf_input = 4.0      # 핀 높이 [mm]
ratio_input = 1.0   # 피치 비율 (s1/s2)

# 실행:
# python cfd_porous_calc.py
```

### 5.2 출력 예시

```
============================================================================================================================
CFD Porous Zone 파라미터 계산기
15°C 공기, 1.5m/s 유속, LTNE 조건 기준
============================================================================================================================
입력값:
  - 핀 간격 (Fs): 2.0 mm
  - 핀 높이 (hf): 4.0 mm
  - 피치 비율 (Ratio = s1/s2): 1.00
----------------------------------------------------------------------------------------------------------------------------
계산 결과 (Fluent 입력값):
  1. Porosity (ε): 0.800000 [-]
  2. Surface Area Density (a_fs): 1234.56 [1/m]
  3. Viscous Resistance (1/K): 1.234567e+06 [1/m²]
  4. Inertial Resistance (C2): 12.345678 [1/m]
  5. Interfacial Heat Transfer Coeff (h_fs): 123.45 [W/m²·K]
----------------------------------------------------------------------------------------------------------------------------
참고값:
  - 종방향 피치 (s2): 55.33 mm
  - 설계점 레이놀즈 수 (Re_Dc): 1234 [-]
  - Nusselt 수 (Nu): 12.34 [-]
  - 수력 직경 (Dh): 2.345 mm
  - 최소 면적 비율 (σ): 0.5678 [-]
============================================================================================================================
```

### 5.3 Fluent 입력 예시

**Porous Zone 설정:**
```
Porosity:
  Direction-1: 0.800000
  Direction-2: 0.800000
  Direction-3: 0.800000

Viscous Resistance:
  Direction-1: 1.234567e+06
  Direction-2: 1.234567e+09
  Direction-3: 1.234567e+09

Inertial Resistance:
  Direction-1: 12.345678
  Direction-2: 12345.678
  Direction-3: 12345.678
```

**LTNE 설정:**
```
Enable Non-Equilibrium Thermal Model: ✓

Surface Area Density: 1234.56 [1/m]

Fluid Solid Heat Transfer Coefficient: 123.45 [W/m²·K]
```

---

## 6. 주의사항 및 제한사항

### 6.1 물성치 고정

- 현재 코드는 **15°C 공기 물성치로 고정**되어 있습니다.
- 다른 온도 조건을 사용할 경우 코드 내부 물성치 값을 수정해야 합니다.

### 6.2 피치 비율 제한

- **0.5 < Ratio < 2.0** 범위를 벗어나면 계산이 중단됩니다.
- 이는 Wang 상관식의 유효 범위에 따른 제한입니다.

### 6.3 기하학적 제약

- **Overlap Check:** $D_o > s_1$ 또는 $D_o > s_2$인 경우 충돌 경고가 발생합니다.
- 설계 변수를 변경할 때는 항상 간섭 여부를 확인하세요.

### 6.4 배열 형태 가정

- 현재 코드는 **Staggered 배열**을 가정하고 있습니다.
- Inline 배열을 사용할 경우 Wang 상관식 계수나 열전달 상관식을 수정해야 할 수 있습니다.

### 6.5 Wang 상관식 적용 범위

- Wang 상관식은 평판 핀(Plain Fin) 기하에 대해 개발되었습니다.
- Annular Fin에 적용 시 근사치로 사용되며, 실험 데이터와의 검증이 권장됩니다.

---

## 7. 참고 문헌 및 상관식

### 7.1 압력강하 상관식

- **Wang, C.C., et al.** (1997). "Heat transfer and friction characteristics of fin-and-tube heat exchangers."
  - 마찰 계수: $f = 0.0267 \cdot Re_{Dc}^{F_1} \cdot (s_1/s_2)^{F_2} \cdot (F_s/D_c)^{F_3}$

### 7.2 열전달 상관식

- **Briggs, D.E., & Young, E.H.** (1963). "Convection heat transfer and pressure drop of air flowing across triangular pitch banks of finned tubes."
  - Nusselt 수: $Nu = 0.134 \cdot Re_D^{0.681} \cdot Pr^{0.33} \cdot (s/l)^{0.2} \cdot (s/t)^{0.113}$

---

## 8. 문제 해결 (Troubleshooting)

### 8.1 "COLLISION!" 오류

**원인:** 핀 외경이 피치보다 큼  
**해결:** $h_f$를 줄이거나 $F_s$를 늘려서 $D_o$를 감소시키세요.

### 8.2 "Invalid pitch_ratio" 오류

**원인:** Ratio가 0.5~2.0 범위를 벗어남  
**해결:** Ratio 값을 0.5와 2.0 사이로 조정하세요.

### 8.3 Fluent에서 압력강하가 예상과 다름

**확인 사항:**
1. Viscous/Inertial Resistance의 Direction 설정 확인
2. Porosity 값 확인
3. Physical Velocity (v_max) 대비 계산된 저항 계수 검증

### 8.4 LTNE 모델에서 열전달이 비현실적

**확인 사항:**
1. $h_{fs}$ 값이 합리적인 범위 내인지 확인 (일반적으로 50~500 W/m²·K)
2. Surface Area Density 값 확인
3. 고체 열전도율 설정 확인 (Fluent Materials)

---

## 9. 연락처 및 지원

코드 오류나 개선 사항이 있으면 이슈를 등록하거나 코드를 수정해주세요.

---

**버전:** 1.0  
**최종 업데이트:** 2024  
**작성자:** CFD Porous Zone Calculator
