#!/usr/bin/env python3
"""
================================================================================
Annular Fin → Porous Media 파라미터 변환 계산기
================================================================================
Nir (1991) 마찰계수 상관식을 사용하여 CFD Porous Media의
점성 저항 계수(1/K)와 관성 저항 계수(C2)를 도출합니다.

Reference:
    Nir, A. (1991). "Heat Transfer and Friction Factor Correlations for
    Crossflow over Staggered Finned Tube Banks", Heat Transfer Engineering,
    Vol.12, No.1, pp.43-58.

Author: Claude (Anthropic)
Date: 2026-01-17
================================================================================
"""

import math
import numpy as np
from dataclasses import dataclass
from typing import Tuple, Dict, List

# =============================================================================
# 1. 물성치 계산 함수
# =============================================================================

def air_properties(T_celsius: float) -> Dict[str, float]:
    """
    주어진 온도에서 공기의 물성치를 계산합니다.

    Parameters:
        T_celsius: 온도 [°C]

    Returns:
        dict: 공기 물성치
            - rho: 밀도 [kg/m³]
            - mu: 동점성계수 [Pa·s]
            - k: 열전도도 [W/(m·K)]
            - Pr: Prandtl 수 [-]
            - nu: 동점성계수 [m²/s]
    """
    T_K = T_celsius + 273.15
    P = 101325  # 대기압 [Pa]
    R = 287.058  # 공기 기체상수 [J/(kg·K)]

    # 밀도 (이상기체 상태방정식)
    rho = P / (R * T_K)

    # 동점성계수 (Sutherland's Law)
    # μ = μ_ref * (T/T_ref)^1.5 * (T_ref + S) / (T + S)
    mu_ref = 1.716e-5  # [Pa·s] at 273.15K
    T_ref = 273.15
    S = 110.4  # Sutherland 상수 [K]
    mu = mu_ref * (T_K / T_ref)**1.5 * (T_ref + S) / (T_K + S)

    # 열전도도 (근사식)
    k = 0.0241 + 7.0e-5 * T_celsius  # 선형 근사

    # Prandtl 수
    cp = 1006  # 정압비열 [J/(kg·K)]
    Pr = mu * cp / k

    # 동점성계수
    nu = mu / rho

    return {
        'rho': rho,
        'mu': mu,
        'k': k,
        'Pr': Pr,
        'nu': nu,
        'T_K': T_K
    }


# =============================================================================
# 2. 기하학적 파라미터 계산
# =============================================================================

@dataclass
class AnnularFinGeometry:
    """Annular Fin 기하학적 파라미터를 저장하는 데이터 클래스"""
    # 입력 파라미터
    Dc: float      # 튜브 외경 [m]
    delta_f: float # 핀 두께 [m]
    s1: float      # 횡방향 피치 (Spanwise) [m]
    Fs: float      # 핀 간격 (Fin Spacing) [m]
    hf: float      # 핀 높이 [m]
    pitch_ratio: float = 1.0  # s1/s2 비율
    N: int = 4     # 튜브 열 수

    def __post_init__(self):
        """종속 파라미터 계산"""
        # 핀 외경
        self.Do = self.Dc + 2 * self.hf

        # 핀 피치 (Center-to-Center)
        self.Fp = self.Fs + self.delta_f

        # 종방향 피치
        self.s2 = self.s1 / self.pitch_ratio

        # 공극률 (Annular Zone 기준)
        self.epsilon = 1 - (self.delta_f / self.Fp)

        # 최소 자유유동 면적비 (sigma)
        # σ = (S1 - Dc - 2*hf*(δf/Fp)) / S1
        self.sigma = (self.s1 - self.Dc - 2*self.hf*(self.delta_f/self.Fp)) / self.s1

        # 면적 계산 (1 피치 구간)
        self._calculate_areas()

    def _calculate_areas(self):
        """1 피치 구간의 면적 계산"""
        # 핀 표면적 (양면 + 팁)
        # A_fin = 2 × π/4 × (Do² - Dc²) + π × Do × δf
        self.A_fin = (2 * (math.pi/4) * (self.Do**2 - self.Dc**2) +
                      math.pi * self.Do * self.delta_f)

        # 튜브 노출 면적 (핀 사이)
        self.A_base = math.pi * self.Dc * self.Fs

        # 총 표면적
        self.A_total = self.A_fin + self.A_base

        # 맨 튜브 표면적 (핀이 없을 때)
        self.A_bare = math.pi * self.Dc * self.Fp

        # 면적비 (Area Ratio) - Nir 상관식의 핵심 인자
        self.area_ratio = self.A_total / self.A_bare

        # Annular Zone 체적 (1 피치 구간)
        self.V_annular = (math.pi/4) * (self.Do**2 - self.Dc**2) * self.Fp

        # 비표면적 [1/m]
        self.a_fs = self.A_total / self.V_annular

    def get_porous_thickness(self) -> float:
        """Porous Zone의 두께 (핀 높이) 반환"""
        return self.hf

    def get_flow_depth(self) -> float:
        """유동 방향 총 깊이 반환"""
        return self.s2 * self.N


# =============================================================================
# 3. Nir (1991) 마찰계수 상관식
# =============================================================================

def nir_friction_factor(Re_Dc: float, s1_Dc: float, area_ratio: float) -> float:
    """
    Nir (1991) 마찰계수 상관식

    f_N = 1.1 × Re^(-0.25) × (S1/Dc)^(-0.4) × (A_tot/A_bare)^(0.15)

    Parameters:
        Re_Dc: 튜브 외경 기준 Reynolds 수 (v_max 기반)
        s1_Dc: 횡방향 피치 / 튜브 외경 [-]
        area_ratio: 면적비 (A_total / A_bare) [-]

    Returns:
        f_N: Nir 마찰계수 [-]

    적용 범위:
        - Re_h = 300 ~ 10,000 (수력직경 기준)
        - Staggered 배열
    """
    term1 = 1.1 * (Re_Dc ** -0.25)
    term2 = (s1_Dc) ** -0.4
    term3 = (area_ratio) ** 0.15

    f_N = term1 * term2 * term3
    return f_N


def calculate_pressure_drop(v_inlet: float, geom: AnnularFinGeometry,
                           air_props: Dict[str, float]) -> Tuple[float, float, float]:
    """
    Nir 상관식을 이용한 압력강하 계산

    ΔP = N × f_N × (ρ × v_max²) / 2

    Parameters:
        v_inlet: 입구 속도 (superficial velocity) [m/s]
        geom: AnnularFinGeometry 객체
        air_props: 공기 물성치 dict

    Returns:
        tuple: (ΔP_total [Pa], ΔP/L [Pa/m], Re_Dc [-])
    """
    rho = air_props['rho']
    mu = air_props['mu']

    # 최대 속도 (최소 단면적 기준)
    v_max = v_inlet / geom.sigma

    # Reynolds 수 (튜브 외경 기준, v_max 사용)
    Re_Dc = (rho * v_max * geom.Dc) / mu

    # Nir 마찰계수
    s1_Dc = geom.s1 / geom.Dc
    f_N = nir_friction_factor(Re_Dc, s1_Dc, geom.area_ratio)

    # 총 압력강하
    dP_total = geom.N * f_N * (rho * v_max**2) / 2

    # 단위 길이당 압력강하
    L = geom.get_flow_depth()
    dP_per_L = dP_total / L

    return dP_total, dP_per_L, Re_Dc


# =============================================================================
# 4. 2-Point Fitting을 통한 Porous 파라미터 도출
# =============================================================================

def fit_darcy_forchheimer(v1: float, v2: float,
                          geom: AnnularFinGeometry,
                          air_props: Dict[str, float]) -> Dict[str, float]:
    """
    2-Point Fitting을 통한 Darcy-Forchheimer 계수 도출

    Nir 상관식: f_N ∝ Re^(-0.25) → ΔP/L ∝ v^1.75 (정확한 Darcy-Forchheimer 형태가 아님)

    따라서 두 속도점에서 계산 후, ΔP/L = A×v + B×v² 형태로 피팅

    수식 유도:
        Y1 = A×v1 + B×v1²
        Y2 = A×v2 + B×v2²

        Z1 = Y1/v1 = A + B×v1
        Z2 = Y2/v2 = A + B×v2

        B = (Z1 - Z2) / (v1 - v2)
        A = Z1 - B×v1

    Parameters:
        v1, v2: 피팅에 사용할 두 속도점 [m/s]
        geom: AnnularFinGeometry 객체
        air_props: 공기 물성치 dict

    Returns:
        dict: Porous 파라미터
            - A: 선형 계수 [Pa·s/m²]
            - B: 2차 계수 [Pa·s²/m³]
            - inv_K: 점성 저항 계수 1/K [1/m²]
            - C2: 관성 저항 계수 [1/m]
            - K: 투과도 [m²]
    """
    rho = air_props['rho']
    mu = air_props['mu']

    # 두 속도점에서 압력강하 계산
    _, Y1, Re1 = calculate_pressure_drop(v1, geom, air_props)
    _, Y2, Re2 = calculate_pressure_drop(v2, geom, air_props)

    # 선형화를 위한 변환
    Z1 = Y1 / v1
    Z2 = Y2 / v2

    # 계수 계산
    B = (Z1 - Z2) / (v1 - v2)
    A = Z1 - B * v1

    # Porous 파라미터 변환
    # ΔP/L = A×v + B×v² = (μ/K)×v + (C2×ρ/2)×v²

    inv_K = A / mu        # 점성 저항 계수 [1/m²]
    K = mu / A            # 투과도 [m²]
    C2 = 2 * B / rho      # 관성 저항 계수 [1/m]

    return {
        'A': A,
        'B': B,
        'inv_K': inv_K,      # Viscous Resistance (Fluent: 1/α)
        'C2': C2,            # Inertial Resistance (Fluent: C2)
        'K': K,              # Permeability
        'Re_low': Re2,
        'Re_high': Re1,
        'v1': v1,
        'v2': v2,
        'Y1': Y1,
        'Y2': Y2
    }


# =============================================================================
# 5. 열전달 계수 (Briggs & Young 상관식)
# =============================================================================

def briggs_young_heat_transfer(v_inlet: float, geom: AnnularFinGeometry,
                                air_props: Dict[str, float]) -> float:
    """
    Briggs & Young (1963) 열전달 상관식

    Nu = 0.134 × Re^0.681 × Pr^(1/3) × (Fs/hf)^0.2 × (Fs/δf)^0.113

    Parameters:
        v_inlet: 입구 속도 [m/s]
        geom: AnnularFinGeometry 객체
        air_props: 공기 물성치 dict

    Returns:
        h_fs: 열전달 계수 [W/(m²·K)]
    """
    rho = air_props['rho']
    mu = air_props['mu']
    k = air_props['k']
    Pr = air_props['Pr']

    v_max = v_inlet / geom.sigma
    Re_Dc = (rho * v_max * geom.Dc) / mu

    Nu = (0.134 * (Re_Dc ** 0.681) * (Pr ** (1/3)) *
          ((geom.Fs / geom.hf) ** 0.2) *
          ((geom.Fs / geom.delta_f) ** 0.113))

    h_fs = Nu * k / geom.Dc

    return h_fs


# =============================================================================
# 6. 메인 계산 함수
# =============================================================================

def calculate_porous_parameters(
    Fs_mm: float,           # 핀 간격 [mm]
    hf_mm: float,           # 핀 높이 [mm]
    T_celsius: float,       # 대기 온도 [°C]
    v_design: float,        # 설계 풍속 [m/s]
    Dc_mm: float = 24.0,    # 튜브 외경 [mm]
    delta_f_mm: float = 0.5,# 핀 두께 [mm]
    s1_mm: float = 55.333,  # 횡방향 피치 [mm]
    pitch_ratio: float = 1.0,# s1/s2 비율
    N: int = 4              # 튜브 열 수
) -> Dict:
    """
    Annular Fin을 Porous Media로 근사하기 위한 모든 파라미터 계산

    Parameters:
        Fs_mm: 핀 간격 [mm]
        hf_mm: 핀 높이 [mm]
        T_celsius: 대기 온도 [°C]
        v_design: 설계 풍속 [m/s]
        Dc_mm: 튜브 외경 [mm] (기본값: 24)
        delta_f_mm: 핀 두께 [mm] (기본값: 0.5)
        s1_mm: 횡방향 피치 [mm] (기본값: 55.333)
        pitch_ratio: s1/s2 비율 (기본값: 1.0)
        N: 튜브 열 수 (기본값: 4)

    Returns:
        dict: 모든 계산 결과
    """
    # 단위 변환 (mm → m)
    Dc = Dc_mm / 1000
    delta_f = delta_f_mm / 1000
    s1 = s1_mm / 1000
    Fs = Fs_mm / 1000
    hf = hf_mm / 1000

    # 공기 물성치
    air_props = air_properties(T_celsius)

    # 기하학적 파라미터
    geom = AnnularFinGeometry(
        Dc=Dc, delta_f=delta_f, s1=s1, Fs=Fs, hf=hf,
        pitch_ratio=pitch_ratio, N=N
    )

    # 2-Point Fitting (설계 속도의 1.0배와 0.3배)
    v_high = v_design
    v_low = v_design * 0.3

    porous_params = fit_darcy_forchheimer(v_high, v_low, geom, air_props)

    # 설계점에서의 압력강하 및 Reynolds 수
    dP_total, dP_per_L, Re_design = calculate_pressure_drop(v_design, geom, air_props)

    # 열전달 계수
    h_fs = briggs_young_heat_transfer(v_design, geom, air_props)

    return {
        # 입력 조건
        'input': {
            'Fs_mm': Fs_mm,
            'hf_mm': hf_mm,
            'T_celsius': T_celsius,
            'v_design': v_design,
            'Dc_mm': Dc_mm,
            'delta_f_mm': delta_f_mm,
            's1_mm': s1_mm,
            'pitch_ratio': pitch_ratio,
            'N': N
        },
        # 공기 물성치
        'air': air_props,
        # 기하학적 파라미터
        'geometry': {
            'Do_mm': geom.Do * 1000,
            'Fp_mm': geom.Fp * 1000,
            's2_mm': geom.s2 * 1000,
            'epsilon': geom.epsilon,
            'sigma': geom.sigma,
            'area_ratio': geom.area_ratio,
            'a_fs': geom.a_fs,
            'porous_thickness_mm': geom.get_porous_thickness() * 1000,
            'flow_depth_mm': geom.get_flow_depth() * 1000
        },
        # Porous 파라미터 (핵심 출력)
        'porous': {
            'inv_K': porous_params['inv_K'],  # [1/m²]
            'C2': porous_params['C2'],         # [1/m]
            'K': porous_params['K'],           # [m²]
            'A': porous_params['A'],           # [Pa·s/m²]
            'B': porous_params['B']            # [Pa·s²/m³]
        },
        # 설계점 결과
        'design_point': {
            'Re_Dc': Re_design,
            'dP_total_Pa': dP_total,
            'dP_per_L_Pa_m': dP_per_L,
            'h_fs_W_m2K': h_fs
        },
        # 피팅 정보
        'fitting': {
            'v_high': v_high,
            'v_low': v_low,
            'Re_high': porous_params['Re_high'],
            'Re_low': porous_params['Re_low']
        }
    }


# =============================================================================
# 7. 결과 출력 함수
# =============================================================================

def print_results(result: Dict):
    """계산 결과를 보기 좋게 출력"""
    print("=" * 70)
    print("       Annular Fin → Porous Media 파라미터 변환 결과")
    print("=" * 70)

    inp = result['input']
    print(f"\n[입력 조건]")
    print(f"  핀 간격 (Fs)        : {inp['Fs_mm']:.1f} mm")
    print(f"  핀 높이 (hf)        : {inp['hf_mm']:.1f} mm")
    print(f"  대기 온도           : {inp['T_celsius']:.2f} °C")
    print(f"  설계 풍속           : {inp['v_design']:.4f} m/s")
    print(f"  튜브 외경 (Dc)      : {inp['Dc_mm']:.1f} mm")
    print(f"  횡방향 피치 (s1)    : {inp['s1_mm']:.3f} mm")
    print(f"  튜브 열 수 (N)      : {inp['N']}")

    air = result['air']
    print(f"\n[공기 물성치 @ {inp['T_celsius']:.2f}°C]")
    print(f"  밀도 (ρ)            : {air['rho']:.4f} kg/m³")
    print(f"  점도 (μ)            : {air['mu']:.4e} Pa·s")
    print(f"  열전도도 (k)        : {air['k']:.4f} W/(m·K)")
    print(f"  Prandtl 수 (Pr)     : {air['Pr']:.3f}")

    geom = result['geometry']
    print(f"\n[기하학적 파라미터]")
    print(f"  핀 외경 (Do)        : {geom['Do_mm']:.1f} mm")
    print(f"  핀 피치 (Fp)        : {geom['Fp_mm']:.2f} mm")
    print(f"  종방향 피치 (s2)    : {geom['s2_mm']:.3f} mm")
    print(f"  공극률 (ε)          : {geom['epsilon']:.4f}")
    print(f"  최소유동면적비 (σ)  : {geom['sigma']:.4f}")
    print(f"  면적비 (Atot/Abare) : {geom['area_ratio']:.3f}")
    print(f"  비표면적 (a_fs)     : {geom['a_fs']:.2f} 1/m")

    porous = result['porous']
    print(f"\n{'='*70}")
    print(f"  ★★★ CFD Porous Media 입력값 (핵심 출력) ★★★")
    print(f"{'='*70}")
    print(f"  점성 저항 (1/K)     : {porous['inv_K']:.4e}  [1/m²]")
    print(f"  관성 저항 (C2)      : {porous['C2']:.4f}  [1/m]")
    print(f"  투과도 (K)          : {porous['K']:.4e}  [m²]")
    print(f"{'='*70}")

    design = result['design_point']
    print(f"\n[설계점 성능]")
    print(f"  Reynolds 수 (Re_Dc) : {design['Re_Dc']:.1f}")
    print(f"  총 압력강하         : {design['dP_total_Pa']:.2f} Pa")
    print(f"  단위길이 압력강하   : {design['dP_per_L_Pa_m']:.2f} Pa/m")
    print(f"  열전달 계수 (h_fs)  : {design['h_fs_W_m2K']:.2f} W/(m²·K)")

    fit = result['fitting']
    print(f"\n[2-Point Fitting 정보]")
    print(f"  고속점: v = {fit['v_high']:.4f} m/s, Re = {fit['Re_high']:.1f}")
    print(f"  저속점: v = {fit['v_low']:.4f} m/s, Re = {fit['Re_low']:.1f}")

    print("\n" + "=" * 70)


def generate_parameter_table(cases: List[Tuple[float, float]],
                             T_celsius: float, v_design: float) -> None:
    """
    여러 케이스에 대한 파라미터 테이블 생성

    Parameters:
        cases: [(Fs_mm, hf_mm), ...] 리스트
        T_celsius: 대기 온도 [°C]
        v_design: 설계 풍속 [m/s]
    """
    print("\n" + "=" * 100)
    print(f"  Porous 파라미터 비교표 (T = {T_celsius:.2f}°C, v = {v_design:.4f} m/s)")
    print("=" * 100)

    header = f"{'Fs':>4} {'hf':>4} | {'ε':>7} {'σ':>7} {'AR':>6} | {'1/K [1/m²]':>12} {'C2 [1/m]':>10} | {'ΔP [Pa]':>9} {'h [W/m²K]':>10}"
    print(header)
    print("-" * 100)

    for Fs_mm, hf_mm in cases:
        result = calculate_porous_parameters(
            Fs_mm=Fs_mm,
            hf_mm=hf_mm,
            T_celsius=T_celsius,
            v_design=v_design
        )

        g = result['geometry']
        p = result['porous']
        d = result['design_point']

        print(f"{Fs_mm:>4.0f} {hf_mm:>4.0f} | "
              f"{g['epsilon']:>7.4f} {g['sigma']:>7.4f} {g['area_ratio']:>6.2f} | "
              f"{p['inv_K']:>12.4e} {p['C2']:>10.4f} | "
              f"{d['dP_total_Pa']:>9.2f} {d['h_fs_W_m2K']:>10.2f}")

    print("=" * 100)
    print("Note: ε=공극률, σ=최소유동면적비, AR=면적비(Atot/Abare)")


# =============================================================================
# 8. 메인 실행
# =============================================================================

if __name__ == "__main__":
    # 사용자 지정 조건
    T_ambient = 14.80177    # 연평균 온도 [°C]
    v_wind = 2.019723       # 연평균 풍속 [m/s]

    print("\n" + "█" * 70)
    print("  Nir (1991) 상관식 기반 Annular Fin Porous 파라미터 계산기")
    print("█" * 70)
    print(f"\n  환경 조건: T = {T_ambient}°C, v = {v_wind} m/s")

    # 단일 케이스 상세 결과
    print("\n\n[Case 1] 상세 결과 - Fs=4mm, hf=4mm")
    result = calculate_porous_parameters(
        Fs_mm=4.0,
        hf_mm=4.0,
        T_celsius=T_ambient,
        v_design=v_wind
    )
    print_results(result)

    # 여러 케이스 비교
    test_cases = [
        (2, 4),   # 조밀한 핀
        (4, 4),   # 기준 케이스
        (6, 4),   # 넓은 간격
        (8, 4),   # 매우 넓은 간격
        (4, 5),   # 핀 높이 증가
        (4, 6),   # 핀 높이 더 증가
        (4, 8),   # 큰 핀
    ]

    generate_parameter_table(test_cases, T_ambient, v_wind)

    # Fluent/COMSOL 입력 가이드
    print("\n" + "=" * 70)
    print("  CFD 소프트웨어 입력 가이드")
    print("=" * 70)
    print("""
  [ANSYS Fluent]
  - Cell Zone Conditions → Porous Zone 체크
  - Viscous Resistance (1/m²): 위 표의 '1/K' 값 입력
  - Inertial Resistance (1/m): 위 표의 'C2' 값 입력
  - 방향: 유동 방향(Streamwise)에 위 값, 횡방향(Transverse)에 1000배 입력

  [SimScale]
  - Advanced Concepts → Porous Media
  - Permeability K (m²): 1/(1/K 값) 계산하여 입력
  - Forchheimer coefficient C_F (1/m): C2 값 직접 입력

  [COMSOL]
  - Porous Media Flow → Darcy-Forchheimer
  - Permeability (m²): K 값 입력
  - Forchheimer coefficient (1/m): C2 값 입력
    """)

    # JSON 출력 (다른 프로그램과 연동용)
    print("\n[JSON 형식 출력 (Fs=4mm, hf=4mm)]")
    import json
    json_output = {
        'T_celsius': T_ambient,
        'v_m_s': v_wind,
        'Fs_mm': 4.0,
        'hf_mm': 4.0,
        'porous_parameters': {
            'inv_K_1_m2': result['porous']['inv_K'],
            'C2_1_m': result['porous']['C2'],
            'K_m2': result['porous']['K'],
            'porosity': result['geometry']['epsilon']
        }
    }
    print(json.dumps(json_output, indent=2))
