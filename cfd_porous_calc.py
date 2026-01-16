"""
CFD Porous Zone 파라미터 계산기
15°C 공기, 1.5m/s 유속, LTNE(비평형) 조건을 반영한 핀-튜브 열교환기 설계

입력 변수: Fs_mm (핀 간격), hf_mm (핀 높이), pitch_ratio (s1/s2 비율)
출력: Porosity, Surface Area Density, Viscous Resistance (1/K), 
      Inertial Resistance (C2), Interfacial Heat Transfer Coefficient (h_fs)
"""

import math


def calculate_advanced_cfd_params(Fs_mm, hf_mm, pitch_ratio):
    """
    15도씨 공기, 1.5m/s 조건, Annular Zone 설계를 위한 CFD 파라미터 계산기
    
    Parameters:
    -----------
    Fs_mm : float
        핀 간격 (Fin Spacing) [mm]
    hf_mm : float
        핀 높이 (Fin Height) [mm]
    pitch_ratio : float
        종방향/횡방향 피치 비율 (Ratio = s1/s2)
        조건: 0.5 < Ratio < 2.0
        s2 = s1 / Ratio로 계산됨
    
    Returns:
    --------
    dict
        계산 결과 딕셔너리:
        - Error: bool - 설계 오류 여부
        - Message: str - 오류 메시지 (오류 시만)
        - Fs_mm, hf_mm, Ratio: 입력값 반환
        - Porosity: 공극률 [-]
        - Area_Density: 비표면적 [1/m]
        - Viscous_Res: 점성 저항 (1/K) [1/m²]
        - Inertial_Res: 관성 저항 (C2) [1/m]
        - h_fs: 계면 열전달 계수 [W/m²·K]
        - Re_Design: 설계점 레이놀즈 수 [-]
    """
    
    # --- 1. 고정 파라미터 및 15도씨 물성치 ---
    # Air properties at 15 degC, 1 atm (명료화된 값)
    rho = 1.225             # 공기 밀도 [kg/m³]
    mu = 1.789e-5           # 공기 점성 [Pa·s]
    k_air = 0.0253          # 공기 열전도율 [W/m·K]
    Pr = 0.71               # Prandtl 수 [-]

    # 기하학적 고정값
    Dc = 0.024              # 튜브 외경 [m]
    delta_f = 0.0005        # 핀 두께 [m]
    s1 = 0.05533            # 횡방향 피치 [m] (Spanwise Pitch)
    N = 4                   # 튜브 열 수 (단위 크기 0.2m / 피치 0.055m ≈ 4열)

    # 입력값 SI 단위 변환
    Fs = Fs_mm / 1000.0     # 핀 간격 [m]
    hf = hf_mm / 1000.0     # 핀 높이 [m]
    
    # 종방향 피치 계산 (Ratio = s1/s2)
    # 0.5 < Ratio < 2.0 조건 확인
    if pitch_ratio <= 0.5 or pitch_ratio >= 2.0:
        return {
            "Error": True,
            "Message": f"Invalid pitch_ratio ({pitch_ratio}). Must be 0.5 < Ratio < 2.0"
        }
    
    s2 = s1 / pitch_ratio   # 종방향 피치 [m] (Streamwise Pitch)

    # 종속 변수 계산
    Fp = Fs + delta_f       # 핀 피치 (Center-to-Center) [m]
    Do = Dc + 2 * hf        # 핀 외경 [m]

    # --- [Safety Check] 물리적 간섭 확인 (Overlap Check) ---
    # 핀끼리 닿거나 튜브가 겹치는지 확인
    transverse_gap = s1 - Do    # 횡방향 간격
    longitudinal_gap = s2 - Do  # 종방향 간격 (In-line 기준)
    
    # Staggered 배열의 경우 대각선 거리도 체크할 수 있으나, 여기서는 단순화
    # 최소 간격이 음수이면 충돌 발생
    if transverse_gap < 0:
        return {
            "Error": True,
            "Message": f"COLLISION! Fin OD ({Do*1000:.1f} mm) > Transverse Pitch ({s1*1000:.2f} mm)"
        }
    
    if longitudinal_gap < 0:
        return {
            "Error": True,
            "Message": f"COLLISION! Fin OD ({Do*1000:.1f} mm) > Longitudinal Pitch ({s2*1000:.2f} mm)"
        }

    # --- 2. 기하학적 파라미터 (Porous Zone Settings) ---
    
    # [결과 1] 공극률 (Porosity, Annular Zone 기준)
    epsilon = 1 - (delta_f / Fp)

    # [결과 2] 비표면적 (Surface Area Density, a_fs)
    # 분자: 핀 양면 + 핀 끝단 + 튜브 노출부
    area_numerator = (0.5 * math.pi * (Do**2 - Dc**2)) + \
                     (math.pi * Do * delta_f) + \
                     (math.pi * Dc * Fs)
    # 분모: 도넛 1피치 부피 (Annular Zone)
    vol_denominator = (0.25 * math.pi * (Do**2 - Dc**2)) * Fp
    
    a_fs = area_numerator / vol_denominator  # [1/m]

    # 수력 직경 (Hydraulic Diameter)
    Dh = 4 * epsilon / a_fs  # [m]
    
    # 최소 유동 단면적 비율 (Sigma)
    # Annular finned tube bank에서 최소 면적은 횡방향 중심선 기준
    sigma = (s1 - Dc - 2*hf*(delta_f/Fp)) / s1

    # --- 3. 유동 저항 (Pressure Drop - Wang Correlation) ---
    
    def get_wang_params(v_inlet_frontal):
        """
        Wang 상관식을 사용한 압력강하 계산
        
        Parameters:
        -----------
        v_inlet_frontal : float
            입구면 유속 (Frontal Velocity) [m/s]
        
        Returns:
        --------
        tuple: (dP_L, Re_Dc, v_max)
            dP_L: 단위 길이당 압력강하 [Pa/m]
            Re_Dc: 튜브 외경 기준 레이놀즈 수 [-]
            v_max: 최대 유속 (Physical Velocity) [m/s]
        """
        # Physical Velocity (Maximum Velocity inside fins)
        v_max = v_inlet_frontal / sigma
        
        # Reynolds Number (튜브 외경 기준)
        Re_Dc = (rho * v_max * Dc) / mu
        
        # Wang Correlation Coefficients
        # F1, F2, F3 계산 (Log base e)
        F1 = -0.764 + 0.739*(s1/s2) + 0.177*(Fp/Dc) - (0.00758/N)
        F2 = -15.689 + (64.021 / math.log(Re_Dc))
        F3 = 1.696 - (15.695 / math.log(Re_Dc))
        
        # 마찰 계수 (Friction Factor)
        # Fs가 반영된 Wang 상관식: f = 0.0267 * Re^F1 * (s1/s2)^F2 * (Fs/Dc)^F3
        f = 0.0267 * (Re_Dc**F1) * ((s1/s2)**F2) * ((Fs/Dc)**F3)
        
        # 압력강하 (단위 길이당) [Pa/m]
        dP_L = (f / Dh) * 0.5 * rho * (v_max**2)
        
        return dP_L, Re_Dc, v_max

    # 2-Point Fitting Method (1.5m/s & 0.5m/s)
    Y1, Re_design, v_max_design = get_wang_params(1.5)  # Design point (1.5 m/s)
    Y2, _, _ = get_wang_params(0.5)                      # Low point (0.5 m/s)
    
    # Physical Velocity 계산
    v1_phy = 1.5 / sigma  # Design point의 Physical Velocity
    v2_phy = 0.5 / sigma  # Low point의 Physical Velocity

    # Darcy-Forchheimer 계수 도출: Y = A*v + B*v²
    # Z = Y/v = A + B*v 형태로 선형 회귀
    Z1 = Y1 / v1_phy
    Z2 = Y2 / v2_phy
    
    B_coeff = (Z1 - Z2) / (v1_phy - v2_phy)  # B 계수
    A_coeff = Z1 - B_coeff * v1_phy          # A 계수

    # [결과 3 & 4] CFD Resistance Coefficients
    inv_K = A_coeff / mu      # Viscous Resistance (1/K) [1/m²]
    C2 = 2 * B_coeff / rho    # Inertial Resistance (C2) [1/m]

    # --- 4. 열전달 계수 (Heat Transfer - Briggs & Young Correlation) ---
    # LTNE 모델을 위해 h_fs가 필요함
    # Briggs and Young Correlation (Annular Finned Tubes 전용)
    # Nu = 0.134 * Re_D^0.681 * Pr^0.33 * (s/l)^0.2 * (s/t)^0.113
    # 여기서:
    #   l = fin height (hf)
    #   t = fin thickness (delta_f)
    #   s = fin spacing (Fs)
    
    Nu = 0.134 * (Re_design**0.681) * (Pr**0.33) * \
         ((Fs/hf)**0.2) * ((Fs/delta_f)**0.113)
         
    # 계면 열전달 계수 (Interfacial Heat Transfer Coefficient) [W/m²·K]
    # h = Nu * k / Dc (Briggs & Young은 Dc를 특성 길이로 사용)
    h_fs = Nu * k_air / Dc

    return {
        "Error": False,
        "Fs_mm": Fs_mm,
        "hf_mm": hf_mm,
        "Ratio": pitch_ratio,
        "s2_mm": s2 * 1000.0,  # 참고용
        "Porosity": epsilon,
        "Area_Density": a_fs,
        "Viscous_Res": inv_K,
        "Inertial_Res": C2,
        "h_fs": h_fs,
        "Re_Design": Re_design,
        "Nu": Nu,  # 참고용
        "Dh_mm": Dh * 1000.0,  # 참고용
        "sigma": sigma  # 참고용
    }


# =================================================================
# 실행 영역 (CLI 입력)
# =================================================================

if __name__ == "__main__":
    print("=" * 120)
    print("CFD Porous Zone 파라미터 계산기")
    print("15°C 공기, 1.5m/s 유속, LTNE 조건 기준")
    print("=" * 120)
    print()
    
    # 사용자 입력 받기
    try:
        Fs_input = float(input("핀 간격 (Fs) [mm]을 입력하세요 (예: 2.0): "))
        hf_input = float(input("핀 높이 (hf) [mm]을 입력하세요 (예: 4.0): "))
        ratio_input = float(input("피치 비율 (Ratio = s1/s2)을 입력하세요 (0.5 < Ratio < 2.0, 예: 1.0): "))
        
        print()
        print("-" * 120)
        print(f"입력값:")
        print(f"  - 핀 간격 (Fs): {Fs_input:.2f} mm")
        print(f"  - 핀 높이 (hf): {hf_input:.2f} mm")
        print(f"  - 피치 비율 (Ratio = s1/s2): {ratio_input:.2f}")
        print("-" * 120)
    except ValueError:
        print("\n[오류] 숫자를 입력해주세요.")
        exit(1)
    except KeyboardInterrupt:
        print("\n\n프로그램이 중단되었습니다.")
        exit(0)
    
    # 계산 실행
    result = calculate_advanced_cfd_params(Fs_input, hf_input, ratio_input)
    
    # 결과 출력
    if result["Error"]:
        print(f"[오류] {result['Message']}")
        print("=" * 120)
    else:
        print("계산 결과 (Fluent 입력값):")
        print(f"  1. Porosity (ε): {result['Porosity']:.6f} [-]")
        print(f"  2. Surface Area Density (a_fs): {result['Area_Density']:.2f} [1/m]")
        print(f"  3. Viscous Resistance (1/K): {result['Viscous_Res']:.6e} [1/m²]")
        print(f"  4. Inertial Resistance (C2): {result['Inertial_Res']:.6f} [1/m]")
        print(f"  5. Interfacial Heat Transfer Coeff (h_fs): {result['h_fs']:.2f} [W/m²·K]")
        print("-" * 120)
        print("참고값:")
        print(f"  - 종방향 피치 (s2): {result['s2_mm']:.2f} mm")
        print(f"  - 설계점 레이놀즈 수 (Re_Dc): {result['Re_Design']:.0f} [-]")
        print(f"  - Nusselt 수 (Nu): {result['Nu']:.2f} [-]")
        print(f"  - 수력 직경 (Dh): {result['Dh_mm']:.3f} mm")
        print(f"  - 최소 면적 비율 (σ): {result['sigma']:.4f} [-]")
        print("=" * 120)
        print("\n[Fluent 입력 가이드는 README_cfd.md를 참조하세요]")
