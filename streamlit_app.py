#!/usr/bin/env python3
"""
================================================================================
Annular Fin → Porous Media Parameter Calculator - Streamlit Web App
================================================================================
Interactive web application for calculating CFD porous media parameters
from annular finned tube heat exchangers using Nir (1991) correlation.

Author: Claude (Anthropic)
Date: 2026-01-17
================================================================================
"""

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import io
import json
from cfd_porous_calc import (
    calculate_porous_parameters,
    air_properties,
    AnnularFinGeometry
)

# Page configuration
st.set_page_config(
    page_title="Porous Media Parameter Calculator",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #555;
        text-align: center;
        margin-bottom: 2rem;
    }
    .metric-card {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 5px solid #1f77b4;
    }
    .key-output {
        background-color: #fff9e6;
        padding: 1.5rem;
        border-radius: 0.5rem;
        border: 3px solid #ffa500;
        margin: 1rem 0;
    }
    .equation-box {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #28a745;
        margin: 0.5rem 0;
    }
</style>
""", unsafe_allow_html=True)

# Title and introduction
st.markdown('<div class="main-header">🔬 Porous Media Parameter Calculator</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Annular Finned Tube Heat Exchanger → CFD Porous Zone Parameters</div>', unsafe_allow_html=True)

st.markdown("""
This tool calculates CFD porous media parameters (viscous resistance **1/K** and inertial resistance **C₂**)
from annular finned tube geometry using the **Nir (1991)** correlation with **multi-point least squares fitting**.
""")

# Sidebar - Input Parameters
st.sidebar.header("📝 Input Parameters")

st.sidebar.subheader("Primary Design Parameters")
Fs_mm = st.sidebar.slider(
    "Fin Spacing (Fs) [mm]",
    min_value=0.5, max_value=30.0, value=4.0, step=0.01,
    help="Distance between adjacent fins (typical: 1-15mm)"
)

hf_mm = st.sidebar.slider(
    "Fin Height (hf) [mm]",
    min_value=0.5, max_value=50.0, value=4.0, step=0.01,
    help="Radial height of the fin (typical: 3-20mm)"
)

T_celsius = st.sidebar.slider(
    "Air Temperature [°C]",
    min_value=-50.0, max_value=100.0, value=14.8, step=0.01,
    help="Ambient air temperature (typical: -40 to 60°C)"
)

v_design = st.sidebar.slider(
    "Design Velocity [m/s]",
    min_value=0.1, max_value=20.0, value=2.0, step=0.01,
    help="Design inlet superficial velocity (typical: 0.5-10 m/s)"
)

st.sidebar.subheader("Advanced Parameters")
with st.sidebar.expander("Tube and Pitch Parameters"):
    Dc_mm = st.number_input(
        "Tube Outer Diameter (Dc) [mm]",
        min_value=5.0, max_value=100.0, value=24.0, step=0.01,
        help="Typical: 10-50mm"
    )

    s1_mm = st.number_input(
        "Transverse Pitch (s1) [mm]",
        min_value=10.0, max_value=300.0, value=55.333, step=0.01,
        help="Must be > tube outer diameter + 2×fin height"
    )

    delta_f_mm = st.number_input(
        "Fin Thickness (δf) [mm]",
        min_value=0.05, max_value=5.0, value=0.5, step=0.01,
        help="Typical: 0.1-2mm"
    )

    pitch_ratio = st.number_input(
        "Pitch Ratio (s1/s2)",
        min_value=0.1, max_value=10.0, value=1.0, step=0.01,
        help="Ratio of transverse to longitudinal pitch"
    )

    N = st.number_input(
        "Number of Tube Rows (N)",
        min_value=1, max_value=50, value=4, step=1,
        help="Typical: 1-20 rows"
    )

with st.sidebar.expander("Fitting Parameters"):
    v_min = st.number_input(
        "Minimum Velocity [m/s]",
        min_value=0.01, max_value=30.0, value=0.5, step=0.01,
        help="Lower bound for curve fitting"
    )

    v_max = st.number_input(
        "Maximum Velocity [m/s]",
        min_value=0.1, max_value=30.0, value=3.5, step=0.01,
        help="Upper bound for curve fitting"
    )

    n_points = st.number_input(
        "Number of Fitting Points",
        min_value=5, max_value=500, value=50, step=1,
        help="More points = higher accuracy but slower"
    )

# Calculate button
calculate_button = st.sidebar.button("🚀 Calculate", type="primary", use_container_width=True)

# Main area
if calculate_button or 'result' not in st.session_state:
    # Perform calculation
    with st.spinner("Calculating porous media parameters..."):
        result = calculate_porous_parameters(
            Fs_mm=Fs_mm,
            hf_mm=hf_mm,
            T_celsius=T_celsius,
            v_design=v_design,
            Dc_mm=Dc_mm,
            delta_f_mm=delta_f_mm,
            s1_mm=s1_mm,
            pitch_ratio=pitch_ratio,
            N=N,
            v_range=(v_min, v_max),
            n_points=n_points
        )
        st.session_state['result'] = result
        st.success("✅ Calculation completed successfully!")
else:
    result = st.session_state['result']

# Display results
if 'result' in st.session_state:
    result = st.session_state['result']

    # Key Equations Section
    with st.expander("📐 **Key Equations**", expanded=False):
        st.markdown("### Nir (1991) Friction Factor Correlation")
        st.latex(r"f_N = 1.1 \times Re_{Dc}^{-0.25} \times \left(\frac{S_1}{D_c}\right)^{-0.4} \times \left(\frac{A_{total}}{A_{bare}}\right)^{0.15}")

        st.markdown("### Pressure Drop Calculation")
        st.latex(r"\Delta P = N \times f_N \times \frac{\rho \times v_{max}^2}{2}")

        st.markdown("### Darcy-Forchheimer Equation")
        st.latex(r"\frac{\Delta P}{L} = \frac{\mu}{K} \cdot v + \frac{C_2 \cdot \rho}{2} \cdot v^2")
        st.latex(r"\frac{\Delta P}{L} = A \cdot v + B \cdot v^2")

        st.markdown("### Multi-Point Least Squares Fitting")
        st.latex(r"\text{Minimize: } \sum_{i=1}^{n} \left(Y_i - A \cdot v_i - B \cdot v_i^2\right)^2")
        st.latex(r"\text{Solution: } (X^T X) \begin{bmatrix} A \\ B \end{bmatrix} = X^T Y")

        st.markdown("### Porous Parameter Conversion")
        st.latex(r"\frac{1}{K} = \frac{A}{\mu} \quad \text{(Viscous Resistance [1/m²])}")
        st.latex(r"C_2 = \frac{2B}{\rho} \quad \text{(Inertial Resistance [1/m])}")

    # Results tabs
    tab1, tab2, tab3, tab4 = st.tabs(["📊 Main Results", "📈 Fitting Analysis", "📋 Detailed Data", "💾 Export"])

    with tab1:
        # CFD Porous Media Parameters (HIGHLIGHTED)
        st.markdown("## ⭐ CFD Porous Media Input Parameters")
        st.markdown('<div class="key-output">', unsafe_allow_html=True)

        col1, col2, col3 = st.columns(3)

        with col1:
            st.metric(
                label="Viscous Resistance (1/K)",
                value=f"{result['porous']['inv_K']:.4e}",
                help="Viscous resistance coefficient [1/m²]"
            )
            st.caption("**Unit:** 1/m²")

        with col2:
            st.metric(
                label="Inertial Resistance (C₂)",
                value=f"{result['porous']['C2']:.4f}",
                help="Inertial resistance coefficient [1/m]"
            )
            st.caption("**Unit:** 1/m")

        with col3:
            st.metric(
                label="Permeability (K)",
                value=f"{result['porous']['K']:.4e}",
                help="Permeability [m²]"
            )
            st.caption("**Unit:** m²")

        col4, col5 = st.columns(2)

        with col4:
            st.metric(
                label="Porosity (ε)",
                value=f"{result['geometry']['epsilon']:.4f}",
                help="Void fraction [-]"
            )
            st.caption("**Unit:** dimensionless")

        with col5:
            st.metric(
                label="Specific Surface Area (a_fs)",
                value=f"{result['geometry']['a_fs']:.2f}",
                help="Surface area per unit volume [1/m]"
            )
            st.caption("**Unit:** 1/m")

        st.markdown('</div>', unsafe_allow_html=True)

        st.markdown("---")

        # Input Conditions and Results
        col_left, col_right = st.columns(2)

        with col_left:
            st.markdown("### 📥 Input Conditions")
            inp = result['input']
            input_data = {
                "Parameter": [
                    "Fin Spacing (Fs)",
                    "Fin Height (hf)",
                    "Air Temperature",
                    "Design Velocity",
                    "Tube Outer Diameter (Dc)",
                    "Transverse Pitch (s1)",
                    "Fin Thickness (δf)",
                    "Number of Tube Rows (N)",
                    "Pitch Ratio (s1/s2)"
                ],
                "Value": [
                    f"{inp['Fs_mm']:.1f} mm",
                    f"{inp['hf_mm']:.1f} mm",
                    f"{inp['T_celsius']:.2f} °C",
                    f"{inp['v_design']:.4f} m/s",
                    f"{inp['Dc_mm']:.1f} mm",
                    f"{inp['s1_mm']:.3f} mm",
                    f"{inp['delta_f_mm']:.2f} mm",
                    f"{inp['N']}",
                    f"{inp['pitch_ratio']:.2f}"
                ]
            }
            st.table(input_data)

            st.markdown("### 🌡️ Air Properties")
            air = result['air']
            air_data = {
                "Property": [
                    "Density (ρ)",
                    "Viscosity (μ)",
                    "Thermal Conductivity (k)",
                    "Prandtl Number (Pr)"
                ],
                "Value": [
                    f"{air['rho']:.4f} kg/m³",
                    f"{air['mu']:.4e} Pa·s",
                    f"{air['k']:.4f} W/(m·K)",
                    f"{air['Pr']:.3f}"
                ]
            }
            st.table(air_data)

        with col_right:
            st.markdown("### 📐 Geometric Parameters")
            geom = result['geometry']
            geom_data = {
                "Parameter": [
                    "Fin Outer Diameter (Do)",
                    "Fin Pitch (Fp)",
                    "Longitudinal Pitch (s2)",
                    "Porosity (ε)",
                    "Min. Flow Area Ratio (σ)",
                    "Area Ratio (A_tot/A_bare)",
                    "Specific Surface Area (a_fs)"
                ],
                "Value": [
                    f"{geom['Do_mm']:.1f} mm",
                    f"{geom['Fp_mm']:.2f} mm",
                    f"{geom['s2_mm']:.3f} mm",
                    f"{geom['epsilon']:.4f}",
                    f"{geom['sigma']:.4f}",
                    f"{geom['area_ratio']:.3f}",
                    f"{geom['a_fs']:.2f} 1/m"
                ]
            }
            st.table(geom_data)

            st.markdown("### 🎯 Design Point Performance")
            design = result['design_point']
            design_data = {
                "Parameter": [
                    "Reynolds Number (Re_Dc)",
                    "Total Pressure Drop",
                    "Pressure Drop per Length",
                    "Heat Transfer Coefficient (h_fs)"
                ],
                "Value": [
                    f"{design['Re_Dc']:.1f}",
                    f"{design['dP_total_Pa']:.2f} Pa",
                    f"{design['dP_per_L_Pa_m']:.2f} Pa/m",
                    f"{design['h_fs_W_m2K']:.2f} W/(m²·K)"
                ]
            }
            st.table(design_data)

    with tab2:
        st.markdown("## 📈 Fitting Analysis")

        # Fitting information
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("Velocity Range", f"{result['input']['v_range'][0]:.1f} - {result['input']['v_range'][1]:.1f} m/s")

        with col2:
            st.metric("Fitting Points", f"{result['input']['n_points']}")

        with col3:
            st.metric("Reynolds Range", f"{result['fitting']['Re_min']:.1f} - {result['fitting']['Re_max']:.1f}")

        with col4:
            r_squared = result['porous']['R_squared']
            st.metric("R² (Fit Quality)", f"{r_squared:.6f}")

        st.markdown("---")

        # Plotting
        v_points = result['fitting']['v_points']
        Y_points = result['fitting']['Y_points']
        Y_fit = result['fitting']['Y_fit']

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Left plot: Pressure drop vs velocity
        ax1.scatter(v_points, Y_points, color='blue', s=60, alpha=0.6,
                   label='Nir Correlation (Calculated)', zorder=3)
        ax1.plot(v_points, Y_fit, 'r-', linewidth=2.5,
                label=f'Darcy-Forchheimer Fit (R²={r_squared:.6f})', zorder=2)
        ax1.set_xlabel('Velocity [m/s]', fontsize=13, fontweight='bold')
        ax1.set_ylabel('Pressure Drop per Length [Pa/m]', fontsize=13, fontweight='bold')
        ax1.set_title('Pressure Drop Fitting', fontsize=15, fontweight='bold')
        ax1.grid(True, alpha=0.3, linestyle='--')
        ax1.legend(fontsize=11, loc='upper left', framealpha=0.9)

        # Add CFD parameters textbox
        inv_K = result['porous']['inv_K']
        C2 = result['porous']['C2']
        K = result['porous']['K']

        textstr = '\n'.join((
            r'$\bf{CFD\ Porous\ Parameters:}$',
            f'1/K = {inv_K:.4e} [1/m²]',
            f'C₂ = {C2:.4f} [1/m]',
            f'K = {K:.4e} [m²]'
        ))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.85)
        ax1.text(0.98, 0.02, textstr, transform=ax1.transAxes, fontsize=10,
                verticalalignment='bottom', horizontalalignment='right', bbox=props)

        # Right plot: Residuals
        residuals = Y_points - Y_fit
        ax2.scatter(v_points, residuals, color='green', s=60, alpha=0.6)
        ax2.axhline(y=0, color='r', linestyle='--', linewidth=2)
        ax2.set_xlabel('Velocity [m/s]', fontsize=13, fontweight='bold')
        ax2.set_ylabel('Residuals [Pa/m]', fontsize=13, fontweight='bold')
        ax2.set_title('Fitting Residuals', fontsize=15, fontweight='bold')
        ax2.grid(True, alpha=0.3, linestyle='--')

        # Residual statistics
        residual_std = np.std(residuals)
        residual_mean = np.mean(residuals)
        textstr2 = '\n'.join((
            r'$\bf{Residual\ Statistics:}$',
            f'Mean = {residual_mean:.4f} Pa/m',
            f'Std Dev = {residual_std:.4f} Pa/m',
            f'Max = {np.max(np.abs(residuals)):.4f} Pa/m'
        ))
        props2 = dict(boxstyle='round', facecolor='lightgreen', alpha=0.85)
        ax2.text(0.98, 0.98, textstr2, transform=ax2.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right', bbox=props2)

        plt.tight_layout()
        st.pyplot(fig)

        # Fitting coefficients
        st.markdown("### 📊 Darcy-Forchheimer Coefficients")
        col1, col2 = st.columns(2)

        with col1:
            st.metric("Coefficient A (Linear)", f"{result['porous']['A']:.4e} Pa·s/m²")
            st.latex(r"A = \frac{\mu}{K}")

        with col2:
            st.metric("Coefficient B (Quadratic)", f"{result['porous']['B']:.4e} Pa·s²/m³")
            st.latex(r"B = \frac{C_2 \cdot \rho}{2}")

    with tab3:
        st.markdown("## 📋 Detailed Calculation Data")

        # All fitting points in a table
        st.markdown("### Fitting Data Points")

        fitting_table_data = {
            "Point": range(1, len(v_points) + 1),
            "Velocity [m/s]": [f"{v:.3f}" for v in v_points],
            "Re_Dc": [f"{re:.1f}" for re in result['fitting']['Re_points']],
            "ΔP/L (Nir) [Pa/m]": [f"{y:.4f}" for y in Y_points],
            "ΔP/L (Fit) [Pa/m]": [f"{y:.4f}" for y in Y_fit],
            "Residual [Pa/m]": [f"{r:.4f}" for r in (Y_points - Y_fit)]
        }

        st.dataframe(fitting_table_data, height=400)

        # Summary statistics
        st.markdown("### 📊 Statistical Summary")
        col1, col2, col3 = st.columns(3)

        with col1:
            st.markdown("**Nir Correlation Results**")
            st.write(f"Mean: {np.mean(Y_points):.4f} Pa/m")
            st.write(f"Std Dev: {np.std(Y_points):.4f} Pa/m")
            st.write(f"Min: {np.min(Y_points):.4f} Pa/m")
            st.write(f"Max: {np.max(Y_points):.4f} Pa/m")

        with col2:
            st.markdown("**Fitted Results**")
            st.write(f"Mean: {np.mean(Y_fit):.4f} Pa/m")
            st.write(f"Std Dev: {np.std(Y_fit):.4f} Pa/m")
            st.write(f"Min: {np.min(Y_fit):.4f} Pa/m")
            st.write(f"Max: {np.max(Y_fit):.4f} Pa/m")

        with col3:
            residuals = Y_points - Y_fit
            st.markdown("**Residuals**")
            st.write(f"Mean: {np.mean(residuals):.4f} Pa/m")
            st.write(f"Std Dev: {np.std(residuals):.4f} Pa/m")
            st.write(f"Min: {np.min(residuals):.4f} Pa/m")
            st.write(f"Max: {np.max(residuals):.4f} Pa/m")

    with tab4:
        st.markdown("## 💾 Export Results")

        st.markdown("Download your calculation results in various formats:")

        col1, col2, col3 = st.columns(3)

        # Generate filename
        inp = result['input']
        filename_base = f"Fs{inp['Fs_mm']:.1f}_hf{inp['hf_mm']:.1f}_T{inp['T_celsius']:.1f}_v{inp['v_design']:.2f}"

        with col1:
            # PNG download
            st.markdown("### 📊 PNG Image")

            # Create plot for download
            fig_download, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

            # Left plot
            ax1.scatter(v_points, Y_points, color='blue', s=50, alpha=0.6,
                       label='Nir Correlation (Calculated)', zorder=3)
            ax1.plot(v_points, Y_fit, 'r-', linewidth=2,
                    label=f'Darcy-Forchheimer Fit (R²={r_squared:.6f})', zorder=2)
            ax1.set_xlabel('Velocity [m/s]', fontsize=13)
            ax1.set_ylabel('Pressure Drop per Length [Pa/m]', fontsize=13)
            ax1.set_title('Pressure Drop Fitting', fontsize=15, fontweight='bold')
            ax1.grid(True, alpha=0.3)
            ax1.legend(fontsize=11, loc='upper left')

            textstr = '\n'.join((
                r'$\bf{CFD\ Porous\ Parameters:}$',
                f'1/K = {inv_K:.4e} [1/m²]',
                f'C₂ = {C2:.4f} [1/m]',
                f'K = {K:.4e} [m²]'
            ))
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
            ax1.text(0.98, 0.02, textstr, transform=ax1.transAxes, fontsize=11,
                    verticalalignment='bottom', horizontalalignment='right', bbox=props)

            # Right plot
            residuals = Y_points - Y_fit
            ax2.scatter(v_points, residuals, color='green', s=50, alpha=0.6)
            ax2.axhline(y=0, color='r', linestyle='--', linewidth=2)
            ax2.set_xlabel('Velocity [m/s]', fontsize=13)
            ax2.set_ylabel('Residuals [Pa/m]', fontsize=13)
            ax2.set_title('Fitting Residuals', fontsize=15, fontweight='bold')
            ax2.grid(True, alpha=0.3)

            residual_std = np.std(residuals)
            residual_mean = np.mean(residuals)
            textstr2 = '\n'.join((
                r'$\bf{Residual\ Statistics:}$',
                f'Mean = {residual_mean:.4f} Pa/m',
                f'Std Dev = {residual_std:.4f} Pa/m'
            ))
            props2 = dict(boxstyle='round', facecolor='lightgreen', alpha=0.8)
            ax2.text(0.98, 0.98, textstr2, transform=ax2.transAxes, fontsize=11,
                    verticalalignment='top', horizontalalignment='right', bbox=props2)

            plt.tight_layout()

            # Save to buffer
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=150, bbox_inches='tight')
            buf.seek(0)

            st.download_button(
                label="📥 Download PNG",
                data=buf,
                file_name=f"{filename_base}.png",
                mime="image/png",
                use_container_width=True
            )

            plt.close(fig_download)

        with col2:
            # TXT download
            st.markdown("### 📄 Text File")

            # Generate text content
            txt_content = f"""================================================================================
       Annular Fin → Porous Media Parameter Conversion Results
       Nir (1991) Correlation with Multi-Point Least Squares Fitting
================================================================================

[INPUT CONDITIONS]
  Fin Spacing (Fs)        : {inp['Fs_mm']:.1f} mm
  Fin Height (hf)         : {inp['hf_mm']:.1f} mm
  Air Temperature         : {inp['T_celsius']:.2f} °C
  Design Velocity         : {inp['v_design']:.4f} m/s
  Tube Outer Diameter (Dc): {inp['Dc_mm']:.1f} mm
  Transverse Pitch (s1)   : {inp['s1_mm']:.3f} mm
  Number of Tube Rows (N) : {inp['N']}

[AIR PROPERTIES @ {inp['T_celsius']:.2f}°C]
  Density (ρ)             : {result['air']['rho']:.4f} kg/m³
  Viscosity (μ)           : {result['air']['mu']:.4e} Pa·s
  Thermal Conductivity (k): {result['air']['k']:.4f} W/(m·K)
  Prandtl Number (Pr)     : {result['air']['Pr']:.3f}

[GEOMETRIC PARAMETERS]
  Fin Outer Diameter (Do) : {result['geometry']['Do_mm']:.1f} mm
  Fin Pitch (Fp)          : {result['geometry']['Fp_mm']:.2f} mm
  Porosity (ε)            : {result['geometry']['epsilon']:.4f}
  Area Ratio (Atot/Abare) : {result['geometry']['area_ratio']:.3f}

================================================================================
  ★★★ CFD POROUS MEDIA INPUT PARAMETERS (KEY OUTPUT) ★★★
================================================================================
  Viscous Resistance (1/K): {result['porous']['inv_K']:.4e}  [1/m²]
  Inertial Resistance (C2): {result['porous']['C2']:.4f}  [1/m]
  Permeability (K)        : {result['porous']['K']:.4e}  [m²]
  Porosity (ε)            : {result['geometry']['epsilon']:.4f}  [-]
  Specific Surface (a_fs) : {result['geometry']['a_fs']:.2f}  [1/m]
================================================================================

[DESIGN POINT PERFORMANCE]
  Reynolds Number (Re_Dc) : {result['design_point']['Re_Dc']:.1f}
  Total Pressure Drop     : {result['design_point']['dP_total_Pa']:.2f} Pa
  Heat Transfer Coef (hfs): {result['design_point']['h_fs_W_m2K']:.2f} W/(m²·K)

[FITTING INFORMATION]
  Velocity Range          : {inp['v_range'][0]:.1f} ~ {inp['v_range'][1]:.1f} m/s
  Number of Fitting Points: {inp['n_points']}
  Coefficient of Determination (R²): {result['porous']['R_squared']:.8f}

Generated by: Nir (1991) Porous Parameter Calculator (Streamlit App)
Date: 2026-01-17
================================================================================
"""

            st.download_button(
                label="📥 Download TXT",
                data=txt_content,
                file_name=f"{filename_base}.txt",
                mime="text/plain",
                use_container_width=True
            )

        with col3:
            # JSON download
            st.markdown("### 🔧 JSON File")

            json_output = {
                'input': result['input'],
                'porous_parameters': {
                    'inv_K_1_m2': result['porous']['inv_K'],
                    'C2_1_m': result['porous']['C2'],
                    'K_m2': result['porous']['K'],
                    'porosity': result['geometry']['epsilon'],
                    'a_fs_1_m': result['geometry']['a_fs'],
                    'R_squared': result['porous']['R_squared']
                },
                'design_point': {
                    'Re_Dc': result['design_point']['Re_Dc'],
                    'dP_total_Pa': result['design_point']['dP_total_Pa'],
                    'h_fs_W_m2K': result['design_point']['h_fs_W_m2K']
                }
            }

            json_str = json.dumps(json_output, indent=2)

            st.download_button(
                label="📥 Download JSON",
                data=json_str,
                file_name=f"{filename_base}.json",
                mime="application/json",
                use_container_width=True
            )

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #888; font-size: 0.9rem;">
    <p><strong>Porous Media Parameter Calculator</strong> | Based on Nir (1991) Correlation</p>
    <p>Developed with ❤️ using Streamlit | © 2026</p>
</div>
""", unsafe_allow_html=True)
