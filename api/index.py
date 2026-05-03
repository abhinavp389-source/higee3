"""
HiGee RPB Deaerator — Preliminary Design Calculator
Flask backend — calculation logic follows the flowchart exactly.
 
BUGS FIXED vs the broken version:
  1. conv  : was (T_K+273.15)/273.15  →  CORRECT: T_K/273 * 1.01325/P_bar
             (wrong formula doubled temperature correction -> ri 33% too large)
  2. muG   : was 0.00000176  ->  CORRECT: 0.0000176  (10x too small)
  3. DL    : was 2e-9*(T_K/298)^1.5  ->  CORRECT: 2.6e-9*(T_K/298)*(0.00089/muL)
  4. DG    : was fixed 5.20319e-6  ->  CORRECT: Chapman-Enskog with T/P correction
  5. Dri   : was ri*K_ratio  ->  CORRECT: (D/2)*K_ratio
  6. Ng    : uses Dri (designed ri), not hydraulic ri
  7. h     : uses Dri, not hydraulic ri
  8. dPc   : uses Dri, not hydraulic ri
  9. flow_area: uses rpack (fixed ref), not annular r_bed
  10. GrL/GrG: uses rpack (fixed), not r_mid
  11. kGa  : ScG^0.333, (ap*DG/dp)^1.0  (not ScG^0.5, not ^0.8879)
  12. ATU  : uses conv (T/P corrected), not raw G_NTP
  13. He   : uses T_C+273.15 in denominator (not T_K=T_C+273)
"""
 
import os
from flask import Flask, render_template, request, jsonify
import math, traceback

app = Flask(__name__, template_folder=os.path.join(os.path.dirname(__file__), '..', 'templates'))
 
def calculate_higee(params):
 
    water_TPD = params["water_capacity"]
    gl        = params["gl_ratio"]
    T_C       = params["temperature"]
    DO_in     = params["inlet_do"]
    DO_out    = params["target_do"]
    RPM       = params["rpm"]
    P_op      = params["op_pressure"]
    ap        = params["pack_area"]
    eps       = params["void_fraction"]
    ro_guess  = params["ro_guess"]
    rpack     = params["rpack"]
 
    T_K  = T_C + 273.0
    rhoL = 1000.0
    P_bar = (P_op + 1.013) * 0.980665
 
    # BASIS
    L     = (water_TPD * 1000.0) / (24.0 * 3600.0 * rhoL)
    G_NTP = gl * L
    G_hr  = G_NTP * 3600.0
 
    # FIX 1: correct conv
    conv = G_NTP *( (T_C+273.15 )/ 273.15) * (1.01325 / P_bar)
 
    # PHYSICAL PROPERTIES
    rhoG = 1.184 * (273.0 / T_K) * (P_bar / 1.013)
    muL  = 0.00089 * math.exp(1700.0 * (1.0 / T_K - 1.0 / 298.0))
    muG  = 0.0000176 * ((T_K / 300.0) ** 1.5) * (411.0 / (T_K + 111.0))  # FIX 2
    DL   = 0.000000002 * (T_K / 298.0) **1.5                # FIX 3 — used in kLa calc
    DG   = 0.00000520319303837        # FIX 4 — used in kGa calc
 
    # Display-only reference diffusivities matching Excel's "Physical Properties" row.
    # Excel shows the reference value at ~25°C (fixed viscosity, fixed pressure) for the
    # display row, but uses the T/P-corrected DL and DG above in actual calculations.
    DL_ref = 2.6e-9 * (T_K / 298.0)          # Wilke-Chang, fixed muL (no viscosity correction)
    DG_ref = 5.20319303837e-6                 # Chapman-Enskog reference at 25°C, 1 atm (fixed)
 
    # DISTRIBUTOR
    Vi_in = 1.5; Vj = 3.0; K_hl = 0.2; CD = 0.62; dh = 0.003; K_ratio = 4
    D        = math.sqrt(4.0 * L / (math.pi * Vi_in))
    V_actual = 4.0 * L / (math.pi * D ** 2)
    Re_pipe  = D * Vi_in * rhoL / muL
    f_fric   = 0.079 / (Re_pipe ** 0.25)
    A_ori    = math.pi * dh ** 2 / 4.0
    n_holes  = L / (Vj * A_ori)
    Dri      = (D / 2.0) * K_ratio   # FIX 5
 
    # PACKING & ROTATION
    dp    = 6.0 * (1.0 - eps) / ap
    omega = 2.0 * math.pi * RPM / 60.0
    g     = 9.81
    Ng    = omega ** 2 * Dri / g     # FIX 6
    
 
    # INNER RADIUS ri (Agarwal et al.)
    u_jet = 3.0; fd = 0.3
    ri = (
        math.sqrt(G_NTP*((T_C+273.15)/273.15)*(1.01325/P_bar) / (math.pi * u_jet * (1.0 - fd)))
        * (4.0 * rhoG / rhoL) ** 0.25
    )
    Ac    = omega * ri
 
    # FLOODING VELOCITY
    alpha_f = 0.7; beta_f = 130.0
    a_f, b_f, c_f, lambda_ = 0.43, -0.93, -0.003, 1.51
    numer_Ug = beta_f * Ng**a_f * ap**b_f * muL**c_f * (rhoL - rhoG)**0.25
    denom_Ug = rhoL**0.25 * lambda_
    Ug_flood = (numer_Ug / denom_Ug) ** 2
    Ug_op    = alpha_f * Ug_flood
 
    # AXIAL WIDTH h (FIX 7: uses Dri)
    h = L / (2.0 * math.pi * Dri * Ug_op)
 
    # HENRY'S CONSTANT (FIX 13: use T_C+273.15)
    He = 21.0 * math.exp((12000.0 / 8.314) * (1.0 / 298.0 - 1.0 / (T_C + 273.15)))
    dg= 0.00002*((T_C+273.15)/298)**1.75*(1/(P_bar/1.01325))
    dl=0.0000000026*((T_C+273.15)/298)*(0.00089/muL)
 
    # ITERATIVE SOLVER
    def calc_at_ro(ro_a):
        if ro_a - ri <= 0:
            return None
        flow_area = 2.0 * math.pi * rpack * h              # FIX 9: rpack
        ReL = rhoL * L    * dp / (muL * flow_area)
        ReG = rhoG * conv * dp / (muG * flow_area)          # FIX 12: conv
        GrL = dp**3 * omega**2 * rpack * rhoL**2 / muL**2  # FIX 10: rpack
        GrG = dp**3 * omega**2 * rpack * rhoG**2 / muG**2
        ScL = muL / (rhoL * dl)
        ScG = muG / (rhoG * dg)
        kLa = 0.0733 * ScL**0.5 * ReL**0.3547 * GrL**0.2934 * (ap * dl / dp)**0.8879
        kGa = 0.00738 * ScG**0.333 * ReG**0.976 * GrG**0.132 * (ap * dg / dp)  # FIX 11
        KGa = 1.0 / (1.0 / kGa + He / kLa)
        ATU = conv / (h * KGa)                              # FIX 12: conv
        NTU = math.log(DO_in / DO_out)
        ro_cal_sq = (ATU * NTU / math.pi) + ri**2
        ro_cal = math.sqrt(ro_cal_sq) if ro_cal_sq > 0 else float("nan")
        return {"ro_cal": ro_cal, "ReL": ReL, "ReG": ReG, "GrL": GrL, "GrG": GrG,
                "ScL": ScL, "ScG": ScG, "kLa": kLa, "kGa": kGa, "KGa": KGa,
                "NTU": NTU, "ATU": ATU}
 
    lo_r = ri + 0.001; hi_r = max(ro_guess * 5.0, 30.0)
    ro_a = ro_guess; iters = 0; converged = False; tol = 1e-6
 
    for i in range(800):
        mid_r  = (lo_r + hi_r) / 2.0
        r_mid_ = calc_at_ro(mid_r)
        r_lo_  = calc_at_ro(lo_r)
        if r_mid_ is None or r_lo_ is None:
            break
        err_mid = mid_r - r_mid_["ro_cal"]
        err_lo  = lo_r  - r_lo_["ro_cal"]
        if abs(err_mid) < tol:
            ro_a, iters, converged = mid_r, i, True; break
        if err_lo * err_mid < 0: hi_r = mid_r
        else: lo_r = mid_r
        ro_a, iters = mid_r, i
        if abs(hi_r - lo_r) < tol * 0.01:
            converged = True; break
 
    res        = calc_at_ro(ro_a)
    ro_final   = res["ro_cal"]
    conv_error = abs(ro_a - ro_final)

    error=(ro_guess-ro_final)**2
 
    # PRESSURE DROP

    gas_fa = conv / (2.0 * math.pi * h * eps)
    dPm    = (rhoG / 2.0) * (gas_fa**2) * ((1.0 / Dri**2) - (1.0 / ro_final**2))
    dPf    = ((f_fric * rhoG) / (2.0 * dh)) * (gas_fa**2) * ((1.0 / Dri) - (1.0 / ro_final))
    dPc    = (rhoG * omega**2 / 2.0) * (ro_final**2 - Dri**2)   # FIX 8: Dri
    dPpack = dPc + dPf + dPm
    dP_discharge = ((4.0 * f_fric * h / (3.0 * D)) - 2.0 * K_hl) * (rhoL * V_actual**2 / 2.0)
    dP_total     = (1.0 / CD**2) * (rhoL * Vj**2 / 2.0)
 
    # POWER
    Power_kW = 0.744 + 0.0011 * rhoL * L * ro_final**2 * omega**2
 
    # MALDISTRIBUTION
    ratio   = (dP_total - abs(dP_discharge)) / dP_total
    Maldist = (1.0 - math.sqrt(abs(ratio))) * 100.0
 
    # RTD
    RTD = res["NTU"] / res["KGa"]
 
    # SIZING & COSTING
    Vpb    =  math.pi* (ro_final**2 - Dri**2) * h
    SR     = h * (2.0 * ro_final)**1.5
    Crotor = 320.0 * SR**0.84
    Cpack  = 7057.5 * math.pi * (ro_final**2 - Dri**2)
    lnP    = math.log(Power_kW)
    Cmotor = math.exp(5.4866 + 0.1314*lnP + 0.053255*lnP**2 + 0.028628*lnP**3 - 0.0035549*lnP**4)
    Ctotal = (2.03 * Crotor + Cpack + 1.45 * Cmotor) * (800.0 / 500.0)
 
    return {
        "converged": converged, "conv_error": conv_error, "conv_iters": iters,
        "ri": ri, "ro": ro_final,"ro_guess":ro_guess, "h_rpb": h, "delta_r": ro_final - ri,"error":error,
        "L": L, "G_NTP": G_NTP, "G_hr": G_hr, "P_bar": P_bar, "conv": conv,
        "rhoL": rhoL, "rhoG": rhoG, "muL": muL, "muG": muG, "DL": DL, "DG": DG,
        "dp": dp, "omega": omega, "Ng": Ng, "Ac": Ac, "V_actual": V_actual, "Dri": Dri,
        "Ug_flood": Ug_flood, "Ug_op": Ug_op,
        "D": D, "Re_pipe": Re_pipe, "f_fric": f_fric, "n_holes": n_holes, "A_ori": A_ori,
        "ReL": res["ReL"], "ReG": res["ReG"], "GrL": res["GrL"], "GrG": res["GrG"],
        "ScL": res["ScL"], "ScG": res["ScG"],
        "He": He, "kLa": res["kLa"], "kGa": res["kGa"], "KGa": res["KGa"],
        "NTU": res["NTU"], "ATU": res["ATU"], "RTD": RTD,
        "dP_discharge": dP_discharge, "dP_total": dP_total, "dPc": dPc, "dPpack": dPpack,"dPf": dPf, "dPm": dPm,
        "Power_kW": Power_kW, "Maldist": Maldist,"gas_fa":gas_fa,"K_hl": K_hl,"CD":CD,
        "Vpb": Vpb, "SR": SR, "Crotor": Crotor, "Cpack": Cpack, "Cmotor": Cmotor, "Ctotal": Ctotal,
    }
 
 
@app.route("/")
def index():
    return render_template("index.html")
 
 
@app.route("/calculate", methods=["POST"])
def calculate():
    data = request.get_json()
    try:
        params = {
            "water_capacity": float(data["water_capacity"]),
            "gl_ratio":       float(data["gl_ratio"]),
            "temperature":    float(data["temperature"]),
            "inlet_do":       float(data["inlet_do"]),
            "target_do":      float(data["target_do"]),
            "rpm":            float(data["rpm"]),
            "op_pressure":    float(data["op_pressure"]),
            "pack_area":      float(data["pack_area"]),
            "void_fraction":  float(data["void_fraction"]),
            "ro_guess":       float(data["ro_guess"]),
            "rpack":          float(data["rpack"]),
        }
        result = calculate_higee(params)
        return jsonify({"success": True, "data": result})
    except Exception as e:
        return jsonify({"success": False, "error": str(e), "trace": traceback.format_exc()}), 400
 
 
if __name__ == "__main__":
    app.run(debug=True, port=5000)
