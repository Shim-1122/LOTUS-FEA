"""LaTeX report generation from a run directory (compiled when pdflatex exists)."""
from __future__ import annotations

import json
import os
import shutil
import subprocess

import numpy as np


def esc(s) -> str:
    s = str(s)
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("$", r"\$"), ("#", r"\#"),
                 ("_", r"\_"), ("{", r"\{"), ("}", r"\}"), ("~", r"\textasciitilde{}"), ("^", r"\^{}")):
        s = s.replace(a, b)
    return s


def _tab(rows, header, spec=None):
    spec = spec or "l" * len(header)
    out = [r"\begin{longtable}{" + spec + r"}\toprule", " & ".join(r"\textbf{" + esc(h) + "}" for h in header) + r"\\\midrule\endhead"]
    for r in rows:
        out.append(" & ".join(esc(x) for x in r) + r"\\")
    out.append(r"\bottomrule\end{longtable}")
    return "\n".join(out)


def _fmt(x, n=4):
    if isinstance(x, (float, np.floating)):
        return f"{x:.{n}g}"
    return str(x)


def _fig(path, caption, width="\\textwidth"):
    return (r"\begin{figure}[H]\centering\includegraphics[width=" + width + "]{" + path + r"}"
            r"\caption{" + caption + r"}\end{figure}")


EQUATIONS = r"""
\paragraph{Strain-energy functions as implemented by FEBio 4.13 (checked against source and against an independent
implementation for homogeneous deformations, Section~\ref{sec:verif}).}
Isochoric invariants use $\tilde{\mathbf C}=J^{-2/3}\mathbf C$.
\begin{align*}
\text{HGO (skin):}\quad & W=\tfrac{c}{2}(\tilde I_1-3)+\frac{k_1}{2k_2}\sum_{a=1,2}\left[e^{k_2\langle E_a\rangle^2}-1\right]
 +\tfrac{k}{2}\left[\tfrac{J^2-1}{2}-\ln J\right],\quad E_a=\kappa(\tilde I_1-3)+(1-3\kappa)(\tilde I_{4a}-1)\\
\text{Ogden (fat, muscle):}\quad & W=\sum_i \frac{c_i}{m_i^2}\left(\tilde\lambda_1^{m_i}+\tilde\lambda_2^{m_i}+\tilde\lambda_3^{m_i}-3\right)+\tfrac{k}{2}(\ln J)^2\\
\text{Mooney--Rivlin (silicone, cushion):}\quad & W=c_1(\tilde I_1-3)+c_2(\tilde I_2-3)+\tfrac{k}{2}(\ln J)^2
\end{align*}
Fibre families lie at $\pm\gamma$ from local axis $\mathbf e_1$ in the $\mathbf e_1$--$\mathbf e_2$ plane and resist only
tension ($E_a\ge 0$). FEBio's $c$ equals ABAQUS $2C_{10}$; FEBio's HGO volumetric term equals ABAQUS GOH with $D=2/k$.
The sources' conventions and the conversions applied are listed per parameter in Table~\ref{tab:mat}.
"""


def build_report(run_path: str, compile_pdf: bool = True) -> dict:
    J = lambda *p: json.load(open(os.path.join(run_path, *p))) if os.path.exists(os.path.join(run_path, *p)) else None
    man, st, cfg = J("manifest.json"), J("status.json"), J("config.resolved.json")
    summ, conv = J("results", "summary.json"), J("solver", "convergence.json")
    lay, sock, load = J("mesh", "layer_report.json"), J("mesh", "socket_report.json"), J("solver", "load_report.json")
    a_src, a_prep, asm, anat = J("geometry", "audit_source.json"), J("geometry", "audit_prepared.json"), J("geometry", "assembly_audit.json"), J("geometry", "anatomy.json")
    figs = J("figures", "index.json") or {}
    blocks = J("mesh", "blocks.json")
    verif = J("verification.json")
    rep = os.path.join(run_path, "report")
    os.makedirs(rep, exist_ok=True)
    state = (st or {}).get("state", "unknown").upper()
    L = []
    L.append(r"""\documentclass[10pt,a4paper]{article}
\usepackage[margin=18mm]{geometry}\usepackage{graphicx,booktabs,longtable,float,amsmath,xcolor,hyperref}
\hypersetup{colorlinks,linkcolor=blue}\setlength{\parindent}{0pt}\setlength{\parskip}{4pt}
\begin{document}""")
    L.append(r"\begin{center}{\LARGE Residual-limb FEA run report}\\[2mm]{\large " + esc(man.get("run_name")) + r"}\\[1mm]"
             + r"run id \texttt{" + esc(man["run_id"]) + r"} \quad created " + esc(man["created_utc"]) + r"\end{center}")
    col = "green!40!black" if state == "CONVERGED" else "red"
    L.append(r"\fcolorbox{" + col + r"}{white}{\parbox{0.97\textwidth}{\textbf{Job state: \textcolor{" + col + "}{" + esc(state) + r"}}"
             + (r" --- last accepted time " + _fmt((st or {}).get("last_converged_time", 0)) if state != "CONVERGED" else "")
             + r". Solver: FEBio " + esc(man["solver"].get("version")) + r". Load model: " + esc(man.get("load_model"))
             + r". \textbf{Numerical verification only: this report does not establish clinical accuracy, socket grip or tissue safety.} "
             + (r"Geometry: \textbf{SYNTHETIC} (no patient data)." if not man.get("patient_geometry") else r"Geometry: supplied private limb model; layers are assumptions, not segmentation.")
             + "}}")
    if man.get("warnings"):
        L.append(r"\section*{Warnings and labelled assumptions}\begin{itemize}" + "".join(r"\item " + esc(w) for w in man["warnings"])
                 + "".join(r"\item " + esc(a) for a in (load or {}).get("assumptions", [])) + r"\end{itemize}")
    # --- provenance
    L.append(r"\section{Inputs and provenance}")
    L.append(_tab([[k, v.get("path", ""), v.get("sha256", "")[:16] + "..."] for k, v in man["inputs"].items()], ["input", "path", "sha256 (prefix)"], "lp{9cm}l"))
    L.append(r"Repository commit: \texttt{" + esc(man["repository"].get("commit")) + r"} (dirty worktree: " + esc(man["repository"].get("dirty_worktree"))
             + r"). Solver input \texttt{solver/model.feb} sha256 \texttt{" + esc(man["solver_input"]["sha256"][:20]) + r"...}. "
             + r"Execution command: \texttt{" + esc(" ".join(man.get("execution", {}).get("command", []))) + "}.")
    L.append(r"\subsection{Geometry audit}")
    rows = []
    for k, a in (a_src or {}).items():
        rows.append([k, a.get("status"), a.get("source_unit", ""), "x".join(_fmt(x, 4) for x in a.get("dimensions_mm", [])),
                     _fmt(a.get("signed_volume_mm3", 0), 6), a.get("watertight"), a.get("connected_components"), a.get("self_intersecting_pairs"), ", ".join(a.get("defects", []))])
    L.append(_tab(rows, ["body", "status", "source unit", "dims mm", "vol mm3", "closed", "comp.", "self-int.", "defects"], "llllllllp{3cm}"))
    if anat:
        L.append(r"\textbf{Anatomy roles (as configured):} " + esc(", ".join(f"{k}: {v}" for k, v in anat["roles"].items())) + ". ")
        if anat.get("missing_anatomy"):
            L.append(r"\textbf{\textcolor{red}{Missing expected anatomy: " + esc(", ".join(anat["missing_anatomy"])) + r"}} --- not created or substituted.")
        for n in anat.get("notes", []):
            L.append(r"\par " + esc(n))
    if asm:
        rows = [[k, _fmt(v.get("fraction_vertices_inside_envelope", "")), _fmt(v.get("min_clearance_to_envelope_mm", "")), v.get("intersects_envelope_tri_pairs", "")]
                for k, v in asm.get("bodies", {}).items()]
        L.append(_tab(rows, ["body", "fraction inside skin", "min clearance mm", "skin intersections"]))
    # --- model
    L.append(r"\section{Model}")
    if lay:
        L.append(esc(lay["assumption"]) + " Skin thickness " + esc(json.dumps({k: round(v, 3) for k, v in lay["thickness_skin_mm"].items()})) +
                 " mm; fat thickness " + esc(json.dumps({k: round(v, 3) for k, v in lay["thickness_fat_mm"].items()})) + " mm. "
                 + f"{lay['vertices_thickness_clamped']} skin vertices had layer thickness reduced near bone. Minimum wedge scaled Jacobian "
                 + _fmt(lay["min_wedge_scaled_jacobian"]) + ".")
    if blocks:
        L.append(_tab([[b["id"], b["name"], b["etype"], b["material"], b["region"], b["formulation"]] for b in blocks],
                      ["id", "block", "element", "material", "region", "formulation"], "rlllp{5cm}l"))
    L.append(EQUATIONS)
    mats = (J("materials.json") or {})
    if mats:
        rows = []
        for k, m in mats.get("materials", {}).items():
            for p, pv_ in m["provenance"].items():
                rows.append([k, m["febio_type"], p, _fmt(m["params"][p]), pv_["source"], pv_["note"][:150]])
        L.append(r"\begin{small}" + _tab(rows, ["material", "FEBio model", "param", "value (MPa/deg/-)", "source", "note"], "llllp{2.2cm}p{6.5cm}") + r"\end{small}\label{tab:mat}")
        refs = mats.get("references", {})
        L.append(r"\textbf{References:} " + " ".join(r"[" + esc(k) + "] " + esc(v["citation"]) + (" doi:" + esc(v["doi"]) if v.get("doi") else "") + "." for k, v in refs.items()))
    if sock:
        L.append(r"\subsection{Socket (provisional parametric model)}" + esc(sock.get("status", "")) + ".")
    # --- loading
    L.append(r"\section{Loading}")
    if load:
        t = load["tension"]
        L.append(f"Applied strand tension {esc(t['applied_N'])} N (basis: {esc(t['basis'])}). Rated tension: {esc(t['rated_N'])}. "
                 f"Output levels: {esc([round(x, 3) for x in t['levels_N']])} N. " + esc(load.get("routing_status", "")))
        if "free_body" in figs:
            L.append(_fig(figs["free_body"], "Lace route, guide forces and petal resultants (free-body check: sums must vanish)."))
    # --- convergence
    L.append(r"\section{Solver and convergence}")
    L.append(f"Tolerances: dtol={esc(cfg['solver']['dtol'])}, etol={esc(cfg['solver']['etol'])}, rtol={esc(cfg['solver']['rtol'])}; "
             f"Newton with {esc(cfg['solver']['qn_method'])} (max\\_ups={esc(cfg['solver']['max_ups'])}); PARDISO; automatic time stepping with must points; "
             f"max retries {esc(cfg['solver']['max_retries'])}. Retries in log: {esc((conv or {}).get('n_retries'))}; negative-Jacobian events: {esc((conv or {}).get('negative_jacobian_events'))}.")
    if conv and conv.get("time_steps"):
        rows = [[s["step"], _fmt(s["time"]), "yes" if s["converged"] else "NO", s["iterations"],
                 _fmt(s["norms"][-1]["residual"][1]) if s["norms"] and "residual" in s["norms"][-1] else "",
                 _fmt(s["norms"][-1]["energy"][1]) if s["norms"] and "energy" in s["norms"][-1] else ""] for s in conv["time_steps"]]
        L.append(r"\begin{small}" + _tab(rows, ["step", "time", "converged", "iterations", "final residual", "final energy norm"]) + r"\end{small}")
    if "convergence" in figs:
        L.append(_fig(figs["convergence"], "Newton convergence history."))
    # --- results
    L.append(r"\section{Results (actual FEBio fields)}")
    for key, cap in (("full_limb", "Full limb, undeformed and deformed at 1$\\times$ (true scale), with sagittal section."),
                     ("full_limb_exaggerated", "\\textbf{EXAGGERATED} deformation (labelled scale) for visual inspection only."),
                     ("sections", "Sections through the deformed limb: displacement, maximum principal log strain, hydrostatic pressure, maximum principal Cauchy stress (cell data, no smoothing). Peak values may be mesh-dependent and are not injury predictions."),
                     ("contact", "Skin contact pressure, tangential traction, derived stick/slip status and gap."),
                     ("regions", "Material regions and skin fibre orientation."), ("quality", "Mesh element quality."),
                     ("force_displacement", "Force--displacement / load--response curves."), ("pullout", "Pull-out: reaction versus relative socket displacement.")):
        if key in figs:
            L.append(_fig(figs[key], cap))
    if summ and summ.get("states"):
        rows = []
        for s in summ["states"]:
            e = s.get("equilibrium", {})
            d = s.get("deformation", {})
            c = s.get("contact", {})
            rows.append([_fmt(s["time"]), _fmt(s["max_tissue_displacement_mm"]), _fmt(d.get("min_J", "")), d.get("inverted_elements", ""),
                         _fmt(c.get("normal_force_N", "")), _fmt(c.get("max_pressure_MPa", 0) * 1000 if c else ""), _fmt(c.get("slip_fraction_of_contact_area", "")),
                         _fmt(e.get("force_imbalance_relative", "")), _fmt(s.get("strain_energy_Nmm", {}).get("total_Nmm", ""))])
        L.append(_tab(rows, ["t", "max |u| tissue mm", "min J", "inverted", "contact N", "max p kPa", "slip area frac.", "force imbalance", "strain energy Nmm"]))
        tight = [x for x in summ["states"] if x["time"] <= 1.0 + 1e-9]
        if tight and tight[-1].get("petal_free_body"):
            fb = tight[-1]["petal_free_body"]
            rows = [[k, *[_fmt(np.linalg.norm(v[q]), 3) for q in ("cable_N", "hinge_on_petal_N", "pad_contact_N", "stabiliser_N", "residual_N")],
                     _fmt(v["residual_over_cable"], 3)] for k, v in sorted(fb.items())]
            L.append(r"\paragraph{Petal free bodies at t = " + _fmt(tight[-1]["time"]) + r" (force magnitudes, N).} "
                     r"Pad contact force is the facet-averaged traction on the skin zone times the deformed facet area (an approximate "
                     r"post-processed integral); the global reaction balance in the table above is the exact check.")
            L.append(_tab(rows, ["petal", "|cable|", "|hinge|", "|pad contact|", "|stabiliser|", "|residual|", "residual/cable"]))
        last = summ["states"][-1]
        if last.get("region_extremes"):
            rows = [[k, *[_fmt(v[x]) for x in ("max_disp_mm", "max_logstrain", "min_logstrain", "max_pressure_MPa", "min_pressure_MPa", "max_von_mises_MPa", "min_J")]]
                    for k, v in last["region_extremes"].items()]
            L.append(r"\begin{small}" + _tab(rows, ["block", "max |u|", "max log e", "min log e", "max p MPa", "min p MPa", "max vM MPa", "min J"]) + r"\end{small}")
    if verif:
        L.append(r"\section{Verification summary}\label{sec:verif}")
        L.append(_tab([[v["name"], v.get("metric", ""), _fmt(v.get("value", "")), _fmt(v.get("tolerance", "")), "PASS" if v.get("pass") else "FAIL"] for v in verif.get("checks", [])],
                      ["check", "metric", "value", "tolerance", "verdict"], "p{6cm}p{4.5cm}lll"))
    L.append(r"\section{Limitations}\begin{itemize}"
             r"\item Numerical verification is separate from experimental or clinical validation; no validation data were supplied."
             r"\item Tissue layers are geometric assumptions; material parameters come from literature on other subjects/sites (see table) or are labelled assumptions."
             r"\item The socket is a provisional parametric model; cable routing is provisional unless the configuration supplies verified routing."
             r"\item Quasi-static, no viscoelasticity, no perfusion, no active muscle, no skin adhesion or liner; friction is Coulomb with a screening coefficient."
             r"\item A converged state does not establish grip; pull-out is evaluated only through the resolved contact tractions and relative motion of this model."
             r"\end{itemize}\end{document}")
    tex = os.path.join(rep, "report.tex")
    open(tex, "w").write("\n".join(L))
    out = {"tex": tex, "pdf": None}
    if compile_pdf and shutil.which("pdflatex"):
        for _ in range(2):
            r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "report.tex"], cwd=rep,
                               capture_output=True, text=True, stdin=subprocess.DEVNULL)
        if os.path.exists(os.path.join(rep, "report.pdf")) and r.returncode == 0:
            out["pdf"] = os.path.join(rep, "report.pdf")
        else:
            out["error"] = r.stdout[-3000:]
    return out
