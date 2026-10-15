"""Independent reference implementations of the strain-energy functions.

These are written from the published equations (not from FEBio source) and
are used to check FEBio's stresses for homogeneous deformations.  Stress is
obtained from W(C) by complex-step differentiation:  S = 2 dW/dC,
sigma = F S F^T / J.
"""
from __future__ import annotations

import numpy as np


def _iso(C):
    J = np.sqrt(np.linalg.det(C))
    Cb = C * J ** (-2.0 / 3.0)
    return J, Cb


def W_mooney_rivlin(C, c1, c2, k):
    J, Cb = _iso(C)
    I1 = np.trace(Cb)
    I2 = 0.5 * (I1 ** 2 - np.trace(Cb @ Cb))
    return c1 * (I1 - 3) + c2 * (I2 - 3) + 0.5 * k * np.log(J) ** 2


def cauchy_ogden(F, cs, ms, k):
    """Exact principal-frame Cauchy stress for FEBio's uncoupled Ogden form
    W = sum_p c_p/m_p^2 sum_i (l~_i^m_p - 1) + k/2 (ln J)^2:
        sigma_i = (1/J) sum_p (c_p/m_p) (l~_i^m_p - 1/3 sum_j l~_j^m_p) + k ln(J)/J
    with principal directions of b = F F^T."""
    F = np.asarray(F, float)
    J = np.linalg.det(F)
    b = F @ F.T
    w, V = np.linalg.eigh(b)
    lam = np.sqrt(w) * J ** (-1.0 / 3.0)
    s = np.zeros(3)
    for c, m in zip(cs, ms):
        lm = lam ** m
        s += (c / m) * (lm - lm.sum() / 3.0)
    s = s / J + k * np.log(J) / J
    return V @ np.diag(s) @ V.T


def W_hgo(C, c, k1, k2, kappa, gamma_deg, k, e1, e2):
    """Gasser-Ogden-Holzapfel (2006), uncoupled, deviatoric invariants,
    two fibre families at +/-gamma from e1 in the e1-e2 plane, U(J) of the
    ABAQUS GOH form k/2((J^2-1)/2 - ln J); fibres active only when E_a > 0."""
    J, Cb = _iso(C)
    I1 = np.trace(Cb)
    g = np.deg2rad(gamma_deg)
    W = 0.5 * c * (I1 - 3)
    for s in (1, -1):
        a = np.cos(g) * e1 + s * np.sin(g) * e2
        I4 = a @ Cb @ a
        E = kappa * (I1 - 3) + (1 - 3 * kappa) * (I4 - 1)
        if E.real > 0:
            W = W + k1 / (2 * k2) * (np.exp(k2 * E ** 2) - 1)
    return W + 0.5 * k * ((J ** 2 - 1) / 2 - np.log(J))


def W_neo_hookean(C, E, v):
    """FEBio coupled neo-Hookean: W = mu/2 (I1-3) - mu ln J + lam/2 (ln J)^2."""
    mu = E / (2 * (1 + v)); lam = E * v / ((1 + v) * (1 - 2 * v))
    J = np.sqrt(np.linalg.det(C))
    return 0.5 * mu * (np.trace(C) - 3) - mu * np.log(J) + 0.5 * lam * np.log(J) ** 2


def cauchy_from_W(Wfun, F, h=1e-20):
    """sigma = (1/J) F S F^T with S = 2 dW/dC (complex step on symmetric C)."""
    F = np.asarray(F, float)
    C = F.T @ F
    S = np.zeros((3, 3))
    for i in range(3):
        for j in range(i, 3):
            dC = np.zeros((3, 3), complex)
            if i == j:
                dC[i, i] = 1j * h
                S[i, i] = 2 * Wfun(C + dC).imag / h
            else:
                dC[i, j] = dC[j, i] = 1j * h
                # W depends on C_ij and C_ji symmetrically: dW/dC_ij = 0.5 * d/dh
                S[i, j] = S[j, i] = Wfun(C + dC).imag / h
    J = np.linalg.det(F)
    return F @ S @ F.T / J
