import math
import torch

torch.set_default_dtype(torch.float64)

_WINDOW_NORM = 0.007029858406609658
_MAIN_LEVELS = [[{'level': 0, 'q_mid': 1.5, 'support_q': 3.3499999999999996}, {'level': 0, 'q_mid': 7.5, 'support_q': 3.3499999999999996}, {'level': 0, 'q_mid': 13.5, 'support_q': 3.3499999999999996}, {'level': 0, 'q_mid': 19.5, 'support_q': 3.3499999999999996}, {'level': 0, 'q_mid': 25.5, 'support_q': 3.3499999999999996}], [{'level': 1, 'q_mid': 4.5, 'support_q': 5.449999999999999}, {'level': 1, 'q_mid': 16.5, 'support_q': 5.449999999999999}], [{'level': 2, 'q_mid': 10.5, 'support_q': 9.649999999999999}], [{'level': 3, 'q_mid': 18.0, 'support_q': 11.75}]]
_LEVEL_PARAMS = [{'lam': -0.35, 'theta': 0.0}, {'lam': -0.4375, 'theta': -0.0}, {'lam': -0.546875, 'theta': 0.0}, {'lam': -0.68359375, 'theta': -0.0}]
_ETA = 0.0
_SUPPORT_R = 1.5
_SUPPORT_P = 2.6
_TRANSLATE = torch.tensor([-14.397002460683188, -0.10265982233563142, -0.3776675923437038, -0.013737259425547538], dtype=torch.float64)


def _bump_unit_interval(u):
    if not (0.0 < u < 1.0):
        return 0.0
    return math.exp(-1.0 / (u * (1.0 - u))) / _WINDOW_NORM


def _stage_weight(t, n_stages):
    if not (0.0 < t < 1.0):
        return -1, 0.0
    s = t * n_stages
    stage = int(math.floor(s))
    if stage >= n_stages:
        stage = n_stages - 1
    u = s - stage
    return stage, n_stages * _bump_unit_interval(u)


def _level_hamiltonian(Q, P, nodes, params):
    q1 = Q[:, 0]
    q2 = Q[:, 1]
    p1 = P[:, 0]
    p2 = P[:, 1]
    lam = params["lam"]
    theta = params["theta"]
    H = 0.0 * (q1 + q2 + p1 + p2)
    for node in nodes:
        z = q1 - node["q_mid"]
        u = z / node["support_q"]
        v = q2 / _SUPPORT_R
        w1 = p1 / _SUPPORT_P
        w2 = p2 / _SUPPORT_P
        cutoff = torch.exp(-(u.pow(4)) - (v.pow(4)) - (w1.pow(4)) - (w2.pow(4)))
        gate = lam * (z * p1 + _ETA * q2 * p2) + theta * (z * p2 - q2 * p1)
        H = H + cutoff * gate
    return H


def Hamiltonian(Q, P, t):
    n_main_stages = len(_MAIN_LEVELS)
    total_stages = n_main_stages + 1
    stage, weight = _stage_weight(t, total_stages)
    H = 0.0 * (Q[:, 0] + P[:, 0])
    if weight == 0.0:
        return H
    if stage < n_main_stages:
        return weight * _level_hamiltonian(Q, P, _MAIN_LEVELS[stage], _LEVEL_PARAMS[stage])
    dq = _TRANSLATE[:2].to(Q.device)
    dp = _TRANSLATE[2:].to(P.device)
    return weight * ((P * dq).sum(dim=1) - (Q * dp).sum(dim=1))
