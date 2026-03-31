import math
import torch

torch.set_default_dtype(torch.float64)

_WINDOW_NORM = 0.007029858406609658
_MAIN_LEVELS = [[{'level': 0, 'q_mid': 1.5, 'support_q': 3.2, 'sign_width': 1.025, 'target_sign_width': 0.97, 'shift_q2': 1.4, 'shift_p2': 0.0, 'collapse_q1': 1.2000000000000002}, {'level': 0, 'q_mid': 7.5, 'support_q': 3.2, 'sign_width': 1.025, 'target_sign_width': 0.97, 'shift_q2': 1.4, 'shift_p2': 0.0, 'collapse_q1': 1.2000000000000002}, {'level': 0, 'q_mid': 13.5, 'support_q': 3.2, 'sign_width': 1.025, 'target_sign_width': 0.97, 'shift_q2': 1.4, 'shift_p2': 0.0, 'collapse_q1': 1.2000000000000002}, {'level': 0, 'q_mid': 19.5, 'support_q': 3.2, 'sign_width': 1.025, 'target_sign_width': 0.97, 'shift_q2': 1.4, 'shift_p2': 0.0, 'collapse_q1': 1.2000000000000002}, {'level': 0, 'q_mid': 25.5, 'support_q': 3.2, 'sign_width': 1.025, 'target_sign_width': 0.97, 'shift_q2': 1.4, 'shift_p2': 0.0, 'collapse_q1': 1.2000000000000002}], [{'level': 1, 'q_mid': 4.5, 'support_q': 5.15, 'sign_width': 1.7000000000000002, 'target_sign_width': 1.3549999999999998, 'shift_q2': 1.4630840896290471, 'shift_p2': 1.5064477908889977, 'collapse_q1': 2.4000000000000004}, {'level': 1, 'q_mid': 16.5, 'support_q': 5.15, 'sign_width': 1.7000000000000002, 'target_sign_width': 1.3549999999999998, 'shift_q2': 1.4630840896290471, 'shift_p2': 1.5064477908889977, 'collapse_q1': 2.4000000000000004}], [{'level': 2, 'q_mid': 10.5, 'support_q': 9.05, 'sign_width': 3.0500000000000003, 'target_sign_width': 1.9325, 'shift_q2': -0.09197849524905977, 'shift_p2': 3.148656849580741, 'collapse_q1': 4.800000000000001}], [{'level': 3, 'q_mid': 18.0, 'support_q': 11.0, 'sign_width': 3.725, 'target_sign_width': 2.79875, 'shift_q2': -3.484185305932386, 'shift_p2': 3.1915635281041865, 'collapse_q1': 6.0}]]
_P1_SCALE = 3.0
_RETURN_FRAC = 0.0
_TRANSLATE = torch.tensor([-13.057398482896478, 1.261112656068954, -0.31975152903014215, -0.5032082308087749], dtype=torch.float64)


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


def _stage_a(Q, P, nodes, shift_scale):
    q1 = Q[:, 0]
    q2 = Q[:, 1]
    p1 = P[:, 0]
    p2 = P[:, 1]
    H = 0.0 * (q1 + q2 + p1 + p2)
    for node in nodes:
        shift_q2 = shift_scale * node["shift_q2"]
        shift_p2 = shift_scale * node["shift_p2"]
        z = q1 - node["q_mid"]
        a = z / node["support_q"]
        b = p1 / _P1_SCALE
        window = torch.exp(-(a.pow(4)) - (b.pow(4)))
        side = torch.tanh(z / node["sign_width"])
        coeff = window * side
        H = H + coeff * (shift_q2 * p2 - shift_p2 * q2)
    return H


def _stage_b(Q, P, nodes):
    q1 = Q[:, 0]
    q2 = Q[:, 1]
    p1 = P[:, 0]
    p2 = P[:, 1]
    H = 0.0 * (q1 + q2 + p1 + p2)
    for node in nodes:
        z = q1 - node["q_mid"]
        a = z / node["support_q"]
        window = torch.exp(-(a.pow(4)))
        shift_norm = math.hypot(node["shift_q2"], node["shift_p2"])
        uq = node["shift_q2"] / shift_norm
        up = node["shift_p2"] / shift_norm
        proj = uq * q2 + up * p2
        side = torch.tanh(proj / node["target_sign_width"])
        H = H - node["collapse_q1"] * window * side * p1
    return H


def Hamiltonian(Q, P, t):
    n_main_stages = 3 * len(_MAIN_LEVELS)
    total_stages = n_main_stages + 1
    stage, weight = _stage_weight(t, total_stages)
    H = 0.0 * (Q[:, 0] + P[:, 0])
    if weight == 0.0:
        return H
    if stage < n_main_stages:
        level = stage // 3
        stage_type = stage % 3
        if stage_type == 0:
            return weight * _stage_a(Q, P, _MAIN_LEVELS[level], 1.0)
        if stage_type == 1:
            return weight * _stage_b(Q, P, _MAIN_LEVELS[level])
        return weight * _stage_a(Q, P, _MAIN_LEVELS[level], -_RETURN_FRAC)
    dq = _TRANSLATE[:2].to(Q.device)
    dp = _TRANSLATE[2:].to(P.device)
    return weight * ((P * dq).sum(dim=1) - (Q * dp).sum(dim=1))
