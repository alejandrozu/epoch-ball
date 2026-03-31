import math
import torch

torch.set_default_dtype(torch.float64)

_WINDOW_NORM = 0.007029858406609658
_MAIN_LEVELS = [[{'level': 0, 'q_mid': 1.5, 'support_q': 2.75, 'sign_width': 1.025, 'target_sign_width': 0.6400000000000001, 'shift_q2': 0.8, 'shift_p2': 0.0, 'collapse_q1': 0.675}, {'level': 0, 'q_mid': 7.5, 'support_q': 2.75, 'sign_width': 1.025, 'target_sign_width': 0.6400000000000001, 'shift_q2': 0.8, 'shift_p2': 0.0, 'collapse_q1': 0.675}, {'level': 0, 'q_mid': 13.5, 'support_q': 2.75, 'sign_width': 1.025, 'target_sign_width': 0.6400000000000001, 'shift_q2': 0.8, 'shift_p2': 0.0, 'collapse_q1': 0.675}, {'level': 0, 'q_mid': 19.5, 'support_q': 2.75, 'sign_width': 1.025, 'target_sign_width': 0.6400000000000001, 'shift_q2': 0.8, 'shift_p2': 0.0, 'collapse_q1': 0.675}, {'level': 0, 'q_mid': 25.5, 'support_q': 2.75, 'sign_width': 1.025, 'target_sign_width': 0.6400000000000001, 'shift_q2': 0.8, 'shift_p2': 0.0, 'collapse_q1': 0.675}], [{'level': 1, 'q_mid': 4.5, 'support_q': 4.25, 'sign_width': 1.7000000000000002, 'target_sign_width': 0.75, 'shift_q2': 0.4535961214255773, 'shift_p2': 0.8912073600614354, 'collapse_q1': 1.35}, {'level': 1, 'q_mid': 16.5, 'support_q': 4.25, 'sign_width': 1.7000000000000002, 'target_sign_width': 0.75, 'shift_q2': 0.4535961214255773, 'shift_p2': 0.8912073600614354, 'collapse_q1': 1.35}], [{'level': 2, 'q_mid': 10.5, 'support_q': 7.25, 'sign_width': 3.0500000000000003, 'target_sign_width': 0.8875, 'shift_q2': -0.7356263965691823, 'shift_p2': 1.0106205047744876, 'collapse_q1': 2.7}], [{'level': 3, 'q_mid': 18.0, 'support_q': 8.75, 'sign_width': 3.725, 'target_sign_width': 1.0593750000000002, 'shift_q2': -1.5429371404826013, 'shift_p2': -0.24647764709882602, 'collapse_q1': 3.375}]]
_P1_SCALE = 3.0
_TRANSLATE = torch.tensor([-13.444575630541275, 0.39837567737931584, -0.005217264936792967, 0.07131437533152989], dtype=torch.float64)


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


def _stage_a(Q, P, nodes):
    q1 = Q[:, 0]
    q2 = Q[:, 1]
    p1 = P[:, 0]
    p2 = P[:, 1]
    H = 0.0 * (q1 + q2 + p1 + p2)
    for node in nodes:
        z = q1 - node["q_mid"]
        a = z / node["support_q"]
        b = p1 / _P1_SCALE
        window = torch.exp(-(a.pow(4)) - (b.pow(4)))
        side = torch.tanh(z / node["sign_width"])
        coeff = window * side
        H = H + coeff * (node["shift_q2"] * p2 - node["shift_p2"] * q2)
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
        H = H + node["collapse_q1"] * window * side * p1
    return H


def Hamiltonian(Q, P, t):
    n_main_stages = 2 * len(_MAIN_LEVELS)
    total_stages = n_main_stages + 1
    stage, weight = _stage_weight(t, total_stages)
    H = 0.0 * (Q[:, 0] + P[:, 0])
    if weight == 0.0:
        return H
    if stage < n_main_stages:
        level = stage // 2
        if stage % 2 == 0:
            return weight * _stage_a(Q, P, _MAIN_LEVELS[level])
        return weight * _stage_b(Q, P, _MAIN_LEVELS[level])
    dq = _TRANSLATE[:2].to(Q.device)
    dp = _TRANSLATE[2:].to(P.device)
    return weight * ((P * dq).sum(dim=1) - (Q * dp).sum(dim=1))
