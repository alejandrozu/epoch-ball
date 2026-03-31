import torch

torch.set_default_dtype(torch.float64)

_LONG_SCALE = 0.2
_TWIST = 0.8
_RADIAL_OFFSET = 1.2
_RADIAL_GROWTH = 0.16
_NORMAL_SCALE = 0.5
_STRENGTH = 1.2
_SUPPORT_X = 18.0
_SUPPORT_Y = 3.5
_X_START = -1.0
_X_CENTER = 13.5


def _vector_field(Q):
    x = Q[:, 0]
    y = Q[:, 1]
    u = _LONG_SCALE * (x - _X_START)
    theta = _TWIST * u
    radius = _RADIAL_OFFSET + _RADIAL_GROWTH * u

    cos_theta = torch.cos(theta)
    sin_theta = torch.sin(theta)

    centerline = torch.stack([radius * cos_theta, radius * sin_theta], dim=1)

    v = torch.stack(
        [
            _RADIAL_GROWTH * cos_theta - radius * _TWIST * sin_theta,
            _RADIAL_GROWTH * sin_theta + radius * _TWIST * cos_theta,
        ],
        dim=1,
    )
    speed = torch.linalg.norm(v, dim=1, keepdim=True).clamp_min(1e-12)
    tangent = v / speed
    normal = torch.stack([-tangent[:, 1], tangent[:, 0]], dim=1)

    accel = torch.stack(
        [
            -2.0 * _RADIAL_GROWTH * _TWIST * sin_theta - radius * (_TWIST ** 2) * cos_theta,
            2.0 * _RADIAL_GROWTH * _TWIST * cos_theta - radius * (_TWIST ** 2) * sin_theta,
        ],
        dim=1,
    )
    tangent_accel = (tangent * accel).sum(dim=1, keepdim=True)
    tangent_u_deriv = (accel - tangent_accel * tangent) / speed
    normal_u_deriv = torch.stack([-tangent_u_deriv[:, 1], tangent_u_deriv[:, 0]], dim=1)

    target = centerline + _NORMAL_SCALE * y[:, None] * normal
    sx = (x - _X_CENTER) / _SUPPORT_X
    sy = y / _SUPPORT_Y
    cutoff = torch.exp(-(sx.pow(8)) - (sy.pow(8)))
    return _STRENGTH * cutoff[:, None] * (target - Q)


def Hamiltonian(Q, P, t):
    return (P * _vector_field(Q)).sum(dim=1)
