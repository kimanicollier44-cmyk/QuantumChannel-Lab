#!/usr/bin/env python3
"""QuantumChannel-Lab - educational single-qubit quantum channel lab.

by Darion Kellar

Pure Python standard library. Zero dependencies. Zero network.
Education only - and the central lesson is stated up front:

  * A quantum channel is the honest bookkeeping of what noise does to
    a qubit: a set of Kraus operators, the affine map it induces on
    the Bloch ball, the average fidelity it leaves behind, and the
    point at which it destroys every entangled state fed through it.
  * This lab computes those four views from the same Kraus operators
    and grades each view against the others and against closed forms.
  * It models ONE qubit and Markovian noise only. It says nothing
    about any real device: no real channel was measured, no real
    qubit was protected, and no capacity of any real link is computed
    or claimed anywhere in this file.

Conventions, stated once and used everywhere:

  * Density matrices are 2x2 complex matrices; rho = (I + r.sigma)/2.
  * Depolarizing means E(rho) = (1 - p) rho + p * I/2 exactly, with
    Kraus weights (1 - 3p/4) on I and (p/4) on each Pauli.
  * Amplitude damping gamma sends |1> to |0>; phase damping lambda
    multiplies coherences by sqrt(1 - lambda).
  * The Choi state is the normalised J = (E x I)(|Phi+><Phi+|),
    trace 1, channel acting on the FIRST qubit.
  * Average gate fidelity is against the identity channel.

Every real-world deployment use is a permanent HOLD gate - with
Gemini Argon 4 Extended Thinking as DECIDER (names-only,
tighten-only), Muse as advisor (never overrides HOLD), and Google
Flow plan-only seats.
"""
import cmath
import hashlib
import json
import math
import sys

# ---------------------------------------------------------------- constants
VERSION = "1.0.0"

# Labelled channel parameters used by the demo book. Teaching values,
# not measurements of any device.
DEMO_PARAMS = {
    "bit_flip": 0.2,
    "phase_flip": 0.2,
    "bit_phase_flip": 0.2,
    "depolarizing": 0.3,
    "amplitude_damping": 0.5,
    "phase_damping": 0.5,
}
PARAM_GRID = (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)

# ---------------------------------------------------------------- HOLD gates
HOLD_GATES = {
    "quantum_network_deployment_claim": "HOLD - no deploying or dimensioning a real quantum network from this single-qubit model",
    "error_correction_product_claim": "HOLD - no error-correction product, code, or fault-tolerance claim is made or certified here",
    "vendor_channel_claim": "HOLD - no vendor channel, link, or device fidelity figure is stated, endorsed, or compared",
    "hardware_specification_claim": "HOLD - no hardware specification (T1, T2, gate error) for any real device is asserted",
    "real_device_characterisation": "HOLD - no characterising a real device's noise from this teaching model",
    "capacity_or_rate_claim": "HOLD - no channel capacity, key rate, or throughput figure for any real link is computed or claimed",
    "qkd_security_claim": "HOLD - no QKD or cryptographic security claim follows from a toy channel model",
    "entanglement_distribution_claim": "HOLD - no claim that entanglement can be distributed over any real distance or hardware",
    "teleportation_service_claim": "HOLD - no teleportation service, fidelity guarantee, or network offering",
    "medical_or_sensing_claim": "HOLD - no medical, sensing, or metrology claim from channel mathematics",
    "investment_claim": "HOLD - no investment claim; a lab model is not a prospectus",
    "certification_claim": "HOLD - no certification that any system meets any standard",
    "live_hardware_control": "HOLD - no control of real qubits, links, or laboratory hardware, ever",
    "model_as_measurement_claim": "HOLD - a computed number is a statement about the model at stated inputs, never a measurement",
}


def gate(action):
    """Return 'HOLD' if the action is gated, else 'PROCEED' (sandbox sim)."""
    return "HOLD" if action in HOLD_GATES else "PROCEED"


# ================================================================ 1. Matrix kit
# 2x2 (and 4x4) complex matrices as plain lists of lists. Everything
# the lab needs - products, daggers, Kronecker products, traces - is
# written out here so every number can be read and diffed.

def mat(rows):
    return [[complex(v) for v in row] for row in rows]


def eye(n):
    return [[complex(1 if i == j else 0) for j in range(n)] for i in range(n)]


def zeros(n):
    return [[0j] * n for _ in range(n)]


def madd(A, B):
    return [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def mscale(A, s):
    return [[a * s for a in row] for row in A]


def mmul(A, B):
    n, m, k = len(A), len(B[0]), len(B)
    out = zeros(n)
    for i in range(n):
        for j in range(m):
            out[i][j] = sum(A[i][t] * B[t][j] for t in range(k))
    return [[out[i][j] for j in range(m)] for i in range(n)]


def dagger(A):
    return [[A[j][i].conjugate() for j in range(len(A))]
            for i in range(len(A[0]))]


def mtrace(A):
    return sum(A[i][i] for i in range(len(A)))


def kron(A, B):
    """Kronecker product: out[(a,b),(c,d)] = A[a][c] * B[b][d]."""
    na, nb = len(A), len(B)
    out = zeros(na * nb)
    for a in range(na):
        for c in range(na):
            for b in range(nb):
                for d in range(nb):
                    out[a * nb + b][c * nb + d] = A[a][c] * B[b][d]
    return out


def frobenius(A):
    return math.sqrt(sum(abs(v) ** 2 for row in A for v in row))


def mat_max_dev(A, B):
    return max(abs(a - b) for ra, rb in zip(A, B) for a, b in zip(ra, rb))


PAULI_I = mat([[1, 0], [0, 1]])
PAULI_X = mat([[0, 1], [1, 0]])
PAULI_Y = mat([[0, -1j], [1j, 0]])
PAULI_Z = mat([[1, 0], [0, -1]])
PAULIS = (PAULI_X, PAULI_Y, PAULI_Z)


def hermitian_eigenvalues(A):
    """Eigenvalues of a Hermitian matrix by cyclic Jacobi rotations.

    Each rotation first phases the pivot element real, then applies
    the stable real Jacobi angle derived in-line below. Sweeps run
    until the off-diagonal Frobenius mass is below 1e-14 or 100
    sweeps pass. Returns eigenvalues sorted ascending.
    """
    n = len(A)
    a = [[A[i][j] for j in range(n)] for i in range(n)]
    for _ in range(100):
        off = math.sqrt(sum(abs(a[p][q]) ** 2
                            for p in range(n) for q in range(p + 1, n)))
        if off < 1e-14:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                b = abs(apq)
                if b < 1e-15:
                    continue
                app, aqq = a[p][p].real, a[q][q].real
                alpha = cmath.phase(apq)
                # Rotation U (identity except the p,q block):
                #   U_pp = c, U_qp = s e^{-i alpha}, U_pq = -s,
                #   U_qq = c e^{-i alpha}
                # makes (U^dagger A U)_pq real and zero when t = s/c
                # solves b t^2 - (aqq - app) t - b = 0 (small root).
                tau = (aqq - app) / (2.0 * b)
                if tau == 0.0:
                    t = -1.0
                else:
                    t = -math.copysign(
                        1.0 / (abs(tau) + math.sqrt(tau * tau + 1.0)), tau)
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                e_neg = cmath.exp(-1j * alpha)
                U = eye(n)
                U[p][p] = c
                U[q][q] = c * e_neg
                U[q][p] = s * e_neg
                U[p][q] = -s
                a = mmul(dagger(U), mmul(a, U))
    return sorted(a[i][i].real for i in range(n))


# ================================================================ 2. Kraus library
# Each builder returns the Kraus list for one channel family at one
# parameter value. Completeness sum K^dagger K = I is asserted for
# every family on the whole parameter grid in CHECKS.

def kraus_bit_flip(p):
    return [mscale(PAULI_I, math.sqrt(1.0 - p)),
            mscale(PAULI_X, math.sqrt(p))]


def kraus_phase_flip(p):
    return [mscale(PAULI_I, math.sqrt(1.0 - p)),
            mscale(PAULI_Z, math.sqrt(p))]


def kraus_bit_phase_flip(p):
    return [mscale(PAULI_I, math.sqrt(1.0 - p)),
            mscale(PAULI_Y, math.sqrt(p))]


def kraus_depolarizing(p):
    # E(rho) = (1 - p) rho + p I/2 exactly (stated convention).
    return [mscale(PAULI_I, math.sqrt(1.0 - 0.75 * p)),
            mscale(PAULI_X, math.sqrt(p / 4.0)),
            mscale(PAULI_Y, math.sqrt(p / 4.0)),
            mscale(PAULI_Z, math.sqrt(p / 4.0))]


def kraus_amplitude_damping(gamma):
    return [mat([[1, 0], [0, math.sqrt(1.0 - gamma)]]),
            mat([[0, math.sqrt(gamma)], [0, 0]])]


def kraus_phase_damping(lam):
    return [mat([[1, 0], [0, math.sqrt(1.0 - lam)]]),
            mat([[0, 0], [0, math.sqrt(lam)]])]


CHANNEL_BUILDERS = {
    "bit_flip": kraus_bit_flip,
    "phase_flip": kraus_phase_flip,
    "bit_phase_flip": kraus_bit_phase_flip,
    "depolarizing": kraus_depolarizing,
    "amplitude_damping": kraus_amplitude_damping,
    "phase_damping": kraus_phase_damping,
}


def channel_kraus(name, param):
    return CHANNEL_BUILDERS[name](param)


def apply_channel(op, kraus):
    """E(op) = sum_k K op K^dagger, for any 2x2 operator op."""
    out = zeros(2)
    for K in kraus:
        out = madd(out, mmul(K, mmul(op, dagger(K))))
    return out


def completeness_deviation(kraus):
    """max |sum K^dagger K - I| over entries."""
    total = zeros(2)
    for K in kraus:
        total = madd(total, mmul(dagger(K), K))
    return mat_max_dev(total, PAULI_I)


def compose_kraus(outer, inner):
    """Kraus list for E_outer after E_inner: {K_o K_i}."""
    return [mmul(Ko, Ki) for Ko in outer for Ki in inner]


# ================================================================ 3. Bloch picture
# rho = (I + r.sigma)/2. A channel acts affinely: r -> M r + t with
#   M[i][j] = (1/2) Tr[ sigma_i E(sigma_j) ],  t[i] = (1/2) Tr[ sigma_i E(I) ]
# computed by applying the Kraus operators to the Pauli operators
# themselves (apply_channel works on any operator, not just states).

def bloch_vector(rho):
    return [mtrace(mmul(s, rho)).real for s in PAULIS]


def rho_from_bloch(r):
    out = mscale(PAULI_I, 0.5)
    for i, s in enumerate(PAULIS):
        out = madd(out, mscale(s, 0.5 * r[i]))
    return out


def bloch_map(kraus):
    """Return (M, t): the 3x3 real matrix and 3-vector of the channel."""
    M = [[0.0] * 3 for _ in range(3)]
    t = [0.0] * 3
    e_i = apply_channel(PAULI_I, kraus)
    for i, s in enumerate(PAULIS):
        t[i] = 0.5 * mtrace(mmul(s, e_i)).real
        for j in range(3):
            e_j = apply_channel(PAULIS[j], kraus)
            M[i][j] = 0.5 * mtrace(mmul(s, e_j)).real
    return M, t


def apply_bloch(M, t, r):
    return [sum(M[i][j] * r[j] for j in range(3)) + t[i] for i in range(3)]


def fibonacci_sphere(n):
    """Deterministic near-uniform unit vectors (golden-angle spiral)."""
    pts = []
    golden = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n):
        z = 1.0 - (2.0 * i + 1.0) / n
        rad = math.sqrt(max(0.0, 1.0 - z * z))
        theta = golden * i
        pts.append((rad * math.cos(theta), rad * math.sin(theta), z))
    return pts


def sphere_grid():
    """Fibonacci points plus the six axis directions."""
    pts = list(fibonacci_sphere(24))
    for axis in range(3):
        for sign in (1.0, -1.0):
            v = [0.0, 0.0, 0.0]
            v[axis] = sign
            pts.append(tuple(v))
    return pts


# ================================================================ 4. Fidelity
# F_avg(E vs identity) = integral over pure states of <psi|E(|psi>)|psi>.
# With r' = M r + t the integrand is (1 + r.M r + r.t)/2; averaging
# over the sphere kills the t term (<r> = 0) and sends <r_i r_j> to
# delta_ij / 3, giving the closed form F_avg = 1/2 + Tr(M) / 6.
# The sampling referee uses a 26-point set (octahedron + cube corners
# + icosahedron vertices) whose first and second moments match the
# sphere EXACTLY, so sampling and closed form agree to machine
# precision rather than to a Monte-Carlo error bar.

def fidelity_point_set():
    pts = []
    for axis in range(3):  # octahedron
        for sign in (1.0, -1.0):
            v = [0.0, 0.0, 0.0]
            v[axis] = sign
            pts.append(tuple(v))
    for sx in (1.0, -1.0):  # cube corners
        for sy in (1.0, -1.0):
            for sz in (1.0, -1.0):
                pts.append(tuple(v / math.sqrt(3.0) for v in (sx, sy, sz)))
    phi = (1.0 + math.sqrt(5.0)) / 2.0  # icosahedron
    raw = []
    for s1 in (1.0, -1.0):
        for s2 in (1.0, -1.0):
            raw.append((0.0, s1, s2 * phi))
            raw.append((s1, s2 * phi, 0.0))
            raw.append((s2 * phi, 0.0, s1))
    for v in raw:
        nrm = math.sqrt(sum(x * x for x in v))
        pts.append(tuple(x / nrm for x in v))
    return pts


def avg_fidelity_closed(M):
    return 0.5 + (M[0][0] + M[1][1] + M[2][2]) / 6.0


def avg_fidelity_sampled(kraus):
    """Direct average of <psi|E(|psi><psi|)|psi> over the point set."""
    total = 0.0
    pts = fidelity_point_set()
    for r in pts:
        rho = rho_from_bloch(r)
        out = apply_channel(rho, kraus)
        total += mtrace(mmul(rho, out)).real
    return total / len(pts)


def fidelity_report():
    rows = {}
    for name, p in DEMO_PARAMS.items():
        kraus = channel_kraus(name, p)
        M, _ = bloch_map(kraus)
        closed = avg_fidelity_closed(M)
        sampled = avg_fidelity_sampled(kraus)
        rows[name] = {
            "param": p,
            "f_avg_closed": round(closed, 6),
            "f_avg_sampled": round(sampled, 6),
            "abs_dev": abs(closed - sampled),
        }
    return rows


# ================================================================ 5. Choi + PT
# J = (E x I)(|Phi+><Phi+|) = (1/2) sum_ij E(|i><j|) x |i><j|.
# A channel is entanglement-breaking iff its Choi state is separable;
# for a qubit channel the Peres-Horodecki test is exact: J is
# separable iff its partial transpose has no negative eigenvalue.

def choi_state(kraus):
    J = zeros(4)
    for i in range(2):
        for j in range(2):
            e_ij = zeros(2)
            e_ij[i][j] = 1.0
            J = madd(J, mscale(kron(apply_channel(e_ij, kraus), e_ij), 0.5))
    return J


def partial_transpose_second(J):
    """J^Gamma[(a,b),(c,d)] = J[(a,d),(c,b)] (transpose qubit 2)."""
    out = zeros(4)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    out[a * 2 + b][c * 2 + d] = J[a * 2 + d][c * 2 + b]
    return out


def choi_pt_min_eigenvalue(kraus):
    return hermitian_eigenvalues(partial_transpose_second(choi_state(kraus)))[0]


def depolarizing_eb_crossing():
    """Bisect p for the depolarizing Choi PT minimum crossing zero."""
    lo, hi = 0.5, 0.9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if choi_pt_min_eigenvalue(kraus_depolarizing(mid)) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def choi_report():
    rows = {}
    for name, p in DEMO_PARAMS.items():
        kraus = channel_kraus(name, p)
        J = choi_state(kraus)
        eigs = hermitian_eigenvalues(J)
        pt_min = choi_pt_min_eigenvalue(kraus)
        rows[name] = {
            "param": p,
            "choi_eigs": [round(v, 6) for v in eigs],
            "choi_trace": round(mtrace(J).real, 12),
            "pt_min_eig": round(pt_min, 6),
            "entanglement_breaking": bool(pt_min >= -1e-9),
        }
    rows["depolarizing_eb_crossing_p"] = round(depolarizing_eb_crossing(), 9)
    return rows


# ================================================================ 6. Decay bookkeeping
# Per application of amplitude damping gamma, populations in |1>
# shrink by (1 - gamma) and coherences by sqrt(1 - gamma): the log
# decay rate of coherence is exactly half the population rate, the
# T2 = 2 T1 statement for pure amplitude damping. Composing with
# phase damping lambda multiplies coherences by sqrt(1 - lambda)
# more, so log rates ADD: rate_T2 = rate_T1 / 2 + rate_dephasing,
# and T2 <= 2 T1 holds across the whole sweep.

def decay_rates(gamma, lam):
    """Log decay rates per step, read off the composed Bloch map."""
    kraus = compose_kraus(kraus_phase_damping(lam),
                          kraus_amplitude_damping(gamma))
    M, _ = bloch_map(kraus)
    pop_factor = M[2][2]
    coh_factor = M[0][0]
    return {
        "gamma": gamma,
        "lambda": lam,
        "population_factor": pop_factor,
        "coherence_factor": coh_factor,
        "population_rate": -math.log(pop_factor) if pop_factor > 0 else None,
        "coherence_rate": -math.log(coh_factor) if coh_factor > 0 else None,
    }


def iterated_decay(gamma, lam, steps):
    """Apply AD then PD `steps` times to |+> and to |1>; report the
    surviving coherence and population, measured from the states."""
    kraus = compose_kraus(kraus_phase_damping(lam),
                          kraus_amplitude_damping(gamma))
    rho_coh = rho_from_bloch((1.0, 0.0, 0.0))
    rho_pop = rho_from_bloch((0.0, 0.0, -1.0))
    for _ in range(steps):
        rho_coh = apply_channel(rho_coh, kraus)
        rho_pop = apply_channel(rho_pop, kraus)
    return {
        "steps": steps,
        "coherence_surviving": abs(rho_coh[0][1]) / 0.5,
        "population_surviving": rho_pop[1][1].real,
    }


# ================================================================ 7. kNN classifier
# Features of a channel: (s_min, s_max) = smallest / largest diagonal
# shrink factor of the Bloch matrix M, |t| = length of the Bloch
# translation, pt_min = smallest partial-transpose eigenvalue of the
# Choi state. The training grid is labelled by the lab's own policy
# (model_label below) - the classifier interpolates the policy, it
# does not replace the computed quantities.

def channel_features(name, param):
    kraus = channel_kraus(name, param)
    M, t = bloch_map(kraus)
    diag = sorted(M[i][i] for i in range(3))
    t_norm = math.sqrt(sum(v * v for v in t))
    pt_min = choi_pt_min_eigenvalue(kraus)
    return (diag[0], diag[2], t_norm, pt_min)


def model_label(s_min, s_max, t_norm, pt_min):
    """The channel-regime policy distilled to four verdicts."""
    if pt_min >= -1e-6:
        return "entanglement-breaking"
    if s_min >= 0.95 and t_norm <= 0.05:
        return "near-identity"
    if t_norm > 0.05:
        return "damping-dominated"
    if s_max - s_min <= 0.05:
        return "depolarizing"
    return "damping-dominated"


def _training_grid():
    rows = []
    for s_min in (0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.97, 1.0):
        for spread in (0.0, 0.03, 0.2, 0.5):
            s_max = min(1.0, s_min + spread)
            for t_norm in (0.0, 0.02, 0.2, 0.6):
                for pt_min in (-0.45, -0.2, -0.02, 0.0, 0.15):
                    feats = (s_min, s_max, t_norm, pt_min)
                    rows.append((feats, model_label(*feats)))
    return rows


TRAINING_GRID = _training_grid()


def classify_features(features, k=5):
    """Distance-weighted kNN over the policy-labelled grid."""
    ranked = []
    for feats, label in TRAINING_GRID:
        dist = math.sqrt(sum((a - b) ** 2
                             for a, b in zip(features, feats)))
        ranked.append((dist, feats, label))
    ranked.sort(key=lambda row: (row[0], row[1]))
    votes = {}
    for dist, _, label in ranked[:k]:
        votes[label] = votes.get(label, 0.0) + 1.0 / (dist + 1e-9)
    priority = {"near-identity": 0, "depolarizing": 1,
                "damping-dominated": 2, "entanglement-breaking": 3}
    return sorted(votes.items(), key=lambda kv: (-kv[1], priority[kv[0]]))[0][0]


CLASSIFIER_SCENARIOS = {
    "identity_channel": ("bit_flip", 0.0),
    "weak_bit_flip_0.02": ("bit_flip", 0.02),
    "depolarizing_0.3": ("depolarizing", 0.3),
    "depolarizing_0.6": ("depolarizing", 0.6),
    "amplitude_damping_0.4": ("amplitude_damping", 0.4),
    "phase_damping_0.75": ("phase_damping", 0.75),
    "depolarizing_0.8": ("depolarizing", 0.8),
    "amplitude_damping_1.0": ("amplitude_damping", 1.0),
    "phase_flip_0.5": ("phase_flip", 0.5),
}


def classifier_report():
    rows = {}
    for name, (chan, p) in CLASSIFIER_SCENARIOS.items():
        feats = channel_features(chan, p)
        rows[name] = {
            "channel": chan,
            "param": p,
            "features": [round(v, 6) for v in feats],
            "policy": model_label(*feats),
            "knn": classify_features(feats),
        }
    return rows


# ================================================================ 8. Advisory seats
def _trace(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def argon4_decide(topic="quantumchannel-demo"):
    """Gemini Argon 4 DECIDER seat: MAP->MODEL->STRESS->DECIDE.

    Names-only keys, live calls 0, tighten-only: can add HOLDs, never lift one.
    """
    stages = ["MAP", "MODEL", "STRESS", "DECIDE"]
    payload = {"seat": "gemini-argon-4", "topic": topic, "stages": stages}
    return {
        "seat": "Gemini Argon 4 Extended Thinking",
        "model_id": "gemini-argon-4",
        "role": "DECIDER",
        "pipeline": stages,
        "tiers": [2048, 8192, 32768],
        "key": "names-only (GEMINI_API_KEY name checked, value never read)",
        "live_calls": 0,
        "decision": "PROCEED - sandbox simulation only, all deployment gates stay HOLD",
        "can_override_hold": False,
        "trace_sha256_16": _trace(payload),
    }


def claude_advisor(topic="quantumchannel-demo"):
    """Muse advisor seat: advises, can never override a HOLD."""
    return {
        "seat": "Muse",
        "role": "ADVISOR",
        "key": "names-only (MODEL_API_KEY name checked, value never read)",
        "live_calls": 0,
        "advice": "Print the convention before the number: depolarizing means (1-p) rho + p I/2 here, the Choi state is normalised, and the entanglement-breaking crossing at p = 2/3 belongs to those conventions - then show the Bloch map, the fidelity, and the PT eigenvalue that were actually computed.",
        "can_override_hold": False,
    }


def google_flow_plan():
    """Google Flow seat: plan-only, nothing rendered, nothing live."""
    return {
        "seat": "Google Flow",
        "role": "PLAN-ONLY",
        "shots": ["a Bloch sphere shrinking uniformly under depolarizing noise",
                  "the sphere sliding toward the north pole as amplitude damping bites",
                  "an ellipsoid flattening as phase damping eats the equator",
                  "the Choi state's partial transpose eigenvalue crossing zero at p = 2/3",
                  "two channels composing while their log decay rates add on a ledger"],
        "rendered": False,
        "live": 0,
    }


def brains_report():
    return {"decider": argon4_decide(), "advisor": claude_advisor(),
            "flow": google_flow_plan()}


# ================================================================ 9. Topic reports
def kraus_report():
    rows = {}
    for name, p in DEMO_PARAMS.items():
        kraus = channel_kraus(name, p)
        rows[name] = {
            "param": p,
            "n_kraus": len(kraus),
            "completeness_dev": completeness_deviation(kraus),
        }
    return rows


def bloch_report():
    rows = {}
    for name, p in DEMO_PARAMS.items():
        kraus = channel_kraus(name, p)
        M, t = bloch_map(kraus)
        rows[name] = {
            "param": p,
            "M_diag": [round(M[i][i], 6) for i in range(3)],
            "M_offdiag_max": round(max(abs(M[i][j]) for i in range(3)
                                       for j in range(3) if i != j), 12),
            "t": [round(v, 6) for v in t],
        }
    return rows


def compose_report():
    rows = {}
    for gamma in (0.1, 0.3, 0.5):
        for lam in (0.0, 0.2, 0.5):
            r = decay_rates(gamma, lam)
            rows[f"gamma{gamma}_lambda{lam}"] = {
                "population_rate": round(r["population_rate"], 6),
                "coherence_rate": round(r["coherence_rate"], 6),
                "t1_over_2_plus_dephasing": round(
                    0.5 * r["population_rate"]
                    + (-0.5 * math.log(1.0 - lam) if lam < 1 else 0.0), 6),
                "t2_le_2t1": bool(r["coherence_rate"]
                                  >= 0.5 * r["population_rate"] - 1e-12),
            }
    rows["iterated_gamma0.3_lambda0.2_40steps"] = iterated_decay(0.3, 0.2, 40)
    return rows


# ================================================================ checks
def run_checks():
    c = []

    def chk(name, ok):
        c.append((name, bool(ok)))

    # 0. Jacobi eigensolver refereed on known spectra first
    chk("jacobi_pauli_y_spectrum",
        all(abs(a - b) < 1e-12 for a, b in
            zip(hermitian_eigenvalues(PAULI_Y), (-1.0, 1.0))))
    chk("jacobi_complex_hermitian_2x2",
        all(abs(a - b) < 1e-12 for a, b in
            zip(hermitian_eigenvalues(mat([[2, 1j], [-1j, 3]])),
                ((5 - math.sqrt(5)) / 2, (5 + math.sqrt(5)) / 2))))
    xx_zz = madd(kron(PAULI_X, PAULI_X), kron(PAULI_Z, PAULI_Z))
    chk("jacobi_xx_plus_zz_spectrum",
        all(abs(a - b) < 1e-12 for a, b in
            zip(hermitian_eigenvalues(xx_zz), (-2.0, 0.0, 0.0, 2.0))))

    # 1. Kraus completeness on the full parameter grid
    chk("kraus_completeness_all_channels_grid",
        all(completeness_deviation(channel_kraus(name, p)) < 1e-12
            for name in CHANNEL_BUILDERS for p in PARAM_GRID))
    chk("kraus_counts",
        len(kraus_bit_flip(0.3)) == 2 and len(kraus_depolarizing(0.3)) == 4
        and len(kraus_amplitude_damping(0.3)) == 2)
    chk("zero_param_is_identity_channel",
        all(mat_max_dev(apply_channel(rho_from_bloch(r),
                                      channel_kraus(name, 0.0)),
                        rho_from_bloch(r)) < 1e-12
            for name in CHANNEL_BUILDERS for r in sphere_grid()))
    chk("depolarizing_p1_is_completely_mixed",
        max(abs(v) for row in bloch_map(kraus_depolarizing(1.0))[0]
            for v in row) < 1e-12)
    chk("amplitude_damping_g1_sends_all_to_ground",
        all(mat_max_dev(apply_channel(rho_from_bloch(r),
                                      kraus_amplitude_damping(1.0)),
                        mat([[1, 0], [0, 0]])) < 1e-12
            for r in sphere_grid()))
    chk("phase_damping_l1_kills_coherence_keeps_population",
        all(abs(apply_channel(rho_from_bloch(r),
                              kraus_phase_damping(1.0))[0][1]) < 1e-12
            and abs(apply_channel(rho_from_bloch(r),
                                  kraus_phase_damping(1.0))[1][1]
                    - rho_from_bloch(r)[1][1]) < 1e-12
            for r in sphere_grid()))

    # 2. Bloch affine map: closed forms + reproduction of direct evolution
    M, t = bloch_map(kraus_depolarizing(0.3))
    chk("bloch_depolarizing_uniform_shrink",
        all(abs(M[i][j] - (0.7 if i == j else 0.0)) < 1e-12
            for i in range(3) for j in range(3))
        and all(abs(v) < 1e-12 for v in t))
    M, t = bloch_map(kraus_amplitude_damping(0.5))
    chk("bloch_amplitude_damping_form",
        abs(M[0][0] - math.sqrt(0.5)) < 1e-12
        and abs(M[1][1] - math.sqrt(0.5)) < 1e-12
        and abs(M[2][2] - 0.5) < 1e-12
        and abs(t[2] - 0.5) < 1e-12 and abs(t[0]) < 1e-12)
    M, _ = bloch_map(kraus_bit_flip(0.2))
    chk("bloch_bit_flip_form",
        abs(M[0][0] - 1.0) < 1e-12 and abs(M[1][1] - 0.6) < 1e-12
        and abs(M[2][2] - 0.6) < 1e-12)
    M, _ = bloch_map(kraus_phase_flip(0.2))
    chk("bloch_phase_flip_form",
        abs(M[0][0] - 0.6) < 1e-12 and abs(M[1][1] - 0.6) < 1e-12
        and abs(M[2][2] - 1.0) < 1e-12)
    chk("bloch_affine_reproduces_direct_all_channels",
        all(max(abs(a - b) for a, b in
                zip(bloch_vector(apply_channel(
                    rho_from_bloch(r), channel_kraus(name, p))),
                    apply_bloch(*bloch_map(channel_kraus(name, p)), r)))
            < 1e-12
            for name, p in DEMO_PARAMS.items() for r in sphere_grid()))
    chk("bloch_unital_channels_contract_ball",
        all(math.sqrt(sum(v * v for v in apply_bloch(
            *bloch_map(channel_kraus(name, p)), r))) <= 1.0 + 1e-12
            for name in ("bit_flip", "phase_flip", "bit_phase_flip",
                         "depolarizing", "phase_damping")
            for p in (0.25, 0.5, 0.9) for r in sphere_grid()))

    # 3. T1 / T2 bookkeeping for amplitude damping
    rates = decay_rates(0.5, 0.0)
    chk("ad_coherence_rate_is_half_population_rate",
        abs(rates["coherence_rate"] - 0.5 * rates["population_rate"])
        < 1e-12)
    chk("ad_t2_equals_2t1_ratio",
        abs(rates["population_rate"] / rates["coherence_rate"] - 2.0)
        < 1e-12)

    # 4. Average gate fidelity: sampling referee vs closed form
    chk("fidelity_sampling_matches_closed_all_grid",
        all(abs(avg_fidelity_sampled(channel_kraus(name, p))
                - avg_fidelity_closed(bloch_map(channel_kraus(name, p))[0]))
            < 1e-9
            for name in CHANNEL_BUILDERS for p in PARAM_GRID))
    chk("fidelity_depolarizing_closed_form",
        all(abs(avg_fidelity_closed(bloch_map(kraus_depolarizing(p))[0])
                - (1.0 - p / 2.0)) < 1e-12 for p in PARAM_GRID))
    chk("fidelity_bit_flip_closed_form",
        all(abs(avg_fidelity_closed(bloch_map(kraus_bit_flip(p))[0])
                - (1.0 - 2.0 * p / 3.0)) < 1e-12 for p in PARAM_GRID))
    chk("fidelity_amplitude_damping_closed_form",
        all(abs(avg_fidelity_closed(
            bloch_map(kraus_amplitude_damping(g))[0])
            - (0.5 + (2.0 * math.sqrt(1.0 - g) + (1.0 - g)) / 6.0)) < 1e-12
            for g in PARAM_GRID))
    chk("fidelity_phase_damping_closed_form",
        all(abs(avg_fidelity_closed(
            bloch_map(kraus_phase_damping(l))[0])
            - (0.5 + (2.0 * math.sqrt(1.0 - l) + 1.0) / 6.0)) < 1e-12
            for l in PARAM_GRID))
    chk("fidelity_identity_is_one",
        abs(avg_fidelity_closed(bloch_map(kraus_bit_flip(0.0))[0]) - 1.0)
        < 1e-12)

    # 5. Choi state: CP, trace, analytic PT spectra, EB crossing
    chk("choi_trace_one_and_hermitian",
        all(abs(mtrace(choi_state(channel_kraus(name, p))) - 1.0) < 1e-12
            and mat_max_dev(choi_state(channel_kraus(name, p)),
                            dagger(choi_state(channel_kraus(name, p))))
            < 1e-12
            for name, p in DEMO_PARAMS.items()))
    chk("choi_eigenvalues_nonnegative_cp",
        all(hermitian_eigenvalues(choi_state(channel_kraus(name, p)))[0]
            > -1e-12
            for name in CHANNEL_BUILDERS for p in PARAM_GRID))
    chk("choi_depolarizing_spectrum",
        all(abs(a - b) < 1e-12 for a, b in
            zip(hermitian_eigenvalues(choi_state(kraus_depolarizing(0.4))),
                (0.1, 0.1, 0.1, 0.7))))
    chk("choi_pt_depolarizing_analytic",
        all(abs(choi_pt_min_eigenvalue(kraus_depolarizing(p))
                - (3.0 * p - 2.0) / 4.0) < 1e-12 for p in PARAM_GRID))
    chk("choi_pt_bit_flip_analytic",
        all(abs(choi_pt_min_eigenvalue(kraus_bit_flip(p))
                + abs(1.0 - 2.0 * p) / 2.0) < 1e-12 for p in PARAM_GRID))
    chk("choi_pt_phase_damping_analytic",
        all(abs(choi_pt_min_eigenvalue(kraus_phase_damping(l))
                + math.sqrt(1.0 - l) / 2.0) < 1e-12 for l in PARAM_GRID))
    chk("choi_pt_amplitude_damping_negative_below_g1",
        all(choi_pt_min_eigenvalue(kraus_amplitude_damping(g)) < -1e-3
            for g in (0.0, 0.25, 0.5, 0.75, 0.9)))
    chk("choi_pt_amplitude_damping_zero_at_g1",
        abs(choi_pt_min_eigenvalue(kraus_amplitude_damping(1.0))) < 1e-9)
    chk("depolarizing_eb_crossing_is_two_thirds",
        abs(depolarizing_eb_crossing() - 2.0 / 3.0) < 1e-9)

    # 6. Composition: Bloch maps multiply, log rates add, T2 <= 2 T1
    chk("compose_kraus_completeness",
        all(completeness_deviation(compose_kraus(
            kraus_phase_damping(l), kraus_amplitude_damping(g))) < 1e-12
            for g in (0.1, 0.5, 0.9) for l in (0.1, 0.5, 0.9)))
    chk("compose_bloch_maps_multiply",
        all(max(abs(a - b) for a, b in
                zip(bloch_vector(apply_channel(
                    rho_from_bloch(r),
                    compose_kraus(kraus_phase_damping(0.2),
                                  kraus_amplitude_damping(0.3)))),
                    apply_bloch(*bloch_map(kraus_phase_damping(0.2)),
                                apply_bloch(
                                    *bloch_map(kraus_amplitude_damping(0.3)),
                                    r)))) < 1e-12
            for r in sphere_grid()))
    chk("compose_log_rates_add",
        all(abs(decay_rates(g, l)["coherence_rate"]
                - (0.5 * decay_rates(g, l)["population_rate"]
                   - 0.5 * math.log(1.0 - l))) < 1e-12
            for g in (0.1, 0.3, 0.5, 0.9) for l in (0.0, 0.2, 0.5, 0.9)))
    chk("t2_le_2t1_across_sweep",
        all(decay_rates(g, l)["coherence_rate"]
            >= 0.5 * decay_rates(g, l)["population_rate"] - 1e-12
            for g in (0.1, 0.3, 0.5, 0.9) for l in (0.0, 0.2, 0.5, 0.9)))
    chk("t2_equals_2t1_when_no_pure_dephasing",
        all(abs(decay_rates(g, 0.0)["coherence_rate"]
                - 0.5 * decay_rates(g, 0.0)["population_rate"]) < 1e-12
            for g in (0.1, 0.3, 0.5, 0.9)))
    it = iterated_decay(0.3, 0.0, 25)
    chk("iterated_ad_population_power_law",
        abs(it["population_surviving"] - 0.7 ** 25) < 1e-12)
    chk("iterated_ad_coherence_power_law",
        abs(it["coherence_surviving"] - math.sqrt(0.7) ** 25) < 1e-12)
    it2 = iterated_decay(0.3, 0.2, 40)
    chk("iterated_composed_coherence_power_law",
        abs(it2["coherence_surviving"]
            - (math.sqrt(0.7) * math.sqrt(0.8)) ** 40) < 1e-9)

    # 7. Classifier: labels come from the lab's own policy
    rep = classifier_report()
    chk("classifier_scenario_labels",
        rep["identity_channel"]["knn"] == "near-identity"
        and rep["weak_bit_flip_0.02"]["knn"] == "near-identity"
        and rep["depolarizing_0.3"]["knn"] == "depolarizing"
        and rep["depolarizing_0.6"]["knn"] == "depolarizing"
        and rep["amplitude_damping_0.4"]["knn"] == "damping-dominated"
        and rep["phase_damping_0.75"]["knn"] == "damping-dominated"
        and rep["depolarizing_0.8"]["knn"] == "entanglement-breaking"
        and rep["amplitude_damping_1.0"]["knn"] == "entanglement-breaking"
        and rep["phase_flip_0.5"]["knn"] == "entanglement-breaking")
    chk("classifier_knn_matches_policy_scenarios",
        all(row["knn"] == row["policy"] for row in rep.values()))
    chk("classifier_knn_matches_policy_grid",
        all(classify_features(feats) == label
            for feats, label in TRAINING_GRID))

    # gates + seats
    chk("hold_gates_14", len(HOLD_GATES) == 14
        and all(v.startswith("HOLD") for v in HOLD_GATES.values()))
    chk("gate_network_deployment_hold",
        gate("quantum_network_deployment_claim") == "HOLD")
    chk("gate_error_correction_product_hold",
        gate("error_correction_product_claim") == "HOLD")
    chk("gate_capacity_claim_hold", gate("capacity_or_rate_claim") == "HOLD")
    chk("gate_vendor_hold", gate("vendor_channel_claim") == "HOLD")
    chk("gate_sandbox_proceed", gate("run_simulation") == "PROCEED")
    brains = brains_report()
    chk("argon4_no_override", brains["decider"]["can_override_hold"] is False
        and brains["decider"]["seat"]
        == "Gemini Argon 4 Extended Thinking"
        and brains["decider"]["model_id"] == "gemini-argon-4")
    chk("claude_no_override", brains["advisor"]["can_override_hold"] is False)
    chk("flow_plan_only", brains["flow"]["rendered"] is False
        and brains["flow"]["live"] == 0)
    return c


# ================================================================ CLI / demo
def demo():
    print("=" * 68)
    print("QuantumChannel-Lab - single-qubit quantum channels, four ways")
    print("by Darion Kellar | MIT | pure stdlib | zero network | EDUCATION ONLY")
    print("=" * 68)
    print("\n[1] THE HONEST CORE - one qubit, Markovian noise, nothing real")
    print("    A channel is bookkeeping: Kraus operators in, and out come")
    print("    the Bloch affine map, the average fidelity, and the point")
    print("    where entanglement dies. All four views are computed from")
    print("    the same Kraus operators and graded against each other.")
    print("    No real device was measured; no capacity is claimed.")
    print("\n[2] Kraus library - completeness sum K^dagger K = I")
    for name, row in kraus_report().items():
        print(f"    {name:20s} param {row['param']:.2f} | "
              f"{row['n_kraus']} Kraus ops | completeness dev "
              f"{row['completeness_dev']:.2e}")
    print("\n[3] Bloch picture - r -> M r + t, read off the Kraus operators")
    for name, row in bloch_report().items():
        print(f"    {name:20s} M diag {row['M_diag']} | "
              f"t {[round(v, 4) for v in row['t']]}")
    print("    -> depolarizing shrinks every axis by (1 - p); amplitude")
    print("       damping shrinks x,y by sqrt(1-gamma), z by (1-gamma),")
    print("       and drags the ball toward |0> with t = (0, 0, gamma).")
    rates = decay_rates(0.5, 0.0)
    print("\n[4] T1 / T2 bookkeeping - pure amplitude damping (gamma = 0.5)")
    print(f"    population log-rate {rates['population_rate']:.6f} | "
          f"coherence log-rate {rates['coherence_rate']:.6f}")
    print(f"    ratio exactly "
          f"{rates['population_rate'] / rates['coherence_rate']:.6f} "
          f"-> T2 = 2 T1 when there is no pure dephasing.")
    print("\n[5] Average gate fidelity vs identity - closed form refereed")
    print("    by a 26-point exact-moment Bloch sampling:")
    for name, row in fidelity_report().items():
        print(f"    {name:20s} F_avg {row['f_avg_closed']:.6f} | "
              f"sampled {row['f_avg_sampled']:.6f} | "
              f"dev {row['abs_dev']:.2e}")
    print("    -> depolarizing: F_avg = 1 - p/2 under the stated")
    print("       convention E(rho) = (1 - p) rho + p I/2.")
    print("\n[6] Choi state + partial transpose - where entanglement dies")
    choi = choi_report()
    for name in DEMO_PARAMS:
        row = choi[name]
        print(f"    {name:20s} PT min eig {row['pt_min_eig']:+.6f} | "
              f"entanglement-breaking {row['entanglement_breaking']}")
    print(f"    depolarizing EB crossing: p* = "
          f"{choi['depolarizing_eb_crossing_p']:.6f} (exact 2/3 = "
          f"{2.0 / 3.0:.6f}); amplitude damping stays non-EB for")
    print("    every gamma < 1 and reaches the boundary only at gamma = 1.")
    print("\n[7] Composition - phase damping after amplitude damping")
    comp = compose_report()
    for key in ("gamma0.3_lambda0.2", "gamma0.3_lambda0.5"):
        row = comp[key]
        print(f"    {key:22s} coherence rate {row['coherence_rate']:.6f} "
              f"= T1/2 + dephasing {row['t1_over_2_plus_dephasing']:.6f} | "
              f"T2 <= 2T1 {row['t2_le_2t1']}")
    it = comp["iterated_gamma0.3_lambda0.2_40steps"]
    print(f"    after 40 composed steps: coherence surviving "
          f"{it['coherence_surviving']:.6f}, population surviving "
          f"{it['population_surviving']:.6f} (measured off the states)")
    print("\n[8] Channel classifier (TinyML kNN, labels from the lab policy)")
    for name, row in classifier_report().items():
        print(f"    {name:28s} -> {row['knn']}")
    print("\n[9] Advisory seats")
    brains = brains_report()
    print(f"    decider: {brains['decider']['decision']} "
          f"trace={brains['decider']['trace_sha256_16']}")
    print(f"    advisor can_override_hold = {brains['advisor']['can_override_hold']}")
    print(f"    flow rendered = {brains['flow']['rendered']} "
          f"live = {brains['flow']['live']}")
    print(f"\n[10] HOLD gates: {len(HOLD_GATES)}")
    for key in HOLD_GATES:
        print(f"    HOLD  {key}")
    checks = run_checks()
    passed = sum(1 for _, ok in checks if ok)
    print(f"\nCHECKS {passed}/{len(checks)}")
    for name, ok in checks:
        if not ok:
            print(f"    FAIL {name}")
    verdict = "PROCEED" if passed == len(checks) else "HOLD"
    print(f"VERDICT {verdict} (sandbox simulation only) | HOLD gates {len(HOLD_GATES)}")
    return 0 if passed == len(checks) else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] == "demo":
        return demo()
    cmd = argv[0]
    if cmd == "checks":
        checks = run_checks()
        for name, ok in checks:
            print(f"{'PASS' if ok else 'FAIL'} {name}")
        print(f"CHECKS {sum(1 for _, ok in checks if ok)}/{len(checks)}")
        return 0 if all(ok for _, ok in checks) else 1
    if cmd == "kraus":
        print(json.dumps(kraus_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "bloch":
        print(json.dumps(bloch_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "fidelity":
        print(json.dumps(fidelity_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "choi":
        print(json.dumps(choi_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "compose":
        print(json.dumps(compose_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "classify":
        print(json.dumps(classifier_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "brains":
        print(json.dumps(brains_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "holds":
        for key, value in HOLD_GATES.items():
            print(f"{key}: {value}")
        return 0
    if cmd == "hold-demo":
        action = argv[1] if len(argv) > 1 else "quantum_network_deployment_claim"
        print(f"{action}: {HOLD_GATES.get(action, 'PROCEED - not gated')}")
        return 2 if gate(action) == "HOLD" else 0
    print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
