"""Tests for QuantumChannel-Lab. Every test grades a computed number
against a closed form or an independent implementation in main.py."""
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import main as qc

TOL = 1e-9


def test_kraus_completeness_grid():
    for name in qc.CHANNEL_BUILDERS:
        for p in qc.PARAM_GRID:
            assert qc.completeness_deviation(qc.channel_kraus(name, p)) < 1e-12


def test_zero_param_identity():
    for name in qc.CHANNEL_BUILDERS:
        kraus = qc.channel_kraus(name, 0.0)
        for r in qc.sphere_grid():
            out = qc.apply_channel(qc.rho_from_bloch(r), kraus)
            assert qc.mat_max_dev(out, qc.rho_from_bloch(r)) < 1e-12


def test_jacobi_known_spectra():
    assert qc.hermitian_eigenvalues(qc.PAULI_Y) == [-1.0, 1.0] or all(
        abs(a - b) < 1e-12 for a, b in
        zip(qc.hermitian_eigenvalues(qc.PAULI_Y), (-1.0, 1.0)))
    xx_zz = qc.madd(qc.kron(qc.PAULI_X, qc.PAULI_X),
                    qc.kron(qc.PAULI_Z, qc.PAULI_Z))
    assert all(abs(a - b) < 1e-12 for a, b in
               zip(qc.hermitian_eigenvalues(xx_zz), (-2.0, 0.0, 0.0, 2.0)))


def test_bloch_depolarizing():
    M, t = qc.bloch_map(qc.kraus_depolarizing(0.3))
    for i in range(3):
        assert abs(M[i][i] - 0.7) < 1e-12
        assert abs(t[i]) < 1e-12


def test_bloch_amplitude_damping():
    M, t = qc.bloch_map(qc.kraus_amplitude_damping(0.5))
    assert abs(M[0][0] - math.sqrt(0.5)) < 1e-12
    assert abs(M[2][2] - 0.5) < 1e-12
    assert abs(t[2] - 0.5) < 1e-12


def test_bloch_flip_forms():
    M, _ = qc.bloch_map(qc.kraus_bit_flip(0.2))
    assert abs(M[0][0] - 1.0) < 1e-12 and abs(M[1][1] - 0.6) < 1e-12
    M, _ = qc.bloch_map(qc.kraus_phase_flip(0.2))
    assert abs(M[2][2] - 1.0) < 1e-12 and abs(M[0][0] - 0.6) < 1e-12
    M, _ = qc.bloch_map(qc.kraus_bit_phase_flip(0.2))
    assert abs(M[1][1] - 1.0) < 1e-12 and abs(M[0][0] - 0.6) < 1e-12


def test_affine_reproduces_direct():
    for name, p in qc.DEMO_PARAMS.items():
        kraus = qc.channel_kraus(name, p)
        M, t = qc.bloch_map(kraus)
        for r in qc.sphere_grid():
            direct = qc.bloch_vector(qc.apply_channel(qc.rho_from_bloch(r),
                                                      kraus))
            affine = qc.apply_bloch(M, t, r)
            assert max(abs(a - b) for a, b in zip(direct, affine)) < 1e-12


def test_t2_equals_2t1_pure_ad():
    rates = qc.decay_rates(0.5, 0.0)
    assert abs(rates["coherence_rate"]
               - 0.5 * rates["population_rate"]) < 1e-12


def test_fidelity_sampling_matches_closed():
    for name in qc.CHANNEL_BUILDERS:
        for p in (0.1, 0.5, 0.9):
            kraus = qc.channel_kraus(name, p)
            closed = qc.avg_fidelity_closed(qc.bloch_map(kraus)[0])
            sampled = qc.avg_fidelity_sampled(kraus)
            assert abs(closed - sampled) < 1e-9


def test_fidelity_closed_forms():
    for p in qc.PARAM_GRID:
        M, _ = qc.bloch_map(qc.kraus_depolarizing(p))
        assert abs(qc.avg_fidelity_closed(M) - (1 - p / 2)) < 1e-12
        M, _ = qc.bloch_map(qc.kraus_bit_flip(p))
        assert abs(qc.avg_fidelity_closed(M) - (1 - 2 * p / 3)) < 1e-12
        M, _ = qc.bloch_map(qc.kraus_amplitude_damping(p))
        expect = 0.5 + (2 * math.sqrt(1 - p) + (1 - p)) / 6
        assert abs(qc.avg_fidelity_closed(M) - expect) < 1e-12


def test_choi_trace_and_cp():
    for name, p in qc.DEMO_PARAMS.items():
        J = qc.choi_state(qc.channel_kraus(name, p))
        assert abs(qc.mtrace(J) - 1.0) < 1e-12
        assert qc.hermitian_eigenvalues(J)[0] > -1e-12


def test_choi_pt_analytics():
    for p in qc.PARAM_GRID:
        assert abs(qc.choi_pt_min_eigenvalue(qc.kraus_depolarizing(p))
                   - (3 * p - 2) / 4) < 1e-12
        assert abs(qc.choi_pt_min_eigenvalue(qc.kraus_bit_flip(p))
                   + abs(1 - 2 * p) / 2) < 1e-12
        assert abs(qc.choi_pt_min_eigenvalue(qc.kraus_phase_damping(p))
                   + math.sqrt(1 - p) / 2) < 1e-12


def test_eb_crossing_two_thirds():
    assert abs(qc.depolarizing_eb_crossing() - 2 / 3) < 1e-9


def test_ad_not_eb_below_gamma1():
    for g in (0.0, 0.25, 0.5, 0.75, 0.9):
        assert qc.choi_pt_min_eigenvalue(qc.kraus_amplitude_damping(g)) < -1e-3
    assert abs(qc.choi_pt_min_eigenvalue(
        qc.kraus_amplitude_damping(1.0))) < 1e-9


def test_compose_rates_add():
    for g in (0.1, 0.3, 0.5):
        for l in (0.0, 0.2, 0.5):
            r = qc.decay_rates(g, l)
            expect = 0.5 * r["population_rate"] - 0.5 * math.log(1 - l)
            assert abs(r["coherence_rate"] - expect) < 1e-12
            assert r["coherence_rate"] >= 0.5 * r["population_rate"] - 1e-12


def test_iterated_power_laws():
    it = qc.iterated_decay(0.3, 0.0, 25)
    assert abs(it["population_surviving"] - 0.7 ** 25) < 1e-12
    assert abs(it["coherence_surviving"] - math.sqrt(0.7) ** 25) < 1e-12


def test_classifier_scenarios():
    rep = qc.classifier_report()
    assert rep["identity_channel"]["knn"] == "near-identity"
    assert rep["depolarizing_0.3"]["knn"] == "depolarizing"
    assert rep["amplitude_damping_0.4"]["knn"] == "damping-dominated"
    assert rep["depolarizing_0.8"]["knn"] == "entanglement-breaking"
    assert rep["phase_flip_0.5"]["knn"] == "entanglement-breaking"
    assert all(row["knn"] == row["policy"] for row in rep.values())


def test_classifier_grid_consistency():
    for feats, label in qc.TRAINING_GRID:
        assert qc.classify_features(feats) == label


def test_hold_gates():
    assert len(qc.HOLD_GATES) == 14
    assert qc.gate("quantum_network_deployment_claim") == "HOLD"
    assert qc.gate("error_correction_product_claim") == "HOLD"
    assert qc.gate("run_simulation") == "PROCEED"


def test_seats_never_override():
    brains = qc.brains_report()
    assert brains["decider"]["can_override_hold"] is False
    assert brains["decider"]["model_id"] == "gemini-argon-4"
    assert brains["advisor"]["can_override_hold"] is False
    assert brains["flow"]["rendered"] is False
    assert brains["flow"]["live"] == 0


def test_checks_all_pass():
    checks = qc.run_checks()
    assert all(ok for _, ok in checks), [n for n, ok in checks if not ok]


def test_cli_hold_demo_exit_2():
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent.parent
                             / "main.py"), "hold-demo"],
        capture_output=True, text=True)
    assert proc.returncode == 2
    assert "HOLD" in proc.stdout


def test_cli_checks_command():
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent.parent
                             / "main.py"), "checks"],
        capture_output=True, text=True)
    assert proc.returncode == 0
    assert "CHECKS 52/52" in proc.stdout


def test_demo_deterministic():
    main_py = str(Path(__file__).resolve().parent.parent / "main.py")
    a = subprocess.run([sys.executable, main_py, "demo"],
                       capture_output=True, text=True)
    b = subprocess.run([sys.executable, main_py, "demo"],
                       capture_output=True, text=True)
    assert a.stdout == b.stdout and a.returncode == 0
