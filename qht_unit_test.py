import numpy as np
import pandas as pd
from qiskit import QuantumRegister, AncillaRegister, QuantumCircuit
from qiskit.quantum_info import Statevector

def controlled_subtract_from_half(circuit, control, target, anc, index=0):
    """
    Perform: |1⟩|y⟩ → |1⟩|~y + 1⟩ = |1⟩|2^n - y⟩ (identity if control=0)
    """
    n = len(target)
    # Step 1: controlled bitwise NOT
    for i in range(n):
        circuit.cx(control, target[i])
    # Step 2: controlled increment (ripple-carry)
    circuit.ccx(control, target[0], anc[0])
    circuit.cx(control, target[0])
    for i in range(1, n):
        circuit.ccx(anc[i-1], target[i], anc[i])
        circuit.cx(anc[i-1], target[i])
    # Step 3: uncompute to restore ancillas and target
    for i in range(n):
        circuit.x(target[i])
    for i in reversed(range(1, n)):
        circuit.ccx(anc[i-1], target[i], anc[i])
    circuit.ccx(control, target[0], anc[0])
    for i in range(n):
        circuit.x(target[i])

def _build_statevector(n, control_val, y_val):
    """Build & run the circuit under test; return its Statevector."""
    ctrl = QuantumRegister(1, 'ctrl')
    tgt  = QuantumRegister(n, 'tgt')
    anc  = AncillaRegister(n, 'anc')
    qc = QuantumCircuit(ctrl, tgt, anc)
    if control_val:
        qc.x(ctrl[0])
    for i in range(n):
        if (y_val >> i) & 1:
            qc.x(tgt[i])
    controlled_subtract_from_half(qc, ctrl[0], tgt, anc)
    return Statevector.from_instruction(qc)

def _build_expected(n, control_val, y_val):
    """Build & run the circuit that prepares the *expected* state."""
    ctrl = QuantumRegister(1, 'ctrl')
    tgt  = QuantumRegister(n, 'tgt')
    anc  = AncillaRegister(n, 'anc')
    qc = QuantumCircuit(ctrl, tgt, anc)
    if control_val:
        qc.x(ctrl[0])
    val = y_val if control_val == 0 else (2**n - y_val) % 2**n
    for i in range(n):
        if (val >> i) & 1:
            qc.x(tgt[i])
    return Statevector.from_instruction(qc)

def test_controlled_subtract_from_half(n):
    """Test all (control, y) pairs and print each result in binary."""
    total = 0
    for control_val in (0, 1):
        for y in range(2**n):
            total += 1
            # get actual and expected statevectors
            sv_act = _build_statevector(n, control_val, y)
            sv_exp = _build_expected(n, control_val, y)
            # compute expected y'
            y_prime = y if control_val == 0 else (2**n - y) % 2**n

            # format binary strings
            y_bin      = format(y,      '0{}b'.format(n))
            y_prime_bin= format(y_prime,'0{}b'.format(n))

            if np.allclose(sv_act.data, sv_exp.data, atol=1e-8):
                print(f"PASS: control={control_val} | input y={y_bin} | output y'={y_prime_bin}")
            else:
                print(f"FAIL: control={control_val} | input y={y_bin} | expected y'={y_prime_bin}")
                print(f"  Actual statevector:   {sv_act.data}")
                print(f"  Expected statevector: {sv_exp.data}")
                raise AssertionError("Mismatch detected, aborting tests.")
    print(f"\nAll {total} test cases passed!")


if __name__ == "__main__":
    test_controlled_subtract_from_half(3)

#-------------------------------------test rotation-------------------------------------------------
def apply_UR(circuit, ctrl, target_y, target_b, N):
    """
    Apply U_R rotation: Ry(-4π·2^j / N), controlled by y[j] and b, on target.
    """
    for j, y_qubit in enumerate(target_y):
        theta = -4 * np.pi * (2 ** j) / N
        circuit.mcry(theta, [y_qubit, target_b[0]], ctrl, [], mode=None)

results = []

if __name__ == "__main__":
    for y_val in range(4):      # y ∈ {00, 01, 10, 11}
        for b_val in [0,1]:    # b ∈ {0, 1}
            qr = QuantumRegister(2, 'q')  # y[0], y[1], b, ctrl
            qrb = QuantumRegister(1, 'b')
            qrc = QuantumRegister(1, 'c')
            qc = QuantumCircuit(qr,qrb,qrc)


            #qc.reset(qr)

            if y_val & 1: qc.x(qr[0])      # y[0]
            #print(y_val & 1)
            if y_val >> 1: qc.x(qr[1])     # y[1]
            #print(y_val >> 1)
            if b_val: qc.x(qrb[0])          # b

            N = 8
            apply_UR(qc, qrc, qr,qrb, N)
            #qc.mcry(-4 * np.pi / N, [0,2], 3, [], mode=None)
            #qc.mcry(-8 * np.pi / N, [1,2], 3, [], mode=None)



            state = Statevector.from_instruction(qc)
            state = state.data.reshape((2, 2, 2, 2))  # 4 qubits
            amp = state[1, b_val, (y_val >> 1) & 1, y_val & 1]
            theory_theta = -2 * np.pi * y_val * b_val / N
            results.append({
                "y": f"{y_val:02b}",
                "b": b_val,
                "theory sin(-θ)": f"{np.sin(theory_theta):.4f}",
                "amp": f"{amp.real:.4f}"
            })

    df = pd.DataFrame(results)
    print(df.to_string(index=False))
