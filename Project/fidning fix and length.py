import numpy as np
import matplotlib.pyplot as plt

def deg2rad(deg):
    return deg * np.pi / 180

def synthesize_and_plot_matlab_style(pps, theta_deg, beta_deg, gamma_deg):
    # Step 1: Convert inputs
    theta1, theta2, theta3 = map(deg2rad, theta_deg)
    beta2, beta3 = map(deg2rad, beta_deg)
    gamma2, gamma3 = map(deg2rad, gamma_deg)

    p1 = complex(*pps[0])
    p2 = complex(*pps[1])
    p3 = complex(*pps[2])

    alpha2 = theta2 - theta1
    alpha3 = theta3 - theta1
    delta2 = p2 - p1
    delta3 = p3 - p1

    # Step 2: Solve vectors
    A = np.array([[delta2, np.exp(1j*alpha2) - 1],
                  [delta3, np.exp(1j*alpha3) - 1]])
    B = np.array([[np.exp(1j*beta2) - 1, delta2],
                  [np.exp(1j*beta3) - 1, delta3]])
    C = np.array([[np.exp(1j*beta2) - 1, np.exp(1j*alpha2) - 1],
                  [np.exp(1j*beta3) - 1, np.exp(1j*alpha3) - 1]])
    AA = A
    BB = np.array([[np.exp(1j*gamma2) - 1, delta2],
                   [np.exp(1j*gamma3) - 1, delta3]])
    CC = np.array([[np.exp(1j*gamma2) - 1, np.exp(1j*alpha2) - 1],
                   [np.exp(1j*gamma3) - 1, np.exp(1j*alpha3) - 1]])

    WA = np.linalg.det(A) / np.linalg.det(C)
    ZA = np.linalg.det(B) / np.linalg.det(C)
    WB = np.linalg.det(AA) / np.linalg.det(CC)
    ZB = np.linalg.det(BB) / np.linalg.det(CC)

    OA = p1 - WA - ZA
    OB = p1 - WB - ZB

    P2 = p1 - ZA
    P3 = p1 - ZB

    # Step 3: Plot just like MATLAB
    plt.figure(figsize=(8, 6))
    plt.title("3PP Synthesized Mechanism (MATLAB Style)")
    plt.grid(True)
    plt.axis('equal')

    # 🔴 WA and ZA vectors
    plt.quiver(OA.real, OA.imag, WA.real, WA.imag, angles='xy', scale_units='xy', scale=1, color='r', label='W_A')
    A = OA + WA
    plt.quiver(A.real, A.imag, ZA.real, ZA.imag, angles='xy', scale_units='xy', scale=1, color='b', label='Z_A')

    # 🟣 WB and ZB vectors
    plt.quiver(OB.real, OB.imag, WB.real, WB.imag, angles='xy', scale_units='xy', scale=1, color='m', label='W_B')
    B = OB + WB
    plt.quiver(B.real, B.imag, ZB.real, ZB.imag, angles='xy', scale_units='xy', scale=1, color='g', label='Z_B')

    # ⚫ Precision Points & Coupler Line
    plt.plot([P2.real, P3.real], [P2.imag, P3.imag], 'k--', label='Coupler (P2→P3)')
    plt.plot(p1.real, p1.imag, 'ko')
    plt.plot(P2.real, P2.imag, 'ks')
    plt.plot(P3.real, P3.imag, 'ks')
    plt.text(p1.real + 1, p1.imag, "PP1")
    plt.text(P2.real + 1, P2.imag, "PP2")
    plt.text(P3.real + 1, P3.imag, "PP3")

    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend()
    plt.show()

    print("\n✅ Ground Pivots:")
    print("OA:", OA)
    print("OB:", OB)
    L1 = abs(WA)
    L2 = abs(WA + ZA)
    L3 = abs(WB)
    L4 = np.linalg.norm([OA.real - OB.real, OA.imag - OB.imag])  # Ground

    lengths = [L1, L2, L3, L4]
    S = min(lengths)
    L = max(lengths)
    P, Q = sorted(lengths)[1:3]

    print("L1 (Crank):", L1)
    print("L2 (Coupler):", L2)
    print("L3 (Rocker):", L3)
    print("L4 (Ground):", L4)
    print("Grashof Check:", S + L, "<=", P + Q, "→", S + L <= P + Q)

    # ✅ Return useful values too
    return {
        "OA": OA, "OB": OB,
        "WA": WA, "ZA": ZA,
        "WB": WB, "ZB": ZB,
        "L1": abs(WA),
        "L2": abs(WA + ZA),
        "L3": abs(WB)
    }

# === USAGE ===
if __name__ == "__main__":
    pps = [(0, 0), (15, 36), (18, 10)]
    thetas = [0, 200, 281]
    beta = [250, 190]
    gamma = [240, 220]

    result = synthesize_and_plot_matlab_style(pps, thetas, beta, gamma)
