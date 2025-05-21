import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def rotation_matrix(theta):
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta),  np.cos(theta)]])

def freudenstein_theta4(OA, OB, A, L2, L3):
    x1, y1 = A
    x2, y2 = OB
    L = np.linalg.norm([x2 - x1, y2 - y1])

    if L > L2 + L3 or L < abs(L2 - L3):
        return None, None

    # Law of cosines
    a = (L2**2 - L3**2 + L**2) / (2 * L)
    h = np.sqrt(abs(L2**2 - a**2))

    P2 = A + a * (OB - A) / L
    offset = h * np.array([-(OB[1] - A[1]) / L, (OB[0] - A[0]) / L])

    B1 = P2 + offset
    B2 = P2 - offset
    return B1, B2

def simulate_physics_based_linkage(OA, OB, L1, L2, L3, steps=180):
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.grid(True)
    ax.set_title("4-Bar Linkage (Freudenstein-Based)")

    ln1, = ax.plot([], [], 'b-', label="Input")
    ln2, = ax.plot([], [], 'r-', label="Coupler")
    ln3, = ax.plot([], [], 'm-', label="Output")
    pivots, = ax.plot([], [], 'ko')
    path, = ax.plot([], [], 'g--', lw=1, alpha=0.6, label="Tip Path")

    tip_trace_x, tip_trace_y = [], []

    def init():
        ax.set_xlim(-40, 40)
        ax.set_ylim(-40, 40)
        return ln1, ln2, ln3, pivots, path

    def update(i):
        angle = np.radians(i * 360 / steps)
        A = OA + rotation_matrix(angle) @ np.array([L1, 0])

        B1, B2 = freudenstein_theta4(OA, OB, A, L2, L3)
        if B1 is None:
            return ln1, ln2, ln3, pivots, path

        B = B1  # Choose the "elbow up" configuration

        ln1.set_data([OA[0], A[0]], [OA[1], A[1]])
        ln2.set_data([A[0], B[0]], [A[1], B[1]])
        ln3.set_data([B[0], OB[0]], [B[1], OB[1]])
        pivots.set_data([OA[0], OB[0], A[0], B[0]], [OA[1], OB[1], A[1], B[1]])

        tip_trace_x.append(B[0])
        tip_trace_y.append(B[1])
        path.set_data(tip_trace_x, tip_trace_y)

        return ln1, ln2, ln3, pivots, path

    ani = animation.FuncAnimation(fig, update, frames=steps, init_func=init, blit=True, interval=30)
    plt.legend()
    plt.show()

# === Run It ===
if __name__ == "__main__":
    OA = np.array([15.989230156977088,11.460117781555518])
    OB = np.array([14.805051905004728,11.343752880475726])
    L1 = 11 # Crank (OA→A)
    L2 = 6  # Coupler (A→B)
    L3 = 12  # Rocker (B→OB)

    simulate_physics_based_linkage(OA, OB, L1, L2, L3)
