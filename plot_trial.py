import matplotlib.pyplot as plt

def plot_spectral_data(filename):
    freqs = []
    col2 = []
    col3 = []

    with open(filename, "r") as f:
        lines = f.readlines()

    data_start = False
    for line in lines:
        line = line.strip()

        if line.startswith("Frequency (MHz)"):
            data_start = True
            continue

        if not data_start or not line:
            continue

        parts = line.split(",")
        if len(parts) < 3:
            continue

        try:
            fval = float(parts[0])
            v2 = float(parts[1])
            v3 = float(parts[2])
        except ValueError:
            continue

        freqs.append(fval)
        col2.append(v2)
        col3.append(v3)

    plt.figure()
    plt.plot(freqs, col2, label="avgps")
    plt.plot(freqs, col3, label="avgps2")
    plt.plot(freqs, [a - b for a, b in zip(col3, col2)], label="avgps2 - avgps")
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("Value")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    plot_spectral_data("elevation_110/trial_2025-11-26_13-02-15.csv")