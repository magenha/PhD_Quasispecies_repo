import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from os import listdir
from os.path import isfile, join


def extract_fName_params(filename):
    """
    Function that, given a file name, extract the parameters
    """
    params_dict = {}
    filename = filename.replace(".csv", "")
    helper = filename.split("_")
    for i in range(1, len(helper)):
        if helper[i] == "mu":
            params_dict["mu"] = float(helper[i + 1])
        elif helper[i] == "p":
            params_dict["p"] = float(helper[i + 1])
        elif helper[i] == "nrw":
            params_dict["nrw"] = int(helper[i + 1])
        elif helper[i] == "L":
            params_dict["L"] = int(helper[i + 1])
        elif helper[i] == "origin":
            params_dict["origin"] = str(helper[i + 1])
        elif helper[i] == "t":
            params_dict["t"] = int(helper[i + 1])

    return params_dict


greek_letterz = [chr(code) for code in range(945, 970)]
# Declaration of variables
DataPath = "./Data/"

onlyfiles = [f for f in listdir(DataPath) if isfile(join(DataPath, f))]
onlyfiles.sort()

import matplotlib

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 14

desired_origins = ["n0", "n10", "n150"]
# desired_origins = ["n10"]

# Plot the complete saturated curve for mu=1, for both euclidean and hamming
if 0 == 1:
    plt.figure()
    run_counter = 0
    for j, fName in enumerate(onlyfiles):
        params = extract_fName_params(fName)

        if (
            params["mu"] == 1.0
            and params["L"] == 1000
            and params["origin"] in desired_origins
        ):
            # Read the file
            df = pd.read_csv(DataPath + fName)
            n0 = int(params["origin"].replace("n", ""))
            # Plot analytical predictions first
            if run_counter == 0:
                x = np.arange(1, max(df["mean_hamming"]))
                plt.plot(x, x, "--", c="blue")
                x = np.arange(1, len(df["mean_hamming"]))
                y = np.array([0.5 * params["L"] * params["mu"] for t in x])
                plt.plot(x, y, "--", c="black")
                # Plot analytical curve
                y = (
                    0.5
                    * params["L"]
                    * params["mu"]
                    * (
                        1
                        - (1 - 2 * n0 / params["L"])
                        * np.e ** (-1.0 * x / (0.5 * params["L"] * params["mu"]))
                    )
                )

                plt.plot(x, y, "--", c="green", linewidth=3)

            # Plot mean Hamming and Euclidean
            plt.plot(
                df.index + 1,
                df["mean_hamming"],
                label=f"{params['origin']} $D_H$",
            )
            plt.plot(
                df.index + 1,
                df["mean_euclidean"] ** 2,
                label=f"{params['origin']} $D_E^2$",
            )

    plt.legend(loc="best")

    plt.ylabel("Average distance <m>")
    plt.xlabel(f"Time (t)")
    plt.title(f"$\mu$={params['mu']}")
    # plt.title(f"{greek_letterz[11]}=1.0")

    logscale = False
    if logscale:
        plt.yscale("log")
        plt.xscale("log")
    plt.locator_params(axis="y", nbins=5)
    plt.locator_params(axis="x", nbins=5)
    # plt.locator_params(numticks=12)
    plt.savefig("./Figures/1.a.1.pdf", dpi=300, format="pdf")
    plt.savefig("./Figures/1.a.1.png", dpi=300, format="png")


# Plot different mu and Ls in a plot, and then collapse them.
if 0 == 1:
    plt.figure()
    run_counter = 0
    for j, fName in enumerate(onlyfiles):
        params = extract_fName_params(fName)

        if params["mu"] != 1.0 or params["L"] != 1000:
            run_counter += 1
            # Read the file
            df = pd.read_csv(DataPath + fName)

            # Plot mean Hamming and Euclidean
            plt.plot(
                df.index + 1,
                df["mean_hamming"],
                label=f"{params['L']}, {params['mu']} Hamming",
            )

    plt.legend(loc="best")

    plt.ylabel("Average distance <m>")
    plt.xlabel(f"Time (t)")
    plt.title(f"$\mu$={params['mu']}")
    # plt.title(f"{greek_letterz[11]}=1.0")

    logscale = False
    if logscale:
        plt.yscale("log")
        plt.xscale("log")
    plt.locator_params(axis="y", nbins=5)
    plt.locator_params(axis="x", nbins=5)
    # plt.locator_params(numticks=12)
    plt.savefig("./Figures/1.a.2.pdf", dpi=300, format="pdf")
    plt.savefig("./Figures/1.a.2.png", dpi=300, format="png")

    plt.figure()
    run_counter = 0
    for j, fName in enumerate(onlyfiles):
        params = extract_fName_params(fName)

        if params["mu"] != 1.0 or params["L"] != 1000:
            run_counter += 1
            # Read the file
            df = pd.read_csv(DataPath + fName)

            # Plot mean Hamming and Euclidean
            x = (df.index + 1) * params["mu"] / (params["L"] / 2)
            y = (df["mean_hamming"]) / (params["L"] / 2)
            plt.plot(
                x,
                y,
                label=f"{params['L']}, {params['mu']} Hamming",
            )

    plt.legend(loc="best")

    plt.ylabel("Collapsed Average distance <m/$\lambda$>")
    plt.xlabel(f"Collapsed Time ($\mu$t/$\lambda$)")
    plt.title(f"$\mu$={params['mu']}")
    # plt.title(f"{greek_letterz[11]}=1.0")

    logscale = False
    if logscale:
        plt.yscale("log")
        plt.xscale("log")
    plt.locator_params(axis="y", nbins=5)
    plt.locator_params(axis="x", nbins=5)
    # plt.locator_params(numticks=12)
    plt.savefig("./Figures/1.a.3.pdf", dpi=300, format="pdf")
    plt.savefig("./Figures/1.a.3.png", dpi=300, format="png")

# Plot Index of dispersion and variance
if 0 == 1:
    plt.figure()
    run_counter = 0
    for j, fName in enumerate(onlyfiles):
        params = extract_fName_params(fName)

        if params["L"] == 1000 and params["origin"] in desired_origins:
            # Read the file
            print(fName)
            df = pd.read_csv(DataPath + fName)
            n0 = int(params["origin"].replace("n", ""))
            x = df.index + 1
            y = df["variance_hamming"]

            plt.plot(x, y, label=f"$D_H$, origin {n0}")
            y = df["variance_euclidean"]

            plt.plot(x, y, label=f"$D_E$, origin {n0}")
            # plt.plot(x, y, label=f"$\mu$={params['mu']}")

    plt.xlabel("Time (t)")
    plt.ylabel("Variance (${\sigma}^2$)")

    plt.legend(loc="best")
    plt.tight_layout()

    plt.savefig("./Figures/1.a.4.svg", dpi=300, format="svg")
    plt.savefig("./Figures/1.a.4.png", dpi=300, format="png")

    plt.figure()
    run_counter = 0
    for j, fName in enumerate(onlyfiles):
        params = extract_fName_params(fName)
        print(fName)
        if params["L"] == 1000 and params["origin"] in desired_origins:
            # Read the file
            df = pd.read_csv(DataPath + fName)
            n0 = int(params["origin"].replace("n", ""))
            x = df.index + 1
            y = df["variance_hamming"] / df["mean_hamming"]

            # plt.plot(x, y, label=f"$\mu$={params['mu']}")
            plt.plot(x, y, label=f"$D_H$, origin {n0}")
            y = df["variance_euclidean"] / df["mean_euclidean"]

            # plt.plot(x, y, label=f"$\mu$={params['mu']}")
            plt.plot(x, y, label=f"$D_E$, origin {n0}")
    y_aster = np.array([1 for a in x])
    plt.plot(x, y_aster, "--", c="black")

    plt.xlabel("Time (t)")
    plt.ylabel("Index of Dispersion (Variance/Average)")

    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig("./Figures/1.a.5.svg", dpi=300, format="svg")
    plt.savefig("./Figures/1.a.5.png", dpi=300, format="png")

# Plot Distribution of distances

if 0 == 1:
    DataPathDistrs = "./Data/Distribs/"

    onlyfiles = [f for f in listdir(DataPathDistrs) if isfile(join(DataPathDistrs, f))]
    onlyfiles.sort()

    for j, fName in enumerate(onlyfiles):
        df = pd.read_csv(DataPathDistrs + fName)
        params = extract_fName_params(fName)
        if params["origin"] in desired_origins:
            data = df["d_hamming"].to_numpy()
            hist, bin_edges = np.histogram(data, bins=10)
            x = np.array(
                [
                    (bin_edges[i + 1] + bin_edges[i + 1]) / 2
                    for i in range(len(bin_edges) - 1)
                ]
            )
            plt.plot(bin_edges[1:], hist, "o-", label=params["t"])

            # Locate mean and variance
            file = DataPath + fName.replace(f"_t_{params['t']}", "")
            df_2 = pd.read_csv(file)
            mean = float(df_2["mean_hamming"][params["t"] - 1])
            std_dev = np.sqrt(float(df_2["variance_hamming"][params["t"] - 1]))
            print("mean value=", mean)

            y_grid = np.linspace(1, 300, 100)
            x_grid = np.repeat(mean, y_grid.shape[0])

            plt.plot(x_grid, y_grid, "--")
            plt.annotate(
                "",
                xy=(mean + std_dev, max(hist) * 1.1),
                xytext=(mean, max(hist) * 1.1),
                arrowprops=dict(
                    facecolor="red", arrowstyle="-|>", shrinkA=0, shrinkB=0
                ),
            )
            plt.annotate(
                "",
                xy=(mean - std_dev, max(hist) * 1.1),
                xytext=(mean, max(hist) * 1.1),
                arrowprops=dict(
                    facecolor="red", arrowstyle="-|>", shrinkA=0, shrinkB=0
                ),
            )
    plt.ylabel("probability density")
    plt.xlabel("Hamming distance to origin (m)")
    plt.tight_layout()
    plt.legend(loc="best")

    plt.savefig(f"./Figures/1.a.6_{params['origin']}.svg", dpi=300, format="svg")
    plt.savefig(f"./Figures/1.a.6_{params['origin']}.png", dpi=300, format="png")


# Search when they collapse, first plot  the graphs
if 0 == 1:
    eps = 0.1

    plt.figure()
    run_counter = 0
    for j, fName in enumerate(onlyfiles):
        params = extract_fName_params(fName)

        if params["origin"] in desired_origins:
            # Read the file
            print(fName)
            df = pd.read_csv(DataPath + fName)
            n0 = int(params["origin"].replace("n", ""))
            x = df.index.to_numpy() + 1
            y = df["mean_hamming"]
            delta = 1 - y / ((1) * params["mu"] * x)
            time = (delta - eps).apply(abs).idxmin() + 1
            print(time)

            # plt.plot(x, y, label=f"$D_H$, origin {n0}")
            plt.plot(x, delta, label=f"top, origin {n0}")
            plt.plot(
                [time for a in range(100)],
                np.linspace(0, delta[time - 1], 100),
                "--",
                c="black",
            )
            plt.plot(
                np.linspace(0, time, 10),
                [delta[time - 1] for a in range(10)],
                "--",
                c="black",
            )

    plt.xlabel("Time (t)")
    plt.ylabel("Relative deviation $(m_L - m(t))/(m_L)$")

    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig("./Figures/dtime.pdf", dpi=300, format="pdf")

    plt.xlim(0, 120)
    plt.ylim(0, 0.125)
    plt.savefig("./Figures/dtime-limited.pdf", dpi=300, format="pdf")


# Search when they collapse, now plot the time for different conditions
if 1 == 1:
    eps = 0.1

    df = pd.read_csv("./Data/deviation_time_ok.csv")
    t_results = df["time"].to_numpy()
    L_results = df["L"].to_numpy()
    mu_results = df["mu"].to_numpy()

    plt.figure()
    t_teo = L_results * eps / mu_results
    log_L_values = np.log10(L_results)

    # Scatter plot with color intensity proportional to log(mu_results)
    scatter = plt.scatter(t_teo, t_results, c=log_L_values, cmap="viridis")
    plt.plot(
        np.linspace(0, max(t_results), 100),
        np.linspace(0, max(t_results), 100),
        "--",
        c="black",
    )

    # Add colorbar and label
    cbar = plt.colorbar(scatter)
    cbar.set_label("Log L", rotation=270, labelpad=15)

    plt.ylabel("simulated")
    plt.xlabel("predicted")
    plt.yscale("log")
    plt.xscale("log")
    plt.tight_layout()
    plt.savefig("./Figures/tau-sim-red_Lscale.pdf", dpi=300, format="pdf")

    plt.figure()
    t_teo = L_results * eps / mu_results
    log_L_values = np.log10(mu_results)

    # Scatter plot with color intensity proportional to log(mu_results)
    scatter = plt.scatter(t_teo, t_results, c=log_L_values, cmap="viridis")
    plt.plot(
        np.linspace(0, max(t_results), 100),
        np.linspace(0, max(t_results), 100),
        "--",
        c="black",
    )

    # Add colorbar and label
    cbar = plt.colorbar(scatter)
    cbar.set_label("Log mu", rotation=270, labelpad=15)

    plt.ylabel("simulated")
    plt.xlabel("predicted")
    plt.yscale("log")
    plt.xscale("log")
    plt.tight_layout()
    plt.savefig("./Figures/tau-sim-red_muscale.pdf", dpi=300, format="pdf")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    # Sample data for illustration
    x = np.linspace(min(L_results), max(L_results), 100)
    y = np.linspace(min(mu_results), max(mu_results), 100)
    X, Y = np.meshgrid(x, y)
    t_teo = X * eps / Y
    # Plotting the theoretical surface with 50% transparency
    surf = ax.plot_surface(X, Y, t_teo, cmap="viridis", alpha=0.5)

    # Scatter plot for actual data points
    ax.scatter(L_results, mu_results, t_results, color="red")

    # Set labels for axes
    ax.set_xlabel("L")
    ax.set_ylabel("mu")
    ax.set_zlabel("Time")
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # ax.set_zscale("log")


plt.show()
