import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import math




# =============================================================================
# DANE JĄDROWE – przykładowy słownik z danymi (λ [1/s])
# =============================================================================
nuclear_data = {
    "C11": 0.0005672978299613248,
    "N13": 0.00115907,
    "O15": 0.00566,
    "F18": 0.000105277,
    "Na22":  0.0000000084387,
    # Dodaj kolejne izotopy w razie potrzeby
}


# =============================================================================
# FUNKCJE POMOCNICZE DOTYCZĄCE MACIERZY A ORAZ KOREKT CZASOWYCH
# =============================================================================
def load_rotation_scheme(filename):
    """
    Wczytuje schemat rotacji z pliku tekstowego.
    Każda linia powinna zawierać: shift (int) oraz czas pomiaru (float).
    Zwraca słownik: {"shifts": shifts, "dt": dt}.
    """
    try:
        data = np.loadtxt(filename)
        shifts = data[:, 0].astype(int).tolist()
        dt = data[:, 1].tolist()
        return {"shifts": shifts, "dt": dt}
    except Exception as e:
        raise Exception(f"Błąd przy wczytywaniu schematu rotacji: {e}")


def load_rotation_times(filename):
    """
    Wczytuje plik z czasami rotacji.
    Jeżeli pierwszy element jest bardzo duży (np. w pikosekundach) to przeskalowujemy.
    Zwraca listę czasów.
    """
    try:
        times = np.loadtxt(filename)
        # Jeśli skala jest bardzo duża, przeskaluj (przykładowa logika)
        if np.any(times > 1e10):
            times = times / 1e12
        return times.tolist() if isinstance(times, np.ndarray) else times
    except Exception as e:
        raise Exception(f"Błąd przy wczytywaniu pliku z czasami rotacji: {e}")



def build_A_matrix(A_single, detector, general_params, lam_vector, n_measurements):
    """
    Buduje macierz A o rozmiarze (n_measurements, n_iso) – dla każdego pomiaru nakładamy
    przesuniętą wersję wektora A_single, skalowaną wektorem czynników zależnym od czasu oraz λ.
    
    Parametry:
      A_single      - lista (o długości n_sources) wartości dla każdego źródła
      detector      - słownik zawierający ścieżki do plików schematu rotacji,
                      np. klucze "rotation_scheme_file" oraz opcjonalnie "rotation_times_file"
      general_params- słownik z parametrami ogólnymi (m.in. "intensity", "time_delay")
      lam_vector    - wektor stałych rozpadu (lista liczb) używanych do budowy A
      n_measurements- maksymalna liczba pomiarów do rozważenia (lub None)
    
    Zwraca:
      A_mat, m, scheme, A_single, cs, dt_list, t_starts
    """
    n_sources = len(A_single)
    intensity = general_params["intensity"]
    time_delay = general_params["time_delay"]

    # Wczytanie schematu rotacji
    scheme = load_rotation_scheme(detector["rotation_scheme_file"])
    num_shifts = len(scheme["shifts"])
    # Upewniamy się, że lista dt ma odpowiednią długość
    if len(scheme["dt"]) < num_shifts:
        scheme["dt"] += [scheme["dt"][-1]] * (num_shifts - len(scheme["dt"]))
    
    m_total = num_shifts

    # Wyznaczenie dt_list – w oparciu o plik z czasami roatacji lub schemat
    if "rotation_times_file" in detector and detector["rotation_times_file"]:
        rotation_times = load_rotation_times(detector["rotation_times_file"])
        dt_list = [rotation_times[2 * i + 1] - rotation_times[2 * i] 
                   for i in range((len(rotation_times) - 1) // 2)]
        m_file = len(dt_list)
        m = n_measurements if n_measurements is not None else m_file
        m = min(m, m_file, m_total)
    else:
        m = n_measurements if n_measurements is not None else m_total
        m = min(m, m_total)
        dt_list = scheme["dt"][:m]

    # Obliczenie przesunięć (rolling index) według schematu
    cs = [0] * m
    cs[0] = scheme["shifts"][0] % n_sources
    for i in range(1, m):
        cs[i] = (cs[i - 1] + scheme["shifts"][i]) % n_sources

    # Wyznaczenie czasów rozpoczęcia pomiarów
    t_starts = [time_delay]
    for i in range(1, m):
        t_starts.append(t_starts[i - 1] + dt_list[i - 1])
    
    # Budowa macierzy A
    n_iso = len(lam_vector)
    A_mat = np.zeros((m, n_sources * n_iso), dtype=np.float64)
    for i in range(m):
        # Przesunięty wektor – wykorzystujemy np.roll z dodatkowym przesunięciem (+1)
        rolled = np.roll(A_single, cs[i] + 1)
        dt_i = dt_list[i]
        for j, lam_val in enumerate(lam_vector):
            # Obliczamy czynnik skalujący dla konkretnej λ
            factor = (1 - np.exp(-lam_val * dt_i)) * np.exp(-lam_val * t_starts[i]) * intensity
            # Produkt element-po-elemencie – skalujemy cały wektor
            A_mat[i, 16 * j:16 * (j + 1)] = np.array(rolled) * factor
    return A_mat, m, scheme, A_single, cs, dt_list, t_starts

def correct_counts(y_raw, y_err, detector, general_params, m, scheme):
    """
    Korekta zliczeń – mnoży przez czas pomiaru.
    Jeśli schemat rotacji jest dostępny, dla i-tego pomiaru używamy scheme["dt"][i],
    w przeciwnym razie dzielimy pomiary na krótkie i długie.
    """
    y_cor = y_raw.copy()
    y_err_cor = y_err.copy()
    for i in range(m):
        dt = scheme["dt"][i]
        y_cor[i] #*= dt
        y_err_cor[i] #*= dt
    return y_cor, y_err_cor


def plot_each_fit(A_single, cs, dt_list, t_starts, y, intensity, lam, tol=1e-6):
    """
    Dla każdego źródła (slotu) dopasowujemy model f(F)=A0*F,
    gdzie czynnik F obliczany jest jako:
      F = (1 - exp(-lam*dt)) * exp(-lam*t_start) * intensity * A_single[slot]
    i rysujemy dane oraz dopasowaną krzywą.
    """
    m = len(dt_list)
    n_sources = len(A_single)
    ncols = 4
    nrows = int(np.ceil(n_sources / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
    axs = axs.flatten() if n_sources > 1 else [axs]

    for s in range(n_sources):
        F_list, y_list = [], []
        for i in range(m):
            rolled = np.roll(A_single[::-1], cs[i] + 1)
            if abs(rolled[s] - A_single[s]) < tol and A_single[s] > 0:
                F_value = (1 - np.exp(-lam * dt_list[i])) * np.exp(-lam * t_starts[i]) * intensity * A_single[s]
                F_list.append(F_value)
                y_list.append(y[i])
        ax = axs[s]
        if F_list:
            F_arr = np.array(F_list)
            y_arr = np.array(y_list)
            # Definicja modelu liniowego
            def model(F, A0): 
                return A0 * F
            try:
                popt, _ = curve_fit(model, F_arr, y_arr)
                A0_val = popt[0]
                F_fit = np.linspace(F_arr.min(), F_arr.max(), 100)
                y_fit = model(F_fit, A0_val)
                ax.plot(F_arr, y_arr, 'bo', label='Dane')
                ax.plot(F_fit, y_fit, 'ro', label=f'Fit: A0={A0_val:.2e}')
            except Exception as e:
                ax.text(0.5, 0.5, f"Błąd fit: {e}", transform=ax.transAxes, ha='center')
            ax.set_title(f"Slot {s+1}")
            ax.set_xlabel("Czynnik (F)")
            ax.set_ylabel("Zliczenia")
            ax.legend(fontsize=8)
            ax.grid(True)
        else:
            ax.text(0.5, 0.5, "Brak danych", transform=ax.transAxes, ha='center')
            ax.set_title(f"Slot {s+1}")
    plt.tight_layout()
    plt.show()


def alternative_exponential_decay_fit(A_single, cs, dt_list, t_starts, y, u_y, lam, intensity=1.99):
    """
    Dla każdego źródła (slotu) dopasowuje model sumy funkcji zaniku promieniotwórczego
    dla n izotopów według wzoru:
      y_i ≈ Σ[j=0..n_iso-1] A0_j * (1 - exp(-λ_j * dt_i)) * exp(-λ_j * t0_i) * max(A_single)
    Przy wykorzystaniu ważonej regresji liniowej (WLS):
      A_est = (FᵀWF)⁻¹ FᵀWy
    Zwraca listy dopasowanych wartości A₀ oraz ich niepewności dla każdego izotopu (lista rozmiaru n_iso, każda zawiera 16 elementów).
    """
    n_sources = len(A_single)
    m = min(len(dt_list), len(cs))
    n_iso = len(lam)

    # Grupowanie pomiarów dla poszczególnych slotów
    source_data = {s: {"t0": [], "dt": [], "y": [], "u_y": []} for s in range(n_sources)}
    A_single_arr = np.array(A_single)
    for i in range(m):
        rolled = np.roll(A_single_arr, cs[i] + 1)
        s_measured = int(np.argmax(rolled))
        source_data[s_measured]["t0"].append(t_starts[i])
        source_data[s_measured]["dt"].append(dt_list[i])
        source_data[s_measured]["y"].append(y[i])
        source_data[s_measured]["u_y"].append(u_y[i])
    
    # Inicjalizacja wynikowych tablic
    fitted_A0 = [[0.0] * n_sources for _ in range(n_iso)]
    fitted_A0_err = [[0.0] * n_sources for _ in range(n_iso)]

    # Przygotowanie wykresów dla każdego slotu
    ncols = 4
    nrows = int(np.ceil(n_sources / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
    axs = axs.flatten()

    # Model sumy funkcji zaniku dla n izotopów
    def model_func(x, *A0_params):
        t0, dt = x
        result = np.zeros_like(t0)
        for j in range(n_iso):
            lambda_j = nuclear_data[lam[j]]
            result += A0_params[j] * (1 - np.exp(-lambda_j * dt)) * np.exp(-lambda_j * t0) * max(A_single) * intensity
        return result

    for s in range(n_sources):
        t0_arr = np.array(source_data[s]["t0"])
        dt_arr = np.array(source_data[s]["dt"])
        y_arr = np.array(source_data[s]["y"])
        uy_arr = np.array(source_data[s]["u_y"])
        ax = axs[s]

        if t0_arr.size == 0:
            for j in range(n_iso):
                fitted_A0[j][s] = 0.0
                fitted_A0_err[j][s] = 0.0
            ax.text(0.5, 0.5, "Brak danych", transform=ax.transAxes, ha="center", va="center")
            ax.set_title(f"Slot {s+1}")
            ax.set_xlabel("t₀ [jedn.]")
            ax.set_ylabel("Zliczenia")
            ax.grid(True)
            continue

        xdata = (t0_arr, dt_arr)
        p0 = [max(y_arr) / n_iso] * n_iso
        try:
            popt, pcov = curve_fit(model_func, xdata, y_arr, p0=p0, sigma=uy_arr, absolute_sigma=True)
            for j in range(n_iso):
                # Skalujemy oszacowania przez odpowiednią stałą λ (nie wariancję, lecz odchylenie)
                fitted_A0[j][s] = popt[j] * nuclear_data[lam[j]]
                var_j = pcov[j, j]
                fitted_A0_err[j][s] = math.sqrt(var_j) * nuclear_data[lam[j]] if var_j > 0 else 0.0

            t0_fine = np.linspace(t0_arr.min(), t0_arr.max(), 300)
            # Interpolujemy dt_fine na podstawie t0_arr i odpowiadających im dt_arr
            dt_fine = np.interp(t0_fine, t0_arr, dt_arr)
            y_model = model_func((t0_fine, dt_fine), *popt)
            ax.errorbar(t0_arr, y_arr, yerr=uy_arr, fmt='o', color='blue', label='Dane')
            ax.plot(t0_fine, y_model, 'r-', label="Fit")
            legend_text = "\n".join([f"Iso{j+1}: A₀={popt[j]:.2e}±{fitted_A0_err[j][s]:.2e}" for j in range(n_iso)])
            ax.text(0.05, 0.95, legend_text, transform=ax.transAxes, va="top", fontsize=7, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.7))
        except Exception as e:
            ax.text(0.5, 0.5, f"Błąd fit: {e}", transform=ax.transAxes, ha="center")
        ax.set_title(f"Slot {s+1}")
        ax.set_xlabel("t₀ [jedn.]")
        ax.set_ylabel("Zliczenia")
        ax.legend(fontsize=8)
        ax.grid(True)
    plt.tight_layout()
    # plt.show()

    print("Wyniki dopasowania A₀ (λ znane) dla kolejnych izotopów i slotów:")
    for j in range(n_iso):
        print(f"\nIzotop {lam[j]} (λ={nuclear_data[lam[j]]:.2e}):")
        for s in range(n_sources):
            print(f"  Slot {s+1}: A₀ = {(fitted_A0[j][s]):.2e} ± {(fitted_A0_err[j][s]):.2e}")
    return fitted_A0, fitted_A0_err