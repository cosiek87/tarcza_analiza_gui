import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import N_A, pi
from scipy.optimize import nnls
from numpy.linalg import inv

np.set_printoptions(precision=15, suppress=False)

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
# FUNKCJE METRYK JAKOŚCI
# =============================================================================
def aFactor(y_exp, y_theor):
    """Metryka aFactor – średnia wartość |y - y_est|/(y+y_est)."""
    eps = 1e-12
    return np.mean(np.abs(y_exp - y_theor) / (y_exp + y_theor + eps))

def compute_rmse(y, y_est):
    return np.sqrt(np.mean((y - y_est) ** 2))

def compute_mae(y, y_est):
    return np.mean(np.abs(y - y_est))

def compute_r2(y, y_est):
    ss_res = np.sum((y - y_est) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1 - ss_res / ss_tot if ss_tot != 0 else 0

def compute_pearson(y, y_est):
    if np.std(y) == 0 or np.std(y_est) == 0:
        return 0
    return np.corrcoef(y, y_est)[0, 1]

# =============================================================================
# FUNKCJE POMOCNICZE DOTYCZĄCE MACIERZY A ORAZ KOREKT CZASOWYCH
# =============================================================================
def load_rotation_scheme(filename):
    """
    Wczytuje schemat rotacji z pliku tekstowego.
    Każda linia powinna zawierać: przesunięcie (int) oraz czas pomiaru (float).
    """
    try:
        data = np.loadtxt(filename, delimiter=None)
        shifts = data[:, 0].astype(int).tolist()
        dt = data[:, 1].tolist()
        return {"shifts": shifts, "dt": dt}
    except Exception as e:
        raise Exception(f"Błąd przy wczytywaniu schematu rotacji: {e}")


def load_rotation_times(filename):
    """
    Wczytuje plik z czasami rotacji.
    Każda linia zawiera jedną wartość czasu.
    Jeśli pierwsza wartość nie wynosi 0, zostanie ona dodana na początku.
    """
    try:
        times = np.loadtxt(filename)
        # Jeśli pierwszy czas nie jest równy 0, dodajemy 0 na początku
        if times[3] > 1e10:
            times = times / 1e12
            
        times = times.tolist() if isinstance(times, np.ndarray) else times
        if times[0] != 0:
            times = [0.0] + times
        return times
    except Exception as e:
        raise Exception(f"Błąd przy wczytywaniu pliku z czasami rotacji: {e}")


def build_A_matrix(A_single, detector, general_params, lam_vector, n_measurements):
    """
    Buduje macierz A o rozmiarze (n_measurements, n_iso).
    
    Jeśli detektor posiada plik schematu rotacji (klucz "rotation_scheme_file"),
    wykorzystujemy go do ustalenia przesunięć tarczy.
    
    Jeśli dodatkowo podany jest plik z czasami rotacji (klucz "rotation_times_file"),
    obliczamy czasy pomiarów (dt) jako:
      - pierwszy pomiar: od 0 do pierwszej wartości (po dodaniu 0, jeśli nie ma)
      - kolejne pomiary: różnica między 2n-1 a 2n wartością z pliku.
    """
    n_sources = len(A_single)
    intensity = general_params["intensity"]
    time_delay = general_params["time_delay"]

    # Wczytanie schematu rotacji
    scheme = load_rotation_scheme(detector["rotation_scheme_file"])
    if len(scheme["dt"]) < len(scheme["shifts"]):
        scheme["dt"] = scheme["dt"] + [scheme["dt"][-1]] * (len(scheme["shifts"]) - len(scheme["dt"]))
    
    m_total = len(scheme["shifts"])

    # Jeśli podano plik z czasami rotacji, wykorzystujemy go do wyznaczenia dt_list
    if "rotation_times_file" in detector and detector["rotation_times_file"]:
        rotation_times = load_rotation_times(detector["rotation_times_file"])
        dt_list = []
        num_pairs = (len(rotation_times) - 1) // 2
        for i in range(num_pairs):
            dt_i = rotation_times[2 * i + 1] - rotation_times[2 * i]
            dt_list.append(dt_i)
        m_file = len(dt_list)
        m = n_measurements if n_measurements is not None else m_file
        m = min(m, m_file, m_total)
    else:
        m = n_measurements if n_measurements is not None else m_total
        m = min(m, m_total)
        dt_list = scheme["dt"][:m]

    # Obliczanie przesunięć (rolling) – na podstawie schematu
    cs = [0] * m
    cs[0] = scheme["shifts"][0] % n_sources
    for i in range(1, m):
        cs[i] = (cs[i - 1] + scheme["shifts"][i]) % n_sources

    # Wyznaczenie czasów rozpoczęcia pomiarów
    t_starts = [0] * m
    t_starts[0] = time_delay  # czas opóźnienia pomiędzy startem a pierwszym pomiarem
    for i in range(1, m):
        t_starts[i] = t_starts[i - 1] + dt_list[i - 1]

    # Budowa macierzy A: dla każdego pomiaru nakładamy przesuniętą wersję wektora A_single,
    # skalowaną czynnikiem zależnym od czasu pomiaru oraz λ.
    A_mat = np.zeros((m, n_sources*len(lam_vector)), dtype=np.float64)
    for i in range(m):
        rolled = np.roll(A_single, cs[i])
        dt_i = dt_list[i]
        for j in range(len(lam_vector)):
            factor = (1 - np.exp(-lam_vector[j] * dt_i)) * np.exp(-lam_vector[j] * t_starts[i]) * intensity
            A_mat[i][16*(j):16*(j+1)] = rolled * factor
    return A_mat, m, scheme

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

# =============================================================================
# FUNKCJA URUCHAMIAJĄCA ANALIZĘ – z obsługą błędów i komunikatów
# =============================================================================
def run_analysis(general_params, detectors):
    big_A_list = []
    big_y_list = []
    big_sigma_list = []
    detector_segments = []      # liczba pomiarów dla każdego detektora
    detector_names_success = [] # nazwy detektorów, dla których udało się wczytać dane
    A_blocks = []               # macierze A dla poszczególnych detektorów
    detector_has_scheme = []    # True, jeśli dla danego detektora dostępny był schemat rotacji
    errors = []                 # komunikaty o brakujących lub nieprawidłowych danych

    for det in detectors:
        try:
            all_data = np.loadtxt(det["A_file"])
        except Exception as e:
            errors.append(f"Detektor {det.get('name','Unnamed')}: Nie można wczytać pliku A ({det.get('A_file','')}) - {e}")
            continue

        # Jeśli podano izotopy, zakładamy, że:
        # - lista izotopów jest definiowana pojedynczo (np. "Cs-137, Co-60, ...")
        # - plik A_file zawiera pojedyńczy wektor (1D) z efektywnościami dla tych izotopów
        if det.get("isotopy", "").strip():
            isotopy = [iso.strip() for iso in det["isotopy"].split(",")]
            n_iso = len(isotopy)
            lam_vector = []
            for iso in isotopy:
                if iso in nuclear_data:
                    lam_vector.append(nuclear_data[iso])
                else:
                    messagebox.showwarning("Ostrzeżenie", f"Detektor {det['name']}: Nie znaleziono danych dla izotopu {iso}. Używam domyślnego lam.")
                    lam_vector.append(general_params["lam"])
            # Upewnij się, że all_data to wektor (1D)
            if all_data.ndim > 1:
                all_data = all_data[0]
        else:
            # Brak definicji izotopów – traktujemy A_file jako pojedyńczy wektor
            if all_data.ndim > 1:
                all_data = all_data[0]
            n_iso = len(all_data)
            lam_vector = [general_params["lam"]] * n_iso

        n_meas = int(det["n_measurements"]) if det.get("n_measurements", "") else n_iso * general_params["a"]

        try:
            A_block, m, scheme = build_A_matrix(all_data, det, general_params, lam_vector, n_meas)
        except Exception as e:
            errors.append(f"Detektor {det.get('name','Unnamed')}: Błąd przy budowaniu macierzy A - {e}")
            continue

        big_A_list.append(A_block)
        detector_segments.append(m)
        A_blocks.append(A_block)
        detector_names_success.append(det["name"])
        detector_has_scheme.append(scheme is not None)

        try:
            yz = np.loadtxt(det["counts_file"])
        except Exception as e:
            errors.append(f"Detektor {det.get('name','Unnamed')}: Nie można wczytać pliku zliczeń ({det.get('counts_file','')}) - {e}")
            continue
        try:
            y_raw = yz[:m, 0]
            y_err = yz[:m, 1]
        except Exception as e:
            errors.append(f"Detektor {det.get('name','Unnamed')}: Dane zliczeń mają nieprawidłowy format - {e}")
            continue
        try:
            y_cor, y_err_cor = correct_counts(y_raw, y_err, det, general_params, m, scheme)
        except Exception as e:
            errors.append(f"Detektor {det.get('name','Unnamed')}: Błąd przy korekcie zliczeń - {e}")
            continue
        big_y_list.append(y_cor)
        big_sigma_list.append(y_err_cor)

    if not big_A_list:
        # Zwracamy przynajmniej puste dane eksperymentalne, aby można było wyrysować wykres
        return {"errors": errors, "y": []}

    A_wysokie = np.vstack(big_A_list)
    y = np.hstack(big_y_list)
    y_sigma = np.hstack(big_sigma_list)
    N = A_wysokie.shape[0]
    k = A_wysokie.shape[1]

    try:
        W_sqrt = np.diag(1.0 / y_sigma)
    except Exception as e:
        return {"errors": errors + [f"Błąd przy obliczaniu wag: {e}"], "y": y.tolist()}
    Aw = W_sqrt @ A_wysokie
    yw = W_sqrt @ y

    try:
        x_nnls, rnorm = nnls(Aw, yw)
    except Exception as e:
        return {"errors": errors + [f"Błąd przy dopasowaniu NNLS: {e}"], "y": y.tolist()}

    residuals = Aw @ x_nnls - yw
    chi2 = np.sum(residuals ** 2)
    res_var = chi2 / (N - k) if N > k else chi2
    M = Aw.T @ Aw
    try:
        # Uncomment the following line if logging is needed for debugging
        # logging.debug(f"A_block: {A_block}")
        M_inv = inv(M)
    except Exception as e:
        return {"errors": errors + [f"Błąd przy inwersji macierzy: {e}"], "y": y.tolist()}
    cov_x = res_var * M_inv
    param_errors = np.sqrt(np.diag(cov_x))
    y_est = A_wysokie @ x_nnls
    # Obliczenie macierzy kowariancji dla y_est
    cov_y_est = A_wysokie @ cov_x @ A_wysokie.T
    # Obliczenie błędów dla każdego punktu y_est
    y_est_errors = np.sqrt(np.diag(cov_y_est))

    metrics_all = []
    index_start = 0
    for m in detector_segments:
        index_end = index_start + m
        y_det = y[index_start:index_end]
        y_est_det = y_est[index_start:index_end]
        met = {
            "aFactor": aFactor(y_det, y_est_det),
            "RMSE": compute_rmse(y_det, y_est_det),
            "MAE": compute_mae(y_det, y_est_det),
            "R2": compute_r2(y_det, y_est_det),
            "Pearson": compute_pearson(y_det, y_est_det)
        }
        metrics_all.append(met)
        index_start = index_end

    results = {
        "x_nnls": x_nnls.tolist(),
        "param_errors": param_errors.tolist(),
        "chi2": chi2,
        "detector_metrics": metrics_all,
        "detector_names": detector_names_success,
        "y": y.tolist(),
        "y_unc": y_sigma.tolist(),
        "y_est": y_est.tolist(),
        "y_est_errors": y_est_errors.tolist(),
        "A_wysokie_shape": A_wysokie.shape,
        "y_shape": y.shape,
        "A_blocks": A_blocks,
        "detector_segments": detector_segments,
        "detector_has_scheme": detector_has_scheme,
        "errors": errors
    }
    return results


# =============================================================================
# INTERFEJS GRAFICZNY – DIALOG DLA DETEKTORA
# =============================================================================
class DetectorDialog(tk.Toplevel):
    def __init__(self, master, detector=None):
        super().__init__(master)
        self.title("Dodaj/Edytuj detektor")
        self.geometry("800x300")
        self.detector = detector
        self.result = None
        self.create_widgets()
        if detector:
            self.populate_fields(detector)
        else:
            default_isotopy = "C11"
            self.entry_isotopy.insert(0, default_isotopy)

    def create_widgets(self):
        r = 0
        tk.Label(self, text="Nazwa:").grid(row=r, column=0, sticky="e")
        self.entry_name = tk.Entry(self, width=100)
        self.entry_name.grid(row=r, column=1, columnspan=2)
        r += 1
        tk.Label(self, text="Plik A (wydajność):").grid(row=r, column=0, sticky="e")
        self.entry_A_file = tk.Entry(self, width=100)
        self.entry_A_file.grid(row=r, column=1)
        tk.Button(self, text="Przeglądaj", command=self.browse_A_file).grid(row=r, column=2)
        r += 1
        tk.Label(self, text="Plik zliczeń:").grid(row=r, column=0, sticky="e")
        self.entry_counts_file = tk.Entry(self, width=100)
        self.entry_counts_file.grid(row=r, column=1)
        tk.Button(self, text="Przeglądaj", command=self.browse_counts_file).grid(row=r, column=2)
        r += 1
        tk.Label(self, text="Liczba pomiarów:").grid(row=r, column=0, sticky="e")
        self.entry_n_meas = tk.Entry(self, width=30)
        self.entry_n_meas.grid(row=r, column=1, columnspan=2)
        r += 1
        tk.Label(self, text="Izotopy (przecinek):").grid(row=r, column=0, sticky="e")
        self.entry_isotopy = tk.Entry(self, width=100)
        self.entry_isotopy.grid(row=r, column=1, columnspan=2)
        # Usunięto niekondycjonalne ustawianie "C11" – teraz wartość domyślna jest ustawiana w __init__
        r += 1
        tk.Label(self, text="Plik schematu rotacji:").grid(row=r, column=0, sticky="e")
        self.entry_rot_scheme = tk.Entry(self, width=100)
        self.entry_rot_scheme.grid(row=r, column=1)
        tk.Button(self, text="Przeglądaj", command=self.browse_rot_scheme).grid(row=r, column=2)
        r += 1
        tk.Label(self, text="Plik kroków rotacji:").grid(row=r, column=0, sticky="e")
        self.entry_steps_scheme = tk.Entry(self, width=100)
        self.entry_steps_scheme.grid(row=r, column=1)
        tk.Button(self, text="Przeglądaj", command=self.browse_steps_scheme).grid(row=r, column=2)
        r += 1
        tk.Button(self, text="OK", command=self.on_ok).grid(row=r, column=0, pady=10)
        tk.Button(self, text="Anuluj", command=self.destroy).grid(row=r, column=1, pady=10)
    
    def browse_A_file(self):
        filename = filedialog.askopenfilename(title="Wybierz plik A", 
                                              filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if filename:
            self.entry_A_file.delete(0, tk.END)
            self.entry_A_file.insert(0, filename)
    
    def browse_counts_file(self):
        filename = filedialog.askopenfilename(title="Wybierz plik zliczeń", 
                                              filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if filename:
            self.entry_counts_file.delete(0, tk.END)
            self.entry_counts_file.insert(0, filename)
    
    def browse_rot_scheme(self):
        filename = filedialog.askopenfilename(title="Wybierz plik schematu rotacji", 
                                              filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if filename:
            self.entry_rot_scheme.delete(0, tk.END)
            self.entry_rot_scheme.insert(0, filename)
    
    def browse_steps_scheme(self):
        filename = filedialog.askopenfilename(title="Wybierz plik kroków rotacji", 
                                              filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if filename:
            self.entry_steps_scheme.delete(0, tk.END)
            self.entry_steps_scheme.insert(0, filename)
    
    def populate_fields(self, detector):
        self.entry_name.insert(0, detector["name"])
        self.entry_A_file.insert(0, detector["A_file"])
        self.entry_counts_file.insert(0, detector["counts_file"])
        self.entry_n_meas.insert(0, str(detector.get("n_measurements", "")))
        self.entry_isotopy.insert(0, detector.get("isotopy", ""))
        self.entry_rot_scheme.insert(0, detector.get("rotation_scheme_file", ""))
        self.entry_steps_scheme.insert(0, detector.get("rotation_times_file", ""))
    
    def on_ok(self):
        try:
            name = self.entry_name.get()
            A_file = self.entry_A_file.get()
            counts_file = self.entry_counts_file.get()
            n_measurements = self.entry_n_meas.get().strip()
            n_measurements = int(n_measurements) if n_measurements else None
            isotopy = self.entry_isotopy.get().strip()
            rotation_scheme_file = self.entry_rot_scheme.get().strip()
            rotation_times_file = self.entry_steps_scheme.get().strip()
            self.result = {
                "name": name,
                "A_file": A_file,
                "counts_file": counts_file,
                "n_measurements": n_measurements,
                "isotopy": isotopy,
                "rotation_scheme_file": rotation_scheme_file,
                "rotation_times_file": rotation_times_file
            }
            self.destroy()
        except Exception as e:
            messagebox.showerror("Błąd", f"Błąd przy wprowadzaniu danych: {e}")

# =============================================================================
# OKNO GŁÓWNE APLIKACJI
# =============================================================================
class MainApplication(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Analiza danych eksperymentalnych")
        self.geometry("1400x1000")  # Duże okno
        self.general_params = {
            "lam": 0.0005672978299613248,
            "time_delay": 479.4,
            "TOB": 332.1,
            "intensity": 1.99534,
            "a": 5
        }
        self.detectors = []  # lista detektorów
        self.create_widgets()
        self.create_menu()
    
    def create_widgets(self):
        frame_params = tk.LabelFrame(self, text="Parametry ogólne")
        frame_params.pack(fill="x", padx=10, pady=5)
        tk.Label(frame_params, text="lam:").grid(row=0, column=0, sticky="e")
        self.entry_lam = tk.Entry(frame_params, width=15)
        self.entry_lam.grid(row=0, column=1, padx=5, pady=2)
        self.entry_lam.insert(0, str(self.general_params["lam"]))
        tk.Label(frame_params, text="Time delay:").grid(row=0, column=2, sticky="e")
        self.entry_time_delay = tk.Entry(frame_params, width=15)
        self.entry_time_delay.grid(row=0, column=3, padx=5, pady=2)
        self.entry_time_delay.insert(0, str(self.general_params["time_delay"]))
        tk.Label(frame_params, text="TOB:").grid(row=0, column=4, sticky="e")
        self.entry_TOB = tk.Entry(frame_params, width=15)
        self.entry_TOB.grid(row=0, column=5, padx=5, pady=2)
        self.entry_TOB.insert(0, str(self.general_params["TOB"]))
        tk.Label(frame_params, text="Intensity:").grid(row=1, column=0, sticky="e")
        self.entry_intensity = tk.Entry(frame_params, width=15)
        self.entry_intensity.grid(row=1, column=1, padx=5, pady=2)
        self.entry_intensity.insert(0, str(self.general_params["intensity"]))
        
        frame_detectors = tk.LabelFrame(self, text="Detektory")
        frame_detectors.pack(fill="both", expand=True, padx=10, pady=5)
        self.detectors_listbox = tk.Listbox(frame_detectors)
        self.detectors_listbox.pack(side="left", fill="both", expand=True, padx=5, pady=5)
        self.detectors_listbox.bind("<Double-Button-1>", self.edit_detector)
        scrollbar = tk.Scrollbar(frame_detectors, orient="vertical", command=self.detectors_listbox.yview)
        scrollbar.pack(side="right", fill="y")
        self.detectors_listbox.config(yscrollcommand=scrollbar.set)
        
        frame_buttons = tk.Frame(self)
        frame_buttons.pack(fill="x", padx=10, pady=5)
        tk.Button(frame_buttons, text="Dodaj detektor", command=self.add_detector).pack(side="left", padx=5)
        tk.Button(frame_buttons, text="Edytuj detektor", command=self.edit_detector).pack(side="left", padx=5)
        tk.Button(frame_buttons, text="Usuń detektor", command=self.remove_detector).pack(side="left", padx=5)
        tk.Button(frame_buttons, text="Uruchom analizę", command=self.run_analysis).pack(side="right", padx=5)
        
        frame_results = tk.LabelFrame(self, text="Wyniki analizy")
        frame_results.pack(fill="both", expand=True, padx=10, pady=5)
        self.text_results = tk.Text(frame_results, wrap="word")
        self.text_results.pack(fill="both", expand=True, padx=5, pady=5)
    
    def create_menu(self):
        menubar = tk.Menu(self)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Zapisz projekt", command=self.save_project)
        filemenu.add_command(label="Wczytaj projekt", command=self.load_project)
        filemenu.add_separator()
        filemenu.add_command(label="Wyjście", command=self.quit)
        menubar.add_cascade(label="Plik", menu=filemenu)
        self.config(menu=menubar)
    
    def add_detector(self):
        dialog = DetectorDialog(self)
        self.wait_window(dialog)
        if dialog.result:
            self.detectors.append(dialog.result)
            self.update_detectors_listbox()
    
    def edit_detector(self, event=None):
        selection = self.detectors_listbox.curselection()
        if not selection:
            messagebox.showinfo("Informacja", "Wybierz detektor do edycji")
            return
        index = selection[0]
        detector = self.detectors[index]
        dialog = DetectorDialog(self, detector)
        self.wait_window(dialog)
        if dialog.result:
            self.detectors[index] = dialog.result
            self.update_detectors_listbox()
    
    def remove_detector(self):
        selection = self.detectors_listbox.curselection()
        if not selection:
            messagebox.showinfo("Informacja", "Wybierz detektor do usunięcia")
            return
        index = selection[0]
        del self.detectors[index]
        self.update_detectors_listbox()
    
    def update_detectors_listbox(self):
        self.detectors_listbox.delete(0, tk.END)
        for det in self.detectors:
            self.detectors_listbox.insert(tk.END, det["name"])
    
    def run_analysis(self):
        try:
            self.general_params["lam"] = float(self.entry_lam.get())
            self.general_params["time_delay"] = float(self.entry_time_delay.get())
            self.general_params["TOB"] = float(self.entry_TOB.get())
            self.general_params["intensity"] = float(self.entry_intensity.get())
        except Exception as e:
            messagebox.showerror("Błąd", f"Błąd przy wczytywaniu parametrów ogólnych: {e}")
            return
        
        if not self.detectors:
            messagebox.showinfo("Informacja", "Dodaj przynajmniej jeden detektor.")
            return
        
        self.text_results.delete("1.0", tk.END)
        self.text_results.insert(tk.END, "Trwa analiza...\n")
        self.update_idletasks()
        
        results = run_analysis(self.general_params, self.detectors)
        if "errors" in results and results["errors"]:
            self.text_results.insert(tk.END, "Wystąpiły następujące błędy:\n")
            for err in results["errors"]:
                self.text_results.insert(tk.END, f"  - {err}\n")
        
        # Jeśli analiza NNLS nie została wykonana, wyświetlamy przynajmniej dane eksperymentalne (jeśli są dostępne)
        if "x_nnls" not in results:
            self.text_results.insert(tk.END, "\nAnaliza nie mogła być wykonana z powodu brakujących lub nieprawidłowych danych.\n")
            if "y" in results and results["y"]:
                y_all = np.array(results["y"])
                plt.figure(figsize=(8,4))
                plt.plot(y_all, 'o-', label="Dane eksperymentalne")
                plt.title("Dane eksperymentalne")
                plt.xlabel("Numer pomiaru")
                plt.ylabel("Zliczenia")
                plt.legend()
                plt.grid(True)
                plt.tight_layout()
                plt.show()
            return
        
        self.text_results.insert(tk.END, f"\nRozmiar A: {results['A_wysokie_shape']}\n")
        self.text_results.insert(tk.END, f"Rozmiar y: {results['y_shape']}\n")
        self.text_results.insert(tk.END, "Parametry (NNLS):\n")
        self.text_results.insert(tk.END, f"{results['x_nnls']}\n")
        self.text_results.insert(tk.END, "Błędy parametrów:\n")
        self.text_results.insert(tk.END, f"{results['param_errors']}\n")
        self.text_results.insert(tk.END, f"Chi2: {results['chi2']}\n\n")
        for i, met in enumerate(results["detector_metrics"]):
            self.text_results.insert(tk.END, f"Detektor {results['detector_names'][i]}:\n")
            self.text_results.insert(tk.END, f"  aFactor:\t\t\t\t\t\t\t{met['aFactor']:.4e}\n")
            self.text_results.insert(tk.END, f"  RMSE -  Root Mean Squared Error:\t\t\t\t\t\t\t{met['RMSE']:.4e}\n")
            self.text_results.insert(tk.END, f"  MAE - Mean Absolute Error:\t\t\t\t\t\t\t{met['MAE']:.4e}\n")
            self.text_results.insert(tk.END, f"  R2 - Współczynnik determinacji:\t\t\t\t\t\t\t{met['R2']:.4e}\n")
            self.text_results.insert(tk.END, f"  Pearson - Współczynnik korelacji Pearsona:\t\t\t\t\t\t\t{met['Pearson']:.4e}\n\n")
        self.text_results.insert(tk.END, "Oszacowane aktywności (dla poszczególnych izotopów):\n")
        isotopy = [iso.strip() for iso in self.detectors[0]["isotopy"].split(",")]
        for j in range(len(isotopy)):
            for k, (act, err) in enumerate(zip(results["x_nnls"][16*j:16*(j+1)], results["param_errors"][16*j:16*(j+1)])):
                self.text_results.insert(tk.END, f"Źródło {k+1} - {isotopy[j]}: {act:.2e} ± {err:.2e}\n")
        
        # Funkcja pomocnicza, która dobiera układ subplotów
        def get_subplot_grid(n):
            if n == 1:
                return (1, 1)
            elif n == 2:
                return (1, 2)
            elif n <= 4:
                return (2, 2)
            else:
                return (2, 3)

        # Przekształcamy dane wejściowe do tablic numpy
        y_all = np.array(results["y"])
        y_all_errors = np.array(results["y_unc"])
        y_est_all = np.array(results["y_est"])
        y_est_all_errors = np.array(results["y_est_errors"])
        index_start = 0
        num_detectors = len(results["detector_names"])

        # Przygotowanie figur i subplotów dla wykresów "dane vs. model"
        rows, cols = get_subplot_grid(num_detectors)
        fig1, axs1 = plt.subplots(rows, cols, figsize=(8 * cols, 4 * rows))
        axs1 = np.array(axs1).flatten()  # spłaszczamy tablicę osi

        # Przygotowanie figur i subplotów dla wykresów scatter
        fig2, axs2 = plt.subplots(rows, cols, figsize=(6 * cols, 6 * rows))
        axs2 = np.array(axs2).flatten()

        # Przygotowanie subplotów dla macierzy wydajności (tylko dla detektorów, które mają schemat)
        detectors_with_scheme = [i for i in range(num_detectors) if results["detector_has_scheme"][i]]
        num_scheme = len(detectors_with_scheme)
        if num_scheme > 0:
            rows_scheme, cols_scheme = get_subplot_grid(num_scheme)
            fig3, axs3 = plt.subplots(rows_scheme, cols_scheme, figsize=(8 * cols_scheme, 4 * rows_scheme))
            axs3 = np.array(axs3).flatten()

        # Indeks dla subplotów macierzy wydajności
        scheme_plot_index = 0

        # Pętla po wszystkich detektorach
        for i in range(num_detectors):
            m = results["detector_segments"][i]
            i1 = index_start
            i2 = index_start + m
            y_det = y_all[i1:i2]
            y_errors = y_all_errors[i1:i2]
            y_est_det = y_est_all[i1:i2]
            y_est_errors = y_est_all_errors[i1:i2]
            index_start = i2

            # Wykres pomiarów vs. model
            ax1 = axs1[i]
            ax1.errorbar(range(m), y_det, yerr=y_errors, fmt='o', color="blue", 
                        label="Pomiar", alpha=0.6)
            ax1.plot(range(m), y_est_det, '-r', label="Model (A·x)", alpha=0.8)
            ax1.fill_between(range(m),
                            np.array(y_est_det) - 3 * np.array(y_est_errors),
                            np.array(y_est_det) + 3 * np.array(y_est_errors),
                            color='red', alpha=0.3, label='Korytarz błędu modelu')
            ax1.set_title(f"Detektor: {results['detector_names'][i]} - dane vs. model")
            ax1.set_xlabel("Numer pomiaru")
            ax1.set_ylabel("Skorygowane zliczenia")
            ax1.legend()
            ax1.grid(True)

            # Wykres scatter: model vs. pomiary
            ax2 = axs2[i]
            ax2.scatter(y_det, y_est_det, c='green', alpha=0.7)
            ax2.plot([min(y_det), max(y_det)], [min(y_det), max(y_det)], 'k--', 
                    label="Idealne dopasowanie")
            ax2.set_title(f"Detektor: {results['detector_names'][i]} - Scatter: pomiary vs. model")
            ax2.set_xlabel("Dane eksperymentalne")
            ax2.set_ylabel("Dane oszacowane")
            ax2.legend()
            ax2.grid(True)

            # Wykres macierzy wydajności, jeśli dostępny
            if results["detector_has_scheme"][i]:
                ax3 = axs3[scheme_plot_index]
                A_block = results["A_blocks"][i]
                im = ax3.imshow(A_block, aspect='auto', interpolation='nearest')
                fig3.colorbar(im, ax=ax3)
                ax3.set_title(f"Macierz wydajności - {results['detector_names'][i]}")
                ax3.set_xlabel("Izotopy")
                ax3.set_ylabel("Pomiary")
                scheme_plot_index += 1

        # Po pętli, dostosowanie układu subplotów dla każdej figury
        fig1.tight_layout()
        fig2.tight_layout()
        if num_scheme > 0:
            fig3.tight_layout()

        plt.show()
    
    def save_project(self):
        project = {
            "general_params": self.general_params,
            "detectors": self.detectors
        }
        filename = filedialog.asksaveasfilename(title="Zapisz projekt", defaultextension=".json",
                                                filetypes=[("JSON", "*.json")])
        if filename:
            try:
                with open(filename, "w") as f:
                    json.dump(project, f, indent=4)
                messagebox.showinfo("Sukces", "Projekt zapisany.")
            except Exception as e:
                messagebox.showerror("Błąd", f"Błąd przy zapisie projektu: {e}")
    
    def load_project(self):
        filename = filedialog.askopenfilename(title="Wczytaj projekt", filetypes=[("JSON", "*.json")])
        if filename:
            try:
                with open(filename, "r") as f:
                    project = json.load(f)
                self.general_params = project.get("general_params", self.general_params)
                self.detectors = project.get("detectors", [])
                self.entry_lam.delete(0, tk.END)
                self.entry_lam.insert(0, str(self.general_params["lam"]))
                self.entry_time_delay.delete(0, tk.END)
                self.entry_time_delay.insert(0, str(self.general_params["time_delay"]))
                self.entry_TOB.delete(0, tk.END)
                self.entry_TOB.insert(0, str(self.general_params["TOB"]))
                self.entry_intensity.delete(0, tk.END)
                self.entry_intensity.insert(0, str(self.general_params["intensity"]))
                self.update_detectors_listbox()
                messagebox.showinfo("Sukces", "Projekt wczytany.")
            except Exception as e:
                messagebox.showerror("Błąd", f"Błąd przy wczytywaniu projektu: {e}")

if __name__ == "__main__":
    app = MainApplication()
    app.mainloop()
