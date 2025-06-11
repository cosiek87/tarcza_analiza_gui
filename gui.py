import tkinter as tk
from tkinter import filedialog, messagebox
import json
import matplotlib.pyplot as plt
import numpy as np
from analysis import run_analysis, nuclear_data


def get_subplot_grid(n):
    """Zwraca układ subplotów (wiersze, kolumny) dla n elementów."""
    if n == 1:
        return (1, 1)
    elif n == 2:
        return (1, 2)
    elif n <= 4:
        return (2, 2)
    else:
        return (2, 3)


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


    def show_alternative_analysis(self, alternatives, x_nnls):

        num_detectors = len(alternatives)
        fig, axs = plt.subplots(num_detectors, 1, figsize=(8, 4 * num_detectors))
        if num_detectors == 1:
            axs = [axs]

        # Lista izotopów – zakładamy, że jest taka sama dla wszystkich detektorów.
        isotopy = [iso.strip() for iso in self.detectors[0]["isotopy"].split(",")]
        n_iso = len(isotopy)

        for i, (det_name, alt_A0, alt_err) in enumerate(alternatives):
            ax = axs[i]
            # Rysujemy alternatywną analizę – mniejsze znaczenie (alpha = 0.6)
            for j, iso in enumerate(isotopy):
                sources = np.arange(1, len(alt_A0[j]) + 1)
                ax.errorbar(sources, alt_A0[j], yerr=alt_err[j], fmt='o', capsize=5,
                            label=f"Alt {iso}", alpha=0.6)
            # Nałożenie wyników NNLS – powinny być główne (gruba linia, większe markery)
            for j, iso in enumerate(isotopy):
                sources = np.arange(1, 17)  # 16 źródeł (slotów)
                nnls_vals = np.array(x_nnls[16 * j:16 * (j + 1)])
                nnls_scaled = nnls_vals * nuclear_data[iso]
                ax.plot(sources, nnls_scaled, 's-', color='red', linewidth=2.5,
                        markersize=8, label=f"NNLS {iso}")
            ax.set_title(f"Alternative Exponential Fit - {det_name}")
            ax.set_xlabel("Źródło (slot)")
            ax.set_ylabel("Rekonstruowana aktywność (Bq)")
            ax.legend(fontsize=9)
            ax.grid(True)

        alt_arr = np.stack([alt_A0 for (_, alt_A0, _) in alternatives], axis=0)
        err_arr = np.stack([alt_err for (_, _, alt_err) in alternatives], axis=0)
        # wagi i średnia ważona
        w = 1.0 / (err_arr**2)
        W_sum = np.sum(w, axis=0)              # kształt (n_iso, n_sources)
        avg_alt = np.sum(alt_arr * w, axis=0) / W_sum
        sigma_avg = np.sqrt(1.0 / W_sum)

        # wykres agregowany – jeden subplot na izotop
        fig_aggr, axs_aggr = plt.subplots(n_iso, 1, figsize=(8, 3*n_iso))
        if n_iso == 1:
            axs_aggr = [axs_aggr]
        for j, iso in enumerate(isotopy):
            ax = axs_aggr[j]
            sources = np.arange(1, avg_alt.shape[1] + 1)
            # alternatywna średnia
            ax.errorbar(sources, avg_alt[j], yerr=sigma_avg[j],
                        fmt='o-', capsize=5, color='blue', label='Alt avg')
            # NNLS (skalowane)
            nnls_vals = np.array(x_nnls[16*j:16*(j+1)]) * nuclear_data[iso]
            ax.plot(sources, nnls_vals, 's--', color='red', 
                    markersize=6, label='NNLS')
            ax.set_title(f"Izotop {iso} – średnia alt vs NNLS")
            ax.set_xlabel("Źródło (slot)")
            ax.set_ylabel("Aktywność [Bq]")
            ax.legend()
            ax.grid(True)
        plt.tight_layout()
        plt.show()
    
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

        # Układ subplotów dla wykresów "dane vs. model"
        num_detectors = len(results["detector_names"])
        rows, cols = get_subplot_grid(num_detectors)
        fig1, axs1 = plt.subplots(rows, cols, figsize=(8 * cols, 4 * rows))
        axs1 = np.array(axs1).flatten()

        # Podobnie dla wykresów scatter: pomiary vs. model
        fig2, axs2 = plt.subplots(rows, cols, figsize=(6 * cols, 6 * rows))
        axs2 = np.array(axs2).flatten()

        # Jeśli dostępne – subploty dla macierzy wydajności
        detectors_with_scheme = [i for i in range(num_detectors) if results["detector_has_scheme"][i]]
        num_scheme = len(detectors_with_scheme)
        if num_scheme > 0:
            rows_scheme, cols_scheme = get_subplot_grid(num_scheme)
            fig3, axs3 = plt.subplots(rows_scheme, cols_scheme, figsize=(8 * cols_scheme, 4 * rows_scheme))
            axs3 = list(np.array(axs3).flatten())
        else:
            fig3 = None

        # Jeśli NNLS nie zostało wykonane, prezentujemy dane eksperymentalne
        if "x_nnls" not in results:
            self.text_results.insert(tk.END, "\nAnaliza nie mogła być wykonana z powodu brakujących lub nieprawidłowych danych.\n")
            if "y" in results and results["y"]:
                plt.figure(figsize=(8, 4))
                plt.plot(results["y"], 'o-', label="Dane eksperymentalne")
                plt.title("Dane eksperymentalne")
                plt.xlabel("Numer pomiaru")
                plt.ylabel("Zliczenia")
                plt.legend()
                plt.grid(True)
                plt.tight_layout()
                plt.show()
            if results["A_blocks"] and results["detector_has_scheme"]:
                scheme_plot_index = 0
                for i in range(num_detectors):
                    if results["detector_has_scheme"][i]:
                        ax3 = axs3[scheme_plot_index]
                        A_block = results["A_blocks"][i]
                        im = ax3.imshow(A_block, aspect='auto', interpolation='nearest')
                        fig3.colorbar(im, ax=ax3)
                        ax3.set_title(f"Macierz wydajności - {results['detector_names'][i]}")
                        ax3.set_xlabel("Izotopy")
                        ax3.set_ylabel("Pomiary")
                        scheme_plot_index += 1
                fig3.tight_layout()
            return

        # Wypisanie wyników NNLS i metryk
        self.text_results.insert(tk.END, f"\nRozmiar A: {results['A_wysokie_shape']}\n")
        self.text_results.insert(tk.END, f"Rozmiar y: {results['y_shape']}\n")
        self.text_results.insert(tk.END, "Parametry (NNLS):\n")
        self.text_results.insert(tk.END, f"{results['x_nnls']}\n")
        self.text_results.insert(tk.END, "Błędy parametrów:\n")
        self.text_results.insert(tk.END, f"{results['param_errors']}\n")
        self.text_results.insert(tk.END, f"Chi2: {results['chi2']}\n\n")
        for i, met in enumerate(results["detector_metrics"]):
            self.text_results.insert(tk.END, f"Detektor {results['detector_names'][i]}:\n")
            self.text_results.insert(tk.END, f"  aFactor:\t\t{met['aFactor']:.4e}\n")
            self.text_results.insert(tk.END, f"  RMSE:\t\t{met['RMSE']:.4e}\n")
            self.text_results.insert(tk.END, f"  MAE:\t\t{met['MAE']:.4e}\n")
            self.text_results.insert(tk.END, f"  R2:\t\t{met['R2']:.4e}\n")
            self.text_results.insert(tk.END, f"  Pearson:\t{met['Pearson']:.4e}\n\n")

        # Wyświetlenie aktywności dla poszczególnych izotopów
        isotopy = [iso.strip() for iso in self.detectors[0]["isotopy"].split(",")]
        self.text_results.insert(tk.END, "Oszacowane aktywności (dla poszczególnych izotopów):\n")
        for j, iso in enumerate(isotopy):
            for k, (act, err) in enumerate(zip(results["x_nnls"][16*j:16*(j+1)], 
                                                results["param_errors"][16*j:16*(j+1)])):
                self.text_results.insert(tk.END, f"Źródło {k+1} - {iso}: {(act * nuclear_data[iso]):.2e} ± {(err * nuclear_data[iso]):.2e}\n")

        self.text_results.insert(tk.END, "\nAlternatywna analiza eksponencjalna (A0 w Bq) dla kolejnych źródeł:\n")
        for (det_name, alt_A0, alt_err) in results["alt_results"]:
            for j, iso in enumerate(isotopy):
                for s in range(16):
                    self.text_results.insert(tk.END, f"  {det_name} - Izotop {iso}, Źródło {s+1}: {alt_A0[j][s]:.2e} ± {alt_err[j][s]:.2e} Bq\n")

        # Wykresy: dane vs. model oraz scatter
        index_start = 0
        for i in range(num_detectors):
            m_det = results["detector_segments"][i]
            i1 = index_start
            i2 = index_start + m_det
            index_start = i2
            y_det = np.array(results["y"])[i1:i2]
            y_err = np.array(results["y_unc"])[i1:i2]
            y_est = np.array(results["y_est"])[i1:i2]
            y_est_err = np.array(results["y_est_errors"])[i1:i2]

            # Wykres "dane vs. model"
            ax1 = axs1[i]
            ax1.errorbar(range(m_det), y_det, yerr=y_err, fmt='o', color="blue", label="Pomiar", alpha=0.6)
            ax1.plot(range(m_det), y_est, '-r', label="Model (A·x)", alpha=0.8)
            ax1.fill_between(range(m_det),
                            y_est - 3 * y_est_err,
                            y_est + 3 * y_est_err,
                            color='red', alpha=0.3, label='Przedział błędu modelu')
            ax1.set_title(f"Detektor: {results['detector_names'][i]} - dane vs. model")
            ax1.set_xlabel("Numer pomiaru")
            ax1.set_ylabel("Skorygowane zliczenia")
            ax1.legend()
            ax1.grid(True)

            # Wykres scatter: pomiary vs. model
            ax2 = axs2[i]
            ax2.scatter(y_det, y_est, c='green', alpha=0.7)
            ax2.plot([y_det.min(), y_det.max()], [y_det.min(), y_det.max()], 'k--', label="Idealne dopasowanie")
            ax2.set_title(f"Detektor: {results['detector_names'][i]} - Scatter")
            ax2.set_xlabel("Dane eksperymentalne")
            ax2.set_ylabel("Dane oszacowane")
            ax2.legend()
            ax2.grid(True)

            # Macierz wydajności, jeśli dostępna
            if results["detector_has_scheme"][i] and fig3:
                ax3 = axs3.pop(0)
                A_block = results["A_blocks"][i]
                im = ax3.imshow(A_block, aspect='auto', interpolation='nearest')
                fig3.colorbar(im, ax=ax3)
                ax3.set_title(f"Macierz wydajności - {results['detector_names'][i]}")
                ax3.set_xlabel("Izotopy")
                ax3.set_ylabel("Pomiary")

        fig1.tight_layout()
        fig2.tight_layout()
        if fig3:
            fig3.tight_layout()

        self.show_alternative_analysis(results["alt_results"], results["x_nnls"])
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
