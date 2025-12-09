from tkinter import Toplevel, Listbox, Text, Menu, END, Tk, BooleanVar, StringVar
from tkinter.ttk import Frame, Label, Entry, Button, Checkbutton, Scrollbar, LabelFrame, Combobox
from os import path, getcwd, makedirs
from re import sub
from datetime import datetime
from tkinter import filedialog, messagebox
from json import dump, load, dumps
from matplotlib.pyplot import figure, subplots, tight_layout, show, plot, title, xlabel, ylabel, legend, grid
from numpy import array, sqrt, arange, stack, sum, ndarray, floating, integer, bool_
from analysis import analysis, nuclear_data


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


# bezpieczne parsowanie numeryczne:
def _safe_float(s, default):
    try:
        return float(str(s).replace(",", "."))
    except Exception:
        return default


def _safe_int(s, default):
    try:
        return int(float(str(s).replace(",", ".")))
    except Exception:
        return default


# =============================================================================
# INTERFEJS GRAFICZNY – DIALOG DLA DETEKTORA
# =============================================================================
class DetectorDialog(Toplevel):
    def __init__(self, master, detector=None):
        super().__init__(master)
        self.title("Dodaj/Edytuj detektor")
        self.geometry("600x350")
        self.resizable(False, False)
        self.detector = detector
        self.result = None
        self.create_widgets()
        if detector:
            self.populate_fields(detector)
        else:
            self.entry_isotopy.insert(0, "C11")

    def create_widgets(self):
        frm = Frame(self, padding=20)
        frm.pack(fill="both", expand=True)

        # Define common options
        label_opts = {"anchor": "e", "padding": (5, 5)}
        entry_opts = {"width": 40}
        btn_opts = {"padding": (2, 2)}

        # Row 0: Nazwa
        Label(frm, text="Nazwa:", **label_opts).grid(row=0, column=0, sticky="e")
        self.entry_name = Entry(frm, **entry_opts)
        self.entry_name.grid(row=0, column=1, columnspan=2, sticky="w", padx=5, pady=5)

        # Row 1: Plik A (wydajność)
        Label(frm, text="Plik A (wydajność):", **label_opts).grid(row=1, column=0, sticky="e")
        self.entry_A_file = Entry(frm, **entry_opts)
        self.entry_A_file.grid(row=1, column=1, sticky="w", padx=5, pady=5)
        Button(frm, text="Przeglądaj", command=self.browse_A_file, **btn_opts).grid(row=1, column=2, padx=5, pady=5)

        # Row 2: Plik zliczeń
        Label(frm, text="Plik zliczeń:", **label_opts).grid(row=2, column=0, sticky="e")
        self.entry_counts_file = Entry(frm, **entry_opts)
        self.entry_counts_file.grid(row=2, column=1, sticky="w", padx=5, pady=5)
        Button(frm, text="Przeglądaj", command=self.browse_counts_file, **btn_opts).grid(row=2, column=2, padx=5, pady=5)

        # Row 3: Liczba pomiarów
        Label(frm, text="Liczba pomiarów:", **label_opts).grid(row=3, column=0, sticky="e")
        self.entry_n_meas = Entry(frm, **entry_opts)
        self.entry_n_meas.grid(row=3, column=1, sticky="w", padx=5, pady=5)

        # Row 4: Izotopy
        Label(frm, text="Izotopy (przecinek):", **label_opts).grid(row=4, column=0, sticky="e")
        self.entry_isotopy = Entry(frm, **entry_opts)
        self.entry_isotopy.grid(row=4, column=1, columnspan=2, sticky="w", padx=5, pady=5)

        # Row 5: Plik schematu rotacji
        Label(frm, text="Plik schematu rotacji:", **label_opts).grid(row=5, column=0, sticky="e")
        self.entry_rot_scheme = Entry(frm, **entry_opts)
        self.entry_rot_scheme.grid(row=5, column=1, sticky="w", padx=5, pady=5)
        Button(frm, text="Przeglądaj", command=self.browse_rot_scheme, **btn_opts).grid(row=5, column=2, padx=5, pady=5)

        # Row 6: Plik kroków rotacji
        Label(frm, text="Plik kroków rotacji:", **label_opts).grid(row=6, column=0, sticky="e")
        self.entry_steps_scheme = Entry(frm, **entry_opts)
        self.entry_steps_scheme.grid(row=6, column=1, sticky="w", padx=5, pady=5)
        Button(frm, text="Przeglądaj", command=self.browse_steps_scheme, **btn_opts).grid(row=6, column=2, padx=5, pady=5)

        # Row 7: Przyciski OK / Anuluj
        btn_frame = Frame(frm)
        btn_frame.grid(row=7, column=0, columnspan=3, pady=15)
        Button(btn_frame, text="OK", command=self.on_ok, width=12).grid(row=0, column=0, padx=10)
        Button(btn_frame, text="Anuluj", command=self.destroy, width=12).grid(row=0, column=1, padx=10)

    def browse_A_file(self):
        filename = filedialog.askopenfilename(title="Wybierz plik A",
                                              filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if filename:
            self.entry_A_file.delete(0, END)
            self.entry_A_file.insert(0, filename)

    def browse_counts_file(self):
        filename = filedialog.askopenfilename(title="Wybierz plik zliczeń",
                                              filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if filename:
            self.entry_counts_file.delete(0, END)
            self.entry_counts_file.insert(0, filename)

    def browse_rot_scheme(self):
        filename = filedialog.askopenfilename(title="Wybierz plik schematu rotacji",
                                              filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if filename:
            self.entry_rot_scheme.delete(0, END)
            self.entry_rot_scheme.insert(0, filename)

    def browse_steps_scheme(self):
        filename = filedialog.askopenfilename(title="Wybierz plik kroków rotacji",
                                              filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if filename:
            self.entry_steps_scheme.delete(0, END)
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


class MainApplication(Tk):
    def __init__(self):
        super().__init__()
        self.title("Analiza danych eksperymentalnych")
        self.geometry("1400x1000")
        self.general_params = {
            "lam": 0.0005672978299613248,
            "time_delay": 479.4,
            "TOB": 332.1,
            "intensity": 1.9983,
            "a": 5
        }
        self.detectors = []  # lista detektorów

        self.fit_method = StringVar(value=self.general_params.get("fit_method", "NNLS"))
        self.mlem_variant = StringVar(value=self.general_params.get("mlem_variant", "CLASSIC"))
        self.support_mask_sources_var = StringVar(value=self.general_params.get("support_mask_sources", ""))

        self.alpha0_var = StringVar(value=str(self.general_params.get("alpha0", 1.0)))
        self.tau_var = StringVar(value=str(self.general_params.get("tau", 25.0)))
        self.subsets_var = StringVar(value=str(self.general_params.get("subsets", 8)))
        self.cov_cap_factor_var = StringVar(value=str(self.general_params.get("cov_cap_factor", 4.0)))


        self.create_widgets()
        self.create_menu()


    def create_widgets(self):
        # Frame for general parameters with a grid layout for consistent alignment
        frame_params = LabelFrame(self, text="Parametry ogólne")
        frame_params.pack(fill="x", padx=10, pady=5)

        Label(frame_params, text="Time delay:").grid(row=0, column=0, sticky="e", padx=5, pady=2)
        self.entry_time_delay = Entry(frame_params, width=15)
        self.entry_time_delay.grid(row=0, column=1, padx=5, pady=2)
        self.entry_time_delay.insert(0, str(self.general_params["time_delay"]))

        Label(frame_params, text="TOB:").grid(row=0, column=2, sticky="e", padx=5, pady=2)
        self.entry_TOB = Entry(frame_params, width=15)
        self.entry_TOB.grid(row=0, column=3, padx=5, pady=2)
        self.entry_TOB.insert(0, str(self.general_params["TOB"]))

        self.exp_analysis = BooleanVar(value=self.general_params.get("exp_analysis", True))
        Checkbutton(frame_params, text="Dodatkowa analiza za pomocą eksponensów",
                        variable=self.exp_analysis).grid(row=0, column=4, padx=5, pady=2)
        
        # Frame for detectors list with scrollbar
        frame_detectors = LabelFrame(self, text="Detektory")
        frame_detectors.pack(fill="both", expand=True, padx=10, pady=5)

        self.detectors_listbox = Listbox(frame_detectors)
        self.detectors_listbox.pack(side="left", fill="both", expand=True, padx=5, pady=5)
        self.detectors_listbox.bind("<Double-Button-1>", self.edit_detector)
        scrollbar = Scrollbar(frame_detectors, orient="vertical", command=self.detectors_listbox.yview)
        scrollbar.pack(side="right", fill="y", padx=5, pady=5)
        self.detectors_listbox.config(yscrollcommand=scrollbar.set)
        
        # Frame for action buttons
        frame_buttons = Frame(self)
        frame_buttons.pack(fill="x", padx=10, pady=5)
        Button(frame_buttons, text="Dodaj detektor", command=self.add_detector).pack(side="left", padx=5)
        Button(frame_buttons, text="Edytuj detektor", command=self.edit_detector).pack(side="left", padx=5)
        Button(frame_buttons, text="Usuń detektor", command=self.remove_detector).pack(side="left", padx=5)
        Button(frame_buttons, text="Uruchom analizę", command=self.run_analysis).pack(side="right", padx=5)
        
        # Frame for displaying analysis results
        frame_results = LabelFrame(self, text="Wyniki analizy")
        frame_results.pack(fill="both", expand=True, padx=10, pady=5)
        scrollbar = Scrollbar(frame_results, orient="vertical")
        scrollbar.pack(side="right", fill="y")
        self.text_results = Text(frame_results, wrap="word", yscrollcommand=scrollbar.set)
        self.text_results.pack(fill="both", expand=True, padx=5, pady=5)
        scrollbar.config(command=self.text_results.yview)

                # --- wybór metody ---
        Label(frame_params, text="Metoda dopasowania:").grid(row=1, column=0, sticky="e", padx=5, pady=2)
        self.combo_fit = Combobox(frame_params, values=["NNLS", "MLEM"], textvariable=self.fit_method, width=12, state="readonly")
        self.combo_fit.grid(row=1, column=1, padx=5, pady=2)

        Label(frame_params, text="Wariant MLEM:").grid(row=1, column=2, sticky="e", padx=5, pady=2)
        self.combo_mlem = Combobox(frame_params, values=["CLASSIC", "R-OS-SPS"], textvariable=self.mlem_variant, width=12, state="readonly")
        self.combo_mlem.grid(row=1, column=3, padx=5, pady=2)

        # --- maska źródeł 0/1 (po slotach) ---
        Label(frame_params, text="Maska źródeł 0/1 (sloty):").grid(row=2, column=0, sticky="e", padx=5, pady=2)
        self.entry_support_mask = Entry(frame_params, width=50, textvariable=self.support_mask_sources_var)
        self.entry_support_mask.grid(row=2, column=1, columnspan=3, sticky="w", padx=5, pady=2)

        # --- parametry R-OS-SPS ---
        Label(frame_params, text="alpha0:").grid(row=3, column=0, sticky="e", padx=5, pady=2)
        self.alpha0_entry = Entry(frame_params, width=10, textvariable=self.alpha0_var)
        self.alpha0_entry.grid(row=3, column=1, sticky="w", padx=5, pady=2)

        Label(frame_params, text="tau:").grid(row=3, column=2, sticky="e", padx=5, pady=2)
        self.tau_entry = Entry(frame_params, width=10, textvariable=self.tau_var)
        self.tau_entry.grid(row=3, column=3, sticky="w", padx=5, pady=2)

        Label(frame_params, text="subsets:").grid(row=4, column=0, sticky="e", padx=5, pady=2)
        self.subsets_entry = Entry(frame_params, width=10, textvariable=self.subsets_var)
        self.subsets_entry.grid(row=4, column=1, sticky="w", padx=5, pady=2)

        # --- cap factor parameter for covariance estimation ---
        cap_factor_label = Label(frame_params, text="Cap factor:")
        cap_factor_label.grid(row=4, column=2, sticky="e", padx=5, pady=2)
    
        # Add StringVar for cap_factor
        self.cov_cap_factor_var = StringVar(value=str(self.general_params.get("cov_cap_factor", 4.0)))
        self.cov_cap_factor_entry = Entry(frame_params, width=10, textvariable=self.cov_cap_factor_var)
        self.cov_cap_factor_entry.grid(row=4, column=3, sticky="w", padx=5, pady=2)
        
        # Add tooltip functionality
        def create_tooltip(widget, text):
            def on_enter(event):
                tooltip = Toplevel()
                tooltip.wm_overrideredirect(True)
                tooltip.wm_geometry(f"+{event.x_root+10}+{event.y_root+10}")
                label = Label(tooltip, text=text, background="lightyellow", relief="solid", borderwidth=1, wraplength=300)
                label.pack()
                widget.tooltip = tooltip

            def on_leave(event):
                if hasattr(widget, 'tooltip'):
                    widget.tooltip.destroy()
                    del widget.tooltip

            widget.bind("<Enter>", on_enter)
            widget.bind("<Leave>", on_leave)

        # Add tooltip to cap factor label
        create_tooltip(cap_factor_label, "Maksymalna dopuszczalna kwadratowa reszta Pearsona - ogranicza wpływ pojedynczych obserwacji na oszacowanie kowariancji w metodach MLEM, zapobiegając niestabilności numerycznej przy dużych resztach.")

        # --- włączanie/wyłączanie kontrolek wariantu i parametrów zależnie od metody ---
        def _toggle_mlem_fields(*_):
            is_mlem = self.fit_method.get().upper() == "MLEM"
            variant = self.mlem_variant.get().upper()
            self.combo_mlem.configure(state=("readonly" if is_mlem else "disabled"))
            # alpha0/tau/subsets tylko dla R-OS-SPS:
            state_params = ("normal" if (is_mlem and variant == "R-OS-SPS") else "disabled")
            self.alpha0_entry.configure(state=state_params)
            self.tau_entry.configure(state=state_params)
            self.subsets_entry.configure(state=state_params)

        _toggle_mlem_fields()
        self.fit_method.trace_add("write", _toggle_mlem_fields)
        self.mlem_variant.trace_add("write", _toggle_mlem_fields)
            
    def create_menu(self):
        menubar = Menu(self)
        filemenu = Menu(menubar, tearoff=0)
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
        self.detectors_listbox.delete(0, END)
        for det in self.detectors:
            self.detectors_listbox.insert(END, det["name"])

    def show_alternative_analysis(self, alternatives, x_nnls, method_label):
        num_detectors = len(alternatives)
        rows, cols = get_subplot_grid(num_detectors)
        fig, axs = subplots(rows, cols, figsize=(8 * cols, 4 * rows))
        axs = array(axs).flatten()
        if num_detectors == 1:
            axs = [axs]

        # Lista izotopów – zakładamy, że jest taka sama dla wszystkich detektorów.
        isotopy = [iso.strip() for iso in self.detectors[0]["isotopy"].split(",")]
        n_iso = len(isotopy)

        for i, (det_name, alt_A0, alt_err) in enumerate(alternatives):
            ax = axs[i]
            # Rysujemy alternatywną analizę – mniejsze znaczenie (alpha = 0.6)
            for j, iso in enumerate(isotopy):
                sources = arange(1, len(alt_A0[j]) + 1)
                ax.errorbar(sources, alt_A0[j], yerr=alt_err[j], fmt='o', capsize=5,
                            label=f"Alt {iso}", alpha=0.6)
            # Nałożenie wyników NNLS – powinny być główne (gruba linia, większe markery)
            for j, iso in enumerate(isotopy):
                sources = arange(1, 17)  # 16 źródeł (slotów)
                nnls_vals = array(x_nnls[16 * j:16 * (j + 1)])
                nnls_scaled = nnls_vals * nuclear_data[iso]
                ax.plot(sources, nnls_scaled, 's-', linewidth=2.5,
                        markersize=8, label=f"{method_label} {iso}")
            ax.set_title(f"Alternative Exponential Fit - {det_name}")
            ax.set_xlabel("Źródło (slot)")
            ax.set_ylabel("Rekonstruowana aktywność (Bq)")
            ax.set_yscale('log')
            ax.legend(fontsize=9)
            ax.grid(True)

        alt_arr = stack([alt_A0 for (_, alt_A0, _) in alternatives], axis=0)
        err_arr = stack([alt_err for (_, _, alt_err) in alternatives], axis=0)
        # wagi i średnia ważona
        w = 1.0 / (err_arr**2)
        W_sum = sum(w, axis=0)              # kształt (n_iso, n_sources)
        avg_alt = sum(alt_arr * w, axis=0) / W_sum
        sigma_avg = sqrt(1.0 / W_sum)

        # wykres agregowany – jeden subplot na izotop
        fig_aggr, axs_aggr = subplots(n_iso, 1, figsize=(8, 3 * n_iso))
        if n_iso == 1:
            axs_aggr = [axs_aggr]
        for j, iso in enumerate(isotopy):
            ax = axs_aggr[j]
            sources = arange(1, avg_alt.shape[1] + 1)
            # alternatywna średnia
            ax.errorbar(sources, avg_alt[j], yerr=sigma_avg[j],
                        fmt='o-', capsize=5, color='blue', label='Alt avg')
            # NNLS (skalowane)
            nnls_vals = array(x_nnls[16 * j:16 * (j + 1)]) * nuclear_data[iso]
            ax.plot(sources, nnls_vals, 's--', color='red', 
                    markersize=6, label=method_label)
            ax.set_title(f"Izotop {iso} – średnia alt vs NNLS")
            ax.set_xlabel("Źródło (slot)")
            ax.set_ylabel("Aktywność [Bq]")
            ax.legend()
            ax.grid(True)
        tight_layout()
        show()

    def run_analysis(self):
        try:
            self.general_params["time_delay"] = float(self.entry_time_delay.get())
            self.general_params["TOB"] = float(self.entry_TOB.get())
            self.general_params["exp_analysis"] = self.exp_analysis.get()
            # zbierz ustawienia z GUI
            self.general_params["fit_method"] = self.fit_method.get()
            self.general_params["mlem_variant"] = self.mlem_variant.get()
            self.general_params["support_mask_sources"] = self.support_mask_sources_var.get().strip()

            self.general_params["alpha0"] = _safe_float(self.alpha0_var.get(), 1.0)
            self.general_params["tau"]     = _safe_float(self.tau_var.get(), 25.0)
            self.general_params["subsets"] = _safe_int(self.subsets_var.get(), 8)
            self.general_params["cov_cap_factor"] = _safe_float(self.cov_cap_factor_var.get(), 4.0)
        except Exception as e:
            messagebox.showerror("Błąd", f"Błąd przy wczytywaniu parametrów ogólnych: {e}")
            return

        if not self.detectors:
            messagebox.showinfo("Informacja", "Dodaj przynajmniej jeden detektor.")
            return

        self.text_results.delete("1.0", END)
        self.text_results.insert(END, "Trwa analiza...\n")
        self.update_idletasks()

        results = analysis(self.general_params, self.detectors)
        if "errors" in results and results["errors"]:
            self.text_results.insert(END, "Wystąpiły następujące błędy:\n")
            for err in results["errors"]:
                self.text_results.insert(END, f"  - {err}\n")

        # Układ subplotów dla wykresów "dane vs. model"
        num_detectors = len(results["detector_names"])
        rows, cols = get_subplot_grid(num_detectors)
        fig1, axs1 = subplots(rows, cols, figsize=(8 * cols, 4 * rows))
        axs1 = array(axs1).flatten()

        # Podobnie dla wykresów scatter: pomiary vs. model
        fig2, axs2 = subplots(rows, cols, figsize=(6 * cols, 6 * rows))
        axs2 = array(axs2).flatten()

        # Jeśli dostępne – subploty dla macierzy wydajności
        detectors_with_scheme = [i for i in range(num_detectors) if results["detector_has_scheme"][i]]
        num_scheme = len(detectors_with_scheme)
        if num_scheme > 0:
            rows_scheme, cols_scheme = get_subplot_grid(num_scheme)
            fig3, axs3 = subplots(rows_scheme, cols_scheme, figsize=(8 * cols_scheme, 4 * rows_scheme))
            axs3 = list(array(axs3).flatten())
        else:
            fig3 = None

        # Jeśli NNLS nie zostało wykonane, prezentujemy dane eksperymentalne
        if "x_est" not in results:
            self.text_results.insert(END, "\nAnaliza nie mogła być wykonana z powodu brakujących lub nieprawidłowych danych.\n")
            if "y" in results and results["y"]:
                figure(figsize=(8, 4))
                plot(results["y"], 'o-', label="Dane eksperymentalne")
                title("Dane eksperymentalne")
                xlabel("Numer pomiaru")
                ylabel("Zliczenia")
                legend()
                grid(True)
                tight_layout()
                show()
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
        self.text_results.insert(END, f"\nRozmiar A: {results['A_wysokie_shape']}\n")
        self.text_results.insert(END, f"Rozmiar y: {results['y_shape']}\n")
        self.text_results.insert(END, f"Metoda: {results.get('method','')}\n")
        self.text_results.insert(END, "Parametry (x):\n")
        self.text_results.insert(END, f"{results.get('x_est', [])}\n")

        param_errs = results.get("param_errors", [])
        if param_errs:
            self.text_results.insert(END, "Błędy parametrów:\n")
            self.text_results.insert(END, f"{param_errs}\n")
        else:
            self.text_results.insert(END, "Błędy parametrów: — (niedostępne dla wybranej metody)\n")

        if "chi2" in results:
            self.text_results.insert(END, f"Chi2: {results['chi2']}\n\n")
        for i, met in enumerate(results["detector_metrics"]):
            self.text_results.insert(END, f"Detektor {results['detector_names'][i]}:\n")
            self.text_results.insert(END, f"  aFactor:\t\t{met['aFactor']:.4e}\n")
            self.text_results.insert(END, f"  RMSE:\t\t{met['RMSE']:.4e}\n")
            self.text_results.insert(END, f"  MAE:\t\t{met['MAE']:.4e}\n")
            self.text_results.insert(END, f"  R2:\t\t{met['R2']:.4e}\n")
            self.text_results.insert(END, f"  Pearson:\t{met['Pearson']:.4e}\n\n")

        # Wyświetlenie aktywności dla poszczególnych izotopów
        isotopy = [iso.strip() for iso in self.detectors[0]["isotopy"].split(",")]
        x_vec = array(results.get("x_est", []))
        if x_vec.size and len(isotopy):
            n_sources = x_vec.size // len(isotopy)
            self.text_results.insert(END, "Oszacowane aktywności (dla poszczególnych izotopów):\n")
            for j, iso in enumerate(isotopy):
                for s in range(n_sources):
                    act = x_vec[n_sources*j + s]
                    err = array(param_errs)[n_sources*j + s] if param_errs else 0.0
                    self.text_results.insert(END, f"Źródło {s+1} - {iso}: {(act * nuclear_data[iso]):.2e} ± {(err * nuclear_data[iso]):.2e}\n")

        if self.general_params.get("exp_analysis", False):
            self.text_results.insert(END, "\nAlternatywna analiza eksponencjalna (A0 w Bq) dla kolejnych źródeł:\n")
            for (det_name, alt_A0, alt_err) in results["alt_results"]:
                for j, iso in enumerate(isotopy):
                    for s in range(16):
                        self.text_results.insert(END, f"  {det_name} - Izotop {iso}, Źródło {s+1}: {alt_A0[j][s]:.2e} ± {alt_err[j][s]:.2e} Bq\n")

        # Zapis raportu analizy do pliku tekstowego
        self._save_analysis_report(results)

        # Wykresy: dane vs. model oraz scatter
        index_start = 0
        for i in range(num_detectors):
            m_det = results["detector_segments"][i]
            i1 = index_start
            i2 = index_start + m_det
            index_start = i2
            y_det = array(results["y"])[i1:i2]
            y_err = array(results["y_unc"])[i1:i2]
            y_est = array(results["y_est"])[i1:i2]
            y_est_err = array(results["y_est_errors"])[i1:i2]

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
        if self.general_params.get("exp_analysis", False):
            self.show_alternative_analysis(results["alt_results"],
                                   results.get("x_est", []),
                                   results.get("method", ""))
        show()


    

    def _json_default(self, o):
        """Konwersja obiektów niewspieranych przez json na typy podstawowe."""
        # NumPy array -> lista
        if isinstance(o, ndarray):
            return o.tolist()
        # NumPy scalary -> natywne typy Pythona
        if isinstance(o, (floating,)):
            return float(o)
        if isinstance(o, (integer,)):
            return int(o)
        if isinstance(o, (bool_,)):
            return bool(o)
        # Path, datetime, inne egzotyczne -> str
        try:
            import pathlib, datetime as _dt
            if isinstance(o, (pathlib.Path,)):
                return str(o)
            if isinstance(o, (_dt.datetime, _dt.date, _dt.time)):
                return o.isoformat()
        except Exception:
            pass
        # fall-back
        return str(o)
    
    def _save_analysis_report(self, results: dict):
        """
        Zapisuje raport tekstowy (UTF-8) z kompletem danych do folderu logs/.
        Nazwa pliku: YYYYMMDD-HHMMSS-mmm_<METHOD>.txt
        """
        try:
            # 1) przygotuj folder
            log_dir = path.join(getcwd(), "logs")
            makedirs(log_dir, exist_ok=True)

            # 2) stempel czasu do milisekund
            now = datetime.now()
            ts = now.strftime("%Y%m%d-%H%M%S") + f"-{int(now.microsecond/1000):03d}"

            # 3) metoda w nazwie pliku (slug)
            method = str(results.get("method", "UNKNOWN"))
            method_slug = sub(r"[^A-Za-z0-9_-]+", "", method.replace(" ", "-"))

            # 4) ścieżka pliku
            fname = f"{ts}_{method_slug}.txt"
            fpath = path.join(log_dir, fname)

            # 5) zbuduj treść raportu (tekst + JSON z kompletem danych)
            header_lines = [
                f"# Raport analizy",
                f"Data/godzina: {now.isoformat(timespec='milliseconds')}",
                f"Metoda: {method}",
                "",
                "# Parametry ogólne (general_params):",
                dumps(self.general_params, ensure_ascii=False, indent=2, default=self._json_default),
                "",
                "# Konfiguracja detektorów (detectors):",
                dumps(self.detectors, ensure_ascii=False, indent=2, default=self._json_default),
                "",
                "# Wyniki (results):",
                dumps(results, ensure_ascii=False, indent=2, default=self._json_default),
                "",
            ]
            content = "\n".join(header_lines)

            # 6) zapis
            with open(fpath, "w", encoding="utf-8") as f:
                f.write(content)

            # 7) sygnał w UI (np. w polu wyników)
            if hasattr(self, "text_results"):
                self.text_results.insert("end", f"\n[LOG] Zapisano raport w {fpath}\n")

        except Exception as e:
            # nie blokuj dalszej pracy GUI
            if hasattr(self, "text_results"):
                self.text_results.insert("end", f"\n[LOG][BŁĄD] Nie udało się zapisać raportu: {e}\n")


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
                    dump(project, f, indent=4)
                messagebox.showinfo("Sukces", "Projekt zapisany.")
            except Exception as e:
                messagebox.showerror("Błąd", f"Błąd przy zapisie projektu: {e}")
    
    def load_project(self):
        filename = filedialog.askopenfilename(title="Wczytaj projekt", filetypes=[("JSON", "*.json")])
        if filename:
            try:
                with open(filename, "r") as f:
                    project = load(f)
                self.general_params = project.get("general_params", self.general_params)
                self.detectors = project.get("detectors", [])
                self.entry_time_delay.delete(0, END)
                self.entry_time_delay.insert(0, str(self.general_params["time_delay"]))
                self.entry_TOB.delete(0, END)
                self.entry_TOB.insert(0, str(self.general_params["TOB"]))
                self.exp_analysis.set(self.general_params.get("exp_analysis", True))
                self.fit_method.set(self.general_params.get("fit_method", "NNLS"))
                self.mlem_variant.set(self.general_params.get("mlem_variant", "CLASSIC"))
                self.support_mask_sources_var.set(self.general_params.get("support_mask_sources", ""))

                self.alpha0_var.set(str(self.general_params.get("alpha0", 1.0)))
                self.tau_var.set(str(self.general_params.get("tau", 25.0)))
                self.subsets_var.set(str(self.general_params.get("subsets", 8)))
                self.cov_cap_factor_var.set(str(self.general_params.get("cov_cap_factor", 4.0)))
                self.update_detectors_listbox()
                messagebox.showinfo("Sukces", "Projekt wczytany.")
            except Exception as e:
                messagebox.showerror("Błąd", f"Błąd przy wczytywaniu projektu: {e}")


if __name__ == "__main__":
    app = MainApplication()
    app.mainloop()
