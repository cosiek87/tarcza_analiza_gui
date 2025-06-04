import numpy as np
from numpy.linalg import inv, pinv
from scipy.optimize import nnls
from utils import build_A_matrix, correct_counts, alternative_exponential_decay_fit, nuclear_data
from metrics import aFactor, compute_rmse, compute_mae, compute_r2, compute_pearson
from tkinter import messagebox

def get_isotopy_and_lam(all_data, det, general_params):
    """
    Ustala listę izotopów oraz odpowiadające im stałe rozpadu.
    Jeżeli definicja izotopów nie jest podana, traktuje A_file jako pojedynczy wektor.
    Zwraca: (isotopy, n_iso, lam_vector, all_data_1d)
    """
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
        # Jeśli all_data jest macierzą, wybieramy pierwszy wiersz
        if all_data.ndim > 1:
            all_data = all_data[0]
    else:
        # Brak definicji – wykorzystujemy cały wektor
        if all_data.ndim > 1:
            all_data = all_data[0]
        n_iso = len(all_data)
        isotopy = []  # Może być pusta lub inne
        lam_vector = [general_params["lam"]] * n_iso
    return isotopy, n_iso, lam_vector, all_data


def get_n_measurements(det, n_iso, general_params):
    """
    Wyznacza liczbę pomiarów: n_meas = det["n_measurements"] lub n_iso * general_params["a"].
    """
    val = det.get("n_measurements", "")
    if isinstance(val, str):
        if val.strip():
            return int(val)
        else:
            return n_iso * general_params["a"]
    elif isinstance(val, int):
        return val
    else:
        return n_iso * general_params["a"]


def process_detector(det, general_params):
    """
    Przetwarza pojedynczy detektor:
      - Wczytuje plik A_file
      - Ustala definicję izotopów i stałych λ
      - Buduje macierz A oraz wyznacza pomocnicze ciągi (cs, dt_list, t_starts)
      - Wczytuje dane zliczeń oraz wykonuje korektę
      - Wykonuje alternatywną analizę eksponencjalną
    Zwraca słownik z wynikami oraz ewentualne komunikaty o błędach.
    """
    errors = []
    try:
        A_data = np.loadtxt(det["A_file"])
    except Exception as e:
        errors.append(f"Detektor {det.get('name','Unnamed')}: Nie można wczytać pliku A ({det.get('A_file','')}) - {e}")
        return None, errors

    isotopy, n_iso, lam_vector, A_data = get_isotopy_and_lam(A_data, det, general_params)
    n_meas = get_n_measurements(det, n_iso, general_params)

    try:
        A_block, m, scheme, A_single_orig, cs, dt_list, t_starts = build_A_matrix(A_data, det, general_params, lam_vector, n_meas)
    except Exception as e:
        errors.append(f"Detektor {det.get('name','Unnamed')}: Błąd przy budowaniu macierzy A - {e}")
        return None, errors

    # Wczytanie danych zliczeń
    try:
        counts_data = np.loadtxt(det["counts_file"])
    except Exception as e:
        errors.append(f"Detektor {det.get('name','Unnamed')}: Nie można wczytać pliku zliczeń ({det.get('counts_file','')}) - {e}")
        return None, errors

    try:
        y_raw = counts_data[:m, 0]
        y_err = counts_data[:m, 1]
    except Exception as e:
        errors.append(f"Detektor {det.get('name','Unnamed')}: Dane zliczeń mają nieprawidłowy format - {e}")
        return None, errors

    try:
        y_cor, y_err_cor = correct_counts(y_raw, y_err, det, general_params, m, scheme)
    except Exception as e:
        errors.append(f"Detektor {det.get('name','Unnamed')}: Błąd przy korekcie zliczeń - {e}")
        return None, errors

    # Alternatywna analiza eksponencjalna
    try:
        alt_A0, alt_err = alternative_exponential_decay_fit(A_single_orig, cs, dt_list, t_starts, y_cor, y_err_cor, isotopy)
    except Exception as e:
        errors.append(f"Detektor {det.get('name','Unnamed')}: Błąd alternatywnej analizy eksponencjalnej - {e}")
        alt_A0, alt_err = None, None

    result = {
        "A_block": A_block,
        "m": m,
        "scheme": scheme,
        "A_single_orig": A_single_orig,
        "cs": cs,
        "dt_list": dt_list,
        "t_starts": t_starts,
        "y_cor": y_cor,
        "y_err_cor": y_err_cor,
        "alt_A0": alt_A0,
        "alt_err": alt_err,
        "detector_name": det["name"],
        "has_scheme": scheme is not None,
        "isotopy": isotopy
    }
    return result, errors

def run_analysis(general_params, detectors):
    big_A_list = []
    big_y_list = []
    big_sigma_list = []
    detector_segments = []      # Liczba pomiarów dla każdego detektora
    successful_detectors = []   # Nazwy detektorów dla których udało się pobrać dane
    A_blocks = []               # Macierze A
    detector_has_scheme = []    # Informacja o schemacie rotacji
    alt_results = []            # Wyniki alternatywnej analizy
    total_errors = []

    for det in detectors:
        det_result, det_errors = process_detector(det, general_params)
        if det_errors:
            total_errors.extend(det_errors)
        if det_result is None:
            continue

        # Zbierzemy wyniki z danego detektora
        A_block = det_result["A_block"]
        m = det_result["m"]
        big_A_list.append(A_block)
        detector_segments.append(m)
        A_blocks.append(A_block)
        successful_detectors.append(det["name"])
        detector_has_scheme.append(det_result["has_scheme"])

        big_y_list.append(det_result["y_cor"])
        big_sigma_list.append(det_result["y_err_cor"])

        # Zapis wyników alternatywnej analizy
        alt_results.append((det["name"], det_result["alt_A0"], det_result["alt_err"]))

    if not big_A_list:
        return {"errors": total_errors, "y": []}

    # Sklejane macierzy i wektorów
    A_wysokie = np.vstack(big_A_list)
    y = np.hstack(big_y_list)
    y_sigma = np.hstack(big_sigma_list)
    N, k = A_wysokie.shape

    try:
        W_sqrt = np.diag(1.0 / y_sigma)
    except Exception as e:
        return {"errors": total_errors + [f"Błąd przy obliczaniu wag: {e}"], "y": y.tolist()}
    Aw = W_sqrt @ A_wysokie
    yw = W_sqrt @ y

    try:
        x_nnls, rnorm = nnls(Aw, yw)
    except Exception as e:
        return {"errors": total_errors + [f"Błąd przy dopasowaniu NNLS: {e}"], "y": y.tolist()}

    residuals = Aw @ x_nnls - yw
    chi2 = np.sum(residuals**2)
    res_var = chi2 / (N - k) if N > k else chi2
    M = Aw.T @ Aw
    try:
        M_inv = inv(M)
    except Exception:
        M_inv = pinv(M)
    cov_x = res_var * M_inv
    param_errors = np.sqrt(np.diag(cov_x))
    y_est = A_wysokie @ x_nnls
    cov_y_est = A_wysokie @ cov_x @ A_wysokie.T
    y_est_errors = np.sqrt(np.diag(cov_y_est))

    # Wyliczenie metryk dla poszczególnych detektorów
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

    return {
        "x_nnls": x_nnls.tolist(),
        "param_errors": param_errors.tolist(),
        "chi2": chi2,
        "detector_metrics": metrics_all,
        "detector_names": successful_detectors,
        "y": y.tolist(),
        "y_unc": y_sigma.tolist(),
        "y_est": y_est.tolist(),
        "y_est_errors": y_est_errors.tolist(),
        "A_wysokie_shape": A_wysokie.shape,
        "y_shape": y.shape,
        "A_blocks": A_blocks,
        "detector_segments": detector_segments,
        "detector_has_scheme": detector_has_scheme,
        "errors": total_errors,
        "alt_results": alt_results
    }