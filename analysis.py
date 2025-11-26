from numpy import loadtxt, vstack, hstack, diag, sum, sqrt, arange, maximum, ones_like, clip, zeros, asarray, ones, log, zeros_like, minimum
from numpy.linalg import inv, pinv, norm
from utils import build_A_matrix, correct_counts, alternative_exponential_decay_fit, nuclear_data
from matplotlib.pyplot import show
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
        A_data = loadtxt(det["A_file"])
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
        counts_data = loadtxt(det["counts_file"])
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
    if general_params.get("exp_analysis", False):
        try:
            alt_A0, alt_err = alternative_exponential_decay_fit(A_single_orig, cs, dt_list, t_starts, y_cor, y_err_cor, isotopy, general_params["intensity"])
        except Exception as e:
            errors.append(f"Detektor {det.get('name','Unnamed')}: Błąd alternatywnej analizy eksponencjalnej - {e}")
            alt_A0, alt_err = None, None
    else:
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


def analysis(general_params, detectors):
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

    show()

    if not big_A_list:
        return {"errors": total_errors, "y": []}

    # Sklejane macierzy i wektorów
    A_wysokie = vstack(big_A_list)
    y = hstack(big_y_list)
    y_sigma = hstack(big_sigma_list)
    N, k = A_wysokie.shape

    # --- ROZWIĄZANIE WYBRANĄ METODĄ (NNLS / MLEM / R-OS-SPS) ---
    x_est, method_used, A_eff, supp = _solve_with_selected_method(A_wysokie, y, y_sigma, general_params)
    y_est = A_wysokie @ x_est

    # NNLS: możemy policzyć błędy parametrów; dla MLEM pomiń (inne założenia stat.)
    param_errors = []
    chi2 = float('nan')
    y_est_errors = zeros_like(y)

    if method_used == "NNLS":
        # odtwórz ważenia jak w solverze NNLS:
        W_sqrt = diag(1.0 / maximum(y_sigma, 1e-12))
        Aw = W_sqrt @ A_wysokie
        yw = W_sqrt @ y
        residuals = Aw @ x_est - yw
        chi2 = sum(residuals**2)
        res_var = chi2 / (N - k) if N > k else chi2
        M = Aw.T @ Aw
        try:
            M_inv = inv(M)
        except Exception:
            M_inv = pinv(M)
        cov_x = res_var * M_inv
        param_errors = sqrt(diag(cov_x)).tolist()
        cov_y_est = A_wysokie @ cov_x @ A_wysokie.T
        y_est_errors = sqrt(diag(cov_y_est))
    else:
        # --- MLEM / R-OS-SPS: ROBUST, ale przycięty ---
        b = zeros_like(y)
        x_red = x_est[supp]
        mu_hat = A_eff @ x_red + b

        ridge = float(general_params.get("cov_ridge", 1e-10))
        cap_factor = float(general_params.get("cov_cap_factor", 4.0))  # 4x Poisson domyślnie
        cov_x_red = _poisson_sandwich_cov_capped(
            A_eff, y, mu_hat,
            ridge=ridge,
            cap_factor=cap_factor
        )

        # overdispersion z deviance – ALE TEŻ PRZYCIĘTY
        dev = _poisson_deviance(y, mu_hat)
        dof = max(1, len(y) - A_eff.shape[1])
        dispersion = max(1.0, dev / dof)
        dispersion_max = float(general_params.get("cov_disp_max", 3.0))  # nie więcej niż 3x
        dispersion = min(dispersion, dispersion_max)

        cov_x_red = cov_x_red * dispersion

        # wstaw w pełny x
        param_errors_full = zeros(A_wysokie.shape[1])
        param_errors_full[supp] = sqrt(diag(cov_x_red))
        param_errors = param_errors_full.tolist()

        # propagacja
        cov_y_est = A_eff @ cov_x_red @ A_eff.T
        y_est_errors = sqrt(diag(cov_y_est))

        # poisson_deviance = float(dev)

    # Wyliczenie metryk per detektor (bez zmian)
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
        "method": method_used,
        "x_est": x_est.tolist(),
        "param_errors": param_errors,
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


def _solve_with_selected_method(A_full, y, y_sigma, general_params):
    """
    Zwraca (x_est, method_used).
    Czyta z general_params:
      - fit_method: 'NNLS' lub 'MLEM'
      - mlem_variant: 'CLASSIC' lub 'R-OS-SPS'
      - support_mask_sources: "0/1,..." (sloty)  [opcjonalnie]
      - support_mask_full:    "0/1,..." (kolumny)[opcjonalnie]
      - alpha0, tau, subsets  [dla R-OS-SPS]
    """
    n_cols = A_full.shape[1]
    # heurystyka liczby slotów: 16; jeśli znasz ją lepiej — nadpisz w general_params["n_sources"]
    n_sources = int(general_params.get("n_sources", 16))
    n_iso = max(1, n_cols // max(1, n_sources))

    def _parse_mask(s):
        return [int(x) != 0 for x in s.replace(";", ",").split(",") if x.strip() != ""]

    mask_sources = None
    mask_full = None
    if str(general_params.get("support_mask_sources", "")).strip():
        mask_sources = _parse_mask(general_params["support_mask_sources"])
    if str(general_params.get("support_mask_full", "")).strip():
        mask_full = _parse_mask(general_params["support_mask_full"])

    support_mask_bool = _make_full_support_mask(
        n_cols, n_sources, n_iso, mask_sources=mask_sources, mask_full=mask_full
    )
    A, supp = _reduce_by_support(A_full, support_mask_bool)

    fit_method = str(general_params.get("fit_method", "NNLS")).upper()
    if fit_method == "MLEM":
        variant = str(general_params.get("mlem_variant", "CLASSIC")).upper()
        b = zeros_like(y)
        if variant == "R-OS-SPS":
            alpha0 = float(general_params.get("alpha0", 1.0))
            tau = float(general_params.get("tau", 25.0))
            subsets = int(general_params.get("subsets", 8))
            x_red = relaxed_os_sps(A, y, x0=None, b=b,
                                   subsets=subsets, max_outer=int(general_params.get("max_iter", 200)),
                                   alpha0=alpha0, tau=tau, backtracking=True,
                                   tol=float(general_params.get("tol", 1e-5)))
            method_used = "MLEM (R-OS-SPS)"
        else:
            x_red = mlem(A, y, x0=init_x_mlem(A, y), b=b,
                         max_iter=int(general_params.get("max_iter", 200)),
                         tol=float(general_params.get("tol", 1e-5)),
                         damping=float(general_params.get("damping", 1.0)),
                         l2=float(general_params.get("l2", 0.0)),
                         l1=float(general_params.get("l1", 0.0)))
            method_used = "MLEM (classic)"
        x = _expand_solution(x_red, supp, n_cols)
        return x, method_used, A, supp

    # domyślnie: NNLS (Gauss, wagi 1/sigma)
    W = diag(1.0 / maximum(y_sigma, 1e-12))
    Aw = W @ A
    yw = W @ y
    from scipy.optimize import nnls
    x_red, _ = nnls(Aw, yw)
    x = _expand_solution(x_red, supp, n_cols)
    return x, "NNLS", A, supp


# Sekcja analizy Poissona #


def _poisson_deviance(y, mu, eps=1e-12):
    """Deviance Poissona: 2 * sum( y*log(y/mu) - (y - mu) ), z definicją 0*log(0)=0."""
    mu = maximum(mu, eps)
    # pierwszy składnik tylko tam, gdzie y>0
    term1 = zeros_like(y, dtype=float)
    mask = y > 0
    term1[mask] = y[mask] * (log(y[mask] / mu[mask]))
    return 2.0 * (sum(term1 - (y - mu)))


def _poisson_sandwich_cov_capped(A_eff, y, mu_hat,
                                 ridge=1e-10, eps=1e-12,
                                 cap_factor=5.0):
    """
    Robust (sandwich) dla Poissona, ale z ograniczeniem wkładu pojedynczych obserwacji.
    Kończy się między:
        I_exp^{-1}   a   cap_factor * I_exp^{-1}
    więc nie wystrzeli jak (y-mu)^2 jest chore.
    """
    # 1) oczekiwany Fisher
    W = diag(1.0 / maximum(mu_hat, eps))
    I_exp = A_eff.T @ W @ A_eff
    if ridge and ridge > 0.0:
        I_exp = I_exp + ridge * diag(ones(I_exp.shape[0]))
    try:
        I_exp_inv = inv(I_exp)
    except Exception:
        I_exp_inv = pinv(I_exp)

    # 2) "middle" – empiryczna kowariancja score'a
    r = y - mu_hat
    # surowy wkład sandwich:
    mid_diag = (r ** 2) / maximum(mu_hat ** 2, eps)

    # 2a) CAP: nie pozwól, żeby jedna obserwacja dała > cap_factor * (1/mu)
    # Poissonowy wkład to 1/mu_hat
    poisson_diag = 1.0 / maximum(mu_hat, eps)
    mid_diag = minimum(mid_diag, cap_factor * poisson_diag)

    midW = diag(mid_diag)
    M = A_eff.T @ midW @ A_eff

    # 3) sandwich
    cov = I_exp_inv @ M @ I_exp_inv
    return cov


def _fisher_cov_poisson(A_eff, mu_hat, ridge=0.0, eps=1e-12):
    """
    Zwraca macierz kowariancji z Fishera: (A^T diag(1/mu_hat) A + ridge*I)^(-1)
    Uwaga: ridge (małe >0) poprawia uwarunkowanie przy bardzo małych mu.
    """
    W = diag(1.0 / maximum(mu_hat, eps))
    H = A_eff.T @ W @ A_eff
    if ridge and ridge > 0.0:
        H = H + ridge * diag(ones(H.shape[0]))
    try:
        H_inv = inv(H)
    except Exception:
        H_inv = pinv(H)
    return H_inv

# === support mask helpers (kolumny A) ===


def _reduce_by_support(A, support_mask_bool):
    supp = asarray(support_mask_bool, dtype=bool)
    if supp.ndim != 1 or supp.shape[0] != A.shape[1]:
        raise ValueError("support_mask ma złą długość względem kolumn A.")
    return A[:, supp], supp


def _expand_solution(x_reduced, supp, n_full):
    x_full = zeros(n_full, dtype=float)
    x_full[supp] = x_reduced
    return x_full


def _make_full_support_mask(n_cols, n_sources, n_iso, mask_sources=None, mask_full=None):
    """Zwraca maskę bool długości n_cols.
    mask_full: pełna maska kolumn (n_cols).
    mask_sources: maska 0/1 po slotach (n_sources), powielana n_iso razy.
    """
    if mask_full is not None and len(mask_full) > 0:
        mf = asarray(mask_full, dtype=bool)
        if mf.size != n_cols:
            raise ValueError(f"mask_full musi mieć długość {n_cols}, a ma {mf.size}.")
        return mf
    if mask_sources is None:
        return ones(n_cols, dtype=bool)
    ms = asarray(mask_sources, dtype=bool)
    if ms.size != n_sources:
        raise ValueError(f"mask_sources musi mieć długość {n_sources}, a ma {ms.size}.")
    return hstack([ms for _ in range(n_iso)])


def _safe_div(a, b, eps=1e-12):
    return a / maximum(b, eps)


# === MLEM init + subsety (jeśli nie masz) ===


def init_x_mlem(A, y):
    return (A.T @ y) / maximum(A.T @ ones_like(y), 1e-12)


def make_subsets(M, S):
    idx = arange(M)
    return [idx[s::S] for s in range(S)]


# === Poisson log-like (do kontroli wzrostu) ===


def _poisson_loglike(A, x, y, b=None, eps=1e-12):
    mu = A @ x + (0.0 if b is None else b)
    return float(sum(y * log(maximum(mu, eps)) - mu))


def mlem(A, y, x0=None, b=None, max_iter=200, tol=1e-5, damping=1.0, l2=0.0, l1=0.0, callback=None):
    """
    Klasyczny MLEM (Richardson–Lucy) dla y ~ Poisson(A x + b), x >= 0.
    damping ∈ (0,1] łagodzi kroki; l2/l1 – delikatna regularizacja.
    """
    m, n = A.shape
    b = zeros(m) if b is None else asarray(b)
    x = clip(init_x_mlem(A, y) if x0 is None else x0, 1e-12, None)

    ones_m = ones_like(y)
    denom_base = A.T @ ones_m  # A^T 1

    for it in range(int(max_iter)):
        mu = A @ x + b
        ratio = A.T @ _safe_div(y, mu)      # A^T (y / (Ax+b))
        denom = denom_base.copy()            # A^T 1

        if l2 > 0.0:
            denom = denom + l2 * x
        if l1 > 0.0:
            denom = denom + l1 * _safe_div(1.0, x)

        mult = _safe_div(ratio, denom)
        x_new = clip(x * (mult ** damping), 0.0, None)

        if callback is not None:
            callback(it, x_new, mu)

        if norm(x_new - x) <= tol * (norm(x) + 1e-12):
            x = x_new
            break
        x = x_new

    return x


# === R-OS-SPS (relaxed Ordered-Subsets SPS) ===


def relaxed_os_sps(A, y, x0=None, b=None,
                   subsets=8, max_outer=200,
                   alpha0=1.0, tau=25.0,
                   backtracking=True, bt_factor=0.5,
                   tol=1e-5, rng=None):
    """
    Zbieżny OS-SPS dla Poissona z kompensacją subsetów i preconditionerem x/(A_s^T 1).
    Powinien numerycznie zbliżać się do MLEM (ta sama skala), a przy tym szybciej konwergować.
    """
    M, N = A.shape
    b = zeros(M) if b is None else asarray(b)

    # start jak w MLEM (blisko ML, dobra skala)
    if x0 is None:
        x = clip((A.T @ y) / maximum(A.T @ ones_like(y), 1e-12), 1e-12, None)
    else:
        x = clip(x0, 1e-12, None)

    S = max(1, int(subsets))
    subs = make_subsets(M, S)

    k = 0
    last_ll = _poisson_loglike(A, x, y, b)

    for it in range(int(max_outer)):
        order = range(S) if rng is None else rng.permutation(S)
        x_old_outer = x.copy()

        for s in order:
            subset = subs[s]
            As, ys, bs = A[subset, :], y[subset], b[subset]
            mu = As @ x + bs

            # gradient subsetu
            grad = As.T @ (ys / maximum(mu, 1e-12) - 1.0)

            # preconditioner zależny od x i subsetu: x / (A_s^T 1)
            denom_s = As.T @ ones_like(ys)
            pre = x / maximum(denom_s, 1e-12)

            # malejąca relaksacja i kompensacja subsetów (≈ 1/S pełnego gradientu)
            alpha = alpha0 / (1.0 + k / tau)
            alpha_eff = alpha * S

            # krok + projekcja na x>=0
            x_cand = clip(x + alpha_eff * (pre * grad), 0.0, None)

            if backtracking:
                ll_new = _poisson_loglike(A, x_cand, y, b)
                # zabezpieczenie przed zbyt dużym krokiem
                while ll_new < last_ll and alpha_eff > 1e-6:
                    alpha_eff *= bt_factor
                    x_cand = clip(x + alpha_eff * (pre * grad), 0.0, None)
                    ll_new = _poisson_loglike(A, x_cand, y, b)
                last_ll = ll_new
            else:
                last_ll = _poisson_loglike(A, x_cand, y, b)

            x = x_cand
            k += 1

        # stop po pełnym przebiegu: względna zmiana x
        rel = norm(x - x_old_outer) / (norm(x_old_outer) + 1e-12)
        if rel < tol:
            break

    return x
