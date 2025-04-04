# Analiza sąsiadów – Kompletny Manual

## Wprowadzenie i Szybki Start

Program **"Analiza sąsiadów"** to narzędzie oparte na interfejsie graficznym (GUI) wykorzystujące bibliotekę Tkinter, które służy do analizy danych z detektorów promieniowania. Umożliwia ono wczytywanie danych z plików, konfigurację parametrów analizy, uruchamianie procesu analizy oraz prezentację wyników w formie tekstowej i graficznej.

### Jak zacząć?
1. **Instalacja i uruchomienie:**
   - Upewnij się, że masz zainstalowane wymagane biblioteki: `tkinter`, `numpy`, `matplotlib`, `scipy`.
   - Pobierz plik `tarcza_analiza_gui.py`.
   - Uruchom program, wpisując w terminalu:
     ```sh
     python tarcza_analiza_gui.py
     ```
2. **Obsługa interfejsu:**
   - W sekcji **Parametry ogólne** ustaw odpowiednie wartości (m.in. opóźnienie czasowe, czas obserwacji, intensywność źródła).
   - Dodaj detektory poprzez przycisk **Dodaj detektor** – podaj nazwę, ścieżki do plików (macierz wydajności, dane zliczeń, opcjonalnie pliki schematu rotacji) oraz określ liczbę pomiarów i listę izotopów.
   - Po dodaniu przynajmniej jednego detektora kliknij **Uruchom analizę**. Wyniki zostaną wyświetlone zarówno w polu tekstowym, jak i w postaci wykresów.

---

## Pełny Opis Działania Programu

### 1. Wczytywanie Danych i Konfiguracja
- **Wczytywanie plików:**  
  Program odczytuje:
  - Plik z macierzą wydajności (plik A) – zawiera efektywności poszczególnych izotopów.
  - Plik zliczeń – dane eksperymentalne z przypisanymi błędami.
  - Opcjonalne pliki ze schematem rotacji i czasami rotacji, które umożliwiają ustalenie przesunięć tarczy oraz określenie czasu trwania pomiarów.
- **Konfiguracja parametrów:**  
  W sekcji **Parametry ogólne** użytkownik ustawia m.in.:
  - **lam:** Stała rozpadu (jeśli izotopy są zdefiniowane, program przypisuje odpowiednie wartości z wbudowanego słownika).
  - **Time delay:** Opóźnienie między startem a pierwszym pomiarem.
  - **TOB:** Czas obserwacji.
  - **Intensity:** Intensywność źródła.

### 2. Dodawanie i Edytowanie Detektorów
- **Dodawanie detektora:**  
  Po kliknięciu przycisku **Dodaj detektor** otwiera się okno dialogowe, w którym należy:
  - Wprowadzić nazwę detektora.
  - Podać ścieżki do plików: macierz wydajności (plik A) oraz plik zliczeń.
  - Określić liczbę pomiarów, listę izotopów (oddzielonych przecinkami) oraz ścieżki do plików ze schematem rotacji.
- **Edytowanie detektora:**  
  Możliwe jest edytowanie wcześniej dodanych detektorów poprzez przycisk **Edytuj detektor** – wprowadzone zmiany są natychmiast widoczne na liście detektorów.

### 3. Proces Analizy Danych
Program realizuje analizę w kilku etapach:

#### a) Przygotowanie danych
- Dane z różnych detektorów są scalane. Jeśli zdefiniowano izotopy, dla każdego przypisywana jest stała rozpadu (λ).
- Na podstawie pliku A oraz danych z schematu rotacji tworzona jest macierz pomiarowa (A) przez funkcję `build_A_matrix`.
- Dane zliczeń są korygowane przy użyciu funkcji `correct_counts`, aby uwzględnić czas trwania pomiaru.

#### b) Dopasowanie modelu
- Z uzyskanej macierzy A oraz wektora zliczeń **y**, wraz z przypisanymi wagami (błędy pomiarowe), rozwiązywany jest problem NNLS (Non-Negative Least Squares). Celem jest znalezienie parametrów $\( x_{\text{nnls}} \)$, tak aby:
  $$
  A \times x \approx y
  $$

#### c) Obliczenie błędów i propagacja niepewności
- Po dopasowaniu modelu obliczana jest macierz kowariancji parametrów.
- Przy użyciu propagacji błędów wyznaczany jest błąd dla każdej wyestymowanej wartości $\( y_{\text{est}} \)$:
  $$
  \text{cov}_{y_{\text{est}}} = A_{\text{wysokie}} \times \text{cov}_x \times A_{\text{wysokie}}^T, \quad y_{\text{est\_errors}} = \sqrt{\text{diag}(\text{cov}_{y_{\text{est}}})}.
  $$

#### d) Obliczanie metryk jakości
Dla każdego detektora wyznaczane są następujące metryki:
- **aFactor:**  
  Średnia względna różnica między danymi eksperymentalnymi a wyestymowanymi wartościami, definiowana jako:
  $$
  \text{aFactor} = \text{mean}\left(\frac{|y_{\text{exp}} - y_{\text{est}}|}{y_{\text{exp}} + y_{\text{est}} + \epsilon}\right)
  $$
- **RMSE (Root Mean Squared Error):**  
  Pierwiastek z średniej kwadratów różnic między danymi a oszacowanymi wartościami.
- **MAE (Mean Absolute Error):**  
  Średnia wartość bezwzględnych różnic między danymi a wynikami modelu.
- **R2 (Współczynnik determinacji):**  
  Określa, jaka część zmienności danych jest wyjaśniona przez model – wartość bliska 1 oznacza bardzo dobre dopasowanie.
- **Pearson:**  
  Współczynnik korelacji Pearsona, który ocenia liniową zależność między danymi a wynikami modelu.

### 4. Prezentacja Wyników
- **Wyniki tekstowe:**  
  Wyniki analizy, w tym oszacowane parametry, błędy, wartość chi-kwadrat oraz metryki jakości, są wyświetlane w polu **Wyniki analizy**. Dla każdego detektora prezentowane są szczegółowe informacje, co ułatwia identyfikację potencjalnych problemów.
- **Wykresy:**  
  Program generuje następujące rodzaje wykresów:
  - **Dane vs. Model:**  
    Dane eksperymentalne są przedstawione jako punkty z errorbarami (odzwierciedlającymi niepewność pomiarową), natomiast model jest przedstawiony jako linia. Wokół linii modelu rysowany jest korytarz błędu (np. $±3·\( y_{\text{est\_errors}} \)$), co pozwala ocenić zakres niepewności.
  - **Scatter Plot (Rozrzut):**  
    Wykres przedstawia zależność między danymi eksperymentalnymi a wartościami wyestymowanymi przez model. Dodatkowo rysowana jest linia idealnego dopasowania $\( y = x \)$, która umożliwia ocenę liniowej korelacji.
  - **Macierz Wydajności:**  
    Jeśli dla detektora dostępny jest schemat rotacji, generowany jest wykres macierzy wydajności, gdzie kolorystyka obrazuje poziom efektywności poszczególnych źródeł. Pozwala to na szczegółową ocenę jakości pomiarów.

### 5. Zapisywanie i Wczytywanie Projektów
- **Zapis projektu:**  
  Użytkownik może zapisać bieżący projekt (konfigurację parametrów oraz listę detektorów) do pliku JSON, wybierając opcję **Plik > Zapisz projekt**.
- **Wczytywanie projektu:**  
  Aby wczytać wcześniej zapisany projekt, należy wybrać **Plik > Wczytaj projekt** i wskazać odpowiedni plik JSON. Dzięki temu konfiguracja zostanie przywrócona, co ułatwia kontynuowanie pracy.