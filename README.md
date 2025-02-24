# Analiza sąsiadów

Dokumentacja użytkownika programu **"Analiza sąsiadów"**.

---

## Wymagania systemowe

- **Python:** 3.x
- **Biblioteki:**
  - tkinter
  - numpy
  - matplotlib
  - scipy

---

## Instalacja

1. Upewnij się, że masz zainstalowanego Pythona 3.x.
2. Zainstaluj wymagane biblioteki (tkinter jest zazwyczaj dołączony do instalacji Pythona):
   ```sh
   pip install numpy matplotlib scipy
   ```

---

## Uruchomienie programu

1. Otwórz terminal lub wiersz poleceń.
2. Przejdź do katalogu, w którym znajduje się plik `tarcza_analiza_gui.py`.
3. Uruchom program:
   ```sh
   python tarcza_analiza_gui.py
   ```

---

## Interfejs użytkownika

Po uruchomieniu programu pojawi się okno GUI z następującymi elementami:

- **Przyciski:**
  - `Wczytaj plik`: Pozwala na wczytanie pliku z danymi pomiarowymi.
  - `Analizuj`: Rozpoczyna analizę wczytanych danych.
  - `Zapisz wyniki`: Zapisuje wyniki analizy do pliku.

- **Pola tekstowe:**
  - Wyświetlają wczytane dane oraz wyniki analizy.

---

## Funkcje programu

- **aFactor(y_exp, y_theor):** Oblicza metrykę *aFactor*, która jest średnią wartością |y - y_est|/(y + y_est).
- **compute_rmse(y, y_est):** Oblicza pierwiastek średniokwadratowy (RMSE) między wartościami rzeczywistymi a estymowanymi.

---

## Przykładowe użycie

1. Kliknij przycisk `Wczytaj plik` i wybierz plik z danymi pomiarowymi.
2. Kliknij przycisk `Analizuj`, aby przeprowadzić analizę danych.
3. Wyniki analizy zostaną wyświetlone w odpowiednich polach tekstowych.
4. Kliknij przycisk `Zapisz wyniki`, aby zapisać wyniki analizy do pliku.

---

## Uwagi

- Upewnij się, że plik z danymi pomiarowymi jest w odpowiednim formacie.
- W razie problemów z uruchomieniem programu sprawdź, czy wszystkie wymagane biblioteki są zainstalowane.

---

## Kontakt

W razie pytań lub problemów skontaktuj się z autorem programu.

---

# Dokumentacja użytkownika programu "Analiza danych eksperymentalnych"

## Wprowadzenie

Program **"Analiza danych eksperymentalnych"** służy do analizy danych z detektorów promieniowania. Umożliwia wczytywanie danych z plików, konfigurowanie parametrów analizy, uruchamianie analizy oraz wyświetlanie wyników w formie tekstowej i graficznej.

---

## Instalacja

1. Upewnij się, że masz zainstalowane wymagane biblioteki:
    - tkinter
    - numpy
    - matplotlib
    - scipy

2. Pobierz plik `tarcza_analiza_gui.py` i uruchom go za pomocą Pythona:
   ```sh
   python tarcza_analiza_gui.py
   ```

---

## Interfejs użytkownika

### Główne okno aplikacji

Po uruchomieniu programu zobaczysz główne okno aplikacji, które składa się z kilku sekcji:

1. **Parametry ogólne:**  
   W tej sekcji możesz ustawić parametry analizy:
   - `lam`: Stała rozpadu izotopu. *(DO USUNIĘCIA)*
   - `Time delay`: Opóźnienie czasowe.
   - `TOB`: Czas obserwacji.
   - `Intensity`: Intensywność źródła. *(DO ZMIANY)*
   - `Liczba obrotów (a)`: Liczba obrotów. *(DO USUNIĘCIA)*

2. **Detektory:**  
   Lista detektorów używanych w analizie. Możesz dodawać, edytować i usuwać detektory.

3. **Wyniki analizy:**  
   Pole tekstowe, w którym wyświetlane są wyniki analizy.

4. **Przyciski:**
   - `Dodaj detektor`: Dodaje nowy detektor.
   - `Edytuj detektor`: Edytuje wybrany detektor.
   - `Usuń detektor`: Usuwa wybrany detektor.
   - `Uruchom analizę`: Uruchamia analizę danych.

### Menu

W menu aplikacji znajdują się opcje:

- **Plik:**
  - `Zapisz projekt`: Zapisuje bieżący projekt do pliku JSON.
  - `Wczytaj projekt`: Wczytuje projekt z pliku JSON.
  - `Wyjście`: Zamyka aplikację.

---

## Dodawanie i edytowanie detektorów

1. Kliknij przycisk `Dodaj detektor` lub `Edytuj detektor`.
2. W oknie dialogowym wprowadź dane detektora:
   - **Nazwa:** Nazwa detektora.
   - **Plik A (wydajność):** Ścieżka do pliku z macierzą wydajności.
   - **Plik zliczeń:** Ścieżka do pliku z danymi zliczeń.
   - **Liczba pomiarów:** Liczba pomiarów (opcjonalnie).
   - **Izotopy (przecinek):** Lista izotopów oddzielona przecinkami.
   - **Plik schematu rotacji:** Ścieżka do pliku z schematem rotacji (opcjonalnie).
3. Kliknij `OK`, aby zapisać detektor, lub `Anuluj`, aby anulować.

---

## Uruchamianie analizy

1. Upewnij się, że wprowadziłeś wszystkie parametry ogólne i dodałeś przynajmniej jeden detektor.
2. Kliknij przycisk `Uruchom analizę`.
3. Wyniki analizy zostaną wyświetlone w polu tekstowym **Wyniki analizy**. Jeśli analiza zakończy się błędami, zostaną one również wyświetlone.

---

## Zapisywanie i wczytywanie projektów

- Aby zapisać projekt, wybierz **Plik > Zapisz projekt**, a następnie podaj nazwę pliku.
- Aby wczytać projekt, wybierz **Plik > Wczytaj projekt**, a następnie wybierz plik projektu.

---

## Wyświetlanie wyników

Wyniki analizy są wyświetlane w formie tekstowej w polu **Wyniki analizy**. Dodatkowo, program generuje wykresy:
- Dane eksperymentalne.
- Model (A·x) vs. pomiary dla każdego detektora.
- Scatter plot: model vs. pomiary.
- Macierz wydajności (jeśli dostępny był schemat rotacji).
