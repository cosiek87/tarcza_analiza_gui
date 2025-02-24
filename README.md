# Analiza sąsiadów

## Wprowadzenie

Program **"Analiza sąsiadów"** służy do analizy danych z detektorów promieniowania. Umożliwia wczytywanie danych z plików, konfigurowanie parametrów analizy, uruchamianie analizy oraz wyświetlanie wyników w formie tekstowej i graficznej.

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
