# Test działania algorytmu – dane syntetyczne typu "piła" (C11 i N13)

## Opis testu

Test sprawdza poprawność działania algorytmu na danych syntetycznych pochodzących z symulacji dwóch zestawów aktywności izotopów C11 i N13. Każdy zestaw składa się z 16 slotów (miejsc), z których 15 zawiera określoną aktywność, a jeden slot pozostaje pusty (nie zawiera żadnej aktywności). 

Zasymulowane źródła emitowały pojedynczą gammę o energii 511 keV. Intensywności linii = 1

## Aktywności źródeł

Aktywności (w Bq) przypisane do kolejnych slotów (licząc od 1 do 16):

- **C11:**  
  3e5, 1e5, 3e5, 1e5, 3e5, 1e5, 3e5, 1e5, 3e5, 1e5, 3e5, 1e5, 3e5, 1e5, 3e5, 0

- **N13:**  
  1e5, 3e5, 1e5, 3e5, 1e5, 3e5, 1e5, 3e5, 1e5, 3e5, 1e5, 3e5, 1e5, 3e5, 1e5, 0

Wartość `0` oznacza pusty slot, czyli brak aktywności w danym miejscu.

## Kryterium zaliczenia testu

Test uznaje się za zaliczony, jeśli wyniki uzyskane z aplikacji dla obu izotopów (C11 i N13) są zgodne z powyższymi wartościami aktywności w granicach niepewności pomiarowej.

## Dodatkowe informacje

- Analiza dotyczy jednoczesnego przetwarzania danych dla obu izotopów.
- Nie ma dodatkowych kryteriów zaliczenia poza zgodnością aktywności.