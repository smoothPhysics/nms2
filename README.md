# Protokoll Abgabe NMS2

## Team
* Philipp Denk, 11714004

## Aufgaben
### propagation.cxx
Aufwand: 6h (Stunden)

## Fragen
* Die Ergebnisse für die beiden Hamiltonmatrize unterscheiden sich im Durchschnitt um 3.48%. 
* Die genaue Dauer der Berechnung beider Werte befindet sich im Outputfile des Programms. Allgemein kann man die Dauer
der Berechnung über die Runtime abschätzen, welche für Punkt a O(n) ist und für Punkt b O(n^3).
Somit kommt man bei N = 65 auf 11ms für Punkt a und auf 36ms bei Punkt b.
* Die Ergebnisse im outputfile eigenvalues.txt zeigen, dass die Genauigkeit mit steigender Anzahl der Gridpunkte wie erwartet steigt. 
Ebenso sieht man, dass der relative Fehler zwischen analytischen Werten aus J. Chem. Phys. 91, 3571 (1989) und berechneten Werten aus dem Programm bei 
Methode b merkbar kleiner ist als bei Methode a. Im Durchschnitt um 2 Größenordnungen. Selbiges gilt für in J. Chem. Phys. 91, 3571 (1989) nummerisch berechneten Werte. 

## Output propagation.cxx
Filename: "eigenvalues.txt"
Inhalt:
    -Angabe der Gridpunkte
    -1. Spalte: Werte aus Methode a
    -2. Spalte: Relative Abweichung zu den nummerischen Werten aus Tabelle II (J. Chem. Phys. 91, 3571 (1989))
    -3. Spalte: Werte aus Methode b
    -4. Spalte: Relative Abweichung zu den nummerischen Werten aus Tabelle II (J. Chem. Phys. 91, 3571 (1989))
    -5. Spalte: Relative Abweichung des Durchscnitts aus a und b von den analytischen Werten aus Tabelle II (J. Chem. Phys. 91, 3571 (1989))
    -6. Spalte: Relative Abweichung zwischen a und b


