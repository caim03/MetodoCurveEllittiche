# MetodoCurveEllittiche
Algoritmo per il calcolo di fattori non banali di n tramite l'utilizzo di curve ellittiche

### Compilazione
Per compilare correttamente il file occorre linkare la libreria "modulus.c":
  - Creare il library object file, responsabile del linking tra il programma e la libreria:
    gcc -o modulus.o -c modulus.c
  - Durante la compilazione del programma aggiungere il file precedente:
    gcc -o curve MetodoCurveEllittiche.c modulus.o
