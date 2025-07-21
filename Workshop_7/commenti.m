% COMMENTI AL WORKBOOK 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Missione ICESat
% il fatto che sia in orbita LEO a 600 km deve essere riscontrato anche dai
% risultati. il fatto che i ricevitore non sia a terra ma abbia una
% altitudine del genere fa capire che alcuni errori possono non essere
% considerati. 
% il gps permette di riferire le misure dell'altimetro laser al sistema
% WGS84, quindi una accurata stima del gps si traduce in una accurata stima
% dello spessore dei ghiacci delle calotte polari. il GLAS lavora in due
% frequenze d'onda, visibile e infrarosso, questo serve per caratterizzare
% l'altezza della copertura nuvolosa nelle calotte polari, mentre l'altra
% fraquenza serve principalmente per la misura dell'altezza delle calotte
% polari. si utilizzano due antenne per ridondanza e non per misure
% differenziali 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Missione CHAMP
% Il magnetometro viene tenuto lontano dai campi magnatici spuri generati
% dalla strumentazione a bordo. vedere sempre la altezza dell'orbita in
% modo da capire quali tipi di errori sulla stima di posizione possono
% diventare rilevanti.
% posizioni dell'antenna del gps (precise orbit determination) stima in
% modo accurato la posizione orbitale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBIETTIVI
% Ricostruire la posizione orbitale dei due satelliti con incertezza per
% almeno 100 minuti in modo da visualizzare almeno un'orbita completa per
% entrambi. Bisogna fare la conversione del tempo in gps time. I data file
% used sono soltanto per il satellite icesat, mentre per il satellite champ
% dobbiamo trovarli noi dal link scritto sulle slide, dobbiamo risalire ai
% dati di navigazione e di osservazione alla data indicata. servono i dati
% rinex osservazione e broadcast di navigazione oppure anche quelli sp3 per
% avere una stima più precisa, che hanno la stima del clock offset dei
% satelliti gps. determianta la posizione dei satelliti bisogna usare
% l'algoritmo per stimare i parametri kepleriani.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITMO PER CALCOLO DEI PARMETRI KEPLERIANI (SCHEMA A BLOCCHI)
% dato che si stima la posizione e non la velocità
% Siccome l'obiettivo è determinare i parametri orbitali dobbiamo calcolare
% il vettore di stato, che non è fornito dal metodo dei minimi quadrati,
% quindi si può usare la soluzione del problema di lambert, che permette di
% calcolare la velocità ad ogni singolo istante di acqusizione. Noto il
% vettore di stato ci sono le routine di astrodinamica per passare da
% coordinate ecef ai classical orbital elements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAMBERT
% dovrebbe essere la solita rountine che è scritta anche nel curtis, vedere
% astrodinamica per il link del libro
% z si trova con il metodo iterativo di newton raphson ad esempio
% bisogna prendere per due punti consecutivi 1-2 poi 2-3 per trovare man
% mano la v1 e poi la v2 poi le altre.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RINEX CONSIDERAZIONI
% versione del rinex è la 2.2 quindi non si usa rinexread ma ci daranno 
% loro una routine, che crea una matrice con gli elementi in ogni colonna.
% Per i file di nav non sono relativi ad icesat, ma da una osservazione
% generale che comunque è valida, la versione è la 3 quindi si può usare
% rinexread che crea una struttura con i vari id della costellazione ecc.
% bisogna scegliere accuratamente un certo istante per poi propagare la
% nostra orbita
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dato che la stima di posizone del gps è affetta da rumore quando si
% considerano le due coppie di punti per lambert si ha la dipendenza dal
% rumore, nel senso che prendere cooppie di punti non consecutivi ma con
% tempi maggiorni fa avere una differenza, diepnde tutto dal rumore. quindi
% si può fare un'interpolazione dell'orbita. le posizioni vanno riportate
% da un ecef a un eci perché altrimenti c'è il contributo relativo alla
% rotazione della terra, allora vado in quello inerziale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% errori posizioni e perturbazioni orbitali fanno in modo che considerando
% il vettore di stato ottenuto per una coppia di punti poi gli elementi
% orbitali non siano uguali a quelli che ottengo un il vettore di stato per
% una coppia di punti diversa