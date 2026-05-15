Il tool del calcolo del rumore tonale di una propeller è contenuto all'interno della cartella noise_computation. 
Sono necessarie analisi BEMT per lo svolgimento del calcolo. Sono presenti due simulazioni già processate.

Sono state implementate due metodologie di calcolo del rumore: formulazioni Gutin e Hanson. 

Script relativi al calcolo *Gutin*:

mic_single -> fornisce il confronto con i reference data del 4412 mid/high-fedelity
mic_array -> fornisce i plot di direttività del rumore relativi a entrambe le simulazioni
thickness_versus_loading -> realizza un confronto tra le singole sorgenti di rumore con i reference data 4412, 
                            calcola i thickness e loading dello 0018
steady_versus_unsteady_loading -> fornisce i risultati di un calcolo con carico instazionario

Codici relativi al calcolo *Hanson*:

hanson_comparison -> realizza un confronto tra i calcoli Hanson e Gutin, in termini di direttività e intensità 
                     del rumore. Il calcolo Hanson viene eseguito dal tool dii J. Goyal, contenuto nella cartella:
                     hanson-model-helicoidal-theory-master, per la configurazione di calcolo fare riferimento al file
                     readme contenuto al suo interno.
hanson_csv_generation -> consente di generare comodamente i file csv di input per il tool Hanson, nel formato richiesto
                         dalla sua modalità di raccolta dati. I file vengon generati nella cartella esterna hanson_csv_repository.


**Nota**: riguardo la sottostima di Gutin rispetto ad Hanson, la invito a prendere visione dei risultati dello script 
"thickness_versus_loading", dove è evidente un off-set al ribasso del rumore tonale rispetto ai risultati del modello 
a fedeltà più alta. Questa sottostima è presente anche nel confronto con Hanson, lanciando "hanson-comparison".
