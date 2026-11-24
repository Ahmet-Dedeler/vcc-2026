Convalida e invio di una sottomissione
======================================

Aggiornato il 2026-09-13. Introdotto da
[CP-0004](checkpoints/0004-primo-trial-locale-e-pacchetti.md), esteso con il percorso
di packaging da [CP-0005](checkpoints/0005-packaging-streaming-trial01.md).

**Una sottomissione è stata inviata e valutata**, il 2026-09-13: entry
`PNn227rxP3bVByS37W41`, punteggio **0,045929**, rango **446 su 920 squadre**. Il
resoconto completo è nella sezione 4 e in
[`submission_PNn227rxP3bVByS37W41.md`](../reports/trial_2026-09-13/submission_PNn227rxP3bVByS37W41.md);
il checkpoint è [CP-0006](checkpoints/0006-prima-sottomissione-e-punteggio.md).
Era `trial-01-transfer`; **`trial-00-controls` non è stato inviato** e non va inviato
(D-017).

> **`trial-01-transfer` è impacchettato; `trial-00-controls` no, e non deve esserlo.**
> Il `.vcc` di trial-01 è stato prodotto il 2026-09-13 con
> `scripts/48_package_prediction.py`, che convalida e scrive **senza caricare la
> matrice**: 0,519 GiB di picco contro i 33,5 che il modello della CLI attribuisce a
> `vcc prep`. Tutte e 24 le convalide attive, payload verificato bit a bit contro
> l'input. Evidenza:
> [`k01pack_packaging.json`](../reports/trial_2026-09-13/k01pack_packaging.json),
> [CP-0005](checkpoints/0005-packaging-streaming-trial01.md).
>
> `vcc prep` continua a non girare qui, ed è documentato in
> [`q01prep_prep_dry_run.log`](../reports/trial_2026-09-12/q01prep_prep_dry_run.log).
> Non serve più: le sue convalide sui metadati sono importate e chiamate dal nostro
> percorso, e la parità sui rifiuti è dimostrata da 46 test su fixture a forma
> ufficiale completa (D-019).
>
> **Il server ha accettato**, il 2026-09-13: la sottomissione è arrivata a `published`
> con `md5_verified: true` ed `error_info: null`. Il percorso di packaging in streaming
> è confermato end-to-end. Resta vero che un `.vcc` valido è un'affermazione sul
> **formato**: il punteggio, misurato, è 0,046.

## 1. Il contratto, verificato oggi

Ogni riga viene dalla pagina ufficiale *Evaluation* di virtualcellchallenge.org
(consultata il 2026-09-12), dalle FAQ ufficiali, o da `vcc prep --help` della CLI
installata (`vcc 0.2.0`). Le prime due sono la fonte; la CLI è ciò che il server usa
per rifiutare.

| Requisito | Valore | Dove è scritto |
|---|---|---|
| Un solo file per tutti e tre i contesti | non uno per contesto | Evaluation, FAQ "Do I submit one file per cell context?" |
| Colonna delle perturbazioni in `.obs` | `target_gene`, simboli (`ADNP`), non id di costrutto | Evaluation |
| Colonna del contesto in `.obs` | `context`, etichette identiche a quelle dei file di controllo | Evaluation |
| Perturbazioni | esattamente quelle di `pert_counts.csv`, per ogni contesto; una in più viene rifiutata | Evaluation |
| Cellule per perturbazione | esattamente 400, in ogni contesto | Evaluation |
| Geni in `.var` | esattamente 18.533, nell'ordine di `gene_names.csv` | Evaluation |
| `.X` | conteggi grezzi: non negativi, interi, finiti. Niente `normalize_total`, niente `log1p` | Evaluation, FAQ |
| Righe di controllo | **vietate**: nessuna riga `non-targeting`. I controlli scaricati sono soltanto input del modello | Evaluation |
| Valori memorizzati in `.X` | al massimo 4.750.000.000 in tutto, circa 13.200 per cellula a 360.000 cellule. Conta ciò che la matrice *memorizza*, quindi anche uno zero salvato esplicitamente | Evaluation, FAQ "rejected as too dense" |
| Conteggi per cellula | nessuna cellula oltre 1.000.000 sommando i geni | Evaluation |
| Forma attesa | 360.000 × 18.533 = 300 perturbazioni × 400 cellule × 3 contesti | Evaluation |
| Sottomissioni al giorno | **due** per squadra, il giorno cambia a mezzanotte UTC; contano solo quelle che arrivano al punteggio; una sola in volo alla volta | Rules, Evaluation |
| Codice e pesi | **non** vanno consegnati. Solo i finalisti devono fornire una descrizione di alto livello che verrà pubblicata (dataset usati, preprocessing, modelli, componenti non appresi, fasi di training; il compute è opzionale) | Rules §"Finalists, Winners, and Prizes"; FAQ "Am I required to submit my code?" |
| Dati esterni | ammessi, se se ne hanno i diritti; Arc non li richiede | FAQ "Can I use other data…" |
| Predizioni | devono essere generate **soltanto** da modelli di machine learning. Vietato inviare risultati sperimentali o di letteratura come parte della predizione (usarli per allenare è invece ammesso) | Rules §"Machine Learning Predictions Only" |

Confermato anche dal bundle locale: `raw/controls/manifest.json` dichiara
`partition: val`, `panel_id: vcc2026-val-1`, `cells_per_pert: 400`, 46 `ntc_id` per
contesto e 18.400 cellule di controllo per contesto.

Verifiche di sola lettura sull'account, eseguite il 2026-09-12 con
`vcc whoami --json`: profilo `default`, endpoint di produzione, token dal keyring,
account attivo e verificato, `can_submit: true`, `blockers: []`, squadra presente.
**Nessun token è stato stampato, copiato o scritto in un file di questo repository.**

## 2. Le due domande da risolvere prima di inviare trial-00

`trial-00-controls` è un ricampionamento delle cellule di controllo reali etichettate
con le perturbazioni richieste. `vcc prep` lo accetta — non contiene righe
`non-targeting` e formalmente è in regola — ma due frasi delle regole ufficiali lo
riguardano direttamente:

1. «The control cells you downloaded are model inputs only — do not copy them into
   your prediction.» (pagina Evaluation)
2. «your predictions must be generated solely by one or more machine learning
   models you use.» (Rules)

Un ricampionamento non è una previsione di un modello. **Non inviare trial-00 senza
aver chiarito questo punto** con help@virtualcellchallenge.org, o senza la decisione
esplicita del proprietario del progetto. Il pacchetto esiste perché serve come
controllo operativo locale e come riferimento di dispersione: è reale, quindi la sua
variabilità fra cellule è giusta per costruzione. Questo non lo rende inviabile.

`trial-01-transfer` non ha questo problema: è la risposta prevista da un modello
calibrato, applicata al profilo basale del contesto.

## 3. Comandi esatti, per Windows

Percorsi di questa macchina; `<run>` è la cartella del run sotto
`C:\Users\ferra\vcc2026-data\artifacts\`.

### Rigenerare tutto da zero

`$RUN` è l'id del run nuovo, e **il passo 3 usa lo stato fittato del passo 2, non
quello di un run precedente.** Un comando che rigenera la calibrazione e poi punta a
`c002` produrrebbe una previsione che non corrisponde alla calibrazione appena
