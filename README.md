# dbsuper Pipeline

## One-shot SLURM driver (`submit_driver.sh`)

The whole pipeline — init → fetch → download → samplesheet → nf-core/chipseq — can be submitted as a single SLURM job. The driver wires together every step below into one re-runnable workflow with per-step logs.

### Usage

```bash
# metadata + URL samplesheet only (no FASTQ download, no nextflow)
sbatch submit_driver.sh --member 1

# same, but with a custom accession subset
sbatch submit_driver.sh --member 1 --accessions ENCSR000AKC,ENCSR000AKJ

# full run: download FASTQs, build local samplesheet, launch nf-core/chipseq
sbatch submit_driver.sh --member 1 --download

# resume-only fast path: skip fetch + samplesheet, go straight to nextflow -resume
sbatch submit_driver.sh --member 1 --download --skip-fetch
```

### Flags

| flag | required | purpose |
|---|---|---|
| `--member <1-5>` | yes | picks the pre-assigned accession list for that member |
| `--accessions A,B,C` | no | override with a custom comma-separated accession subset |
| `--download` | no | also pull FASTQ files + build local samplesheet + launch nextflow (without this the driver stops after fetching manifest/URL samplesheet) |
| `--skip-fetch` | no | skip steps 3–4, go directly to `nextflow run -resume`. Requires `raw/encode_results/nfcore_chipseq_samplesheet.local.csv` to already exist |

### Execution order inside the driver

1. **`init.sh`** — guarded by `.init_done` marker; runs only once.
2. **`module load apptainer/system nextflow/25.10.4 encodefetch/0.5.0`** — Lmod modules on Karakoram. `nextflow/25.10.4` bundles its own OpenJDK 23 and `encodefetch/0.5.0` bundles Python 3.11 with pandas (used by `make_nfcode_samplesheet.py`). `apptainer/system` must load first so its `/usr/bin` PATH entry doesn't shadow the encodefetch env's `python3`.
3. **Unrename pass** — renames any `*.fastq.gz` back to `*.fastq` so encodefetch's built-in skip-existing check works on re-runs (encodefetch keys on the `.fastq` suffix; `make_nfcode_samplesheet.py` later renames them back to `.gz`).
4. **`fetch_member_accessions.sh`** — runs encodefetch for the selected accessions.
5. **Gate** — if `--download` was not given, the driver exits here.
6. **`make_nfcode_samplesheet.py`** — renames `.fastq` → `.fastq.gz`, rewrites URL samplesheet with local paths into `nfcore_chipseq_samplesheet.local.csv`, sanitises/dedupes rows.
7. **`nextflow run nf-core/chipseq -r 2.1.0 -resume`** — with SLURM executor; child sbatch jobs dispatched per `nfcore_chipseq.config`.

### Logs

Each step writes to its own timestamped file under `logs/`:

```
logs/
├── driver/        driver-<jobname>-<jobid>.log|err   (SLURM's -o/-e)
├── init/          init-<TS>.log
├── fetch/         fetch-member<N>-<TS>.log
├── samplesheet/   samplesheet-<TS>.log
├── nextflow/      .nextflow.log (+ rotated .log.1, .log.2, ...)
└── reports/       nfcore_trace.txt, nfcore_report.html, nfcore_timeline.html, nfcore_dag.png
```

### Re-run / resume behavior

- Safe to re-submit with the exact same args — `init.sh` skips, encodefetch skips already-downloaded files, `make_nfcode_samplesheet.py` is idempotent, Nextflow `-resume` reuses task cache from `workDir`.
- To force re-init: `rm .init_done`.
- To force a clean nextflow run: `rm -rf .nextflow/` (leaves `workDir/` cached tasks intact; those get re-hashed on next run).
- To wipe Nextflow task cache: `rm -rf <workDir>` (path in `nfcore_chipseq.config`).

### Cancelling

```bash
scancel <driver-jobid>                                # kills driver; Nextflow's shutdown hook scancels children
squeue -u $USER -h -o '%A' | xargs -r scancel         # nuclear: kill every job you own
```

---

## Setup

### Load required modules

The pipeline uses Karakoram's Lmod modules instead of a project-specific conda
env. `submit_driver.sh` does this automatically; for interactive shells:

```bash
module load apptainer/system nextflow/25.10.4 encodefetch/0.5.0
```

Module set:

| Module | Provides |
|---|---|
| `apptainer/system` | singularity/apptainer container runtime |
| `nextflow/25.10.4` | Nextflow CLI + bundled OpenJDK 23 |
| `encodefetch/0.5.0` | `encodefetch` CLI + Python 3.11 with pandas |

Load order matters — `apptainer/system` prepends `/usr/bin` to `PATH`, which
would shadow the encodefetch env's `python3`. Load it first.

The legacy `environment.yml` is kept for reference only.

### Initialize directory structure

If you don't already have the directory structure, run:

```bash
bash init.sh
```

This creates: `out`, `raw`, `ref`,`logs` and `Enhancerflow` directories.

## Download ENCODE data and run nf-core/chipseq

### Get manifest and samplesheet

Adjust `--threads` based on your available CPUs. By default, encodefetch 0.5.0+
downloads FASTQs. Use `--metadata-only` to fetch just the manifest and
samplesheet without downloading files:

```bash
cd raw
encodefetch --assay-title "Histone ChIP-seq" --target-label H3K27ac --organism "Homo sapiens" --file-type fastq --status released --threads 20 --progress --nfcore --metadata-only
```

Drop `--metadata-only` to also download the FASTQs.

### Fetching your assigned accessions (split across 5 members)

Use `fetch_member_accessions.sh` — all accession lists are pre-loaded. Pass `--member` with your number (1–5).

**Step 1 — fetch manifest and samplesheet only (no files downloaded yet):**

```bash
# Fetch metadata and generate the nf-core samplesheet for your assigned accessions
bash fetch_member_accessions.sh --member 2

# Or override with a specific subset of accessions
bash fetch_member_accessions.sh --member 2 --accessions ENCSR000AKC,ENCSR000AKJ,ENCSR000AKP
```

**Step 2 — re-run with `--download` to also pull the FASTQ files:**

```bash
# Same command, add --download to fetch the actual FASTQ files
bash fetch_member_accessions.sh --member 2 --download

# Or for a specific subset
bash fetch_member_accessions.sh --member 2 --accessions ENCSR000AKC,ENCSR000AKJ,ENCSR000AKP --download
```

The script automatically `cd`s into `raw/` before running and prints which mode it is in.

### Work distribution (913 accessions across 5 members)

| Member 1 | Member 2 | Member 3 | Member 4 | Member 5 |
|----------|----------|----------|----------|----------|
| ENCSR000AKC | ENCSR033XKQ | ENCSR023HXE | ENCSR001UFN | ENCSR001KCX |
| ENCSR000AKJ | ENCSR043UWF | ENCSR055YMD | ENCSR011GWK | ENCSR049LCN |
| ENCSR000AKP | ENCSR071DPR | ENCSR058XPK | ENCSR018KWK | ENCSR059VTG |
| ENCSR000AKY | ENCSR096QCF | ENCSR059QYS | ENCSR038JBA | ENCSR071ILJ |
| ENCSR000ALB | ENCSR103IXF | ENCSR061EXX | ENCSR049JKQ | ENCSR072QNY |
| ENCSR000ALG | ENCSR113XPE | ENCSR061VJM | ENCSR064RLF | ENCSR089ZLR |
| ENCSR000ALK | ENCSR116QUS | ENCSR063AUO | ENCSR079HOU | ENCSR099RQB |
| ENCSR000ALS | ENCSR126UEY | ENCSR074WNF | ENCSR085WWL | ENCSR109IWL |
| ENCSR000ALW | ENCSR130HAG | ENCSR080WMI | ENCSR099LGH | ENCSR120ITZ |
| ENCSR000AMI | ENCSR139WMA | ENCSR086NVH | ENCSR107RQX | ENCSR123WDC |
| ENCSR000AMN | ENCSR143DQF | ENCSR107PGG | ENCSR113CMI | ENCSR124GAS |
| ENCSR000AMO | ENCSR149EOW | ENCSR111ZAD | ENCSR113QDZ | ENCSR131SEB |
| ENCSR000AMR | ENCSR151NEZ | ENCSR119AVG | ENCSR118FQN | ENCSR145RXI |
| ENCSR000AMZ | ENCSR153LYY | ENCSR122GIL | ENCSR124KQG | ENCSR193TCX |
| ENCSR000ANF | ENCSR156MYA | ENCSR147KFH | ENCSR128HTH | ENCSR212SJH |
| ENCSR000ANN | ENCSR156XNC | ENCSR170MJE | ENCSR138WRC | ENCSR214MNL |
| ENCSR000ANP | ENCSR163UEW | ENCSR185OHN | ENCSR154SJE | ENCSR232ABT |
| ENCSR000ANT | ENCSR163XFK | ENCSR200QAA | ENCSR155GOF | ENCSR276SJT |
| ENCSR000ANV | ENCSR164POT | ENCSR201KGX | ENCSR161JTT | ENCSR278ACD |
| ENCSR000AOB | ENCSR165CDY | ENCSR201VAH | ENCSR178FVP | ENCSR293XMP |
| ENCSR000AOC | ENCSR170MAJ | ENCSR255SDP | ENCSR178MWD | ENCSR300BKE |
| ENCSR000AOP | ENCSR174GMQ | ENCSR261PLD | ENCSR189UPU | ENCSR309LQG |
| ENCSR000AOQ | ENCSR175ABH | ENCSR278VUZ | ENCSR195BZZ | ENCSR309MYO |
| ENCSR000APH | ENCSR175HRD | ENCSR317VER | ENCSR196CXJ | ENCSR314KMO |
| ENCSR000APN | ENCSR175XPY | ENCSR329HFV | ENCSR198WIH | ENCSR324VMZ |
| ENCSR000APT | ENCSR177QFY | ENCSR344CSQ | ENCSR236NLU | ENCSR327PYZ |
| ENCSR000APU | ENCSR179YLS | ENCSR347VUG | ENCSR237ZMF | ENCSR328AHC |
| ENCSR000AQW | ENCSR180YHA | ENCSR354HOM | ENCSR243LHN | ENCSR339WBT |
| ENCSR000ASA | ENCSR181OHF | ENCSR380KOO | ENCSR247FUH | ENCSR361GKY |
| ENCSR000ASJ | ENCSR189JIH | ENCSR386CKJ | ENCSR255BOF | ENCSR376YUP |
| ENCSR000ASR | ENCSR189QAD | ENCSR391EQV | ENCSR255MIG | ENCSR384KIB |
| ENCSR000ASS | ENCSR191ZQT | ENCSR391NPE | ENCSR264HOB | ENCSR385ODY |
| ENCSR000AUI | ENCSR194KGO | ENCSR397NQK | ENCSR295FXG | ENCSR394DWD |
| ENCSR000AUP | ENCSR195CFR | ENCSR399QWL | ENCSR306YLG | ENCSR394ZWY |
| ENCSR000AUU | ENCSR195JWZ | ENCSR400XSW | ENCSR315NIO | ENCSR399JUL |
| ENCSR000AUW | ENCSR200ETW | ENCSR401KZW | ENCSR316ZPQ | ENCSR403VPZ |
| ENCSR000AVF | ENCSR200JVJ | ENCSR402HFW | ENCSR343GIV | ENCSR405BDN |
| ENCSR000DPL | ENCSR203KCB | ENCSR402JWL | ENCSR344ERJ | ENCSR452GUC |
| ENCSR000DPN | ENCSR203KEU | ENCSR403PZH | ENCSR347KLR | ENCSR454FJM |
| ENCSR000EUT | ENCSR204OJS | ENCSR405ESP | ENCSR355RJN | ENCSR459LGY |
| ENCSR000EUX | ENCSR204TAU | ENCSR405FZE | ENCSR360NJD | ENCSR470SYZ |
| ENCSR000EVA | ENCSR207ABA | ENCSR413QLR | ENCSR363KHN | ENCSR490NCK |
| ENCSR000EWR | ENCSR207FYD | ENCSR413YCI | ENCSR367UDO | ENCSR494YJW |
| ENCSR000EWW | ENCSR208QRN | ENCSR419BNY | ENCSR368NLE | ENCSR524NCH |
| ENCSR000EXK | ENCSR209QGZ | ENCSR424IJM | ENCSR404IRS | ENCSR526EXI |
| ENCSR000EXM | ENCSR210ZPC | ENCSR425FUS | ENCSR405VDU | ENCSR529SWW |
| ENCSR000FCH | ENCSR212FST | ENCSR425PQI | ENCSR420XKX | ENCSR538MSN |
| ENCSR000FCU | ENCSR213SMK | ENCSR426VHO | ENCSR422TGJ | ENCSR541DBO |
| ENCSR000FDA | ENCSR214UZE | ENCSR427DIF | ENCSR437DCH | ENCSR545WFH |
| ENCSR000JGO | ENCSR217SDT | ENCSR429YAE | ENCSR449RSE | ENCSR554RPI |
| ENCSR000NPF | ENCSR220TRW | ENCSR430RVP | ENCSR456XCT | ENCSR559IPI |
| ENCSR001BSB | ENCSR222QLW | ENCSR432DPY | ENCSR461VPV | ENCSR562HZE |
| ENCSR001SHB | ENCSR223UPC | ENCSR435JKM | ENCSR475NOQ | ENCSR563SVE |
| ENCSR002YRE | ENCSR227FYJ | ENCSR436JNB | ENCSR479HSY | ENCSR569NEE |
| ENCSR004EKY | ENCSR230BWN | ENCSR437QMD | ENCSR503BIB | ENCSR572GRM |
| ENCSR004HIE | ENCSR230IMS | ENCSR438KQU | ENCSR528SMB | ENCSR583XAM |
| ENCSR004YQD | ENCSR233RWF | ENCSR438SPO | ENCSR556OGN | ENCSR586IQI |
| ENCSR006ANF | ENCSR235KRX | ENCSR440PMP | ENCSR569IIN | ENCSR589NSH |
| ENCSR007HLH | ENCSR235ZBF | ENCSR443UYU | ENCSR585UEE | ENCSR590YGY |
| ENCSR007WYC | ENCSR239LAP | ENCSR447OHF | ENCSR593INW | ENCSR652EPR |
| ENCSR007YOT | ENCSR240GDT | ENCSR447ZGY | ENCSR593KDJ | ENCSR659WYW |
| ENCSR010OVI | ENCSR242AHB | ENCSR449AUD | ENCSR596PFU | ENCSR688CIB |
| ENCSR011BHU | ENCSR242TBH | ENCSR449AXO | ENCSR597BWL | ENCSR698YQL |
| ENCSR011MGQ | ENCSR245GEV | ENCSR450JZC | ENCSR597RXN | ENCSR707TMM |
| ENCSR012PII | ENCSR247NUF | ENCSR451WXC | ENCSR597UDW | ENCSR720GSE |
| ENCSR013KEC | ENCSR248IFK | ENCSR453MUW | ENCSR597ULV | ENCSR732BRB |
| ENCSR014TDK | ENCSR249IKQ | ENCSR454VRA | ENCSR600TOW | ENCSR745AVM |
| ENCSR015GFK | ENCSR249INE | ENCSR455AFK | ENCSR601VHO | ENCSR750BQT |
| ENCSR016XBE | ENCSR250EHC | ENCSR456AMX | ENCSR603KDJ | ENCSR774MLK |
| ENCSR020OIW | ENCSR250NHD | ENCSR458RRZ | ENCSR604JDV | ENCSR775IMI |
| ENCSR025KPY | ENCSR251PPB | ENCSR461GYN | ENCSR606WJA | ENCSR778NQS |
| ENCSR027BPE | ENCSR262JRZ | ENCSR461UOU | ENCSR607IHY | ENCSR778QHG |
| ENCSR028NXO | ENCSR266XMB | ENCSR463NAD | ENCSR608XIG | ENCSR779ISJ |
| ENCSR028QEA | ENCSR267YXV | ENCSR465UAX | ENCSR612BWE | ENCSR783SNV |
| ENCSR029SIG | ENCSR268BDW | ENCSR473PNT | ENCSR614NPG | ENCSR786VUJ |
| ENCSR034RQV | ENCSR268JQE | ENCSR480NNC | ENCSR615HXA | ENCSR791ISZ |
| ENCSR034ZKE | ENCSR268ZCF | ENCSR480OHP | ENCSR619POC | ENCSR792VLP |
| ENCSR039AHR | ENCSR274VPI | ENCSR489LNU | ENCSR619TOU | ENCSR798RTU |
| ENCSR041UZZ | ENCSR278JAH | ENCSR492PXH | ENCSR620AZM | ENCSR799SRL |
| ENCSR051QLZ | ENCSR279KIX | ENCSR494MDB | ENCSR620XWM | ENCSR799VTP |
| ENCSR054BKO | ENCSR285NFQ | ENCSR494WCX | ENCSR627RBH | ENCSR801IPH |
| ENCSR055XHN | ENCSR291VCJ | ENCSR495HIT | ENCSR632MPN | ENCSR803UDX |
| ENCSR059MVB | ENCSR299EQJ | ENCSR498NGT | ENCSR637ISS | ENCSR804MAP |
| ENCSR061IAQ | ENCSR300ZUS | ENCSR498ZRC | ENCSR638CSE | ENCSR807FIZ |
| ENCSR065SKV | ENCSR303IKJ | ENCSR499IIR | ENCSR639HUV | ENCSR807TBS |
| ENCSR066GUY | ENCSR305ISQ | ENCSR500YBS | ENCSR640XRV | ENCSR807XUB |
| ENCSR067BGS | ENCSR307DQT | ENCSR505OPZ | ENCSR641QPH | ENCSR810EPZ |
| ENCSR067CHO | ENCSR310TYO | ENCSR505YFA | ENCSR641SDI | ENCSR814KSK |
| ENCSR069EGE | ENCSR312HLG | ENCSR506FPQ | ENCSR641YLG | ENCSR814QIV |
| ENCSR069UMW | ENCSR313CEH | ENCSR507SRD | ENCSR642HHF | ENCSR814XCZ |
| ENCSR074ECR | ENCSR314BEX | ENCSR507UDH | ENCSR642LNB | ENCSR820DOR |
| ENCSR078LIZ | ENCSR315IRO | ENCSR510RPC | ENCSR643TBI | ENCSR821DYI |
| ENCSR081OTO | ENCSR315UUP | ENCSR516LQO | ENCSR645MXO | ENCSR822YWX |
| ENCSR082SHT | ENCSR316LNO | ENCSR517NSQ | ENCSR645SYH | ENCSR822ZIG |
| ENCSR086BNR | ENCSR316UPM | ENCSR518COE | ENCSR650OZF | ENCSR826UTD |
| ENCSR086QZY | ENCSR318HUC | ENCSR519CFV | ENCSR655XLM | ENCSR830EXC |
| ENCSR086XCT | ENCSR321KDV | ENCSR520BIM | ENCSR656PGJ | ENCSR835MMN |
| ENCSR091KXI | ENCSR321LKT | ENCSR520KVD | ENCSR656ZEQ | ENCSR835OJV |
| ENCSR094VJC | ENCSR322TJD | ENCSR522MTS | ENCSR659RHV | ENCSR836KPL |
| ENCSR095YMD | ENCSR324JDC | ENCSR524MPL | ENCSR660IQS | ENCSR837CSL |
| ENCSR102GGG | ENCSR325PXS | ENCSR526RGE | ENCSR661KMA | ENCSR837DVF |
| ENCSR102XUM | ENCSR325VCV | ENCSR531AUV | ENCSR664MLM | ENCSR837SGJ |
| ENCSR105CRB | ENCSR326YRW | ENCSR531MLI | ENCSR666TFS | ENCSR838ZAU |
| ENCSR105EMQ | ENCSR327XTS | ENCSR531VWS | ENCSR668DIS | ENCSR839GXY |
| ENCSR105VFR | ENCSR329FXI | ENCSR532SRK | ENCSR668EVA | ENCSR840FCT |
| ENCSR108NVQ | ENCSR330SVG | ENCSR534WQA | ENCSR668GBL | ENCSR841AJO |
| ENCSR111QCU | ENCSR334DRN | ENCSR535GFO | ENCSR672UXO | ENCSR842NGQ |
| ENCSR116MKX | ENCSR339XMR | ENCSR535VLV | ENCSR674VPA | ENCSR842TXG |
| ENCSR119XNK | ENCSR340NAL | ENCSR540ADS | ENCSR677LBH | ENCSR846OTK |
| ENCSR120WKZ | ENCSR340OJR | ENCSR540KQC | ENCSR678LND | ENCSR847AIA |
| ENCSR123HEE | ENCSR340ZTB | ENCSR543CPW | ENCSR679OVD | ENCSR854OXF |
| ENCSR124VOE | ENCSR342KXD | ENCSR543ZVZ | ENCSR680IWU | ENCSR855NCG |
| ENCSR125KVG | ENCSR343DFX | ENCSR546SDM | ENCSR681ELN | ENCSR857GMX |
| ENCSR126ZPV | ENCSR344PHP | ENCSR549SVK | ENCSR681RQP | ENCSR860TEJ |
| ENCSR128HTH | ENCSR346FVK | ENCSR549ZDH | ENCSR683RXF | ENCSR863BVD |
| ENCSR131DVD | ENCSR349VAW | ENCSR550WUX | ENCSR685XVI | ENCSR864KVZ |
| ENCSR133NBJ | ENCSR350EFV | ENCSR554HDT | ENCSR687ZCM | ENCSR864OOO |
| ENCSR136KIM | ENCSR355GNZ | ENCSR555HJT | ENCSR693VHX | ENCSR868ZOR |
| ENCSR136ZQZ | ENCSR355UYP | ENCSR557DFM | ENCSR694RCH | ENCSR869TDT |
| ENCSR137ZID | ENCSR367WYJ | ENCSR560BEL | ENCSR697MUJ | ENCSR872YGQ |
| ENCSR138DOM | ENCSR368CYW | ENCSR561KOM | ENCSR700LWI | ENCSR875OQJ |
| ENCSR138TQF | ENCSR370OQJ | ENCSR561VNL | ENCSR702OVJ | ENCSR875QDS |
| ENCSR143XNJ | ENCSR379DNM | ENCSR561YSH | ENCSR705BTW | ENCSR876RGF |
| ENCSR150QXE | ENCSR379WXM | ENCSR562HZE | ENCSR709ABP | ENCSR880SUY |
| ENCSR152HLZ | ENCSR386UHY | ENCSR563EFJ | ENCSR711URW | ENCSR884SIF |
| ENCSR152MAZ | ENCSR405NKC | ENCSR564IGJ | ENCSR714CKG | ENCSR886ATT |
| ENCSR156WTZ | ENCSR416ESX | ENCSR564WJA | ENCSR714TJD | ENCSR889TEU |
| ENCSR158PLZ | ENCSR417AQU | ENCSR569IBY | ENCSR716BTZ | ENCSR891BTJ |
| ENCSR163VJS | ENCSR425NYG | ENCSR574WQE | ENCSR716FUH | ENCSR891KSP |
| ENCSR168STG | ENCSR438NOO | ENCSR577BCL | ENCSR716XDB | ENCSR891XGQ |
| ENCSR204GMP | ENCSR481QDJ | ENCSR577GVS | ENCSR717HIA | ENCSR892HPQ |
| ENCSR214UPJ | ENCSR481SJM | ENCSR579YLO | ENCSR718BTD | ENCSR892XFG |
| ENCSR215WNN | ENCSR514EZP | ENCSR580WGM | ENCSR719FEJ | ENCSR894MOX |
| ENCSR234FYG | ENCSR519UXW | ENCSR581KPP | ENCSR723FET | ENCSR897RRA |
| ENCSR247BVB | ENCSR529OPX | ENCSR582UTE | ENCSR723RUM | ENCSR899XXQ |
| ENCSR247XFQ | ENCSR534CZJ | ENCSR603SVD | ENCSR724GUS | ENCSR905TYC |
| ENCSR258GSA | ENCSR548SPY | ENCSR613JAK | ENCSR726HTS | ENCSR909UAG |
| ENCSR264DYL | ENCSR557HEX | ENCSR628GIZ | ENCSR726WVB | ENCSR909ZLE |
| ENCSR293ZPH | ENCSR557RDB | ENCSR629GZZ | ENCSR729ENO | ENCSR910PDW |
| ENCSR308JFE | ENCSR564GJK | ENCSR642TGI | ENCSR729GQT | ENCSR912TVO |
| ENCSR373OML | ENCSR568JFZ | ENCSR644XSS | ENCSR734FLK | ENCSR914HOJ |
| ENCSR444GJM | ENCSR583BXC | ENCSR648VOG | ENCSR735SLW | ENCSR917QEH |
| ENCSR462ZQL | ENCSR585FQO | ENCSR650AIN | ENCSR736ALU | ENCSR917XOL |
| ENCSR474WLZ | ENCSR595SMV | ENCSR663GIS | ENCSR736ZEG | ENCSR918ESG |
| ENCSR479BTD | ENCSR597PTH | ENCSR675ZDM | ENCSR738SXD | ENCSR919BBI |
| ENCSR491TAD | ENCSR607UWA | ENCSR687YOP | ENCSR740ROX | ENCSR919WLM |
| ENCSR510OLN | ENCSR608IVH | ENCSR694AEQ | ENCSR741HEF | ENCSR927IOS |
| ENCSR534KYQ | ENCSR617ONQ | ENCSR696HNV | ENCSR741STU | ENCSR928HSI |
| ENCSR539HNK | ENCSR621OKB | ENCSR700CAX | ENCSR742WWB | ENCSR931WLE |
| ENCSR557RJE | ENCSR625BDY | ENCSR704JMY | ENCSR743DDX | ENCSR932QRC |
| ENCSR563OLJ | ENCSR642QDW | ENCSR704QBF | ENCSR747HAM | ENCSR937EVN |
| ENCSR564VBC | ENCSR650EVQ | ENCSR708OCZ | ENCSR748OJA | ENCSR937OPV |
| ENCSR566AUD | ENCSR652ZLE | ENCSR714GKF | ENCSR748TFF | ENCSR938FQG |
| ENCSR570YPW | ENCSR687HYO | ENCSR726FHN | ENCSR751BHO | ENCSR944KAZ |
| ENCSR576VXL | ENCSR687ZYY | ENCSR739QEI | ENCSR752UOD | ENCSR945LPX |
| ENCSR583DLA | ENCSR698TWO | ENCSR744SMV | ENCSR754DRC | ENCSR946LTR |
| ENCSR583MIW | ENCSR728EIP | ENCSR746QZL | ENCSR758KRK | ENCSR948TOS |
| ENCSR592KRF | ENCSR729MRJ | ENCSR749PIF | ENCSR758OEC | ENCSR948YYZ |
| ENCSR617MEB | ENCSR749YQP | ENCSR750DLX | ENCSR763IDK | ENCSR953ULK |
| ENCSR633OJL | ENCSR770VVH | ENCSR770VVH | ENCSR764HZU | ENCSR954IGQ |
| ENCSR634QZP | ENCSR773UBR | ENCSR772LYM | ENCSR768LHG | ENCSR954JMZ |
| ENCSR644CBZ | ENCSR780UWK | ENCSR785LQK | ENCSR769CSE | ENCSR955IXZ |
| ENCSR646FGM | ENCSR793ISY | ENCSR810RRR | ENCSR769FOC | ENCSR955XFL |
| ENCSR651ISB | ENCSR818VNL | ENCSR815BDB | ENCSR770OIC | ENCSR958BKG |
| ENCSR652ULP | ENCSR828HXJ | ENCSR819JOI | ENCSR770XRC | ENCSR958CPV |
| ENCSR683VJG | ENCSR828XQV | ENCSR831CKP | ENCSR771YJT | ENCSR960VLF |
| ENCSR704GTT | ENCSR835ARG | ENCSR831HVP | ENCSR772ABL | ENCSR963LVX |
| ENCSR708TMV | ENCSR839VEU | ENCSR845BPQ | ENCSR773IYZ | ENCSR966IGL |
| ENCSR754WVA | ENCSR844YME | ENCSR869KBA | ENCSR775FTU | ENCSR966RWA |
| ENCSR771IOP | ENCSR845IKY | ENCSR872BBM | ENCSR776KLS | ENCSR969HBF |
| ENCSR792FSD | ENCSR846VGG | ENCSR874WOQ | ENCSR776TJT | ENCSR971ETA |
| ENCSR828DDP | ENCSR849NEM | ENCSR888RBQ | ENCSR777AJC | ENCSR972GGW |
| ENCSR836XKP | ENCSR857YQI | ENCSR889RIF | ENCSR777CQW | ENCSR981UJA |
| ENCSR841CWN | ENCSR859LXC | ENCSR890MJO | ENCSR777VQP | ENCSR982QIF |
| ENCSR843OBM | ENCSR860QVB | ENCSR933UGQ | ENCSR799TWV | ENCSR985CHJ |
| ENCSR854GEV | ENCSR887SZF | ENCSR936CZX | ENCSR803JPO | ENCSR986ZZK |
| ENCSR867SJE | ENCSR891ZHA | ENCSR954ZQC | ENCSR833HBT | ENCSR987FAM |
| ENCSR871SEB | ENCSR895KPB | ENCSR955OYW | ENCSR837RXS | ENCSR987PNT |
| ENCSR913DKN | ENCSR930RIA | ENCSR975PED | ENCSR858ROP | ENCSR989PTS |
| ENCSR922FLS | ENCSR962EJE | ENCSR995BPJ | ENCSR889ZOH | ENCSR993YSN |
| ENCSR961QMI | ENCSR968MKZ |  | ENCSR891YGO |  |
| ENCSR979YKY | ENCSR968XKG |  | ENCSR895CJW |  |
| ENCSR980TBV | ENCSR978XYH |  | ENCSR898WPF |  |
|  | ENCSR979XBX |  | ENCSR906OSH |  |
|  | ENCSR988CHU |  | ENCSR927MXE |  |
|  |  |  | ENCSR932WVK |  |
|  |  |  | ENCSR933UGQ |  |
|  |  |  | ENCSR948DNK |  |
|  |  |  | ENCSR982AYN |  |

### Run nf-core/chipseq pipeline

Modify the `nfcore_chipseq.config` file as needed, then run:

```bash
nextflow run nf-core/chipseq -r 2.1.0 -profile apptainer \
  -c nfcore_chipseq.config \
  --input raw/encode_results/nfcore_chipseq_samplesheet.csv \
  --outdir out/nfcore_chipseq \
  --genome GRCh38 \
  --narrow_peak \
  -resume
```

## Download raw data and reference genomes locally

If you prefer to download the raw data and reference genomes locally:

### Generate local samplesheet only if you want to use local reference genome

```bash
python3 make_nfcode_samplesheet.py
```

### Run nf-core/chipseq with local data

```bash
nextflow run nf-core/chipseq -r 2.1.0 -profile apptainer \
  -c nfcore_chipseq.config \
  --igenomes_base /scratch/$USER/work/dbsuper_nf_new/ref/igenomes \
  --input raw/encode_results/nfcore_chipseq_samplesheet.local.csv \
  --outdir out/nfcore_chipseq \
  --genome GRCh38 \
  --narrow_peak \
  -resume
```

## Run nf-core/enhancerflow

Change into the Enhancerflow directory before launching the pipeline.

Because Nextflow caches assets under `NXF_HOME`, export it to a scratch path first:

```bash
export NXF_HOME=/scratch/$USER/.nextflow

cd /scratch/$USER/work/dbsuper_nf_new/Enhancerflow

nextflow run khan-lab/enhancerflow \
  -r main \
  -c nextflow.config \
  --input enhancerflow_samplesheet.csv \
  --genome hg38 \
  --fasta /scratch/$USER/work/dbsuper_nf_new/ref/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa \
  --gtf /scratch/$USER/work/dbsuper_nf_new/ref/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf \
  --outdir results \
  --skip_motifs true \
  --skip_comparison true \
  -resume
```
