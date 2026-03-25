# MOLECULAR_DOCKING

Pipeline de docking molecular con DOCK6. Produce outputs compatibles con **dock2profile**.

## Requisitos

| Dependencia | Instalación | Notas |
|---|---|---|
| Python 3.9-3.12 | conda | |
| RDKit, OpenBabel, AmberTools | conda | Via `environment.yaml` |
| PDB2PQR | pip | Protonación pH-aware del receptor |
| **DOCK6** | Manual | Licencia académica gratuita de [UCSF](https://dock.compbio.ucsf.edu/DOCK_6/index.htm) |
| **ChimeraX** | Manual | [Descargar](https://www.cgl.ucsf.edu/chimerax/download.html) — necesario para preparar el receptor |

## Instalación

```bash
# 1. Clonar
git clone https://github.com/benjaminprieto/reference_docking.git
cd reference_docking

# 2. Crear entorno conda
conda env create -f environment.yaml
conda activate reference_docking_env

# 3. Instalar paquete en modo editable
pip install -e ".[dev]"

# 4. Verificar dependencias
bash check_dependencies.sh
```

### DOCK6

DOCK6 requiere licencia académica gratuita. Después de obtenerla:

```bash
# Compilar e instalar en /opt/dock6/
tar -xzf dock.6.X.tar.gz
cd dock6
./configure gnu
make

# Agregar al PATH (en ~/.bashrc)
export PATH=/opt/dock6/bin:$PATH
```

### ChimeraX

Descargar desde https://www.cgl.ucsf.edu/chimerax/download.html e instalar. El pipeline busca el binario en `/usr/bin/chimerax-daily`, luego `chimerax` en PATH.

## Estructura del proyecto

```
reference_docking/
├── 01_src/reference_docking/       Core modules (lógica, sin CLI)
│   ├── m00_preparation/
│   │   ├── molecule_parser.py          00a — parseo de moléculas
│   │   ├── receptor_preparation.py     00b — receptor → mol2 DOCK6-ready
│   │   ├── ionization_profiling.py     00c — protonación de ligandos al pH
│   │   └── binding_site_definition.py  00d — recorte del receptor (opcional)
│   ├── m01_docking/                    DOCK6 engine
│   │   ├── antechamber_preparation.py  01a — mol2 con cargas AM1-BCC
│   │   ├── grid_generation.py          01b — DMS → spheres → grids
│   │   ├── dock6_runner.py             01c — dock6 por molécula
│   │   └── score_collector.py          01d — parseo de scores → Excel
│   └── m02_vina/                       Vina engine
│       ├── vina_preparation.py         02a — receptor/ligandos → PDBQT
│       ├── vina_runner.py              02b — Vina/Vina-GPU docking
│       └── vina_score_collector.py     02c — scores → CSV/Excel
├── 02_scripts/                     CLI scripts (argparse + YAML → core)
├── 03_configs/                     YAML por módulo (parámetros algorítmicos)
├── 04_data/campaigns/              Campañas (receptor + moléculas + grids)
│   └── example_campaign/
│       ├── campaign_config.yaml        Fuente de verdad de la campaña
│       ├── receptor/                   PDB del receptor
│       └── molecules/                  Moléculas de entrada
├── 05_results/                     Outputs por campaña/módulo
├── environment.yaml
├── pyproject.toml
└── check_dependencies.sh
```

## Flujo del pipeline

```
Shared preparation:
00a molecule_parser       → unique_molecules.csv + .sdf
00b receptor_preparation  → rec_charged.mol2 + rec_noH.pdb
00c ionization_profiling  → SDF protonados por pH
00d binding_site_def      → rec_noH_site.pdb (opcional)

DOCK6 engine:
01a antechamber           → mol2 con AM1-BCC charges
01b grid_generation       → DMS, spheres, box, grid.nrg/bmp
01c dock6_run             → scored mol2 por molécula
01d score_collection      → Excel compatible dock2profile

Vina engine:
02a vina_preparation      → receptor.pdbqt + ligandos PDBQT
02b vina_runner           → Vina/Vina-GPU docking
02c vina_score_collector  → scores CSV/Excel
```

## Uso

### 1. Crear campaña

```bash
cp -r 04_data/campaigns/example_campaign 04_data/campaigns/mi_campana
```

Editar `campaign_config.yaml` con el receptor, moléculas, y pH de docking.

### 2. Correr pipeline

Desde la raíz del proyecto (o como Run Configurations en PyCharm):

```bash
# === Shared preparation ===

# 00a — Parsear moléculas
python 02_scripts/00a_molecule_parser.py --config 03_configs/00a_molecule_parser.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 00b — Preparar receptor (ChimeraX + PDB2PQR)
python 02_scripts/00b_receptor_preparation.py --config 03_configs/00b_receptor_preparation.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 00c — Protonar ligandos
python 02_scripts/00c_ionization_profiling.py --config 03_configs/00c_ionization_profiling.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 00d — Binding site definition (recorte del receptor)
python 02_scripts/00d_binding_site_definition.py --config 03_configs/00d_binding_site_definition.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# === DOCK6 engine ===

# 01a — Antechamber (AM1-BCC charges)
python 02_scripts/01a_antechamber_preparation.py --config 03_configs/01a_antechamber_preparation.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 01b — Generar grids
python 02_scripts/01b_grid_generation.py --config 03_configs/01b_grid_generation.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 01c — DOCK6 docking
python 02_scripts/01c_dock6_run.py --config 03_configs/01c_dock6_run.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 01d — DOCK6 score collection
python 02_scripts/01d_score_collection.py --config 03_configs/01d_score_collection.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# === Vina engine ===

# 02a — Vina preparation (PDBQT conversion + binding box)
python 02_scripts/02a_vina_preparation.py --config 03_configs/02a_vina_preparation.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 02b — Vina docking
python 02_scripts/02b_vina_runner.py --config 03_configs/02b_vina_runner.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# 02c — Vina score collection
python 02_scripts/02c_vina_score_collector.py --config 03_configs/02c_vina_score_collector.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml

# === Or run everything at once ===
bash run_pipeline.sh mi_campana              # both engines
bash run_pipeline.sh mi_campana dock6        # DOCK6 only
bash run_pipeline.sh mi_campana vina         # Vina only
```

### PyCharm Run Configuration

Para cada módulo crear un Run Configuration:

| Campo | Valor |
|---|---|
| Script | `02_scripts/00b_receptor_preparation.py` |
| Parameters | `--config 03_configs/00b_receptor_preparation.yaml --campaign 04_data/campaigns/mi_campana/campaign_config.yaml` |
| Working directory | Raíz del proyecto |
| Python interpreter | `reference_docking_env` |

## Configuración

Cada campaña se define en `campaign_config.yaml`. Los parámetros clave:

```yaml
campaign_id: "mi_campana"
docking_ph: 7.2

receptor:
  pdb: "receptor/mi_receptor.pdb"
  protonation:
    enabled: true
    tool: "pdb2pqr"       # pdb2pqr | chimerax | obabel

molecules:
  input_file: "molecules/"
  protonation:
    tool: "obabel"        # obabel | dimorphite_dl

grids:
  generate: true
  binding_site:
    method: "reference_ligand"
    reference_mol2: "molecules/ligando_cristalografico.mol2"
    radius: 10.0
```

## Notas técnicas

**DOCK6 y el límite de 80 caracteres.** Los programas Fortran de DOCK6 (sphgen, showbox) truncan paths a ~80 chars. El pipeline usa symlinks y filenames cortos automáticamente — no necesitas preocuparte de esto.

**Receptor preparation.** El módulo 00b usa ChimeraX para generar el mol2 del receptor con Sybyl atom types y AMBER ff14SB charges. Si usas la estrategia `pdb2pqr`, PDB2PQR+PROPKA predicen pKa por residuo antes de que ChimeraX asigne las cargas.

**Coordenadas cristalográficas.** Si tienes un ligando co-cristalizado, asegúrate de usar las coordenadas extraídas del PDB (no regeneradas desde SMILES) para `binding_site.reference_mol2`.