"""Centralized paths and physical hardware constants for the optimization task.

Example:
    >>> import polykernels
    >>> _P = polykernels.paths_and_constants.Paths()
    >>> _C = polykernels.paths_and_constants.Constants()
    >>> print(_P.data_dir)
"""

from pathlib import Path


class Paths:
    """Centralized project directories and file paths."""

    root: Path = Path(__file__).resolve().parents[1]
    data_dir: Path = root / "data"
    models_dir: Path = root / "models"
    results_dir: Path = root / "results"
    data_csv: Path = data_dir / "processed" / "gs-lf_combined.csv"

    #DATASETS OF INTEREST--------------
    dataset1_csv: Path = data_dir / "processed" / "Beyond molecular structure_seifrid_imputed.csv"
    dataset2_1_csv: Path = data_dir / "processed" / "flux_data_imputed.csv"
    dataset2_2_csv: Path = data_dir / "processed" / "separation_data_imputed.csv"
    dataset3_csv: Path = data_dir / "processed" / "cleaned_dataset_Ultrafiltration Membrane.csv"
    dataset4_csv: Path = data_dir / "processed" / "cleaned_dataset_pervaporation_membranes_wang.csv"
    dataset5_csv: Path = data_dir / "processed" / "Rg data with clusters aging imputed.csv"
    dataset6_csv: Path = data_dir / "processed" / "cleaned_suzuki_synthesis.csv"
    #---------------------------------

    features_npz: Path = data_dir / "processed" / "features.npz"
    splits_tvt_npz: Path = data_dir / "processed" / "splits_tvt.npz"
    splits_5cv_npz: Path = data_dir / "processed" / "splits_5cv.npz"

    # Standard model training artifact filenames
    predictions_file: str = "predictions.npz"
    hparams_file: str = "hparams.yaml"
    stats_file: str = "stats.json"
    history_file: str = "history.csv"
    model_file: str = "model.safetensors"
    metrics_file: str = "metrics.csv"
    production_model_file: str = "production_model.safetensors"
    calibration_file: str = "calibration_weights.safetensors"
    production_metrics_file: str = "production_metrics.csv"

    # Centralized filenames for Optuna studies (duplicated/shared across files)
    optuna_db_file: str = "optuna_study.db"
    study_summary_file: str = "study_summary.csv"


class Constants:
    """Global constants used across the project."""

    smiles_col: str = "IsomericSMILES" #obsolete
    descriptors_col: str = "descriptors" #obsolete
    tvt_split_key: str = "tvt_split"


    # Dataset-specific SMILES column names
    dataset1_donor_smiles: str = "Donor SMILES"
    dataset1_acceptor_smiles: str = "Acceptor SMILES"
    dataset2_smiles : str = "smile"
    dataset3_smiles : str = "SMILES"
    dataset4_smiles : str = "SMILES"
    dataset5_smiles : str = "Monomer SMILES"
    # Reaction miniaturization is a small-molecule HTE dataset, not a polymer one: the
    # structural entity is the reaction substrate rather than a repeat unit.
    dataset6_smiles : str = "Substrate SMILES"
   

    # Dataset-specific feature and target column names
    x_cols = {
    1: ['HOMO_D (eV)', 'LUMO_D (eV)', 'Eg_D (eV)', 'Ehl_D (eV)',
            'HOMO_A (eV)', 'LUMO_A (eV)', 'Eg_A (eV)', 'Ehl_A (eV)',
            "D:A ratio (m/m)", "solvent additive conc. (% v/v)",
            "temperature of thermal annealing",
            "HTL energy level (eV)", "ETL energy level (eV)"],
    2: [
        'Selective_layer_thickness',
        'Filler_size',
        'Filler_concentration',
        'Crosslinker_concentration',
        'Organic_Hildebrand_solubility_parameter',
        'Organic_molecular_volume',
        'Effective_area',
        'Mass_ratio(A/B)',
        'Experimental_temperature',
        'Downstream_pressure',
        ],
    3: [
        'P_MW',
        'surface tension (mN/m)',
        'P_wt%',
        'pore maker molecular weight (Da)',
        'pore maker wt% (pore maker to total weight)',
        'additive wt% (additive to polymer)',
        'solubility parameter (MPa1/2)',
        'TMP (bar)',
        'oraganic compound concentration (mg/L)',
        'organic compound size (Da)',
        'foulant concentration (mg/L)',
        ],
    4:
    [
    'Contact angle θ （°）',
    'Thickness ℓ (um)',
    'xw (wt%)',
    'Temperature (°C)',
    'Permeate pressure (mbar)',
    'Solvent solubility parameter （MPa1/2）',
    ],
    5:
    [
    "Xn", "Mw (g/mol)", "PDI", "Concentration (mg/ml)",
    "Temperature SANS/SLS/DLS/SEC (K)",
    "polymer dP", "polymer dD", "polymer dH",
    "solvent dP", "solvent dD", "solvent dH"
    ],
    6:
    [
    "Equiv AcOH", "Equiv TTIP", "Solvent_Encoded", "Reaction Conc (mM)"
    ]
    }

    y_cols = {
        1: ['calculated PCE (%)'],
        2: ['log (Total flux)', 'log (Separation factor)'], #there are 2 prediction targets
        3: [
            'water permeability (LMH/bar)',
            'organic compound removal (%)',
            'flux decline ratio (%)',
            'flux recovery ratio (%)',
            'reversible fouling ratio (%)',
            'irreversible fouling ratio(%)',
            ],
        4:  ['log (Total flux)', 'log (Separation factor)'],
        5:  ["log Rg (nm)"],
        6:  ["Approx Conv (%)"]
    }


    

    # fmt: off
    labels: list[str] = [
        "alcoholic", "aldehydic", "alliaceous", "almond", "amber", "animal",
        "anisic", "apple", "apricot", "aromatic", "balsamic", "banana", "beefy",
        "bergamot", "berry", "bitter", "black currant", "brandy", "burnt",
        "buttery", "cabbage", "camphoreous", "caramellic", "cedar", "celery",
        "chamomile", "cheesy", "cherry", "chocolate", "cinnamon", "citrus",
        "clean", "clove", "cocoa", "coconut", "coffee", "cognac", "cooked",
        "cooling", "cortex", "coumarinic", "creamy", "cucumber", "dairy", "dry",
        "earthy", "ethereal", "fatty", "fermented", "fishy", "floral", "fresh",
        "fruit skin", "fruity", "garlic", "gassy", "geranium", "grape",
        "grapefruit", "grassy", "green", "hawthorn", "hay", "hazelnut",
        "herbal", "honey", "hyacinth", "jasmin", "juicy", "ketonic", "lactonic",
        "lavender", "leafy", "leathery", "lemon", "lily", "malty", "meaty",
        "medicinal", "melon", "metallic", "milky", "mint", "muguet", "mushroom",
        "musk", "musty", "natural", "nutty", "odorless", "oily", "onion",
        "orange", "orangeflower", "orris", "ozone", "peach", "pear", "phenolic",
        "pine", "pineapple", "plum", "popcorn", "potato", "powdery", "pungent",
        "radish", "raspberry", "ripe", "roasted", "rose", "rummy", "sandalwood",
        "savory", "sharp", "smoky", "soapy", "solvent", "sour", "spicy",
        "strawberry", "sulfurous", "sweaty", "sweet", "tea", "terpenic",
        "tobacco", "tomato", "tropical", "vanilla", "vegetable", "vetiver",
        "violet", "warm", "waxy", "weedy", "winey", "woody",
    ]
    # fmt: on
    num_labels: int = len(labels)  # 138
