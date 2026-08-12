from pathlib import Path

from all_factories import imputer_factory
from filter_data import _get_dataset_features


HERE = Path(__file__).resolve().parent
DATASETS = HERE.parent.parent / "datasets" / "Validation datasets"

PAPER = "Understanding and Designing a High-Performance Ultrafiltration Membrane Using Machine Learning"
DATASET_NAME = "cleaned_dataset_Ultrafiltration Membrane"
OUTPUT_NAME = "cleaned_dataset_Ultrafiltration Membrane_imputed"

columns_to_impute = [
    "P_MW",
    "surface tension (mN/m)",
    "pore maker molecular weight (Da)",
    "organic compound size (Da)",
    "solubility parameter (MPa1/2)",
]


def impute_and_save_dataset():
    dataset, _, _, _ = _get_dataset_features(DATASETS, PAPER, DATASET_NAME)

    mean_imputer = imputer_factory["mean"]
    dataset[columns_to_impute] = mean_imputer.fit_transform(
        dataset[columns_to_impute]
    )

    output_folder = DATASETS / PAPER
    dataset.to_csv(output_folder / f"{OUTPUT_NAME}.csv", index=False)
    dataset.to_pickle(output_folder / f"{OUTPUT_NAME}.pkl")

    return dataset


if __name__ == "__main__":
    imputed_dataset = impute_and_save_dataset()
    print(f"Saved {OUTPUT_NAME}.csv and {OUTPUT_NAME}.pkl")
