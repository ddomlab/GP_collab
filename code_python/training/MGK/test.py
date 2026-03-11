### help
import pandas as pd
from data_handler.data import Dataset

df = pd.DataFrame({
    'smiles': ['CCO', 'CCC', 'CCCC', 'CCCCC', 'CCCCCC'],
    'target': [1.0, 2.0, 3.0, 4.0, 5.0]
})

dataset = Dataset.from_df(
    df,
    smiles_columns=['smiles'],
    targets_columns=['target']
)
dataset.set_status(graph_kernel_type='graph')
dataset.create_graphs(n_jobs=4)
