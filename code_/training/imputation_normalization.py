from typing import Callable, Optional, Union, List
import numpy as np
import pandas as pd
from sklearn.preprocessing import FunctionTransformer

from sklearn.compose import ColumnTransformer
# from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import Pipeline
from all_factories import (imputer_factory, 
                           representation_scaling_factory,
                           transforms)



def drop_columns(df:pd.DataFrame, columns_to_drop: List):
    return df.drop(columns=columns_to_drop)


def calculate_mw(df: pd.DataFrame,  # df containing Mw, Mn & PDI columns only
                 mw: str  , mn: str  , pdi: str   # Column names for Mw, Mn & PDI
                 ) -> pd.DataFrame:

  df.loc[df[mw].isna(), mw] = df[pdi] * df[mn]  # df.loc may give you some issues
  return df

def preprocessing_workflow(imputer: Optional[str]=None,
                           feat_to_impute: Optional[list[str]] = None,
                           representation: Optional[str] = None,
                           numerical_feat: Optional[list] = None,
                           structural_feat: Optional[list] = None,
                           special_column: Optional[str] = None,
                           scaler: str = None,
                           ) -> Pipeline:
    # structure_feats and scalar_feats for scaling
    # all_columns = list(set(columns_to_impute)|set(special_column))
    # imputation of columns
    steps = []
    # imputation
    if imputer:
        steps = [
            ("Impute feats", ColumnTransformer(
                transformers=[
                    (f'imputer_{imputer}', imputer_factory[imputer], feat_to_impute),
                ],
                remainder="passthrough", verbose_feature_names_out=False
            )),
            
                ]
    
    if special_column:
        steps.extend([
            ("Calculate Mw", ColumnTransformer(
                transformers=[
                    (f'calculator_{special_column}',
                     (FunctionTransformer(calculate_mw,
                                          kw_args={'mw': special_column, 'mn': 'Mn (g/mol)', 'pdi': 'PDI'},
                                          validate=False)), ['Mn (g/mol)', special_column, 'PDI'])
                ],
                remainder="passthrough", verbose_feature_names_out=False
            )),

            ("Impute the rest of Mw values", ColumnTransformer(
                transformers=[
                    (f"imputer_{special_column}", imputer_factory[imputer], [special_column])
                ],
                remainder="passthrough", verbose_feature_names_out=False
            )),
            ("drop Mn", ColumnTransformer(
                transformers=[
                    ("Drop Mn column",
                     FunctionTransformer(drop_columns, kw_args={'columns_to_drop': ['Mn (g/mol)']}, validate=False),
                     ['Mn (g/mol)'])
                ], remainder="passthrough", verbose_feature_names_out=False))

        ])
    # Normalization
    if numerical_feat is not None and 'Mn (g/mol)' in numerical_feat:
        numerical_feat.remove('Mn (g/mol)')
        
    if scaler:
        transformers = []
        if numerical_feat:

            transformers.append(
                ("numerical_scaling", transforms[scaler], numerical_feat)
                )
        # elif representation_scaling_factory[representation]['callable']:
        elif representation:
            transformers.append(
                ("structural_scaling", transforms[scaler], structural_feat)
                )
            
        scaling = ("scaling features",
                ColumnTransformer(transformers=[*transformers], remainder="passthrough", verbose_feature_names_out=False)
                )
        steps.append(scaling)
    return Pipeline(steps)