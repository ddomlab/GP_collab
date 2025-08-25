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


def preprocessing_workflow(
        numerical_feat: Optional[list] = None,
        structural_feat: Optional[list] = None,
        scaler: str = None,
    ) -> Pipeline:
    
    steps = []

    # Remove Mn from numerical feats if present
    if numerical_feat is not None and 'Mn (g/mol)' in numerical_feat:
        numerical_feat.remove('Mn (g/mol)')

    # Apply scaling only if NO representation
    if scaler:
        transformers = []
        if numerical_feat:
            transformers.append(
                ("numerical_scaling", transforms[scaler], numerical_feat)
            )
        
        if structural_feat:
            transformers.append(
                ("structural_passthrough", "passthrough", structural_feat)
            )

        scaling = (
            "scaling features",
            ColumnTransformer(transformers=transformers, remainder="passthrough",
                              verbose_feature_names_out=False)
        )
        steps.append(scaling)

    return Pipeline(steps)
