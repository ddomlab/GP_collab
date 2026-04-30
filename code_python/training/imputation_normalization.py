from typing import Callable, Optional, Union, List
import numpy as np
import pandas as pd
from sklearn.preprocessing import FunctionTransformer

from sklearn.compose import ColumnTransformer, make_column_selector
# from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import Pipeline
from all_factories import (imputer_factory, 
                           representation_scaling_factory,
                           transforms)

import numbers


class make_array_column_selector:
    """
    Like sklearn's make_column_selector, but for 2D numpy object arrays.
    Selects columns whose first cell is an instance of any type in `dtype_include`
    and not an instance of any type in `dtype_exclude`.

    Parameters
    ----------
    dtype_include : type or tuple of types, default=None
    dtype_exclude : type or tuple of types, default=None
    """

    def __init__(self, *, dtype_include=None, dtype_exclude=None):
        self.dtype_include = dtype_include
        self.dtype_exclude = dtype_exclude

    @staticmethod
    def _as_tuple(t):
        if t is None:
            return None
        return tuple(t) if isinstance(t, (list, tuple)) else (t,)

    def __call__(self, X):
        X = np.asarray(X)
        if X.ndim != 2:
            raise ValueError(f"Expected 2D array, got {X.ndim}D.")

        include = self._as_tuple(self.dtype_include)
        exclude = self._as_tuple(self.dtype_exclude)

        selected = []
        for j in range(X.shape[1]):
            sample = X[0, j]            # probe first row
            if include is not None and not isinstance(sample, include):
                continue
            if exclude is not None and isinstance(sample, exclude):
                continue
            selected.append(j)
        return selected


def drop_columns(df:pd.DataFrame, columns_to_drop: List):
    return df.drop(columns=columns_to_drop)


def calculate_mw(df: pd.DataFrame,  # df containing Mw, Mn & PDI columns only
                 mw: str  , mn: str  , pdi: str   # Column names for Mw, Mn & PDI
                 ) -> pd.DataFrame:

  df.loc[df[mw].isna(), mw] = df[pdi] * df[mn]  # df.loc may give you some issues
  return df


def preprocessing_workflow(imputer: Optional[str]=None,
                           regressor_type: Optional[str]=None,
                           feat_to_impute: Optional[list[str]] = None,
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
    if "mgk" in regressor_type.lower():
        from mgktools.graph.hashgraph import HashGraph
        if scaler:
            scaling = ("scaling features",
                    ColumnTransformer(
                        transformers=[
                            ("pass", "passthrough", make_array_column_selector(dtype_include=HashGraph)),
                            ("scaling features", transforms[scaler], make_array_column_selector(dtype_include=numbers.Real))
                        ],
                        # remainder="passthrough",
                        verbose_feature_names_out=False)
                    )
            steps.append(scaling)
            print("MGK regressor, only scaling numerical features")
              
    else:
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
            if structural_feat:
                transformers.append(
                ("structural_scaling", "passthrough", structural_feat)  # graphs stay at front
                )
                
            if numerical_feat:
                transformers.append(
                    ("numerical_scaling", transforms[scaler], numerical_feat)
                )
            # elif representation_scaling_factory[representation]['callable']:

                
            scaling = ("scaling features",
                    ColumnTransformer(transformers=[*transformers],
                                    remainder="passthrough",
                                    verbose_feature_names_out=False)
                    )
            steps.append(scaling)
    return Pipeline(steps)