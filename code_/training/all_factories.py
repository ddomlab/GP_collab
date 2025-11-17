from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor,RandomForestClassifier
from xgboost import XGBRegressor, XGBClassifier
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.svm import SVR
from ngboost import NGBRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import Lasso
# from GPR_model import GPRegressor
from GPytorch_kernel_mix import GPMixRegressor
from sklearn.neural_network import MLPRegressor
# from sklearn.multioutput import MultiOutputRegressor
from sklearn.preprocessing import FunctionTransformer
import numpy as np
from sklearn.preprocessing import (StandardScaler,
                                   QuantileTransformer,
                                   MinMaxScaler,
                                #    FunctionTransformer,
                                   RobustScaler)
from sklearn.base import TransformerMixin
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer

from typing import Callable, Optional, Union, Dict
from types import NoneType
from skopt.space import Integer, Real, Categorical

from sklearn.gaussian_process import GaussianProcessRegressor
# from sklearn.gaussian_process.kernels import PairwiseKernel
from sklearn.gaussian_process.kernels import RBF, Matern, RationalQuadratic, ConstantKernel as C
from Sklearn_kernel_mix import (JaccardKernel, custom_RBF,
                            custom_Matern, AdditiveRBF,
                            AdditiveMatern, ProductKernel,
                            AddKernel, weighted_jaccard)

cutoffs = {
            # "Rh1 (nm)":1000,
            # "Rg1 (nm)":1000,
            # "Lc (nm)":1000,
            # "Lp (nm)":(None,100),
            # "Concentration (mg/ml)":(None,50)
        }


unrolling_feature_factory: dict[str, list[str]] = {
                                                "polysize": ['Mw (g/mol)', 'PDI'],
                                                "HSPs": ['polymer dP', 'polymer dD' , 'polymer dH', 'solvent dP', 'solvent dD', 'solvent dH'],
                                                "Ra": ["Ra"],
                                                "solvent properties": ['Concentration (mg/ml)', 'Temperature SANS/SLS/DLS/SEC (K)'],
                                                 }


def generate_acronym_string(feats):
    acronym_list = []
    for key, values in unrolling_feature_factory.items():
        if any(feat in feats for feat in values):
            acronym_list.append(key)
    return "_".join(acronym_list)


imputer_factory: Dict[str, TransformerMixin] = {
    "mean": SimpleImputer(strategy="mean"),
    "median": SimpleImputer(strategy="median"),
    "most-frequent": SimpleImputer(strategy="most_frequent"),
    "uniform KNN": KNNImputer(weights="uniform"),
    "distance KNN": KNNImputer(weights="distance"),
    "iterative": IterativeImputer(sample_posterior=True),
}

def inverse_log_transform(X):
    return np.power(10, X)

transforms: dict[str, Callable] = {
    None:                None,
    "Log":                  FunctionTransformer(np.log10, inverse_func=inverse_log_transform,
                                                    check_inverse=True, validate=False),
    "MinMax":               MinMaxScaler(),
    "Standard":             StandardScaler(),
    "RobustScaler":        RobustScaler(),
    "Uniform Quantile":     QuantileTransformer(),
}


radius_to_bits: dict[int, int] = {3: 512, 4: 1024, 5: 2048, 6: 4096}



representation_scaling_factory: dict[str, dict[str, Union[Callable, str]]] = {
    "ECFP":                {"callable": None, "type": None},
    "MACCS":                {"callable": None, "type": None},
    "Mordred":             {"callable": StandardScaler(),
                            "type":     "Standard"},
    # "graph embeddings":    {"callable": MinMaxScaler,
    #                         "type": "MinMax"},
    # "BRICS":               {"callable": MinMaxScaler, "type": "MinMax"},
    # "SELFIES":             {"callable": MinMaxScaler, "type": "MinMax"},
    # "SMILES":              {"callable": MinMaxScaler, "type": "MinMax"},
    # "OHE":                 {"callable": MinMaxScaler, "type": "MinMax"},
    # "GNN":     {"callable": MinMaxScaler, "type": "MinMax"},
}


regressor_factory: dict[str, type]={
    "MLR": LinearRegression(),
    "KNN": KNeighborsRegressor(),
    "SVR": SVR(),
    "XGBR": XGBRegressor(),
    "RF": RandomForestRegressor(),
    "Lasso": Lasso(),
    "DT": DecisionTreeRegressor(),
    "NGB": NGBRegressor(),
    "XGBC":XGBClassifier(),
    "RFC": RandomForestClassifier(),
    # "GPR": GPRegressor,
    "sklearn-GPR":GaussianProcessRegressor(),
    "MLP": MLPRegressor(),
    "HGBR": HistGradientBoostingRegressor(),
    "GPMixR": GPMixRegressor(),
}

def optimized_models(model_name:str,random_state:int=42, **kwargs):
    if 'NGB'==model_name:
        return NGBRegressor(n_estimators=500, learning_rate=0.01, tol=1e-4,
                             random_state=None, verbose=False,**kwargs)
    if 'XGBR'==model_name:
        return  XGBRegressor(eval_metric="rmse", 
                               random_state=None, n_jobs=-1,**kwargs)
    
    if 'RF'==model_name:
        return RandomForestRegressor(n_estimators=100, max_depth=None, 
                                        random_state=None, n_jobs=-1,
                                        max_features="sqrt"
                                        )
    if 'HGBR'==model_name:
        return HistGradientBoostingRegressor(max_iter=2000, max_depth=None, 
                                            #  min_samples_leaf=20, max_leaf_nodes=1000,
                                             learning_rate=0.01, l2_regularization=1e-6,
                                             scoring='neg_root_mean_squared_error', random_state=None,**kwargs)
    if 'MLP'==model_name:
        return MLPRegressor(max_iter=200)

    if 'sklearn-GPR'==model_name:
        # kernel = kwargs.get('kernel')
        return GaussianProcessRegressor(random_state=random_state, **kwargs)
    
    if "GPMixR" == model_name:
        return GPMixRegressor(**kwargs)
    
    if "GPMixMCMC" == model_name:
        return GPMixMCMCRegressor(**kwargs)
    
    return None


kernel_space = {
    'matern': custom_Matern,
    'rbf': custom_RBF,
    'additive_rbf': AdditiveRBF, #just for count
    'additive_matern': AdditiveMatern, #just for count
    'product': ProductKernel,
    'add': AddKernel,
    'jaccard': JaccardKernel
}


def _parse_kernel(kernel_name: str):
    """Return (base_name, kwargs) — e.g. ('matern', {'nu': 2.5})"""
    if kernel_name is None:
        return None, {}
    # Detect Matern with nu suffix
    if kernel_name.lower().startswith("matern"):
        # Try to extract numeric part after 'matern'
        try:
            nu_val = float(kernel_name.lower().replace("matern", ""))
            return "matern", {"nu": nu_val}
        except ValueError:
            return "matern", {}
    return kernel_name, {}


def _get_gp_kernel(gp_info: dict, idx: Optional[dict[str, list[int]]] = None):
    """
    Build GP kernel based on gp_info and feature index dictionary.
    """

    # --- Mixing is OFF ---
    if gp_info['mixing'] is None:
        active_key = next((k for k, v in gp_info['kernel'].items() if v is not None), None)
        if active_key is None:
            raise ValueError("No active kernel found when mixing is off.")

        if active_key == 'fp':
            length_scale = np.full(len(idx['fp']), 1.0) if idx and idx.get('fp') else 1.0
            active_kernel_name, extra_args = _parse_kernel(gp_info['kernel']['fP'])
            kernel_class = kernel_space[active_kernel_name]
            kernel = kernel_class(length_scale=length_scale, metric=weighted_jaccard, **extra_args)
        else:
            length_scale = np.full(len(idx['count']), 1.0) if idx and idx.get('count') else 1.0
            active_kernel_name, extra_args = _parse_kernel(gp_info['kernel']['count'])
            kernel_class = kernel_space[active_kernel_name]
            kernel = kernel_class(length_scale=length_scale, **extra_args)

    # --- Additive Mixing ---
    elif gp_info['mixing'] == 'additive':
        length_scale_fp = np.full(len(idx['fp']), 1.0) if idx and idx.get('fp') else 1.0
        length_scale_count = np.full(len(idx['count']), 1.0) if idx and idx.get('count') else 1.0

        fp_name, fp_args = _parse_kernel(gp_info['kernel']['fP'])
        count_name, count_args = _parse_kernel(gp_info['kernel']['count'])

        kernel = AddKernel(
            kernel_cont=kernel_space[count_name](length_scale=length_scale_count, **count_args),
            kernel_fp=kernel_space[fp_name](length_scale=length_scale_fp, metric=weighted_jaccard, **fp_args),
            cont_idx=idx['count'],
            fp_idx=idx['fp']
        )

    # --- Product Mixing ---
    elif gp_info['mixing'] == 'product':
        length_scale_fp = np.full(len(idx['fp']), 1.0) if idx and idx.get('fp') else 1.0
        length_scale_count = np.full(len(idx['count']), 1.0) if idx and idx.get('count') else 1.0

        fp_name, fp_args = _parse_kernel(gp_info['kernel']['fP'])
        count_name, count_args = _parse_kernel(gp_info['kernel']['count'])

        kernel = ProductKernel(
            kernel_cont=kernel_space[count_name](length_scale=length_scale_count, **count_args),
            kernel_fp=kernel_space[fp_name](length_scale=length_scale_fp, metric=weighted_jaccard, **fp_args),
            cont_idx=idx['count'],
            fp_idx=idx['fp']
        )

    else:
        raise ValueError(f"Unknown mixing method: {gp_info['mixing']}")

    return kernel






def get_regressor_search_space(algortihm:str, kernel:str=None) -> Dict :
    if algortihm == "MLR":
        return {
            "regressor__regressor__fit_intercept": [True, False]
    }
    
    if algortihm == "Lasso":
        return {
        "regressor__alpha": Real(1e-3, 1e3, prior="log-uniform"),
        "regressor__fit_intercept": [True, False],
        "regressor__selection": Categorical(["cyclic", "random"]),
    }

    if algortihm == "KNN":
        return {
        "regressor__n_neighbors": Integer(1, 50),
        "regressor__weights": Categorical(["uniform", "distance"]),
        "regressor__algorithm": Categorical(["ball_tree", "kd_tree", "brute"]),
        "regressor__leaf_size": Integer(1, 100),
    }

    if algortihm == "SVR":
        return {
        "regressor__kernel": Categorical(["linear", "rbf"]),
    }

    if algortihm == "MLP":
        return {
                    # "regressor__regressor__hidden_layer_sizes": Categorical([
                    # (64), (128), (256), (512)
                    # ]),
                    "regressor__regressor__learning_rate_init": Real(1e-5, 1e-1, prior="log-uniform"),
                    "regressor__regressor__alpha": Real(1e-8, 1e-6, prior="log-uniform"),  # This is weight decay
                    "regressor__regressor__batch_size": Categorical([32, 64, 128, 256]),
                    "regressor__regressor__max_iter": [200],
                }

    if algortihm == "RF":
        return {
        "regressor__regressor__n_estimators": Integer(10, 2000, prior="log-uniform"),
        "regressor__regressor__max_depth": [None],
        "regressor__regressor__min_samples_split": Real(0.001, 0.99),
        "regressor__regressor__min_samples_leaf": Real(0.001, 0.99),
        "regressor__regressor__max_features": Categorical(["sqrt", "log2", None]),
    }

    if algortihm == "XGBR":
        return {
        "regressor__regressor__n_estimators": Integer(50, 2000, prior="log-uniform"),
        "regressor__regressor__max_depth": Integer(10, 10000, prior="log-uniform"),
        # "regressor__grow_policy": Categorical(["depthwise", "lossguide"]),
        "regressor__regressor__n_jobs": [-2],
        "regressor__regressor__learning_rate": Real(1e-4, 1e-1, prior="log-uniform"),
    }


    if algortihm == "XGBC":
        return {
        "regressor__n_estimators": Integer(50, 2000, prior="log-uniform"),
        "regressor__max_depth": Integer(3, 10000, prior="log-uniform"),
        "regressor__learning_rate": Real(1e-3, 1e-1, prior="log-uniform"),
        # "classifier__classifier__subsample": Real(0.5, 1.0, prior="uniform"),
        # "classifier__classifier__colsample_bytree": Real(0.5, 1.0, prior="uniform"),
        # "classifier__classifier__gamma": Real(1e-8, 1.0, prior="log-uniform"),
        # "classifier__classifier__min_child_weight": Integer(1, 10),
        # "classifier__classifier__reg_alpha": Real(1e-5, 1.0, prior="log-uniform"),
        # "classifier__classifier__reg_lambda": Real(1e-5, 1.0, prior="log-uniform"),
        "regressor__n_jobs": [-2],
        "regressor__objective":['binary:logistic']
    }


    if algortihm == "RFC":
        return {
        "regressor__n_estimators": Integer(10, 2000, prior="log-uniform"),
        "regressor__max_depth": Integer(2, 10000, prior="log-uniform"),
    }


    if algortihm == "DT":
        return {
        "regressor__regressor__min_samples_split": Real(0.05, 0.99),
        "regressor__regressor__min_samples_leaf": Real(0.05, 0.99),
        "regressor__regressor__max_features": Categorical([None,"sqrt", "log2"]),
        "regressor__regressor__max_depth": [None],
        "regressor__regressor__ccp_alpha": Real(0.05, 0.99),
    }

    if algortihm == "NGB":
        return {
        "regressor__regressor__n_estimators": Integer(50, 2000, prior="log-uniform"),
        "regressor__regressor__learning_rate": Real(1e-4, 1e-1, prior="log-uniform"),
        # "regressor__regressor__minibatch_frac": [1],
        # "regressor__regressor__minibatch_size":   Integer(1, 100),
        #  "regressor__regressor__Base":             Categorical(["DecisionTreeRegressor", "Ridge", "Lasso",
        #                                            "KernelRidge", "SVR"]),
        "regressor__regressor__natural_gradient": [True],
        "regressor__regressor__verbose": [False],
        # "regressor__regressor__min_samples_split": Real(0.05, 0.99),
        # "regressor__regressor__min_samples_leaf": Real(0.05, 0.99),
        "regressor__regressor__tol": Real(1e-4, 1e-2, prior="log-uniform"),
    }


    if algortihm == "GPR":
        if kernel == 'rbf':
            return {
        "regressor__regressor__lr": [1e-2], 
        "regressor__regressor__n_epoch": [100],
        "regressor__regressor__lengthscale": Real(0.05, 3.0), 
            }
        if kernel == 'matern':
            return {
        "regressor__regressor__lr": [1e-2], 
        "regressor__regressor__n_epoch": [100],
        "regressor__regressor__lengthscale": Real(0.05, 3.0), 
        "regressor__regressor__nu": Real(0.5, 2.5),

            }
        if kernel == 'tanimoto':
            return {
        "regressor__regressor__lr": [1e-2], 
        "regressor__regressor__n_epoch": [100],
        }


    if algortihm == "sklearn-GPR":
        if kernel == 'rbf':
            return {
            "regressor__regressor__kernel__length_scale": Real(0.05, 3.0),
            }

        if kernel == 'matern':
            return {
            "regressor__regressor__n_restarts_optimizer": [25],
            "regressor__regressor__kernel__nu": Real(0.5, 2.5),
            "regressor__regressor__kernel__length_scale": Real(0.05, 3.0),
            }
    
    else:
        return None



results = {
    "model": [],
    "best_params":[],
    "imputer": [],
    "cv_score":[],
    "r2_train":[],
    "rmse_train":[],
    "r2_test":[],
    "rmse_test":[]
}


def construct_kernel(algorithm:str, kernel:str=None):
    if algorithm == "GPR":
        return kernel
    elif algorithm == "sklearn-GPR":
        if kernel == "rbf":
            return RBF()
        elif kernel == "matern":
            return Matern()    
    else:
        return None

