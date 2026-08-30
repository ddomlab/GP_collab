from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor,RandomForestClassifier
from xgboost import XGBRegressor, XGBClassifier, callback
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.svm import SVR
from ngboost import NGBRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import Lasso
from GPytorch_kernel_mix import GPytorchMCMCRegressor, GPytorchMAPsklearnRegressor
from GPyro_kernel_mix import GpyroHMCsklearnRegressor

from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import FunctionTransformer
import numpy as np
from sklearn.preprocessing import (StandardScaler,
                                   QuantileTransformer,
                                   MinMaxScaler,
                                   RobustScaler)
from sklearn.base import TransformerMixin
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer

from typing import Callable, Optional, Union, Dict
from types import NoneType
from skopt.space import Integer, Real, Categorical

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, Matern, RationalQuadratic, ConstantKernel as C


unrolling_feature_factory: dict[str, list[str]] = {
                                                "polysize": ['Mw (g/mol)', 'PDI'],
                                                "HSPs": ['polymer dP', 'polymer dD' , 'polymer dH', 'solvent dP', 'solvent dD', 'solvent dH'],
                                                "Ra": ["Ra"],
                                                "solvent properties": ['Concentration (mg/ml)', 'Temperature SANS/SLS/DLS/SEC (K)'],
                                                 }


imputer_factory: Dict[str, TransformerMixin] = {
    "mean": SimpleImputer(strategy="mean"),
    "median": SimpleImputer(strategy="median"),
    "most-frequent": SimpleImputer(strategy="most_frequent"),
    "uniform KNN": KNNImputer(weights="uniform"),
    "distance KNN": KNNImputer(weights="distance"),
    "iterative": IterativeImputer(sample_posterior=True),
}

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
}


def optimized_models(
                    model_name:str,
                    kernel_parameters:dict=None,
                    feat_group:dict[str, list[int]]=None,
                    random_state:int=42, 
                    kernel_type:Optional[str]=None,
                    kernel_mixing_method:Optional[str]=None,
                    graph_kernel_config:Optional[str]=None,
                    alpha:Optional[float]=0.01,
                    target_transformer:str|None=None,
                     **kwargs):
    normalize_y = bool(target_transformer)
    if 'NGB'==model_name:
        return NGBRegressor(n_estimators=500, learning_rate=0.01, tol=1e-4,
                             random_state=random_state, verbose=False,
                            #  early_stopping_rounds=50,
                            #  **kwargs,
                             
                             )
    if 'XGBR'==model_name:
        # params = callback.EarlyStopping(rounds=50, metric_name="rmse")
        return  XGBRegressor(
                            eval_metric="rmse", 
                            random_state=random_state, 
                            n_jobs=-1,
                            n_estimators=100,
                            verbose=False,
                            # callbacks=[params]
                            # early_stopping_rounds=3
                            #    **kwargs,
                        )
    
    if 'RF'==model_name:
        return RandomForestRegressor(n_estimators=100, max_depth=None, 
                                        random_state=random_state, n_jobs=-1,
                                        max_features="sqrt"
                                        )
    if 'HGBR'==model_name:
        return HistGradientBoostingRegressor(max_iter=2000, max_depth=None, 
                                            #  min_samples_leaf=20, max_leaf_nodes=1000,
                                             learning_rate=0.01, l2_regularization=1e-6,
                                             scoring='neg_root_mean_squared_error', random_state=random_state,**kwargs)
    if 'MLP'==model_name:
        return MLPRegressor(max_iter=200)

    if 'sklearn-GPR'==model_name:
        return GaussianProcessRegressor(random_state=random_state, **kwargs)
    
    if "GPytorchMAP"==model_name:
        return GPytorchMAPsklearnRegressor(
                                    feat_group=feat_group,
                                    kernel_type=kernel_type,
                                    kernel_mixing_method=kernel_mixing_method,
                                    ssk_parameters=kernel_parameters,
                                    n_epoch=kwargs.get('n_epochs', 400),
                                    lr=kwargs.get('lr', 0.01),
                                    progbar=kwargs.get('progbar', False),
                                    prior=kwargs.get('prior', True),
                                    random_state=random_state,
                                    normalize_y=normalize_y,
                                    use_cuda=kwargs.get('use_cuda', True),
                                    )
    if "GPytorchHMC" == model_name:
        return GPytorchMCMCRegressor(**kwargs)
    
    if "GpyroHMC" == model_name:
        return GpyroHMCsklearnRegressor(
                                    feat_group=feat_group,
                                    num_samples=kwargs.get('num_samples', 200),
                                    warmup_steps=kwargs.get('warmup_steps', 200),
                                    num_chains=kwargs.get('num_chains', 1),
                                    use_cuda=kwargs.get('use_cuda', False),
                                    random_state=random_state,
                                    kernel_mixing_method=kernel_mixing_method,
                                    kernel_type=kernel_type,
                                    normalize_y=normalize_y,
                                    )
    

    if "MGK" == model_name:
        from mgktools.models.regression.gpr.gpr import MGKRegressorSklearn

        return MGKRegressorSklearn(
                kernel=graph_kernel_config.kernel,
                optimizer="L-BFGS-B",
                alpha=alpha,
                normalize_y=normalize_y,
                loss='likelihood', 
                repeat=1,
                verbose=False
                )
    
    else:
        raise ValueError(f"Model {model_name} not recognized in optimized_models factory.")



def get_regressor_search_space(algortihm:str, kernel:str=None) -> Dict :
    SEARCH_SPACES: Dict[str, Dict] = {
        "MLR": {
            "regressor__regressor__fit_intercept": [True, False],
        },

        "Lasso": {
            "regressor__alpha": Real(1e-3, 1e3, prior="log-uniform"),
            "regressor__fit_intercept": [True, False],
            "regressor__selection": Categorical(["cyclic", "random"]),
        },

        "KNN": {
            "regressor__n_neighbors": Integer(1, 50),
            "regressor__weights": Categorical(["uniform", "distance"]),
            "regressor__algorithm": Categorical(["ball_tree", "kd_tree", "brute"]),
            "regressor__leaf_size": Integer(1, 100),
        },

        "SVR": {
            "regressor__kernel": Categorical(["linear", "rbf"]),
        },

        "MLP": {
            "regressor__regressor__learning_rate_init": Real(1e-5, 1e-1, prior="log-uniform"),
            "regressor__regressor__alpha": Real(1e-8, 1e-6, prior="log-uniform"),
            "regressor__regressor__batch_size": Categorical([32, 64, 128, 256]),
            "regressor__regressor__max_iter": [200],
        },

        "RF": {
            "regressor__regressor__n_estimators": Integer(10, 2000, prior="log-uniform"),
            "regressor__regressor__max_depth": [None],
            "regressor__regressor__min_samples_split": Real(0.001, 0.99),
            "regressor__regressor__min_samples_leaf": Real(0.001, 0.99),
            "regressor__regressor__max_features": Categorical(["sqrt", "log2", None]),
        },

        "XGBR": {
            "regressor__regressor__n_estimators": Integer(50, 2000, prior="log-uniform"),
            "regressor__regressor__max_depth": Integer(10, 10000, prior="log-uniform"),
            "regressor__regressor__n_jobs": [-2],
            "regressor__regressor__learning_rate": Real(1e-4, 1e-1, prior="log-uniform"),
        },

        "DT": {
            "regressor__regressor__min_samples_split": Real(0.05, 0.99),
            "regressor__regressor__min_samples_leaf": Real(0.05, 0.99),
            "regressor__regressor__max_features": Categorical([None, "sqrt", "log2"]),
            "regressor__regressor__max_depth": [None],
            "regressor__regressor__ccp_alpha": Real(0.05, 0.99),
        },

        "NGB": {
            "regressor__regressor__n_estimators": Integer(50, 2000, prior="log-uniform"),
            "regressor__regressor__learning_rate": Real(1e-4, 1e-1, prior="log-uniform"),
            "regressor__regressor__natural_gradient": [True],
            "regressor__regressor__verbose": [False],
            "regressor__regressor__tol": Real(1e-4, 1e-2, prior="log-uniform"),
        },
    }
    return SEARCH_SPACES.get(algortihm, {})


