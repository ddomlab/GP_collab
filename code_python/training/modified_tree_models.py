import numpy as np
import shap
from sklearn.ensemble import RandomForestRegressor
from joblib import Parallel, delayed
from xgboost import XGBRegressor as _XGBRegressor
from ngboost import NGBRegressor as _NGBRegressor
from ngboost.distns import Normal
from ngboost.scores import LogScore
from ngboost.learners import default_tree_learner


def _named_feature_importances(estimator, importances):
    """Pair fitted feature importances with their training-column names.

    Scikit-learn records ``feature_names_in_`` when ``fit`` receives a
    dataframe whose column names are strings. Array inputs have no feature
    metadata, so use scikit-learn's conventional ``x0``, ``x1``, ... names in
    that case.
    """
    importances = np.asarray(importances, dtype=float).reshape(-1)
    feature_names = getattr(estimator, "feature_names_in_", None)

    if feature_names is None:
        feature_names = np.asarray(
            [f"x{i}" for i in range(importances.size)],
            dtype=object,
        )
    else:
        feature_names = np.asarray(feature_names, dtype=object).reshape(-1)

    if feature_names.size != importances.size:
        raise RuntimeError(
            "The number of stored feature names does not match the number "
            "of fitted feature importances."
        )

    if len(set(feature_names.tolist())) != feature_names.size:
        raise ValueError(
            "Feature names must be unique to create a named importance mapping."
        )

    return {
        str(name): float(importance)
        for name, importance in zip(feature_names, importances)
    }


def _shap_feature_importances(estimator, *, model_output, target_scale=1.0):
    """Explain the X saved by predict using the X saved by fit."""
    if not hasattr(estimator, "X_train"):
        raise AttributeError(
            "Training data do not exist. Call fit() before _get_SHAP()."
        )
    if not hasattr(estimator, "X_test"):
        raise AttributeError(
            "Prediction data do not exist. Call predict() before _get_SHAP()."
        )

    background = shap.maskers.Independent(
        estimator.X_train,
        max_samples=len(estimator.X_train),
    )
    explainer = shap.TreeExplainer(
        model=estimator,
        data=background,
        model_output=model_output,
        feature_perturbation="interventional",
    )
    shap_values = explainer(estimator.X_test, check_additivity=False)
    values = np.asarray(shap_values.values) * np.asarray(target_scale)

    global_importance = np.abs(values).mean(axis=0)
    return _named_feature_importances(estimator, global_importance)


class RFRegressor(RandomForestRegressor):
    """Random forest with optional target scaling and ensemble uncertainty.

    All inherited constructor defaults match ``sklearn``.  ``normalize_y`` is
    an extension and defaults to ``False`` so constructing this class with no
    arguments retains the package model's target-scale behaviour.
    """

    def __init__(
        self,
        n_estimators=100,
        *,
        criterion="squared_error",
        max_depth=None,
        min_samples_split=2,
        min_samples_leaf=1,
        min_weight_fraction_leaf=0.0,
        max_features=1.0,
        max_leaf_nodes=None,
        min_impurity_decrease=0.0,
        bootstrap=True,
        oob_score=False,
        n_jobs=None,
        random_state=None,
        verbose=0,
        warm_start=False,
        ccp_alpha=0.0,
        max_samples=None,
        monotonic_cst=None,
        normalize_y: bool = False,
    ):
        super().__init__(
            n_estimators=n_estimators,
            criterion=criterion,
            max_depth=max_depth,
            min_samples_split=min_samples_split,
            min_samples_leaf=min_samples_leaf,
            min_weight_fraction_leaf=min_weight_fraction_leaf,
            max_features=max_features,
            max_leaf_nodes=max_leaf_nodes,
            min_impurity_decrease=min_impurity_decrease,
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
            warm_start=warm_start,
            ccp_alpha=ccp_alpha,
            max_samples=max_samples,
            monotonic_cst=monotonic_cst,
        )

        self.normalize_y = normalize_y

    @property
    def y(self):
        if not hasattr(self, "_y_train"):
            raise AttributeError(
                "Training targets do not exist. Call fit() first."
            )

        return self._y_train * self.y_std_ + self.y_mean_

    @y.setter
    def y(self, y):
        y = np.asarray(y, dtype=float)
        if y.ndim not in (1, 2):
            raise ValueError("y must be one- or two-dimensional.")

        if self.normalize_y:
            self.y_mean_ = np.mean(y, axis=0)
            self.y_std_ = np.std(y, axis=0)
            self.y_std_ = np.where(self.y_std_ == 0, 1.0, self.y_std_)

            y_scaled = (y - self.y_mean_) / self.y_std_
        else:
            # Scalars broadcast correctly for both single- and multi-output y.
            self.y_mean_ = 0.0
            self.y_std_ = 1.0
            y_scaled = y

        self._y_train = y_scaled

    def fit(self, X, y, sample_weight=None):
        self.X_train = X
        self.y = y

        super().fit(
            X,
            self._y_train,
            sample_weight=sample_weight,
        )

        return self

    def predict(self, X, return_std=False):
        self.X_test = X

        # RF mean prediction
        y_pred_scaled = super().predict(X)

        # Transform back to original scale
        y_pred = y_pred_scaled * self.y_std_ + self.y_mean_

        if return_std:
            # The forest validates dataframe feature names, but its individual
            # trees are fitted internally on unnamed arrays. Pass unnamed data
            # to those trees to avoid scikit-learn's feature-name warning.
            X_tree = X.to_numpy() if hasattr(X, "to_numpy") else X

            # Predictions from individual trees in parallel
            tree_predictions = np.asarray(
                Parallel(n_jobs=self.n_jobs, prefer="threads")(
                    delayed(tree.predict)(X_tree)
                    for tree in self.estimators_
                )
            )

            # Tree-to-tree standard deviation
            y_std_scaled = np.std(
                tree_predictions,
                axis=0,
            )

            # Transform std back to original scale
            y_std = y_std_scaled * self.y_std_

        else:
            y_std = None

        return {
            "y_pred": np.asarray(y_pred).ravel(),
            "y_std": (
                np.asarray(y_std).ravel()
                if y_std is not None
                else None
            ),
        }

    def _get_MDI(self):
        """Return ``{feature_name: MDI importance}`` for the fitted forest."""
        return _named_feature_importances(self, self.feature_importances_)

    def _get_SHAP(self):
        """Return global SHAP importance for the last prediction X."""
        return _shap_feature_importances(
            self,
            model_output="raw",
            target_scale=self.y_std_,
        )


class XGBRegressor(_XGBRegressor):
    """XGBoost regressor with optional target scaling.

    XGBoost does not expose predictive uncertainty through this estimator, so
    ``y_std`` is always ``None``.
    """

    def __init__(
        self,
        *,
        objective="reg:squarederror",
        normalize_y: bool = False,
        verbosity: int = 0,
        **kwargs,
    ):
        super().__init__(
            objective=objective,
            verbosity=verbosity,
            **kwargs,
        )
        self.normalize_y = normalize_y

    @property
    def y(self):
        """Return training targets on their original scale."""
        if not hasattr(self, "_y_train"):
            raise AttributeError(
                "Training targets do not exist. Call fit() first."
            )

        return self._y_train * self.y_std_ + self.y_mean_

    @y.setter
    def y(self, y):
        y = np.asarray(y, dtype=float)
        if y.ndim not in (1, 2):
            raise ValueError("y must be one- or two-dimensional.")

        if self.normalize_y:
            self.y_mean_ = np.mean(y, axis=0)
            self.y_std_ = np.std(y, axis=0)
            self.y_std_ = np.where(
                np.isclose(self.y_std_, 0.0),
                1.0,
                self.y_std_,
            )
            y_scaled = (y - self.y_mean_) / self.y_std_
        else:
            self.y_mean_ = 0.0
            self.y_std_ = 1.0
            y_scaled = y

        self._y_train = y_scaled

    def get_xgb_params(self):
        """Return native booster parameters, excluding wrapper-only options."""
        params = super().get_xgb_params()
        params.pop("normalize_y", None)
        return params

    def fit(self, X, y, *args, **kwargs):
        self.X_train = X
        self.y = y

        # Early-stopping targets must be on the same scale as training y.
        if kwargs.get("eval_set") is not None:
            kwargs["eval_set"] = [
                (
                    X_eval,
                    (np.asarray(y_eval) - self.y_mean_) / self.y_std_,
                )
                for X_eval, y_eval in kwargs["eval_set"]
            ]

        super().fit(X, self._y_train, *args, **kwargs)
        return self

    def predict(self, X, *args, return_std=False, **kwargs):
        """Return XGBoost predictions and ``None`` for uncertainty.

        Positional and keyword prediction options are forwarded unchanged to
        :class:`xgboost.XGBRegressor`. ``return_std`` is accepted only for API
        compatibility with the other project regressors.
        """
        self.X_test = X
        del return_std
        y_pred_scaled = super().predict(X, *args, **kwargs)
        y_pred = y_pred_scaled * self.y_std_ + self.y_mean_

        return {
            "y_pred": np.asarray(y_pred).ravel(),
            "y_std": None,
        }

    def _get_MDI(self):
        """Return ``{feature_name: importance}`` for the fitted booster."""
        return _named_feature_importances(self, self.feature_importances_)

    def _get_SHAP(self):
        """Return global SHAP importance for the last prediction X."""
        return _shap_feature_importances(
            self,
            model_output="raw",
            target_scale=self.y_std_,
        )



class NGBRegressor(_NGBRegressor):
    """NGBoost regressor with optional target scaling and predictive std.

    All inherited constructor defaults match ``ngboost``.  ``normalize_y`` is
    an extension and defaults to ``False`` so the no-argument model preserves
    NGBoost's package behaviour.
    """

    def __init__(
        self,
        Dist=Normal,
        Score=LogScore,
        Base=default_tree_learner,
        natural_gradient=True,
        n_estimators=500,
        learning_rate=0.01,
        minibatch_frac=1.0,
        col_sample=1.0,
        verbose=True,
        verbose_eval=100,
        tol=1e-4,
        random_state=None,
        validation_fraction=0.1,
        early_stopping_rounds=None,
        normalize_y: bool = False,
    ):
        super().__init__(
            Dist=Dist,
            Score=Score,
            Base=Base,
            natural_gradient=natural_gradient,
            n_estimators=n_estimators,
            learning_rate=learning_rate,
            minibatch_frac=minibatch_frac,
            col_sample=col_sample,
            verbose=verbose,
            verbose_eval=verbose_eval,
            tol=tol,
            random_state=random_state,
            validation_fraction=validation_fraction,
            early_stopping_rounds=early_stopping_rounds,
        )

        self.normalize_y = normalize_y

    @property
    def y(self):
        """
        Return training targets on the original scale.
        """
        if not hasattr(self, "_y_train"):
            raise AttributeError(
                "Training targets do not exist. Call fit() first."
            )

        return self._y_train * self.y_std_ + self.y_mean_

    @y.setter
    def y(self, y):
        """
        Store normalized y internally if normalize_y=True.
        """
        y = np.asarray(y, dtype=float).reshape(-1)

        if self.normalize_y:
            self.y_mean_ = float(np.mean(y))
            self.y_std_ = float(np.std(y))

            if self.y_std_ == 0:
                self.y_std_ = 1.0

            y_scaled = (y - self.y_mean_) / self.y_std_

        else:
            self.y_mean_ = 0.0
            self.y_std_ = 1.0
            y_scaled = y

        self._y_train = y_scaled

    def fit(
        self,
        X,
        y,
        X_val=None,
        Y_val=None,
        sample_weight=None,
        val_sample_weight=None,
        train_loss_monitor=None,
        val_loss_monitor=None,
        early_stopping_rounds=None,
    ):
        self.X_train = X
        # Normalize training target
        self.y = y

        # NGBoost converts X to an ndarray internally and therefore does not
        # retain dataframe column names itself. Store them before delegation,
        # following scikit-learn's feature_names_in_ convention.
        columns = getattr(X, "columns", None)
        has_string_columns = columns is not None and all(
            isinstance(name, str) for name in columns
        )
        if has_string_columns:
            self.feature_names_in_ = np.asarray(columns, dtype=object)
        else:
            self.__dict__.pop("feature_names_in_", None)

        # Validation target must use the SAME normalization
        if Y_val is not None:
            Y_val = (
                np.asarray(Y_val).reshape(-1) - self.y_mean_
            ) / self.y_std_

        super().fit(
            X,
            self._y_train,
            X_val=X_val,
            Y_val=Y_val,
            sample_weight=sample_weight,
            val_sample_weight=val_sample_weight,
            train_loss_monitor=train_loss_monitor,
            val_loss_monitor=val_loss_monitor,
            early_stopping_rounds=early_stopping_rounds,
        )

        return self

    def predict(self, X, max_iter=None, return_std=False):
        """Predict a mean and, optionally, a standard deviation.

        ``max_iter`` intentionally remains the second positional parameter,
        matching :class:`ngboost.NGBRegressor`.  ``return_std`` is the added
        project-specific option.
        """
        self.X_test = X

        pred_dist = self.pred_dist(
            X,
            max_iter=max_iter,
        )

        y_pred_scaled = np.asarray(
            pred_dist.predict()
        ).ravel()

        y_pred = (
            y_pred_scaled * self.y_std_
            + self.y_mean_
        )

        if return_std:
            y_std_scaled = np.asarray(
                pred_dist.std()
            ).ravel()

            y_std = y_std_scaled * self.y_std_
        else:
            y_std = None

        return {
            "y_pred": np.asarray(y_pred).ravel(),
            "y_std": (
                np.asarray(y_std).ravel()
                if y_std is not None
                else None
            ),
        }

    def _get_MDI(self):
        """Return named MDI importances for NGBoost's location parameter."""
        return _named_feature_importances(self, self.feature_importances_[0])

    def _get_SHAP(self):
        """Return global SHAP importance for the last prediction X's mean."""
        return _shap_feature_importances(
            self,
            model_output=0,
            target_scale=self.y_std_,
        )

    
__all__ = [
    "RFRegressor",
    "XGBRegressor",
    "NGBRegressor",
]
