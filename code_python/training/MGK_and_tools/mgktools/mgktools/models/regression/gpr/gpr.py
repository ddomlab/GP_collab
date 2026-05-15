#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from graphdot.model.gaussian_process.nystrom import *
from graphdot.model.gaussian_process.gpr import GaussianProcessRegressor as GPR
from sklearn.base import BaseEstimator, RegressorMixin
import numpy as np


class MarginalizedGaussianProcessRegressor(GPR):
    def predict(self, Z, return_std=False, return_cov=False):
        """Predict using the trained GPR model.

        Parameters
        ----------
        Z: list of objects or feature vectors.
            Input values of the unknown data.
        return_std: boolean
            If True, the standard-deviations of the predictions at the query
            points are returned along with the mean.
        return_cov: boolean
            If True, the covariance of the predictions at the query points are
            returned along with the mean.

        Returns
        -------
        ymean: 1D array
            Mean of the predictive distribution at query points.
        std: 1D array
            Standard deviation of the predictive distribution at query points.
        cov: 2D matrix
            Covariance of the predictive distribution at query points.
        """
        if not hasattr(self, 'Kinv'):
            raise RuntimeError('Model not trained.')
        Ks = self._gramian(None, Z, self._X)[:, self._y_mask]
        ymean = (Ks @ self.Ky) * self._ystd + self._ymean
        ymean = np.asarray(ymean).ravel()
        if return_std is True:
            Kss = self._gramian(self.alpha, Z, diag=True)
            std = np.sqrt(
                np.maximum(0, Kss - (Ks @ (self.Kinv @ Ks.T)).diagonal())
            )
            std = np.asarray(std).ravel()
            return ymean, std
        elif return_cov is True:
            Kss = self._gramian(self.alpha, Z)
            cov = np.maximum(0, Kss - Ks @ (self.Kinv @ Ks.T))
            return ymean, cov
        else:
            return ymean

    def predict_loocv(self, Z, z, return_std=False):
        """Compute the leave-one-out cross validation prediction of the given
        data.

        Parameters
        ----------
        Z: list of objects or feature vectors.
            Input values of the unknown data.
        z: 1D array
            Target values of the training data.
        return_std: boolean
            If True, the standard-deviations of the predictions at the query
            points are returned along with the mean.

        Returns
        -------
        ymean: 1D array
            Leave-one-out mean of the predictive distribution at query points.
        std: 1D array
            Leave-one-out standard deviation of the predictive distribution at
            query points.
        """
        z_mask, z_masked = self.mask(z)
        if self.normalize_y is True:
            z_mean, z_std = np.mean(z_masked), np.std(z_masked)
            z = (z_masked - z_mean) / z_std
        else:
            z_mean, z_std = 0, 1
            z = z_masked

        K = self._gramian(self.alpha, Z)[z_mask, :][:, z_mask]
        Kinv, _ = self._invert(K, rcond=self.beta)
        Kinv_diag = Kinv.diagonal()
        ymean = (z - Kinv @ z / Kinv_diag) * z_std + z_mean
        ymean = np.asarray(ymean).ravel()
        if return_std is True:
            std = np.sqrt(1 / np.maximum(Kinv_diag, 1e-14))
            std = np.asarray(std).ravel()
            return ymean, std
        else:
            return ymean

    def predict_interpretable(self, Z):
        assert self.normalize_y is False, "y normalization must be False for molecular attribution."
        if not hasattr(self, 'Kinv'):
            raise RuntimeError('Model not trained.')
        Ks = self._gramian(None, Z, self.X)
        Linv = np.linalg.inv(self.Kinv.L)
        K_left = Ks @ Linv.T @ Linv
        return K_left, np.einsum('ij,j->ij', Ks @ Linv.T @ Linv, self.y * self._ystd + self._ymean)

    def predict_nodal(self, Z):
        assert self.normalize_y is False, "y normalization must be False for atomic attribution."
        if not hasattr(self, 'Kinv'):
            raise RuntimeError('Model not trained.')
        Ks = self.kernel(Z, self.X, nodal_X=True)
        ymean = (Ks @ self.Ky) * self._ystd + self._ymean
        return ymean



class MGKRegressorSklearn(BaseEstimator, RegressorMixin):
    def __init__(self, kernel, alpha, optimizer, normalize_y,
                 loss, repeat, verbose
                 ):
        
        self.kernel = kernel
        self.alpha = alpha
        self.optimizer = optimizer
        self.normalize_y = normalize_y
        self.loss = loss
        self.repeat = repeat
        self.verbose = verbose

    def fit(self, X, y):
        # X is a list/array of graphs — DO NOT call check_array on it.
        # Lazy import so the wrapper module is cheap to import.
        self._estimator = MarginalizedGaussianProcessRegressor(
            kernel=self.kernel,
            alpha=self.alpha,
            optimizer=self.optimizer,
            normalize_y=self.normalize_y,
        )
        # graphdot expects array-like; np.asarray on graphs gives object array
        # for i in range(X.shape[1]):
        #     print(f"col {i}: {type(X[0, i]).__name__}")
        self._estimator.fit(
            X, y,
            loss=self.loss,
            repeat=self.repeat,
            verbose=self.verbose,
        )

        self.is_fitted_ = True
        return self

    def predict(self, X, return_std=False, return_cov=False):
        if not getattr(self, "is_fitted_", False):
            raise RuntimeError("Call fit before predict.")
        
        if return_std:
            y_pred, y_std = self._estimator.predict(
                X, return_std=return_std, return_cov=return_cov
            )

        else:
            y_pred = self._estimator.predict(
                X, return_std=return_std, return_cov=return_cov
            )
            y_std = None

        return {
            "y_pred": y_pred,
            "y_std": y_std,
        }