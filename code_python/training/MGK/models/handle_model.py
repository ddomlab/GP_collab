from models.gpr import GaussianProcessRegressor



def set_model(model_type: str, # 'gpr', 'gpr-nystrom', 'gpr-nle', 'svr', 'gpc', 'svc'
              kernel,
              optimizer = None,
              alpha: float = None,
              n_jobs: int = 1):
    if model_type == 'gpr':
        assert alpha is not None
        model = GaussianProcessRegressor(
            kernel=kernel,
            optimizer=optimizer,
            alpha=alpha,
            normalize_y=False,
        )

    else:
        raise RuntimeError(f'unsupported model_type: {model_type}')
    return model