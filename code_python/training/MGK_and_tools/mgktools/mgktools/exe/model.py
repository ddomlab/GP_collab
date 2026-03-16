#!/usr/bin/env python
# -*- coding: utf-8 -*-
from mgktools.models import set_model
from mgktools.exe.args import TrainArgs


def set_model_from_args(args: TrainArgs, kernel):
    if args.model_type == 'gpr-nle':
        N = args.n_local
    elif args.model_type == 'gpr-nystrom':
        N = args.n_core
    else:
        N = args.n_samples_per_model
    if hasattr(args, 'optimizer'):
        optimizer = args.optimizer
    else:
        optimizer = None
    return set_model(
        model_type=args.model_type,
        graph_kernel_type=args.graph_kernel_type,
        kernel=kernel,
        optimizer=optimizer,
        alpha=args.alpha_,
        C=args.C_,
        n_estimators=args.n_estimators,
        n_samples_per_model=N,
        ensemble_rule=args.ensemble_rule,
        n_jobs=args.n_jobs
    )
