def pytest_addoption(parser):
    """Register custom CLI options for the test suite."""
    parser.addoption(
        "--internet-tests",
        action="store_true",
        default=False,
        help="Run tests that retrieve stuff from the internet. This increases test time.",
    )
    parser.addoption(
        "--multigpu-tests",
        action="store_true",
        default=False,
        help="Run tests that are designed for multiGPU.",
    )
    parser.addoption(
        "--autotune-tests",
        action="store_true",
        default=False,
        help="Run tests that are designed for Ray Autotune.",
    )
    parser.addoption(
        "--mlflow-tests",
        action="store_true",
        default=False,
        help="Run tests that are designed for MLFlow.",
    )
    parser.addoption(
        "--custom-dataloader-tests",
        action="store_true",
        default=False,
        help="Run tests that deal with custom dataloaders. This increases test time.",
    )
    parser.addoption(
        "--optional",
        action="store_true",
        default=False,
        help="Run tests that are optional.",
    )
    parser.addoption(
        "--jax",
        action="store_true",
        default=False,
        help="Run tests that are Jax adopted.",
    )
    parser.addoption(
        "--accelerator",
        action="store",
        default="cpu",
        help="Option to specify which accelerator to use for tests.",
    )
    parser.addoption(
        "--devices",
        action="store",
        default="auto",
        help="Option to specify which devices to use for tests.",
    )
    parser.addoption(
        "--seed",
        action="store",
        default=0,
        help="Option to specify which scvi-tools seed to use for tests.",
    )
    parser.addoption(
        "--private",
        action="store_true",
        default=False,
        help="Run tests that are private.",
    )
