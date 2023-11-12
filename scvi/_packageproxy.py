class OptionalPackageProxy:
    def __init__(self, package_name, dependency_name):
        self.__dict__["package_name"] = package_name
        self.__dict__["dependency_name"] = dependency_name
        self.__dict__["_module"] = None

    def _import_module(self):
        # Split the package name to traverse submodules
        package_parts = self.package_name.split(".")
        try:
            # Import the top-level module
            module = __import__(package_parts[0])
            # Traverse to the desired submodule
            for part in package_parts[1:]:
                module = getattr(module, part)
            return module
        except ImportError:
            raise ImportError(
                f"The package `{package_parts[0]}` is now part of the `{self.dependency_name}` dependency. "
                f"Please install it using `pip install {self.dependency_name}`."
            )

    def __getattr__(self, item):
        if self.__dict__["_module"] is None:
            self.__dict__["_module"] = self._import_module()
        return getattr(self.__dict__["_module"], item)

    def __setattr__(self, key, value):
        if key in ["_module", "package_name"]:
            self.__dict__[key] = value
        else:
            if self.__dict__["_module"] is None:
                self.__dict__["_module"] = self._import_module()
            setattr(self.__dict__["_module"], key, value)


# Create proxy instances for each optional package
jax = OptionalPackageProxy("jax", "scvi-tools[jax]")
jnp = OptionalPackageProxy("jax.numpy", "scvi-tools[jax]")
random = OptionalPackageProxy("jax.random", "scvi-tools[jax]")
jaxlib = OptionalPackageProxy("jaxlib", "scvi-tools[jax]")
optax = OptionalPackageProxy("optax", "scvi-tools[jax]")
flax = OptionalPackageProxy("flax", "scvi-tools[jax]")
nn = OptionalPackageProxy("flax.linen", "scvi-tools[jax]")
chex = OptionalPackageProxy("chex", "scvi-tools[jax]")
numpyro = OptionalPackageProxy("numpyro", "scvi-tools[jax]")
dist = OptionalPackageProxy("numpyro.distributions", "scvi-tools[jax]")
