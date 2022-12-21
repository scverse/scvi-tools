template = """\
---
{card_data}
---

# Description

{description}

# Model properties

Many model properties are in the model tags. Some more are listed below.

**model_init_params**:
```json
{model_init_params}
```

**model_setup_anndata_args**:
```json
{model_setup_anndata_args}
```

**model_summary_stats**:
{model_summary_stats}

**model_data_registry**:
{model_data_registry}

**model_parent_module**: {model_parent_module}

**data_is_latent**: {data_is_latent}

# Training data

This is an optional link to where the training data is stored if it is too large
to host on the huggingface Model hub.

<!-- This field is required for models that haven't been minified by converting to latent
mode. See the scvi-tools documentation for more details. -->

Training data url: {training_data_url}

# Training code

This is an optional link to the code used to train the model.

Training code url: {training_code_url}

# References

{references}\
"""
