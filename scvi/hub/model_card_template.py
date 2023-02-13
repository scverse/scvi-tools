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

**data_is_minified**: {data_is_minified}

# Training data

This is an optional link to where the training data is stored if it is too large
to host on the huggingface Model hub.

<!-- If your model is not uploaded with any data (e.g., minified data) on the Model Hub, then make
sure to provide this field if you want users to be able to access your training data. See the scvi-tools
documentation for details. -->

Training data url: {training_data_url}

# Training code

This is an optional link to the code used to train the model.

Training code url: {training_code_url}

# References

{references}\
"""
