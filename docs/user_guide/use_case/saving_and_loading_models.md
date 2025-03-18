# Saving and loading SCVI models

In scvi-tools (Single-Cell Variational Inference), saving and loading models is straightforward and allows you to store your trained models for later use or share them with others. To save a model, you use the save method, which saves the model's state to a file in the .pt format and optionally includes the data used for training. You can then load the saved model using the load() method, which reloads the model's state from the saved file, allowing you to continue training or perform inference without retraining from scratch. Here's an example:

```python
# Saving a model
model.save("my_model.pt")

# Loading a model
model = scvi.model.SCVI.load("my_model.pt")
```

There are several common use cases that require us to save and load a model, besides the general case:
1. Saving and loading a model to/from scvi-hub, perhaps after doing minification, which reduces the disk space required by the training data.
2. Saving a model and use it as a reference mapping for transfer learning, see our SCVI Arches tutorial {doc}`/tutorials/notebooks/scrna/scarches_scvi_tools`
3. Directly create models from a stored scvi-tools models, like in SCANVI, where we init a model with weights from a pretrained SCVI model.
