import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE

plt.switch_backend('agg')


def show_t_sne(latent, labels, n_samples=1000):
    idx_t_sne = np.random.permutation(len(latent))[:n_samples]
    if latent.shape[1] != 2:
        latent = TSNE().fit_transform(latent[idx_t_sne])
    plt.figure(figsize=(10, 10))
    plt.scatter(latent[:, 0], latent[:, 1], c=(np.array(labels)[idx_t_sne]).ravel(), edgecolors='none')
    plt.axis("off")
    plt.tight_layout()
    print("saving tsne figure as tsne.png")
    plt.savefig("tsne.png")
