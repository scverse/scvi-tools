import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

plt.switch_backend('agg')


def show_t_sne(latent, labels, title):
    if latent.shape[1] != 2:
        latent = TSNE().fit_transform(latent)
    plt.figure(figsize=(10, 10))
    plt.scatter(latent[:, 0], latent[:, 1], c=labels, edgecolors='none')
    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    print("saving tsne figure as tsne.png")
    plt.savefig("tsne.png")
