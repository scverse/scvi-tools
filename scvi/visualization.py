import matplotlib.pyplot as plt
from sklearn.manifold import TSNE


def show_t_sne(latent, labels, title, return_t_sne=False):
    if latent.shape[1] != 2:
        latent = TSNE().fit_transform(latent)
    plt.figure(figsize=(10, 10))
    plt.scatter(latent[:, 0], latent[:, 1], c=labels, edgecolors='none')  # cmap=cmap,cmap=plt.get_cmap("tab10", 7)
    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    #plt.show()

    if return_t_sne:
        return latent
