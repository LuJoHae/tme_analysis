import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def plot_loss(loss_history, output_path):
    plt.figure(figsize=(8, 5))
    plt.plot(loss_history, label='Training Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('VAE Training Loss')
    plt.legend()
    plt.savefig(output_path / 'loss_curves.png')
    plt.close()

def plot_latent_space(z, labels, output_path, filename='latent_space.png'):
    plt.figure(figsize=(14, 12))
    df = pd.DataFrame({'z1': z[:, 0], 'z2': z[:, 1], 'Cancer_Type': labels})
    # Use a large palette since TCGA has ~33 cancer types
    num_types = len(np.unique(labels))
    palette = sns.color_palette("husl", num_types)
    sns.scatterplot(data=df, x='z1', y='z2', hue='Cancer_Type', s=15, alpha=0.8, palette=palette)
    plt.title('TCGA VAE Latent Space')
    # Put legend outside with smaller font and multiple columns if needed
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small', ncol=2)
    plt.tight_layout()
    if not str(filename).endswith('.png'):
        filename = filename + '.png'
    plt.savefig(output_path / filename)
    plt.close()

def plot_reconstruction_error(reconstruction_errors, output_path):
    plt.figure(figsize=(8, 5))
    sns.histplot(reconstruction_errors, bins=50, kde=True)
    plt.xlabel('Reconstruction Error (MSE)')
    plt.ylabel('Frequency')
    plt.title('Reconstruction Error Distribution')
    plt.savefig(output_path / 'reconstruction_error.png')
    plt.close()

def plot_latent_density(z, output_path):
    plt.figure(figsize=(8, 8))
    sns.kdeplot(x=z[:, 0], y=z[:, 1], cmap="Blues", fill=True, bw_adjust=.5)
    plt.xlabel('z1')
    plt.ylabel('z2')
    plt.title('Latent Space Density')
    plt.savefig(output_path / 'latent_density.png')
    plt.close()
