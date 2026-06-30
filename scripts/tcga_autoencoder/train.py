import torch
from .model import vae_loss
from .plot import plot_latent_space

def train_vae(model, dataloader, optimizer, epochs, device, eval_data=None, labels=None, output_dir=None):
    # Plot epoch 0
    if eval_data is not None and labels is not None and output_dir is not None:
        model.eval()
        with torch.no_grad():
            _, mu, _ = model(eval_data.to(device))
            plot_latent_space(mu.cpu().numpy(), labels, output_dir, filename=f"latent_space_epoch_0.png")
            
    model.train()
    loss_history = []
    
    for epoch in range(epochs):
        train_loss = 0
        for batch_idx, (data, batch_labels) in enumerate(dataloader):
            data = data.to(device)
            batch_labels = batch_labels.to(device)
            optimizer.zero_grad()
            
            recon_batch, mu, logvar = model(data)
            loss, recon_loss, kl, triplet = vae_loss(recon_batch, data, mu, logvar, batch_labels, lambda_triplet=100.0)
            
            loss.backward()
            train_loss += loss.item()
            optimizer.step()
            
        avg_loss = train_loss / len(dataloader.dataset)
        loss_history.append(avg_loss)
        print(f"Epoch {epoch+1}/{epochs}, Loss: {avg_loss:.4f}")
        
        if eval_data is not None and labels is not None and output_dir is not None:
            if (epoch + 1) % 25 == 0 or (epoch + 1) == epochs:
                model.eval()
                with torch.no_grad():
                    _, mu, _ = model(eval_data.to(device))
                    plot_latent_space(mu.cpu().numpy(), labels, output_dir, filename=f"latent_space_epoch_{epoch+1}.png")
                model.train()
                
    return loss_history
