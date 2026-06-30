import torch
import torch.nn as nn
import torch.nn.functional as F

class TCGAVAE(nn.Module):
    def __init__(self, input_dim, hidden_dim1=2048, hidden_dim2=1024, hidden_dim3=512, latent_dim=2):
        super(TCGAVAE, self).__init__()
        
        # Encoder
        self.fc1 = nn.Linear(input_dim, hidden_dim1)
        self.fc2 = nn.Linear(hidden_dim1, hidden_dim2)
        self.fc3_enc = nn.Linear(hidden_dim2, hidden_dim3)
        
        self.fc_mu = nn.Linear(hidden_dim3, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim3, latent_dim)
        
        # Decoder
        self.fc4 = nn.Linear(latent_dim, hidden_dim3)
        self.fc5 = nn.Linear(hidden_dim3, hidden_dim2)
        self.fc6 = nn.Linear(hidden_dim2, hidden_dim1)
        self.fc7 = nn.Linear(hidden_dim1, input_dim)
        
        self.dropout = nn.Dropout(0.2)

    def encode(self, x):
        h1 = F.relu(self.fc1(x))
        h1 = self.dropout(h1)
        h2 = F.relu(self.fc2(h1))
        h2 = self.dropout(h2)
        h3 = F.relu(self.fc3_enc(h2))
        h3 = self.dropout(h3)
        return self.fc_mu(h3), self.fc_logvar(h3)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z):
        h4 = F.relu(self.fc4(z))
        h4 = self.dropout(h4)
        h5 = F.relu(self.fc5(h4))
        h5 = self.dropout(h5)
        h6 = F.relu(self.fc6(h5))
        h6 = self.dropout(h6)
        return self.fc7(h6)

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        return self.decode(z), mu, logvar

def batch_hard_triplet_loss(embeddings, labels, margin=1.0):
    pairwise_dist = torch.cdist(embeddings, embeddings, p=2)
    labels = labels.view(-1, 1)
    mask_pos = torch.eq(labels, labels.T).float()
    mask_neg = 1.0 - mask_pos
    mask_pos_no_self = mask_pos - torch.eye(labels.size(0), device=labels.device)
    
    # Hardest positive (max dist)
    hardest_pos = torch.max(pairwise_dist * mask_pos_no_self, dim=1)[0]
    
    # Hardest negative (min dist of negatives)
    max_dist = torch.max(pairwise_dist).item()
    hardest_neg = torch.min(pairwise_dist + mask_pos * (max_dist + 1e-5), dim=1)[0]
    
    loss = F.relu(hardest_pos - hardest_neg + margin)
    return loss.mean()

def vae_loss(recon_x, x, mu, logvar, labels=None, lambda_triplet=1.0):
    recon_loss = F.mse_loss(recon_x, x, reduction='sum')
    kl_divergence = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    
    total_loss = recon_loss + kl_divergence
    
    triplet = torch.tensor(0.0, device=x.device)
    if labels is not None:
        triplet = batch_hard_triplet_loss(mu, labels)
        # Weight triplet loss by dataset size roughly or just a scalar
        # Usually requires tuning, we'll just add it scaled to match magnitude roughly, 
        # but since recon is sum, we'll scale triplet by batch_size * lambda_triplet
        total_loss += triplet * x.size(0) * lambda_triplet
        
    return total_loss, recon_loss, kl_divergence, triplet
