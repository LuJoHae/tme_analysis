= 3. Results

== 3.1 Single-Cell Reference & Selective Inference
Leiden clustering successfully mapped the 40,002 cells into 6 major lineages and 22 cell types (@fig-sc-ref). Analysis of cell quality metrics verified that cells retained high sequencing depth with low mitochondrial proportions across all clusters (@fig-supp-s1). Evaluating Leiden clustering granularities confirmed that resolution 0.5 achieved optimal segregation between lymphoid, myeloid, stromal, and malignant compartments without over-fragmentation (@fig-supp-s2). Dual-method tumor cell detection reliably distinguished malignant epithelial cells via inferred copy number aberrations and marker expression (@fig-supp-s3).

#figure(
  image("../figures/single-cell/umap_cell_type_annotations.png", width: 85%),
  caption: [Hierarchical cell type annotation of the single-cell reference. UMAP projection of 40,002 single cells displaying 22 transcriptionally distinct cell types across 6 major lineages.]
) <fig-sc-ref>

Within these 22 cell types, K-means sub-clustering yielded 58 detailed functional sub-clusters. Truncated normal selective inference p-values for all sub-clusters were highly significant ($p < 10^(-16)$, @fig-supp-s4), providing mathematical confirmation of robust centroid separation beyond sample noise.
