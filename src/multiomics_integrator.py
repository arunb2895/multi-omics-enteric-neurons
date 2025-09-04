"""
multiomics_integrator.py
========================

This module provides a simple, working prototype for integrating multiple
omics datasets — metabolomics, bulk RNA‑seq, spatial transcriptomics and
single‑cell RNA‑seq — into a common low‑dimensional representation.

The implementation uses scikit‑learn to perform principal component
analysis (PCA) on each dataset separately (per‑view dimensionality
reduction), concatenates the resulting latent representations and then
applies a second PCA to learn a joint embedding across all views.  This
two‑step approach is inspired by intermediate integration strategies
highlighted in recent multi‑omics surveys【503067821426615†L585-L599】.  It serves
as a template for more sophisticated models (e.g. deep autoencoders or
multi‑view CCA) that could be substituted later without changing the
interface.

Usage example:

    from multiomics_integrator import MultiOmicsIntegrator
    import numpy as np

    # Synthetic data: 100 samples for each omics layer
    metabolomics   = np.random.randn(100, 50)
    bulk_rna       = np.random.randn(100, 1000)
    spatial_trans  = np.random.randn(100, 200)
    sc_rna         = np.random.randn(100, 500)

    integrator = MultiOmicsIntegrator(n_components=10)
    integrated_latent = integrator.fit_transform(
        metabolomics, bulk_rna, spatial_trans, sc_rna
    )
    print(integrated_latent.shape)  # (100, 10)

Note that in a real application the datasets should be pre‑processed
appropriately (normalization, log‑transformation, batch correction, etc.)
before being passed to the integrator.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple
import numpy as np
from sklearn.decomposition import PCA


@dataclass
class MultiOmicsIntegrator:
    """Simple integrator for multiple omics datasets.

    Parameters
    ----------
    n_components : int
        Number of principal components to retain for each view and in the
        final integrated latent space.
    random_state : Optional[int], default=None
        Random state for reproducible PCA initialisation.
    """

    n_components: int = 10
    random_state: Optional[int] = None

    def __post_init__(self) -> None:
        # PCA models for each omics layer and the final integration step
        self.pca_metabolomics = PCA(n_components=self.n_components, random_state=self.random_state)
        self.pca_bulk_rna = PCA(n_components=self.n_components, random_state=self.random_state)
        self.pca_spatial = PCA(n_components=self.n_components, random_state=self.random_state)
        self.pca_sc_rna = PCA(n_components=self.n_components, random_state=self.random_state)
        self.pca_joint = PCA(n_components=self.n_components, random_state=self.random_state)

    def fit_transform(
        self,
        metabolomics: np.ndarray,
        bulk_rna: np.ndarray,
        spatial_trans: np.ndarray,
        sc_rna: np.ndarray,
    ) -> np.ndarray:
        """Fit the integrator on the provided datasets and return the joint latent representation.

        Each dataset must have the same number of samples (rows).  Features
        (columns) can differ across datasets.

        Parameters
        ----------
        metabolomics : np.ndarray
            Array of shape (n_samples, n_metabolite_features).
        bulk_rna : np.ndarray
            Array of shape (n_samples, n_bulk_genes).
        spatial_trans : np.ndarray
            Array of shape (n_samples, n_spatial_features).  Typically
            summarised spot‑level counts or aggregated per sample.
        sc_rna : np.ndarray
            Array of shape (n_samples, n_sc_features).  This could be
            pseudo‑bulk gene expression aggregated across single cells for
            each sample.

        Returns
        -------
        integrated_latent : np.ndarray
            Low‑dimensional representation of shape (n_samples, n_components).
        """
        n_samples = metabolomics.shape[0]
        # sanity checks
        if not (bulk_rna.shape[0] == n_samples and spatial_trans.shape[0] == n_samples and sc_rna.shape[0] == n_samples):
            raise ValueError("All input datasets must have the same number of samples (rows).")

        # First level PCA per view
        latent_metabolomics = self.pca_metabolomics.fit_transform(metabolomics)
        latent_bulk = self.pca_bulk_rna.fit_transform(bulk_rna)
        latent_spatial = self.pca_spatial.fit_transform(spatial_trans)
        latent_sc = self.pca_sc_rna.fit_transform(sc_rna)

        # Concatenate latent spaces (early integration)
        concatenated = np.hstack([
            latent_metabolomics,
            latent_bulk,
            latent_spatial,
            latent_sc,
        ])

        # Final PCA to learn joint representation (intermediate integration)
        integrated_latent = self.pca_joint.fit_transform(concatenated)
        return integrated_latent



def demo(n_samples: int = 100, random_state: Optional[int] = None) -> Tuple[np.ndarray, MultiOmicsIntegrator]:
    """Demonstrate the integrator on synthetic data.

    Generates random datasets for each omics layer, fits the integrator and
    returns the learned latent representation along with the integrator instance.

    Parameters
    ----------
    n_samples : int, default=100
        Number of synthetic samples to generate.
    random_state : Optional[int]
        Seed for reproducibility.

    Returns
    -------
    integrated_latent : np.ndarray
        Latent representation of shape (n_samples, n_components).
    integrator : MultiOmicsIntegrator
        The fitted integrator.
    """
    rng = np.random.default_rng(seed=random_state)
    # Simulate synthetic data with different feature sizes
    metabolomics = rng.standard_normal((n_samples, 50))
    bulk_rna = rng.standard_normal((n_samples, 1000))
    spatial_trans = rng.standard_normal((n_samples, 200))
    sc_rna = rng.standard_normal((n_samples, 500))

    integrator = MultiOmicsIntegrator(n_components=10, random_state=random_state)
    integrated_latent = integrator.fit_transform(metabolomics, bulk_rna, spatial_trans, sc_rna)
    return integrated_latent, integrator


if __name__ == "__main__":
    # Run demonstration when executed as a script
    latent, integrator = demo(n_samples=50, random_state=42)
    print("Integrated latent representation shape:", latent.shape)
