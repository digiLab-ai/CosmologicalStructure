# Cosmological structure

<!-- <p align="center">
    <img src="./resources/images/density.png" width="800" height="200" />
</p> -->
![image](./resources/images/density.png)

![digiLab](./resources/images/digiLab_badge.svg)
![twinLab](./resources/images/twinLab_badge.svg)
![slack](https://img.shields.io/badge/slack-@digilabglobal-purple.svg?logo=slack)](https://digilabglobal.slack.com)

Here, we use `twinLab` to create an emulator for the statistical properties of the distribution of structure across the Universe.

The large-scale distribution of billions of galaxies contains a wealth of information about the origin, expansion, and contents of the Universe. For example, the galaxy distribution is sensitive to the amounts of dark energy and dark matter, the current and historial expansion speed, and the process of inflation that occurred soon after the big bang and which seeded the diverse array of cosmological structure, including all galaxies, stars and planets, that we see today. Galaxies exist in dense clumps, called groups, clusters, super clusters, depending on the number, at the nodes of the density distribution. The underlying density is governed by the dynamics and properties of dark matter, which defines a skeleton along which galaxies flow and eventually cluster.

To extract cosmological information from the galaxy distribution requires precise models of the statistical properties of the distribution as a function of the underlying parameters. Analytical theories, developed over the last 30 years, work at early times and on extremely large scales, where perturbations to the mean density are small. However, on the (comparatively small) scale of galaxies, the perturbations are huge and modelling their distribution can only be accurately achieved using expensive $n$-body simulations (the image above shows a slice of density through one such simulation). High-fidelity simulations can take up to a month to run distributed across tens of thousands of cores on the top super computers in the world. It is impractical to run accurate simulations at all points in parameter space, especially since the parameter-space of models under investigation is always expanding. In modern cosmology this includes the space of exotic dark energy models, beyond-Einstein gravity theories, and non-standard particle physics models for dark matter and neutrinos.

In this example, we use `twinLab` to create an emulator for the matter power spectrum, a statistical quantity that contains a subset of the possible information from the clustering distribution of galaxies. The power spectrum can both be computed via simulation and measured in observational datasets. The training of the `twinLab` emulator is performed online (in the cloud) and is completed in a matter of minutes. Once trained, the emulator can be used for extremely rapid power-spectrum evaluation across parameter space in a way that interpolates and extrapolates reasonably. The major benefit of using `twinLab` is that we get an accurate estimate of our model uncertainty for free, so that we know exactly how much we should trust our trained surrogate. The model is trained on approximate simulation data that occupy a Latin-hypercube distribution across five parameters of interest to cosmologists. The model can be rapidly retrained if necessary, and additional parameters can be incorporated if desired.

## Quick start

Clone the repository and change directory to the project root:
```shell
git clone https://github.com/digiLab-ai/CosmologicalStructure.git
cd CosmologicalStructure
```

Install the dependencies:
```shell
poetry install
```

Copy the `.env.example` file to `.env` 
```shell
cp .env.example .env
```
and fill out your `twinLab` login details in `.env`

Run the [demo notebook](./demo.ipynb):
```shell
poetry run jupyter notebook demo.ipynb
```
## App

You can visualise the results of the emulator here: [https://alexander-mead.github.io./universe.html](https://alexander-mead.github.io./universe.html)

![image](https://github.com/digiLab-ai/CosmologicalStructure/assets/9140961/71fed753-8a2a-43a9-9642-839557517c05)
