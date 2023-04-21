# Standard imports
import os

# Third-party imports
import numpy as np
from scipy.stats.qmc import scale, LatinHypercube
import camb
import hmcode
import pandas as pd

### Parameters ###

# Directories
datasets_dir = "resources/datasets"
campaign_dir = "resources/campaigns/cosmology"

grid_base = "grid"
train_base = "cosmo"
eval_base = "eval"

# Wavenumbers
kmin, kmax = 1e-3, 1e1  # Minimum and maximum wavenumber [h/Mpc]
nk = 100  # Number of wavenumbers
zs = [0.]  # Redshifts

# CAMB parameters
kmax_CAMB = 10.*kmax

# Random numbers
seed = 42

# Parameter ranges
Omega_c_min, Omega_c_max, Omega_c_def = 0.2, 0.3, 0.255
Omega_b_min, Omega_b_max, Omega_b_def = 0.035, 0.055, 0.045
Omega_k_min, Omega_k_max, Omega_k_def = -0.05, 0.05, 0.
h_min, h_max, h_def = 0.5, 0.9, 0.7
ns_min, ns_max, ns_def = 0.91, 1.01, 0.96
sigma_8_min, sigma_8_max, sigma_8_def = 0.7, 0.9, 0.8
w0_min, w0_max, w0_def = -1.3, -0.7, -1.
wa_min, wa_max, wa_def = -1.73, 1.28, 0.
m_nu_min, m_nu_max, m_nu_def = 0., 1., 0.

# Vary these parameters?
vary_Omega_c = True
vary_Omega_b = True
vary_Omega_k = False
vary_h = True
vary_ns = True
vary_sigma_8 = True
vary_w0 = False
vary_wa = False
vary_m_nu = False

# Number of cosmologies
ntrain = 100
neval = 0
Latin_sampling = True

# Learning choices
power_ratio = False
power_log = False

### ###

###  Calculations ###

ratio_bit = '_ratio' if power_ratio else ''
log_bit = '_log' if power_log else ''
latin_bit = '_latin' if Latin_sampling else ''

# File names
grid_file = grid_base+'.csv'
train_file = train_base+latin_bit+ratio_bit+log_bit+'.csv'
eval_file = eval_base+ratio_bit+log_bit+'.csv'

# File paths
grid_file = os.path.join(campaign_dir, grid_file)
train_file = os.path.join(datasets_dir, train_file)
eval_file = os.path.join(campaign_dir, eval_file)

# Seed random number generator
rng = np.random.default_rng(seed=seed)

# Ranges
k = np.logspace(np.log10(kmin), np.log10(kmax), nk)  # Wavenumbers [h/Mpc]
k_dict = {f"k{ik}": _k for ik, _k in enumerate(k)}
df = pd.DataFrame(k_dict, index=[0])
df.to_csv(grid_file, index=False)

# Create dictionary for dataframe
column_names = [
    "z",
    "Omega_c",
    "Omega_b",
    "Omega_k",
    "h",
    "ns",
    "sigma_8",
    "w0",
    "wa",
    "m_nu",
]

### ###

# Loop over train and eval data
for ival, (n, file) in enumerate(zip([ntrain, neval], [train_file, eval_file])):

    # Initialize dataframe
    df_columns = column_names+list(k_dict.keys())
    df = {name: [] for name in df_columns}

    # Generate random samples
    if Latin_sampling and ival == 0:
        sampler = LatinHypercube(d=9, seed=seed, scramble=False)
        samples = sampler.random(n)
        lower_bounds = [Omega_c_min, Omega_b_min, Omega_k_min,
                        h_min, ns_min, sigma_8_min, w0_min, wa_min, m_nu_min]
        upper_bounds = [Omega_c_max, Omega_b_max, Omega_k_max,
                        h_max, ns_max, sigma_8_max, w0_max, wa_max, m_nu_max]
        samples = scale(samples, lower_bounds, upper_bounds)
        Omega_cs = samples[:, 0]
        Omega_bs = samples[:, 1]
        Omega_ks = samples[:, 2]
        hs = samples[:, 3]
        nss = samples[:, 4]
        sigma_8s = samples[:, 5]
        w0s = samples[:, 6]
        was = samples[:, 7]
        m_nus = samples[:, 8]
    else:
        Omega_cs = rng.uniform(Omega_c_min, Omega_c_max, n)
        Omega_bs = rng.uniform(Omega_b_min, Omega_b_max, n)
        Omega_ks = rng.uniform(Omega_k_min, Omega_k_max, n)
        hs = rng.uniform(h_min, h_max, n)
        nss = rng.uniform(ns_min, ns_max, n)
        sigma_8s = rng.uniform(sigma_8_min, sigma_8_max, n)
        w0s = rng.uniform(w0_min, w0_max, n)
        was = rng.uniform(wa_min, wa_max, n)
        m_nus = rng.uniform(m_nu_min, m_nu_max, n)

    # Loop over cosmologies and generate random samples
    for icos in range(n):
        Omega_c = Omega_cs[icos] if vary_Omega_c else Omega_c_def
        Omega_b = Omega_bs[icos] if vary_Omega_b else Omega_b_def
        Omega_k = Omega_ks[icos] if vary_Omega_k else Omega_k_def
        h = hs[icos] if vary_h else h_def
        ns = nss[icos] if vary_ns else ns_def
        sigma_8 = sigma_8s[icos] if vary_sigma_8 else sigma_8_def
        As = 2e-9  # Needs to be set and corrected later
        w0 = w0s[icos] if vary_w0 else w0_def
        wa = was[icos] if vary_wa else wa_def
        # while True:  # Ensure that dark energy does not dominate the early universe
        #     wa = rng.uniform(wa_min, wa_max) if vary_wa else wa_def
        #     if w0+wa < 0.:
        #         break
        m_nu = m_nus[icos] if vary_m_nu else m_nu_def

        # Add cosmological parameters to dataframe
        for z in zs:
            params = [z, Omega_c, Omega_b, Omega_k,
                      h, ns, sigma_8, w0, wa, m_nu]
            for parameter_name, parameter_value in zip(column_names, params):
                df[parameter_name].append(parameter_value)

        # Write cosmological parameters to screen
        print(f"Cosmology: {icos:d}; (Om_c, Om_b, Om_k, h, ns, sig8, w0, wa, m_nu) = \
({Omega_c:.2}, {Omega_b:.2}, {Omega_k:.2}, {h:.2}, {ns:.2}, {sigma_8:.2}, {w0:.2}, {wa:.2}, {m_nu:.2})")

        # CAMB parameters
        H0 = 100.*h           # Hubble constant [km/s/Mpc]
        ombh2 = Omega_b*h**2  # Physical baryon density
        omch2 = Omega_c*h**2  # Physical CDM density

        # Run CAMB to get linear spectrum
        parameters = camb.CAMBparams(WantCls=False)
        parameters.set_cosmology(
            H0=H0, ombh2=ombh2, omch2=omch2, omk=Omega_k, mnu=m_nu)
        parameters.set_dark_energy(w=w0, wa=wa, dark_energy_model='ppf')
        parameters.InitPower.set_params(As=As, ns=ns)
        parameters.set_matter_power(redshifts=zs, kmax=kmax_CAMB)
        results = camb.get_results(parameters)
        sigma_8_init = results.get_sigma8_0()
        As *= (sigma_8/sigma_8_init)**2
        parameters.InitPower.set_params(As=As, ns=ns)
        results = camb.get_results(parameters)

        # Non-linear power spectrum (from HMcode-2020)
        Pk = hmcode.power(k, zs, results, verbose=False)
        if power_ratio:
            for iz, z in enumerate(zs):  # Linear matter power spectrum
                Pk_lin = results.get_matter_power_interpolator(
                    nonlinear=False).P(z, k)
                Pk[iz, :] /= Pk_lin
        if power_log:
            Pk = np.log(Pk)
        for iz, _ in enumerate(zs):
            for ik, _ in enumerate(k):
                df[f"k{ik}"].append(Pk[iz, ik])

    # Construct dataframes and save to file
    df = pd.DataFrame(df)
    df.to_csv(file, index=False)
    print()
