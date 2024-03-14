import pandas as pd
import matplotlib.pyplot as plt

df_lo = pd.read_csv("build/WilsonCoefficients_LO.csv")
df_nlo = pd.read_csv("build/WilsonCoefficients_NLO.csv")
df_nnlo = pd.read_csv("build/WilsonCoefficients_NNLO.csv")

alpha_s_nlo = df_nlo['alpha_s'] / (2 * 3.141592653589793)
alpha_s_nnlo = (df_nnlo['alpha_s'] / (2 * 3.141592653589793))**2

coefficients = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']

for coef in coefficients:
    plt.figure(figsize=(18, 6))

    plt.subplot(1, 2, 1)
    plt.plot(df_lo['Q'], df_lo[f'{coef}_real'], label=f'{coef} LO Réel', color='blue')
    plt.plot(df_nlo['Q'], df_nlo[f'{coef}_real'] * alpha_s_nlo, label=f'{coef} NLO Réel (ajusté)', color='green')
    plt.plot(df_nnlo['Q'], df_nnlo[f'{coef}_real'] * alpha_s_nnlo, label=f'{coef} NNLO Réel (ajusté)', color='red')
    plt.xlabel('Q (GeV)')
    plt.ylabel('Valeur Réelle')
    plt.title(f'Évolution de la partie réelle de {coef} avec l\'échelle Q')
    plt.legend()
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(df_lo['Q'], df_lo[f'{coef}_imag'], label=f'{coef} LO Imaginaire', color='blue', linestyle='--')
    plt.plot(df_nlo['Q'], df_nlo[f'{coef}_imag'] * alpha_s_nlo, label=f'{coef} NLO Imaginaire (ajusté)', color='green', linestyle='--')
    plt.plot(df_nnlo['Q'], df_nnlo[f'{coef}_imag'] * alpha_s_nnlo, label=f'{coef} NNLO Imaginaire (ajusté)', color='red', linestyle='--')
    plt.xlabel('Q (GeV)')
    plt.ylabel('Valeur Imaginaire')
    plt.title(f'Évolution de la partie imaginaire de {coef} avec l\'échelle Q')
    plt.legend()
    plt.grid(True)

    plt.savefig(f"plot/wilson/{coef}_running_plot.pdf")
    # plt.show()
