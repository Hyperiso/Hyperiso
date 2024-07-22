import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def compare_coefficients(path1, path2, title, model = "sm"):
    df1 = pd.read_csv(path1)
    df2 = pd.read_csv(path2)
    print(title, model)
    print(df1)
    print(df2)
    real_columns = [col for col in df1.columns if 'real' in col]
    imag_columns = [col for col in df1.columns if 'imag' in col]

    diff_real = (df1[real_columns].iloc[0] - df2[real_columns].iloc[0]) / df2[real_columns].iloc[0] * 100
    diff_imag = (df1[imag_columns].iloc[0] - df2[imag_columns].iloc[0]) / df2[imag_columns].iloc[0] * 100

    width = 0.35
    x = np.arange(len(real_columns))

    fig, ax = plt.subplots(2, 2, figsize=(16, 12))

    ax[0, 0].bar(x - width/2, df1[real_columns].iloc[0], width, label='Hyperiso', color='blue', alpha=0.6)
    ax[0, 0].bar(x + width/2, df2[real_columns].iloc[0], width, label='Superiso', color='red', alpha=0.6)
    ax[0, 0].set_title(rf'Real part {title} Wilson Coefficient ($\mu = 81$ GeV) Comparison')
    ax[0, 0].set_xticks(x)
    ax[0, 0].set_xticklabels(real_columns, rotation=90)
    ax[0, 0].legend()

    ax[0, 1].bar(x, diff_real, width, color='green')
    ax[0, 1].set_title(f'Percentage Difference (Real Part) for {title}')
    ax[0, 1].set_xticks(x)
    ax[0, 1].set_xticklabels(real_columns, rotation=90)
    ax[0, 1].axhline(0, color='black', linewidth=0.8)

    ax[1, 0].bar(x - width/2, df1[imag_columns].iloc[0], width, label='Hyperiso', color='blue', alpha=0.6)
    ax[1, 0].bar(x + width/2, df2[imag_columns].iloc[0], width, label='Superiso', color='red', alpha=0.6)
    ax[1, 0].set_title(rf'imaginary part {title} Wilson Coefficient ($\mu = 81$ GeV) Comparison')
    ax[1, 0].set_xticks(x)
    ax[1, 0].set_xticklabels(imag_columns, rotation=90)

    ax[1, 1].bar(x, diff_imag, width, color='green')
    ax[1, 1].set_title(f'Percentage Difference (Imaginary Part) for {title}')
    ax[1, 1].set_xticks(x)
    ax[1, 1].set_xticklabels(imag_columns, rotation=90)
    ax[1, 1].axhline(0, color='black', linewidth=0.8)

    plt.tight_layout()
    plt.savefig(f"plot/wilson/Wilson_comparison_{title}_{model}.png")



compare_coefficients('csv/susy/WilsonCoefficients_LO.csv', 'csv/superiso/susy/WilsonCoefficients_LO.csv', 'LO', model = "susy")
compare_coefficients('csv/susy/WilsonCoefficients_NLO.csv', 'csv/superiso/susy/WilsonCoefficients_NLO.csv', 'NLO', model = "susy")
compare_coefficients('csv/susy/WilsonCoefficients_NNLO.csv', 'csv/superiso/susy/WilsonCoefficients_NNLO.csv', 'NNLO', model = "susy")

compare_coefficients('csv/sm/WilsonCoefficients_LO.csv', 'csv/superiso/sm/WilsonCoefficients_LO.csv', 'LO', model = "sm")
compare_coefficients('csv/sm/WilsonCoefficients_NLO.csv', 'csv/superiso/sm/WilsonCoefficients_NLO.csv', 'NLO', model = "sm")
compare_coefficients('csv/sm/WilsonCoefficients_NNLO.csv', 'csv/superiso/sm/WilsonCoefficients_NNLO.csv', 'NNLO', model = "sm")

compare_coefficients('csv/thdm/WilsonCoefficients_LO.csv', 'csv/superiso/thdm/WilsonCoefficients_LO.csv', 'LO', model = "thdm")
compare_coefficients('csv/thdm/WilsonCoefficients_NLO.csv', 'csv/superiso/thdm/WilsonCoefficients_NLO.csv', 'NLO', model = "thdm")
compare_coefficients('csv/thdm/WilsonCoefficients_NNLO.csv', 'csv/superiso/thdm/WilsonCoefficients_NNLO.csv', 'NNLO', model = "thdm")


compare_coefficients('csv/sm/WilsonCoefficients_PRIMECQ_LO.csv', 'csv/superiso/sm/WilsonCoefficients_PRIMECQ_LO.csv', 'LO prime and scalar', model = "sm")
compare_coefficients('csv/thdm/WilsonCoefficients_PRIMECQ_LO.csv', 'csv/superiso/thdm/WilsonCoefficients_PRIMECQ_LO.csv', 'LO prime and scalar', model = "thdm")
compare_coefficients('csv/susy/WilsonCoefficients_PRIMECQ_LO.csv', 'csv/superiso/susy/WilsonCoefficients_PRIMECQ_LO.csv', 'LO prime and scalar', model = "susy")
compare_coefficients('csv/susy/WilsonCoefficients_PRIMECQ_NLO.csv', 'csv/superiso/susy/WilsonCoefficients_PRIMECQ_NLO.csv', 'NLO prime and scalar', model = "susy")