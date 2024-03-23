import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("build/quark_mass_running.csv")

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

axes[0].plot(data['Q'], data['RunningMassTop'], label='Running Mass Top', color='blue')
axes[0].set_xlabel('Q (GeV)')
axes[0].set_ylabel('Running Mass (GeV)')
axes[0].set_title('Running of Top Quark Mass with Scale Q')
axes[0].legend()
axes[0].grid(True)

axes[1].plot(data['Q'], data['RunningMassBottom'], label='Running Mass Bottom', color='red')
axes[1].set_xlabel('Q (GeV)')
axes[1].set_ylabel('Running Mass (GeV)')
axes[1].set_title('Running of Bottom Quark Mass with Scale Q')
axes[1].legend()
axes[1].grid(True)

axes[2].plot(data['Q'], data['RunningMassCharm'], label='Running Mass Charm', color='green')
axes[2].set_xlabel('Q (GeV)')
axes[2].set_ylabel('Running Mass (GeV)')
axes[2].set_title('Running of Charm Quark Mass with Scale Q')
axes[2].legend()
axes[2].grid(True)

plt.tight_layout()

plt.savefig("plot/mass/quark_mass_running_plot.pdf")

# plt.show()
