# Third-party notices

HyperIso includes and links third-party software. The combined source and binary
distribution is provided under **GPL-3.0-or-later**. HyperIso-authored portions
that were originally published under the MIT License remain available under the
MIT terms in `LICENSES/MIT.txt`.

## MinuitCpp

- Component: `Hyperiso/Hyperiso/core/src/ExternalIntegration/third_party/minuit-cpp`
- Upstream: `doughague/minuit-cpp`
- Copyright: Doug Hague and the LCG ROOT Math team, as described upstream
- License: GNU General Public License, version 3 or (at your option) any later version
- License text: `LICENSES/GPL-3.0-or-later.txt`

## 2HDMC

- Component: `Hyperiso/Hyperiso/core/src/ExternalIntegration/third_party/2HDMC`
- Vendored version: 1.8.0 (from `THDM::version`)
- Authors: David Eriksson, Johan Rathsman and Oscar Stål
- Program reference: *2HDMC — Two-Higgs-Doublet Model Calculator*, Computer
  Physics Communications 181 (2010) 189–205 and the accompanying CPC program
  record
- License: GNU General Public License, as stated by the CPC program record. The
  vendored snapshot does not include a standalone upstream license file or a more
  specific SPDX identifier.
- Combined HyperIso distribution: GPL-3.0-or-later, reflecting the GPLv3-or-later
  MinuitCpp requirement. Confirm the exact 2HDMC licence variant with the upstream
  authors or CPC archive before the final journal deposit.
- License text distributed with HyperIso: `LICENSES/GPL-3.0-or-later.txt`

The vendored snapshots must retain their original copyright and attribution
notices. When updating either component, update this file and record the exact
upstream version or commit in the release notes.
