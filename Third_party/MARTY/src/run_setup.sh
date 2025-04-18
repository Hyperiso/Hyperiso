#!/bin/bash
cd /home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY
script="/home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/MARTY_INSTALL/lib/libmarty.so"
if [ ! -f "$script" ]; then
echo " Installation of marty started "
mkdir MARTY_INSTALL
cd src/MARTY
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/MARTY_INSTALL
sed -i 's/#include "abstract.h"/#include <algorithm>\n#include "abstract.h"/g' /home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/src/MARTY/src/csl/abreviation.h
sed -i 's/#include "vector.h"/#include <algorithm>\n#include "vector.h"/g' /home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/src/MARTY/src/csl/diagonalization.h
sed -i 's/#pragma once/#include <algorithm>\n#pragma once/g' /home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/src/MARTY/src/csl/hardFactor.h
sed -i 's/#include <map>/#include <algorithm>\n#include <map>/g' /home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/src/MARTY/src/csl/index.h
sed -i 's/#include "abstract.h"/#include <type_traits>\n#include "parent.h"\n#include "abstract.h"/g' /home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/src/MARTY/src/csl/replace.h
sed -i 's/\([^:]\)\bfind_if\b/\1std::find_if/g' /home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/src/MARTY/src/csl/variableParent.cpp
sed -i 's/#include "abstract.h"/#include <algorithm>\n#include "abstract.h"/g' /home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/src/MARTY/src/csl/variableParent.h
sed -i 's/#include <map>/#include <algorithm>\n#include <map>/g' /home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/src/MARTY/src/grafed/core/latexLink.h
sed -i 's/#include <map>/#include <algorithm>\n#include <map>/g' /home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/src/MARTY/src/grafed/gui/latexLink.h
make
make install
else
echo "Installation of Marty already done" 
fi
bashrc_file="$HOME/.bashrc"
if ! grep -q "export PATH=\$PATH:/home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/MARTY_INSTALL/bin" "$bashrc_file"; then
        echo "export PATH=\$PATH:/home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/MARTY_INSTALL/bin" >> "$bashrc_file"
        fi

        if ! grep -q "export CPATH=\$CPATH:/home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/MARTY_INSTALL/include" "$bashrc_file"; then
            echo "export CPATH=\$CPATH:/home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/MARTY_INSTALL/include" >> "$bashrc_file"
        fi

        if ! grep -q "export C_INCLUDE_PATH=\$C_INCLUDE_PATH:/home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/MARTY_INSTALL/include" "$bashrc_file"; then
            echo "export C_INCLUDE_PATH=\$C_INCLUDE_PATH:/home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/MARTY_INSTALL/include" >> "$bashrc_file"
        fi

        if ! grep -q "export LIBRARY_PATH=\$LIBRARY_PATH:/home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/MARTY_INSTALL/lib" "$bashrc_file"; then
            echo "export LIBRARY_PATH=\$LIBRARY_PATH:/home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/MARTY_INSTALL/lib" >> "$bashrc_file"
        fi

        if ! grep -q "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/MARTY_INSTALL/lib" "$bashrc_file"; then
            echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/home/nfardeau/Hyperiso/Hyperiso/Hyperiso/Hyperiso/core/src/ExternalIntegration/../../../../../Third_party/MARTY/MARTY_INSTALL/lib" >> "$bashrc_file"
        fi

        source "$bashrc_file"

        echo "Environment variable are uptodate in the .bashrc file."
    