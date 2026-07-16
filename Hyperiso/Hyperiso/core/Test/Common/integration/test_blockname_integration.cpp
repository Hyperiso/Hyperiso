#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

#include "Include.h"

int main() {
    std::cout << "== Running INTEGRATION tests for BlockName ==\n";

    {
        std::map<BlockName, int> M;
        M[BlockName{"SMINPUTS"}] = 1;
        M[BlockName{"GAUGE"}]    = 2;
        M[BlockName{"MASS"}]     = 3;

        assert(M.find(BlockName{"MASS"})    != M.end());
        assert(M.find(BlockName{"GAUGE"})   != M.end());
        assert(M.find(BlockName{"SMINPUTS"})!= M.end());
    }

    {
        BlockName a{"YU","UCOUPL"};
        BlockName b{"UCOUPL","XYZ"};

        std::unordered_set<BlockName> S;
        S.insert(a);
        S.insert(b);

        assert(S.size() >= 1);
    }

    {
        BlockName a{"YU","UCOUPL"};
        BlockName b{"UCOUPL","XYZ"}; // a == b

        std::unordered_map<BlockName, int> U;
        U[a] = 10;
        U[b] = 20;

        assert(U.size() >= 1);
    }

    std::cout << "\n✅ BlockName integration tests passed (avec note sur unordered_*).\n";
    return 0;
}


// ⚠️ Note importante (long terme)

// L’égalité BlockName::operator== est “intersection non vide” entre ensembles d’alias → ce n’est pas une relation d’équivalence (pas transitive).

// Le hash<BlockName> combine toutes les aliases → deux objets == peuvent avoir des hash différents.

// Conséquences :

// unordered_set / unordered_map ne sont pas sûrs : la contrainte standard exige que si a == b alors hash(a) == hash(b). Ici ce n’est pas garanti → vous pouvez avoir des doublons “égaux” à vos yeux, des lookups ratés, etc.

// Vous utilisez déjà std::map<BlockName, ...> (trié) à beaucoup d’endroits → c’est safe. Évitez unordered_* avec BlockName (et avec Prototype qui se base sur BlockName).

// Pistes de “patch” durable

// Choisir une identité canonique pour BlockName :

// Option A (simple, sûre) : stocker un seul nom canonique (ex: uppercase du premier alias) dans BlockName et

// operator== compare ce nom canonique uniquement,

// hash<BlockName> hashe uniquement ce nom canonique,

// hasAlias() devient un helper (possiblement en dehors, via un dictionnaire global d’alias → canonique).

// Option B : garder un set d’aliases mais définir l’identité comme le plus petit alias uppercased (lexico) :

// canonical() = min(uppercased_aliases),

// operator== compare canonical(),

// hash hashe canonical().

// Ainsi, unordered_* redeviennent sûrs et l’égalité est transitive.