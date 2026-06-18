# Hyperiso C++ examples

Ces exemples C++ suivent la structure et les commentaires des exemples Python fournis.

## Compilation

```bash
mkdir -p build
cd build
cmake ..
cmake --build .
```

Les exemples `Statistic/*` nécessitent que l'installation Hyperiso exporte aussi la cible `Hyperiso::StatisticLib`.
Si cette cible existe mais pointe vers un include directory absent, CMake ignore automatiquement les targets Statistic au lieu de bloquer tous les autres exemples.

## Organisation

- `Core/` : initialisation, lecture/écriture de paramètres, blocs dépendants, QCD runner.
- `Wilson/` : coefficients de Wilson, matching/running, exemple MARTY, scans développeur.
- `Observable/` : observables physiques, configurations de désintégration, scan paramétrique.
- `Statistic/` : incertitudes, MLE, contours et scan de vraisemblance.
- `*/dev_examples/` : exemples orientés développement/diagnostic pour tester rapidement des APIs internes.

## Correctif C++ important

Les appels à `WilsonBuildConfig({...}, ...)` sont volontairement typés avec `std::unordered_set<WGroup>{...}` ou `std::unordered_set<WGroupId>{...}`.
Sans ce type explicite, GCC peut hésiter entre les deux constructeurs suivants :

```cpp
WilsonBuildConfig(std::unordered_set<WGroup>, ...);
WilsonBuildConfig(std::unordered_set<WGroupId>, ...);
```

Donc il faut écrire par exemple :

```cpp
WilsonBuildConfig config(
    std::unordered_set<WGroup>{WGroup::B, WGroup::BScalar},
    81.0,
    2.0,
    QCDOrder::LO
);
```

## Notes sur les exemples dynamiques

- `Wilson/dev_examples/wilson_dynamic_ids_example.cpp` montre comment ajouter des groupes Wilson à partir de noms (`GroupMapper::id_of`) et comment requêter des coefficients builtin à partir de noms (`WCoefMapper::id_of`).
- `Observable/dev_examples/observable_dynamic_ids_example.cpp` montre comment ajouter des observables builtin à partir de noms (`ObservableMapper::id_of`) et comment enregistrer un `ObservableId` custom.

Les exemples restent compatibles avec l'API installée actuelle : ils évitent d'appeler directement les overloads dynamiques qui nécessitent un patch du code Hyperiso.

## Patch proposé pour un vrai support dynamique

L'archive contient aussi :

```text
patches/hyperiso_dynamic_enum_support.patch
```

Ce patch ajoute les overloads et le routage nécessaires pour mieux supporter les enums dynamiques :

- `WilsonRequest` stocke `WGroupId` / `WCoefId`, avec un constructeur de compatibilité `WGroup` / `WCoef`.
- `WilsonInterface` ajoute des overloads `WGroupId` / `WCoefId` pour les getters simples : `getM`, `getFM`, `getR`, `getFR`, `getSM`, `getSR`, etc.
- `ObsManager::add_obs(ObservableId)` et `ObsManager::add_obs(BinnedObservableId)` routent par `DecayMapper::get_decay_id_or_throw(...)` au lieu de repasser par `ObservableMapper::enum_of(...).value()`.
- `ObsManager::select_decay(ObservableId)` sélectionne aussi le decay par `DecayId`, ce qui permet aux observables custom enregistrées dans le mapper d'être routées proprement.

Application typique depuis la racine du repo Hyperiso :

```bash
patch -p1 < examples_cpp/patches/hyperiso_dynamic_enum_support.patch
```



## Nouveaux exemples dynamiques par lambdas

- `dev_wilson_lambda_custom_group` montre comment déclarer un `WGroupId` custom, des `WCoefId` custom, puis fournir des lambdas de matching et de running via `CustomWilsonGroupConfig`.
- `dev_observable_lambda_decay` montre comment déclarer un nouveau decay et une nouvelle observable via `LambdaDecayConfig`, puis calculer l'observable avec une lambda utilisateur. La lambda peut lire les Wilson dynamiques via `ctx.W().getFM(group, coeff, ...)`.

Ces deux exemples nécessitent l'archive source `hyperiso_dynamic_lambda_update.zip`, pas seulement l'ancienne installation Hyperiso.
