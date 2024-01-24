#include "Wilson.h"

void WilsonInitializer::init() {
    
}

WilsonManager *WilsonManager::GetInstance(double mu_match)
{
    if (!WilsonManager::instance) {
        if (mu_match == 0.0) {
            // Log an error
        } else {
            WilsonManager::instance = new WilsonManager(mu_match);
        }
    }
    return WilsonManager::instance;
}
