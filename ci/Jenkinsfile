pipeline {
    agent any

    stages {
        stage('Checkout') {
            steps {
                // Récupère le code source du système de contrôle de version
                checkout scm
            }
        }

        stage('Build') {
            steps {
                // Commandes pour construire avec cmake et tout
                sh 'mkdir -p build && cd build && cmake .. && make'
            }
        }

        stage('Test') {
            steps {
                // Exécute les tests unitaires et statistiques
                sh './build/runTests' 
            }
        }

        stage('Results') {
            steps {
                // Collecte et affiche les résultats des tests
                junit '**/test-results/*.xml' 
            }
        }
    }

    post {
        always {
            // Actions à effectuer après le pipeline, nettoyage...
            echo 'Le pipeline est terminé.'
        }
    }
}
