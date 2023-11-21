class Logger {
public:
    static Logger& getInstance() {
        static Logger instance;
        return instance;
    }

    void log(const std::string& message) {
        // Implémentation de la logique de log
    }

private:
    Logger() {}  // Constructeur privé
};
