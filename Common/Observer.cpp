#include <iostream>
#include <string>
#include <list>

//this will be used one day
class IObserver {
public:
    virtual ~IObserver() {}
    virtual void update(const std::string &messageFromSubject) = 0;
};

class ISubject {
public:
    virtual ~ISubject() {}
    virtual void attach(IObserver *observer) = 0;
    virtual void detach(IObserver *observer) = 0;
    virtual void notify() = 0;
};

class ConcreteSubject : public ISubject {
    std::list<IObserver *> observers;
    std::string message;

public:
    virtual ~ConcreteSubject() {}

    void attach(IObserver *observer) override {
        observers.push_back(observer);
    }

    void detach(IObserver *observer) override {
        observers.remove(observer);
    }

    void notify() override {
        for (IObserver *observer : observers) {
            observer->update(message);
        }
    }

    void createMessage(std::string message = "Empty") {
        this->message = message;
        notify();
    }
};

class ConcreteObserver : public IObserver {
    std::string messageFromSubject;
    ConcreteSubject &subject;
    static int static_number;
    int number;

public:
    ConcreteObserver(ConcreteSubject &subject) : subject(subject) {
        this->subject.attach(this);
        number = ++static_number;
    }

    virtual ~ConcreteObserver() {
        subject.detach(this);
    }

    void update(const std::string &messageFromSubject) override {
       this->messageFromSubject = messageFromSubject;
        printInfo();
    }

    void printInfo() {
        std::cout << "Observer " << number << ": a new message is available --> " << messageFromSubject << std::endl;
    }
};
