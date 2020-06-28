#include "Spinlock.h"

Spinlock::Spinlock() {
    this->value.store(true);
}

void Spinlock::lock() {
    bool expected = true;
    bool changed;

    do{
        changed = this->value.compare_exchange_strong(expected, false);
        expected = true;
    } while(!changed);

}

void Spinlock::unlock() {
    bool expected = false;
    if(!this->value.compare_exchange_strong(expected, true))
        throw "Something's not right";
}
