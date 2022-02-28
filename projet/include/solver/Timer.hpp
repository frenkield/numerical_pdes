#ifndef SOLVER_TIMER_H
#define SOLVER_TIMER_H

#include <chrono>

using namespace std;
using namespace std::chrono;

class Timer {

private:
    high_resolution_clock::time_point startTime = high_resolution_clock::now();

public:

    void reset() {
        startTime = high_resolution_clock::now();
    }

    void logElapsedTime(string message) {

        auto elapsedTime =
            (duration_cast<duration<double>>)(high_resolution_clock::now() - startTime);

        cout << message << " => " << elapsedTime.count() << " secs" << endl;
    }
};

#endif