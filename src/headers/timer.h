#ifndef TIMER_H_
#define TIMER_H_

#include <boost/chrono.hpp>

void printDuration(const char* title, boost::chrono::milliseconds duration) {
	//std::cout << title << " Duration: " << duration.count() << std::endl;
	std::cout << duration.count() << std::endl;
}

#endif  // TIMER_H_