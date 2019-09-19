#ifndef __UTILITY_TIMER_H__
#define __UTILITY_TIMER_H__

#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>

using namespace std;
using namespace chrono;

namespace Utility {
class Timer {
 public:
	decltype(high_resolution_clock::now()) now_time;

	void reset()
	{
		now_time = high_resolution_clock::now();
	}

	double stop()
	{
		const auto last_time = now_time;
		now_time = high_resolution_clock::now();

		const auto elapsed_read = duration_cast<duration<double>>(now_time - last_time).count();

		return elapsed_read;
	}

	double passed() const
	{
		const auto current_time = high_resolution_clock::now();

		const auto elapsed_read = duration_cast<duration<double>>(current_time - now_time).count();

		return elapsed_read;
	}

	static auto now() -> decltype(high_resolution_clock::now())
	{
		return high_resolution_clock::now();
	}

	static double elapsed(decltype(high_resolution_clock::now()) start)
	{
		const auto current_time = high_resolution_clock::now();

		const auto elapsed_read = duration_cast<duration<double>>(current_time - start).count();

		return elapsed_read;
	}

	static std::string NowTime()
	{
		const auto now = chrono::system_clock::now();
		const auto now_c = chrono::system_clock::to_time_t(now);

#if defined(__GNUC__)
		std::stringstream ss;
		ss << std::put_time(std::localtime(&now_c), "%F %T");
		return ss.str();
#else
		return "not supposed";
#endif
	}

	Timer()
	{
		reset();
	}
};

class CPUTimer
{
public:
  clock_t now_time;

  void reset()
  {
    now_time = std::clock();
  }

  double stop()
  {
    const auto last_time = now_time;
    now_time = std::clock();

    const auto elapsed_read = (static_cast<double>(now_time) - last_time) / CLOCKS_PER_SEC;

    return elapsed_read;
  }

  double passed() const
  {
    const auto current_time = std::clock();

    const auto elapsed_read = (static_cast<double>(current_time) - now_time) / CLOCKS_PER_SEC;

    return elapsed_read;
  }

  static clock_t now()
  {
    return std::clock();
  }

  CPUTimer()
  {
    reset();
  }
};

}  // namespace Utility
#endif
