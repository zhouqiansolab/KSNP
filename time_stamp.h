//
// Created by ixiaohu on 2022/8/19.
//

#ifndef KSNP_TIME_STAMP_H
#define KSNP_TIME_STAMP_H

#include <sys/resource.h>
#include <sys/time.h>
#include <map>
#include <stack>
#include <string>

struct Time1 {
	double cpu, real;
	Time1():cpu(0), real(0) {}
	Time1(double c, double r):cpu(c), real(r) {}
	Time1 operator + (const Time1 &t) const {
		return {this->cpu + t.cpu, this->real + t.real};
	}
	Time1 operator - (const Time1 &t) const {
		return {this->cpu - t.cpu, this->real - t.real};
	}
	void operator += (const Time1 &t) {
		this->cpu += t.cpu;
		this->real += t.real;
	}
};

class TimeStamp {
private:
	std::stack<Time1> time_stack;
	std::map<std::string, Time1> time_map;

	static inline double cputime() {
		struct rusage r{};
		getrusage(RUSAGE_SELF, &r);
		return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
	}

	static inline double realtime() {
		struct timeval tp{};
		gettimeofday(&tp, nullptr);
		return tp.tv_sec + tp.tv_usec * 1e-6;
	}

public:
	TimeStamp() = default;

	static inline Time1 get_time() { return {cputime(), realtime()}; }

	void push() { time_stack.push(get_time()); }

	void record(const std::string &key) {
		auto now = get_time();
		if (time_stack.empty()) std::abort();
		auto top = time_stack.top(); time_stack.pop();
		if (time_map.find(key) == time_map.end()) {
			time_map[key] = now - top;
		} else {
			time_map[key] += (now - top);
		}
	}

	void output() {
		for (const auto &a : time_map) {
			const auto &stage = a.first;
			const auto &time = a.second;
			fprintf(stderr, "[Stage,CPU,Real]\t%s\t%.2f\t%.2f\n", stage.c_str(), time.cpu, time.real);
		}
	}
};

/** A time monitor at instruction level. rdtsc() / CPU frequency = seconds.
 * Class TimeStamp has a large overhead, do not use it in frequently calling. */
#if defined(__GNUC__) && !defined(__clang__)
	#if defined(__i386__)
		static inline unsigned long long __rdtsc(void) {
		    unsigned long long int x;
		    __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
		    return x;
		}
	#elif defined(__x86_64__)
		static inline unsigned long long __rdtsc(void) {
			unsigned hi, lo;
			__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
			return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
		}
	#endif
#endif

#endif //KSNP_TIME_STAMP_H

