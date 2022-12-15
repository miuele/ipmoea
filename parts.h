#ifndef IPMOEA_PARTS_H_INCLUDED
#define IPMOEA_PARTS_H_INCLUDED

namespace ipmoea {

#include "zsequence.h"

namespace parts {

template <class T>
struct integrator {
	constexpr integrator(T time_step)
		: prev_value_{}, integ_{}, time_step_(time_step)
	{}

	constexpr T operator()(T value) {
		float integ =
			integ_ + time_step_ * (value + prev_value_) / 2;

		prev_value_ = value;
		integ_ = integ;

		return integ;
	}

private:
	T prev_value_;
	T integ_;
	T time_step_;
};

template <class T>
struct hpf_first_order {
	constexpr hpf_first_order(T time_step, T cutoff)
		: prev_value_{}, prev_output_{}, crit_(cutoff * time_step)
	{}

	constexpr T operator()(T value) {
		const float output =
			((2 - crit_)*prev_output_ + 2*(value - prev_value_)) / (2 + crit_);

		prev_value_ = value;
		prev_output_ = output;

		return output;
	}

private:
	T prev_value_;
	T prev_output_;
	T crit_;
};

template <class T>
struct lpf_first_order {
	constexpr lpf_first_order(T time_step, T cutoff)
		: prev_value_{}, prev_output_{}, crit_(cutoff * time_step)
	{}

	constexpr T operator()(T value) {
		const float output =
			((2 - crit_)*prev_output_ + crit_*(value + prev_value_)) / (2 + crit_);

		prev_value_ = value;
		prev_output_ = output;

		return output;
	}

private:
	T prev_value_;
	T prev_output_;
	T crit_;
};

}

}

#endif
