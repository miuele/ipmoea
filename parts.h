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
	constexpr hpf_first_order(T time_step, T cutoff, T init = T{})
		: prev_value_{init}, prev_output_{init}, crit_(cutoff * time_step)
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
	constexpr lpf_first_order(T time_step, T cutoff, T init = T{})
		: prev_value_{init}, prev_output_{init}, crit_(cutoff * time_step)
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

template <class T>
struct differentiator {
	constexpr differentiator(T time_step, T init = T{})
		: time_step_(time_step), prev_value_(init)
	{}

	constexpr T operator()(T value) {
		const float output = (value - prev_value_) / time_step_;

		prev_value_ = value;

		return output;
	}

private:
	T time_step_;
	T prev_value_;
};

template <class T>
struct rate_limiter {
	constexpr rate_limiter(T time_step, T lower_limit, T upper_limit, T init = T{})
		: time_step_(time_step)
		, prev_value_(init)
		, lower_limit_(lower_limit), upper_limit_(upper_limit)
	{}

	constexpr T operator()(T value) {
		const float rate = (value - prev_value_) / time_step_;

		if (rate < lower_limit_) {
			value = prev_value_ + lower_limit_ * time_step_;
		} else if (rate > upper_limit_) {
			value = prev_value_ + upper_limit_ * time_step_;
		}

		prev_value_ = value;

		return value;
	}

private:
	T time_step_;
	T prev_value_;
	T lower_limit_;
	T upper_limit_;
};

}

}

#endif
