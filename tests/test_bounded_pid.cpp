#include <utility>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <cmath>

#include "zsequence.h"
#include "parts.h"
#include "utils.h"

using namespace ipmoea::utils;
using namespace ipmoea::parts;

template <class T>
struct second_order_system_op {

	second_order_system_op(T time_step, T omega_n, T zeta)
		: crit_(time_step * omega_n), zeta_(zeta)
	{}

	auto operator()(iir_input<T, 2> u, iir_output<T, 2> y) {
		T crit_sq = crit_ * crit_;
		return y.template z<0>() = 
			((2*(4 - crit_sq)*y.template z<-1>() - (4 - 4*zeta_*crit_ + crit_sq)*y.template z<-2>())
				+ crit_sq*(u.template z<0>() + 2*u.template z<-1>() + u.template z<-2>()))
			/ (4 + 4*zeta_*crit_ + crit_sq);
	}

private:
	T crit_, zeta_;
};

template <class T>
using second_order_system = iir_filter_from_op<T, 2, 2, second_order_system_op<T>>;

struct simple_pid {

	simple_pid(float time_step, float k_p, float k_i, float k_d, float tf)
		: k_p_(k_p), k_i_(k_i), k_d_(k_d), tf_(tf), integ_(time_step), hpf_(time_step, 1.f / tf),
		input_0_(std::numeric_limits<float>::quiet_NaN())
	{}

	float operator()(float input) {
		if (std::isnan(input_0_)) {
			input_0_ = input;
		}

		return k_p_ * input + k_i_ * integ_(input) + (k_d_ / tf_) * hpf_(input - input_0_);
	}

private:
	float k_p_, k_i_, k_d_;
	float tf_;
	integrator<float> integ_;
	hpf_first_order<float> hpf_;
	float input_0_;
};

struct simple_bounded_pid {

	simple_bounded_pid(float time_step, float lb, float ub, float k_p, float k_i, float k_d, float tf)
		: k_p_(k_p), k_i_(k_i), k_d_(k_d), tf_(tf), lb_(lb), ub_(ub), integ_(time_step), hpf_(time_step, 1.f / tf),
		input_0_(std::numeric_limits<float>::quiet_NaN())
	{}

	float operator()(float input) {

		if (std::isnan(input_0_)) {
			input_0_ = input;
		}

		auto integ = integ_;
		float unbounded = k_p_ * input + k_i_ * integ(input) + (k_d_ / tf_) * hpf_(input - input_0_);

		if (unbounded >= ub_) {
			return ub_;
		} else if (unbounded <= lb_) {
			return lb_;
		} else {
			integ_ = integ;
			return unbounded;
		}
	}

private:
	float k_p_, k_i_, k_d_;
	float tf_;
	float lb_, ub_;
	integrator<float> integ_;
	hpf_first_order<float> hpf_;
	float input_0_;
};

struct simple_bounded_pid_dp {

	simple_bounded_pid_dp(float time_step, float lb, float ub, float k_p, float k_i, float k_d, float tf)
		: k_p_(k_p), k_i_(k_i), k_d_(k_d), tf_(tf), lb_(lb), ub_(ub), integ_(time_step), hpf_(time_step, 1.f / tf),
		pv_0_(std::numeric_limits<float>::quiet_NaN())
	{}

	float operator()(float sp, float pv) {

		float e = sp - pv;

		if (std::isnan(pv_0_)) {
			pv_0_ = pv;
		}

		auto integ = integ_;
		float unbounded = k_p_ * e + k_i_ * integ(e) + (k_d_ / tf_) * hpf_(pv_0_ - pv);

		if (unbounded >= ub_) {
			return ub_;
		} else if (unbounded <= lb_) {
			return lb_;
		} else {
			integ_ = integ;
			return unbounded;
		}
	}

private:
	float k_p_, k_i_, k_d_;
	float tf_;
	float lb_, ub_;
	integrator<float> integ_;
	hpf_first_order<float> hpf_;
	float pv_0_;
};

int main() {

    std::setvbuf(stdin, NULL, _IOLBF, 0);
    std::setvbuf(stdout, NULL, _IOLBF, 0);

	constexpr float h = 0.02f;
	constexpr float w = 1.2f;

	second_order_system<float> sys1(h, w, 0.45f);
	second_order_system<float> sys2(h, w, 0.45f);

	constexpr float ub = 5.f;
	constexpr float lb = -5.f;

	constexpr float kp = 10.f;
	constexpr float ki = 6.f;
	constexpr float kd = 5.f;
	constexpr float tc = 0.2f;

	constexpr int cont_h_mul = 1;

	simple_pid pid(cont_h_mul * h, kp, ki, kd, tc);
	simple_bounded_pid_dp bpid(cont_h_mul * h, lb, ub, kp, ki, kd, tc);

	float out1 = 0.f;
	float out2 = 0.f;

	float in1 = 0.f;
	float in2 = 0.f;
	unsigned int counter = 0;
    for (;;) {
		float v;
		std::scanf("%f", &v);

		if (counter % cont_h_mul == 0) {
			in1 = std::clamp(pid(v - out1), lb, ub);
			in2 = bpid(v, out2);
		}
		std::printf("%f, %f, %f, %f\n", in1 / 20.f, out1, in2 / 20.f, out2);
		out1 = 0.2*sys1(in1);
		out2 = 0.2*sys2(in2);

		++counter;
    }
}
