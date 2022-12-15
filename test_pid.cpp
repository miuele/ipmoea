#include <utility>
#include <cstdio>

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

template <class T>
struct simple_pid_controller {

	simple_pid_controller(float time_step, float k_p, float k_i, float k_d, float cutoff)
		: k_p_(k_p), k_i_(k_i), k_d_(k_d), integ_(time_step), diff_(time_step, cutoff)
	{}

	auto operator()(T input) {
		return k_p_ * input + k_i_ * integ_(input) + k_d_ * diff_(input);
	}

private:
	T k_p_, k_i_, k_d_;
	integrator<float> integ_;
	hpf_first_order<float> diff_;
};

int main() {

    std::setvbuf(stdin, NULL, _IOLBF, 0);
    std::setvbuf(stdout, NULL, _IOLBF, 0);

	constexpr float h = 0.05f;
	constexpr float w = 1.f;

	second_order_system<float> sys(h, w, 0.75f);
	simple_pid_controller<float> pid(h, 120.f, 50.f, 0.1f, 10.f);

	float out = 0.f;
    for (;;) {
		float v;
		std::scanf("%f", &v);

		float in = pid(v - out);
		std::printf("%f\n", out);
		out = 0.2*sys(in);
    }
}
