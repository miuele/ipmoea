#include <utility>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <limits>
#include <cmath>

#include "zsequence.h"
#include "parts.h"
#include "utils.h"
#include "compositor.h"

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

int main() {

    std::setvbuf(stdin, NULL, _IOLBF, 0);
    std::setvbuf(stdout, NULL, _IOLBF, 0);

	constexpr float h = 0.02f;
	constexpr float w = 1.2f;

	second_order_system<float> sys(h, w, 0.45f);

	using namespace std::placeholders;
	auto clamper = std::bind(std::clamp<float>, _1, 0.f, 2.f);
	auto limiter = rate_limiter<float>(h, -0.1f, 0.5f);
	differentiator<float> diff(h);

	auto clampedsys = compose<float>(clamper, std::move(sys));

    for (;;) {
		float v;
		std::scanf("%f", &v);

		float out1 = clampedsys(v);
		float out2 = limiter(out1);

		printf("%f, %f\n", out1, out2);
    }
}
