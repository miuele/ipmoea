#ifndef PRAC_CNTL
#define PRAC_CNTL

#include "zsequence.h"
#include "parts.h"
#include "utils.h"

#include <limits>
#include <cmath>

namespace ipmoea {

namespace prac_cntl {

namespace params {

struct time_step { float delta_t; };

struct clamp { float lb, ub; };

struct p_gain { float k_p; };
struct i_gain { float k_i; };
struct d_gain { float k_d; };

struct d_tf { float tf; };

}

struct simple_pid_der_pv {

	constexpr simple_pid_der_pv(params::time_step dt, params::clamp bounds, params::p_gain k_p, params::i_gain k_i)
		: simple_pid_der_pv(dt, bounds, k_p, k_i, params::d_gain{0.f}, params::d_tf{1.f})
	{}

	constexpr simple_pid_der_pv(params::time_step dt, params::clamp bounds, params::p_gain k_p, params::i_gain k_i, params::d_gain k_d, params::d_tf tf)
		: k_p_(k_p.k_p), k_i_(k_i.k_i), k_d_(k_d.k_d), tf_(tf.tf), lb_(bounds.lb), ub_(bounds.ub), integ_(dt.delta_t), hpf_(dt.delta_t, 1.f / tf_),
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
		}
		if (unbounded <= lb_) {
			return lb_;
		}

		integ_ = std::move(integ);

		return unbounded;
	}

private:
	float k_p_, k_i_, k_d_;
	float tf_;
	float lb_, ub_;
	ipmoea::parts::integrator<float> integ_;
	ipmoea::parts::hpf_first_order<float> hpf_;
	float pv_0_;
};

}

}

#endif
