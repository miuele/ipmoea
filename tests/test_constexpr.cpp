#include <utility>
#include <cstdio>

#include "zsequence.h"
#include "parts.h"
#include "utils.h"

using namespace ipmoea::utils;

auto mylpf(float h, float w) {
	return make_iir_filter<float, 1, 1>(
		[h, w](iir_input<float, 1> in, iir_output<float, 1> out) {
			return out.z<0>() = ((2 - h*w)*out.z<-1>() + h*w*(in.z<0>() + in.z<-1>())) / (2 + h*w);
		});
}

template <std::size_t N>
constexpr std::array<float, N> simulate(float h, float w) {
	decltype(simulate<N>(0.f, 0.f)) result = {};
	ipmoea::parts::lpf_first_order<float> lpf(h, w);
	for (std::size_t i = 0; i < result.size(); ++i) {
		const_cast<float &>(static_cast<const decltype(result) &>(result)[i]) = lpf(1.f);
	}
	return result;
}

constexpr auto expo = simulate<16>(0.05f, 1.f);

int main() {

    std::setvbuf(stdin, NULL, _IOLBF, 0);
    std::setvbuf(stdout, NULL, _IOLBF, 0);

    using namespace ipmoea::parts;

	constexpr float h = 0.05f;
	constexpr float w = 5.f;

	lpf_first_order<float> lpf(h, w);
	
    for (std::size_t i = 0; ; ++i) {
        float v;
        std::scanf("%f", &v);
        std::printf("%f", lpf(v));

		if (i < expo.size()) {
			std::printf(", %f", expo[i]);
		}
		std::printf("\n");
    }
}
