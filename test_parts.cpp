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

int main() {

    std::setvbuf(stdin, NULL, _IOLBF, 0);
    std::setvbuf(stdout, NULL, _IOLBF, 0);

    using namespace ipmoea::parts;

	constexpr float h = 0.05f;
	constexpr float w = 5.f;

	lpf_first_order<float> lpf(h, w);
	auto my_lpf = mylpf(h, w);

	auto av = make_iir_filter<float, 2, 0>([](iir_input<float, 2> in, iir_output<float, 0> out) {
			return out.z<0>() = (in.z<0>() + in.z<-1>() + in.z<-2>()) / 3;
		});

	auto av2 = make_iir_filter<float, 2, 0>([](iir_input<float, 2> in) {
			return (in.z<0>() + in.z<-1>() + in.z<-2>()) / 3;
		});
	
    for (;;) {
        float v;
        std::scanf("%f", &v);
        std::printf("%f, %f, %f, %f\n", lpf(v), my_lpf(v), av(v), av2(v));
    }
}
